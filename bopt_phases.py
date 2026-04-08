#!/usr/bin/env python3
"""
-----------------------------------------------------------------------
Bound-state optimizer with ESEP scan + PHASED correlation optimization.

If no ALPHA_DIR is provided, the run starts from the ctrl-default deck
(no alpha seeding).

Per ESEP scale:
  Phase 0: optimize MAIN (non-SS) correlations (ESEP + global corr)
  Phase 1: optimize FIRST SS block only (SS[0])
  Phase 2: optimize remaining SS blocks (SS[1..]) with MAIN lightly
  Phase 3: full optimize (MAIN + all SS)

-----------------------------------------------------------------------
"""
import os
import sys
import json
import shutil
from datetime import datetime
from pathlib import Path
import re
import copy
import numpy as np
from contextlib import redirect_stdout
from io import StringIO

sys.path.append(os.path.join(os.path.dirname(__file__), "src"))

from utility import utility_t
from deck import deck_t, GenerateOptFile
from wavefunction import wavefunction_t, InitNShellBoundWF
from parameters import parameters_t


def is_finite(x):
    try:
        return x is not None and np.isfinite(float(x))
    except Exception:
        return False

def strip_quotes(s):
    s = str(s).strip()
    if len(s) >= 2 and s[0] == s[-1] and s[0] in ("'", '"'):
        return s[1:-1]
    return s

def quote_path(path):
    return f"'{path}'"

def ensure_dir(p):
    if p:
        Path(p).mkdir(parents=True, exist_ok=True)

def create_pair_name(pot2b, pot3b):
    return f"{pot2b.split('.')[0]}.{pot3b.split('.')[0]}"

def scale_tag(scale):
    return f"s{int(round(scale * 100)):03d}"

def set_esep_scale(instructions, esep_scale):
    out = []
    for inst in instructions:
        d = dict(inst)
        # Only modify ESEP entries that have non-zero scale (preserve flat-only entries like idx 2)
        if d.get("key") == "ESEP" and d.get("scale", 0.0) != 0.0:
            d["scale"] = float(esep_scale)
        out.append(d)
    return out

def get_ss_list(dk):
    ss = getattr(dk, "SS", None)
    if ss is None:
        return []
    if isinstance(ss, (list, tuple)):
        return list(ss)
    try:
        return list(ss)
    except Exception:
        return []


def load_dk_unquoted(target, path):
    p = strip_quotes(path)
    target.DK.FILE_NAME = p
    target.DK.Read(target.PARAMS)
    target.CTRL.INPUT_BRA.DECK_FILE = quote_path(p)

def evaluate_current(target, log_prefix):
    try:
        e, v = target.Evaluate(True, log_prefix)
    except Exception:
        return None, None
    if not is_finite(e):
        return None, None
    return float(e), float(v) if is_finite(v) else float("nan")



def make_ctrl_paths_relative(ctrl_in, nqmcc_dir, ctrl_out):
    nqmcc_dir = os.path.normpath(strip_quotes(nqmcc_dir))
    ctrl_in = strip_quotes(ctrl_in)
    ctrl_out = strip_quotes(ctrl_out)

    with open(ctrl_in, "r") as f:
        text = f.read()

    def rewrite_token(token):
        tok = token.strip()
        quote = None
        if len(tok) >= 2 and tok[0] == tok[-1] and tok[0] in ("'", '"'):
            quote = tok[0]
            tok = tok[1:-1]

        norm = os.path.normpath(tok)

        if os.path.isabs(norm) and os.path.normpath(nqmcc_dir) in norm:
            try:
                norm = os.path.relpath(norm, os.path.normpath(nqmcc_dir))
            except Exception:
                pass

        if quote:
            return f"{quote}{norm}{quote}"
        return norm

    for key in ["CONST_FILE", "L2BP_FILE", "L3BP_FILE", "WALK_FILE", "POT_FILE"]:
        pat = re.compile(rf"({key}\s*=\s*)([^\n#]+)")
        def repl(m):
            return m.group(1) + rewrite_token(m.group(2))
        text = pat.sub(repl, text)

    with open(ctrl_out, "w") as f:
        f.write(text)


def setup_main_directory(util):
    try:
        os.mkdir(util.WORKING_DIR)
    except FileExistsError:
        util.WORKING_DIR = f"{util.WORKING_DIR}{util.NAME}-{datetime.now().strftime('%Y-%m-%d_%H-%M')}/"
        os.mkdir(util.WORKING_DIR)

    ensure_dir(f"{util.NQMCC_DIR}walks")

def setup_potential_pair_directory(util, pair_name):
    pot_dir = f"{util.WORKING_DIR}{pair_name}/"
    ensure_dir(pot_dir)
    ensure_dir(f"{pot_dir}ctrl")
    ensure_dir(f"{pot_dir}logs")
    ensure_dir(f"{pot_dir}dk")
    ensure_dir(f"{pot_dir}opt")
    ensure_dir(f"{util.NQMCC_DIR}walks")
    return pot_dir

def build_target(util, pot_dir, pot_index):
    base_ctrl_in = strip_quotes(util.CTRL_FILE)
    patched_ctrl = f"{pot_dir}ctrl/base.ctrl"
    make_ctrl_paths_relative(base_ctrl_in, util.NQMCC_DIR, patched_ctrl)

    target = wavefunction_t(patched_ctrl, util.NQMCC_DIR, util.BIN_DIR, util.RUN_CMD)

    target.CTRL.FILE_NAME = f"{pot_dir}ctrl/target.ctrl"
    target.CTRL.NUM_BLOCKS = util.NUM_BLOCKS
    target.CTRL.BLOCK_SIZE = util.BLOCK_SIZE
    target.CTRL.WALKERS_PER_NODE = util.WALKERS_PER_NODE
    target.CTRL.NUM_OPT_EVALUATIONS = util.NUM_OPT_EVALUATIONS

    target.CTRL.CONST_FILE = quote_path(f"{util.NQMCC_DIR}constants/{util.CONSTANTS_FILES[pot_index]}")
    target.CTRL.L2BP_FILE  = quote_path(f"{util.NQMCC_DIR}pots/{util.TWO_BODY_FILES[pot_index]}")
    target.CTRL.L3BP_FILE  = quote_path(f"{util.NQMCC_DIR}pots/{util.THREE_BODY_FILES[pot_index]}")
    return target


def find_alpha_core(alpha_cores_dir, pot_pair_name):
    deck_name = f"{pot_pair_name}_optimized.dk"
    path = os.path.join(alpha_cores_dir, pot_pair_name, "dk", deck_name)
    return path if os.path.exists(path) else None

def generate_alpha_seeded_deck(target, pot_pair_name, output_path, alpha_cores_dir, he4_params_path, copy_esep=True):
    alpha_path = find_alpha_core(alpha_cores_dir, pot_pair_name)
    he4_params = parameters_t(he4_params_path)
    alpha_deck = deck_t(he4_params, alpha_path)

    class MinimalAlpha:
        def __init__(self, dk):
            self.DK = dk

    alpha = MinimalAlpha(alpha_deck)
    ensure_dir(os.path.dirname(output_path))
    InitNShellBoundWF(target, alpha, output_path, copy_esep=copy_esep)
    return output_path


def try_optimize(target, opt_file, dk_out, log_prefix, instructions, e_ref, require_improvement=True):
    opt_file = strip_quotes(opt_file)
    dk_out   = strip_quotes(dk_out)
    log_prefix = strip_quotes(log_prefix)

    opt = GenerateOptFile(target.PARAMS, target.DK, quote_path(opt_file), instructions)
    opt.UpdateFloats(target.PARAMS, 6)

    log_txt = f"{log_prefix}.txt"
    buf = StringIO()
    try:
        with redirect_stdout(buf):
            e, v = target.Optimize(opt, quote_path(dk_out), True, log_prefix)

    except Exception as ex:
        with open(log_txt, "w") as f:
            f.write(buf.getvalue())
            f.write(f"\nEXCEPTION: {ex}\n")
        return None, None, False, f"Optimize exception: {ex}"

    with open(log_txt, "w") as f:
        f.write(buf.getvalue())

    if not is_finite(e):
        return None, None, False, "Non-finite energy"
    if not os.path.exists(dk_out):
        return None, None, False, f"No DK output: {dk_out}"

    e = float(e)
    v = float(v) if is_finite(v) else float("nan")

    if require_improvement:
        dE = e - float(e_ref)
        if dE >= 0.0:
            return None, None, False, f"No improvement: ΔE={dE:.6f}"

    return e, v, True, ""


PHASE_SCALES = {
    "p0_main"      : 0.50,
    "p1_ss0"       : 1.00,
    "p2_ss_rest"   : 0.50,
    "p3_full"      : 1.00,
}

SS_KEYS = [
    "SPU", "SPV", "SPR", "SPA","SPB","SPC","SPK","SPL",
    "PPU", "PPV", "PPR", "PPA", "PPB", "PPC", "PPK", "PPL",
    "WSE", "WSV", "WSR", "WSA", "WBRHO", "WBALPH"
]

def build_instructions_base(esep_scale, opt_scale, opt_3b, eta_flat=0.0, target_dk=None):
    inst = [
        {"ss": False, "key": "ESEP", "idx": 0, "scale": esep_scale, "flat": 0.0},
        {"ss": False, "key": "ESEP", "idx": 1, "scale": esep_scale, "flat": 0.0},
        {"ss": False, "key": "ESEP", "idx": 2, "scale": esep_scale, "flat": 0.0},
        {"ss": False, "key": "ESEP", "idx": 3, "scale": esep_scale, "flat": 0.0},

        {"ss": False, "key": "ETA",  "all": True, "scale": 0.0, "flat": float(eta_flat)},
        {"ss": False, "key": "ZETA",  "all": True, "scale": opt_scale, "flat": 0.0},
        {"ss": False, "key": "FSCAL", "all": True, "scale": opt_scale, "flat": 0.0},
        {"ss": False, "key": "ALPHA", "all": True, "scale": opt_scale, "flat": 0.0},
        {"ss": False, "key": "BETA",  "all": True, "scale": opt_scale, "flat": 0.0},

        {"ss": False, "key": "AC",  "all": True, "scale": opt_scale, "flat": 0.0},
        {"ss": False, "key": "AA",  "all": True, "scale": opt_scale, "flat": 0.0},
        {"ss": False, "key": "AR",  "all": True, "scale": opt_scale, "flat": 0.0},

        {"ss": False, "key": "DELTA",   "scale": 0, "flat": opt_3b},
        {"ss": False, "key": "EPSILON", "scale": 0, "flat": opt_3b},
        {"ss": False, "key": "THETA",   "scale": 0, "flat": opt_3b},
        {"ss": False, "key": "UPSILON", "scale": 0, "flat": opt_3b},
        {"ss": False, "key": "RSCAL",   "scale": 0, "flat": opt_3b},
        {"ss": False, "key": "USCAL",   "scale": 0, "flat": opt_3b},
    ]

    q_params = [
        "QPS1", "QPS2", "QSSS1", "QSSS2",
        "QSSP1", "QSSP2",
        "QSPP1", "QSPP2",
        "QPPP1", "QPPP2",
        "QSPD1", "QSPD2", "QSSD1", "QSSD2", "QPPD1", "QPPD2",
        "QSDD1", "QSDD2", "QPDD1", "QPDD2",
        "QDDD1", "QDDD2",
    ]
    
    for q_key in q_params:
        if target_dk is None or hasattr(target_dk, q_key):
            inst.append({"ss": False, "key": q_key, "scale": opt_scale, "flat": 0.0})

    return inst

def build_instructions_ss(ss_indices, flat_scale):
    inst = []
    for ssi in ss_indices:
        inst += [
            {"ss": True, "ss_idx": ssi, "key": "SPU", "scale": 0.0, "flat": flat_scale},
            {"ss": True, "ss_idx": ssi, "key": "SPV", "scale": 0.0, "flat": flat_scale},
            {"ss": True, "ss_idx": ssi, "key": "SPR", "scale": 0.0, "flat": flat_scale},
            {"ss": True, "ss_idx": ssi, "key": "SPA", "scale": 0.0, "flat": flat_scale},
            {"ss": True, "ss_idx": ssi, "key": "SPB", "scale": 0.0, "flat": flat_scale},
            {"ss": True, "ss_idx": ssi, "key": "SPC", "scale": 0.0, "flat": flat_scale},
            {"ss": True, "ss_idx": ssi, "key": "SPK", "scale": 0.0, "flat": flat_scale},
            {"ss": True, "ss_idx": ssi, "key": "SPL", "scale": 0.0, "flat": flat_scale},

            {"ss": True, "ss_idx": ssi, "key": "PPU", "scale": 0.0, "flat": flat_scale},
            {"ss": True, "ss_idx": ssi, "key": "PPV", "scale": 0.0, "flat": flat_scale},
            {"ss": True, "ss_idx": ssi, "key": "PPR", "scale": 0.0, "flat": flat_scale},
            {"ss": True, "ss_idx": ssi, "key": "PPA", "scale": 0.0, "flat": flat_scale},
            {"ss": True, "ss_idx": ssi, "key": "PPB", "scale": 0.0, "flat": flat_scale},
            {"ss": True, "ss_idx": ssi, "key": "PPC", "scale": 0.0, "flat": flat_scale},
            {"ss": True, "ss_idx": ssi, "key": "PPK", "scale": 0.0, "flat": flat_scale},
            {"ss": True, "ss_idx": ssi, "key": "PPL", "scale": 0.0, "flat": flat_scale},

            {"ss": True, "ss_idx": ssi, "key": "WSE", "scale": 0.0, "flat": flat_scale},
            {"ss": True, "ss_idx": ssi, "key": "WSV", "scale": 0.0, "flat": flat_scale},
            {"ss": True, "ss_idx": ssi, "key": "WSR", "scale": 0.0, "flat": flat_scale},
            {"ss": True, "ss_idx": ssi, "key": "WSA", "scale": 0.0, "flat": flat_scale},
            {"ss": True, "ss_idx": ssi, "key": "WBRHO", "scale": 0.0, "flat": flat_scale},
            {"ss": True, "ss_idx": ssi, "key": "WBALPH", "scale": 0.0, "flat": flat_scale},
        ]
    return inst

def copy_ss0_to_rest(target):
    ss_list = get_ss_list(target.DK)
    if len(ss_list) <= 1:
        return False

    src = ss_list[0]
    for ssi in range(1, len(ss_list)):
        dst = ss_list[ssi]
        for k in SS_KEYS:
            if hasattr(src, k) and hasattr(dst, k):
                try:
                    setattr(dst, k, copy.deepcopy(getattr(src, k)))
                except Exception:
                    setattr(dst, k, getattr(src, k))
    return True


def BoundStateOptimizePhases(util, pair_name, pot_dir, pot_index, step=0.2, start_scale=0.0):
    print("=" * 72)
    print(f"{pair_name}: building target + optional alpha-seeding + initial Evaluate")
    print("=" * 72)

    target0 = build_target(util, pot_dir, pot_index)

    starting_deck_label = "ctrl-default"
    if getattr(util, "ALPHA_DIR", None):
        he4_params_path = f"{util.NQMCC_DIR}nuclei/params/he4.params"
        seeded_deck_path = f"{pot_dir}dk/{pair_name}_alpha_seed.dk"
        try:
            generate_alpha_seeded_deck(target0, pair_name, seeded_deck_path, util.ALPHA_DIR, he4_params_path, copy_esep=True)
            target0.DK = deck_t(target0.PARAMS, seeded_deck_path)
            target0.CTRL.INPUT_BRA.DECK_FILE = quote_path(seeded_deck_path)
            starting_deck_label = "alpha-seeded"
            print(f"Alpha seeding successful: {seeded_deck_path}")
        except Exception as ex:
            print(f"Alpha seeding failed ({ex}), continuing with unseeded version")
            starting_deck_label = "ctrl-default"

    log_init = f"{pot_dir}logs/{pair_name}.initial"
    e0, v0 = target0.Evaluate(True, log_init)
    if not is_finite(e0):
        raise RuntimeError("Initial evaluation failed")
    e0 = float(e0)
    v0 = float(v0) if is_finite(v0) else float("nan")
    print(f"INITIAL ({starting_deck_label}): E={e0:.6f}  V={v0:.6f}")

    dk_backup = f"{pot_dir}dk/{pair_name}_original.dk"
    target0.DK.Write(target0.PARAMS, quote_path(dk_backup))
    print(f"Saved ORIGINAL DK: {dk_backup}")

    esep_scale_max = float(util.ESEP_SCALE)
    opt_scale_base = float(util.OPT_SCALE)
    eta_flat = float(util.ETA_FLAT)
    three_body_flat = float(util.THREE_BODY_FLAT)
    ss_flat = float(util.SS_FLAT)

    n_ss0 = len(get_ss_list(target0.DK))
    if n_ss0 > 0:
        print(f"Spatial symmetries detected: SS blocks = {n_ss0}")
    else:
        print("No SS blocks detected: running MAIN-only phases")

    scales = []
    x = float(start_scale)
    while x <= esep_scale_max + 1e-12:
        scales.append(round(x, 10))
        x += float(step)
    if scales and scales[-1] < esep_scale_max - 1e-12:
        scales.append(esep_scale_max)

    print(f"ESEP scan scales: {scales}")

    attempts = []
    best = None

    for scale in scales:
        tag = scale_tag(scale)
        print(f"\n{pair_name}: TRY scale={scale} (tag={tag})")

        target = build_target(util, pot_dir, pot_index)
        load_dk_unquoted(target, dk_backup)

        n_ss = len(get_ss_list(target.DK))

        best_e_scale = e0
        best_v_scale = v0
        best_dk_scale = dk_backup
        best_opt_scale = None

        def run_phase(phase_name, instructions):
            nonlocal best_e_scale, best_v_scale, best_dk_scale, best_opt_scale
            opt_file = f"{pot_dir}opt/{pair_name}_opt_{tag}_{phase_name}.opt"
            dk_out   = f"{pot_dir}dk/{pair_name}_opt_{tag}_{phase_name}.dk"
            log_pref = f"{pot_dir}logs/{pair_name}.opt_{tag}_{phase_name}"

            e, v, ok, msg = try_optimize(
                target,
                opt_file,
                dk_out,
                log_pref,
                set_esep_scale(instructions, scale),
                best_e_scale,
                require_improvement=True
            )

            attempts.append({
                "scale": scale, "tag": tag, "phase": phase_name, "ok": ok,
                "energy": e if ok else None, "variance": v if ok else None,
                "dk": dk_out, "opt": opt_file, "msg": msg if not ok else "",
            })

            if ok:
                best_e_scale = float(e)
                best_v_scale = float(v) if is_finite(v) else float("nan")
                best_dk_scale = dk_out
                best_opt_scale = opt_file
                load_dk_unquoted(target, best_dk_scale)
                print(f"  ✅ {phase_name}: E={best_e_scale:.6f}  ΔE={(best_e_scale-e0):.6f}")
                return True
            else:
                print(f"  ⚠️ {phase_name} failed: {msg}")
                load_dk_unquoted(target, best_dk_scale)
                return False

        opt_scale_p0 = opt_scale_base * PHASE_SCALES["p0_main"]
        three_body_flat_p0 = three_body_flat * PHASE_SCALES["p0_main"]
        inst_p0 = build_instructions_base(esep_scale_max, opt_scale_p0, three_body_flat_p0, eta_flat, target.DK)
        ok0 = run_phase("p0_main", inst_p0)

        if n_ss <= 0:
            if best is None or best_e_scale < best["energy"]:
                best = {"scale": scale, "tag": tag, "energy": best_e_scale, "variance": best_v_scale, "dk": best_dk_scale, "opt": best_opt_scale}
            continue

        s0_flat_p1 = ss_flat * PHASE_SCALES["p1_ss0"]
        inst_ss0 = build_instructions_ss([0], flat_scale=s0_flat_p1)
        ok1 = run_phase("p1_ss0", inst_ss0)
        
        
        rest = list(range(1, n_ss))
        if rest:
            ss_flat_p2 = ss_flat * PHASE_SCALES["p2_ss_rest"]
            opt_scale_p2 = opt_scale_base * PHASE_SCALES["p2_ss_rest"]
            esep_scale_p2 = esep_scale_max * PHASE_SCALES["p2_ss_rest"]
            three_body_flat_p2 = three_body_flat * PHASE_SCALES["p2_ss_rest"]
            inst_base_p2 = build_instructions_base(esep_scale_p2, opt_scale_p2, three_body_flat_p2, eta_flat, target.DK)
            inst_ss_rest = build_instructions_ss(rest, flat_scale=ss_flat_p2)
            ok2 = run_phase("p2_ss_rest", inst_base_p2 + inst_ss_rest)

        opt_scale_p3 = opt_scale_base * PHASE_SCALES["p3_full"]
        ss_flat_p3 = ss_flat * PHASE_SCALES["p3_full"]
        three_body_flat_p3 = three_body_flat * PHASE_SCALES["p3_full"]
        inst_base_p3 = build_instructions_base(esep_scale_max, opt_scale_p3, three_body_flat_p3, eta_flat, target.DK)
        inst_ss_all = build_instructions_ss(list(range(n_ss)), flat_scale=ss_flat_p3)
        ok3 = run_phase("p3_full", inst_base_p3 + inst_ss_all)

        if best is None or best_e_scale < best["energy"]:
            best = {"scale": scale, "tag": tag, "energy": best_e_scale, "variance": best_v_scale, "dk": best_dk_scale, "opt": best_opt_scale}

    final_dk  = f"{pot_dir}dk/{pair_name}_optimized.dk"
    final_opt = f"{pot_dir}opt/{pair_name}_opt.opt"

    if best is None:
        print("No successful improvements; keeping ORIGINAL as optimized")
        shutil.copy2(dk_backup, final_dk)
        best_e, best_v, best_scale = e0, v0, 0.0
    else:
        shutil.copy2(best["dk"], final_dk)
        if best.get("opt") and os.path.exists(best["opt"]):
            shutil.copy2(best["opt"], final_opt)
        best_e, best_v, best_scale = best["energy"], best["variance"], best["scale"]

    result = {
        "POTENTIAL_PAIR": pair_name,
        "STARTING_DECK_MODE": starting_deck_label,
        "ALPHA_CORES_DIR": getattr(util, "ALPHA_DIR", None),
        "E_INITIAL": e0,
        "V_INITIAL": v0,
        "E_OPTIMIZED": best_e,
        "V_OPTIMIZED": best_v,
        "DELTA_E": best_e - e0,
        "BEST_ESEP_SCALE": best_scale,
        "FINAL_DK": final_dk,
        "FINAL_OPT": final_opt,
        "DK_BACKUP": dk_backup,
        "ATTEMPTS": attempts,
    }
    with open(f"{pot_dir}{pair_name}_results.json", "w") as f:
        json.dump(result, f, indent=2)

    return result


def run(util, step=0.2, start_scale=0.0):
    setup_main_directory(util)

    # Extract radius settings from control file
    pair_name_temp = create_pair_name(util.TWO_BODY_FILES[0], util.THREE_BODY_FILES[0])
    pot_dir_temp = setup_potential_pair_directory(util, pair_name_temp)
    target_temp = build_target(util, pot_dir_temp, 0)
    
    limit_radii_str = str(getattr(target_temp.CTRL, "LIMIT_CHARGE_RADII", ".false.")).strip()
    limit_radii = limit_radii_str.lower() == ".true."
    
    radius_settings = {}
    if limit_radii:
        radius_settings = {
            "limit_charge_radii": True,
            "neutron_radius_limit": float(getattr(target_temp.CTRL, "NEUTRON_RADIUS_LIMIT", None)),
            "proton_radius_limit": float(getattr(target_temp.CTRL, "PROTON_RADIUS_LIMIT", None))
        }
        print(f"LIMIT_CHARGE_RADII: True")
        print(f"NEUTRON_RADIUS_LIMIT: {radius_settings['neutron_radius_limit']}")
        print(f"PROTON_RADIUS_LIMIT: {radius_settings['proton_radius_limit']}")
    else:
        radius_settings = {
            "limit_charge_radii": False
        }
        print(f"LIMIT_CHARGE_RADII: False")

    metadata = {
        "timestamp": datetime.now().isoformat(),
        "optimization_settings": {
            "opt_scale": util.OPT_SCALE,
            "num_opt_evaluations": util.NUM_OPT_EVALUATIONS,
            "esep_scale": float(util.ESEP_SCALE),
            "step_size": step,
            "start_scale": start_scale,
            "phase_scales": PHASE_SCALES,
            "eta_flat": util.ETA_FLAT,
            "three_body_flat": util.THREE_BODY_FLAT,
            "ss_flat": util.SS_FLAT,
        },
        "radius_settings": radius_settings
    }

    # Setup output capture
    opt_output_lines = []
    
    class OutputCapture:
        def __init__(self):
            self.original_print = print
            
        def captured_print(self, *args, **kwargs):
            # Capture to our list
            line = ' '.join(str(arg) for arg in args)
            opt_output_lines.append(line)
            # Still print to console using original print
            self.original_print(*args, **kwargs)
    
    # Create output capture instance
    output_capture = OutputCapture()
    
    # Monkey patch print in builtins and this module
    import builtins
    original_builtins_print = builtins.print
    builtins.print = output_capture.captured_print
    
    # Also patch print in global scope of this module
    globals()['print'] = output_capture.captured_print

    try:
        results = []
        for i in range(len(util.TWO_BODY_FILES)):
            pair_name = create_pair_name(util.TWO_BODY_FILES[i], util.THREE_BODY_FILES[i])
            pot_dir = setup_potential_pair_directory(util, pair_name)

            print("\n" + "=" * 72)
            print(f"PROCESSING {i+1}/{len(util.TWO_BODY_FILES)}: {pair_name}")
            print("=" * 72)

            try:
                res = BoundStateOptimizePhases(util, pair_name, pot_dir, i, step=step, start_scale=start_scale)
                results.append(res)
            except Exception as ex:
                print(f"ERROR optimizing {pair_name}: {ex}")
                results.append({"POTENTIAL_PAIR": pair_name, "ERROR": str(ex)})

    finally:
        # Restore original print functions
        builtins.print = original_builtins_print
        globals()['print'] = original_builtins_print

    # Prepare combined results with metadata (excluding attempts from individual results)
    processed_results = []
    for result in results:
        # Create a copy without attempts for combined results
        clean_result = {k: v for k, v in result.items() if k != "ATTEMPTS"}
        processed_results.append(clean_result)

    combined_data = {
        "run_metadata": metadata,
        "results": processed_results
    }

    out = f"{util.WORKING_DIR}combined_results.json"
    with open(out, "w") as f:
        json.dump(combined_data, f, indent=2)
    print("\nSaved:", out)

    # Save all captured output to .out file
    out_file = f"{util.WORKING_DIR}optimization_log.out"
    with open(out_file, "w") as f:
        f.write("# Optimization Log\n")
        f.write(f"# Generated: {datetime.now().isoformat()}\n")
        f.write(f"# OPT_SCALE: {util.OPT_SCALE}\n")
        f.write(f"# NUM_OPT_EVALUATIONS: {util.NUM_OPT_EVALUATIONS}\n")
        f.write(f"# ESEP_SCALE: {util.ESEP_SCALE}\n")
        f.write(f"# STEP_SIZE: {step}\n")
        f.write(f"# START_SCALE: {start_scale}\n")
        f.write(f"# PHASE_SCALES: {PHASE_SCALES}\n")
        f.write(f"# ETA_FLAT: {util.ETA_FLAT}\n")
        f.write(f"# THREE_BODY_FLAT: {util.THREE_BODY_FLAT}\n")
        f.write(f"# SS_FLAT: {util.SS_FLAT}\n")
        f.write(f"# LIMIT_CHARGE_RADII: {radius_settings['limit_charge_radii']}\n")
        if radius_settings['limit_charge_radii']:
            f.write(f"# NEUTRON_RADIUS_LIMIT: {radius_settings['neutron_radius_limit']}\n")
            f.write(f"# PROTON_RADIUS_LIMIT: {radius_settings['proton_radius_limit']}\n")
        f.write("#" + "="*70 + "\n\n")
        for line in opt_output_lines:
            f.write(line + "\n")
    print("Saved optimization log:", out_file)

    return results


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--utility", required=True)
    p.add_argument("--step", type=float, default=0.2)
    p.add_argument("--start", type=float, default=0.0)
    args = p.parse_args()

    util = utility_t(args.utility)
    run(util, step=args.step, start_scale=args.start)
