#!/usr/bin/env python3
"""
Two-stage bound-state optimizer with adaptive ESEP scan.

Stage 1 — ESEP basin search:
  Scan over ESEP scale values running Phase 0 (main correlations only).
  Each scale is retried with shrinking step sizes if it fails,
  so a good basin isn't discarded just because the initial opt_scale was too large.
  The single best ESEP scale is carried forward to Stage 2.

Stage 2 — structural fine-tuning at the winning scale:
  Phase 1: first SS block only
  Phase 2: remaining SS blocks + main correlations (if >1 SS block)
  Phase 3: full optimization (main + all SS)

  Each phase runs a grid search over multipliers [0.5×, 1.0×, 1.5×, 2.0×]
  applied to opt_scale and flat_scale simultaneously, keeping the best result.
  The winning multiplier from Phase 1 rescales effective_ss_flat so Phase 2/3
  start from a calibrated baseline.
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

sys.path.append(os.path.dirname(__file__))

from utility import utility_t
from deck import deck_t, GenerateOptFile
from wavefunction import wavefunction_t, InitNShellBoundWF
from parameters import parameters_t


# ---------------------------------------------------------------------------
# Shared utility helpers
# ---------------------------------------------------------------------------

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
        if d.get("key") == "ESEP" and d.get("scale", 0.0) != 0.0:
            d["scale"] = float(esep_scale)
        out.append(d)
    return out

def get_ss_list(dk):
    ss = getattr(dk, "SS", None)
    if ss is None:
        return []
    try:
        return list(ss)
    except Exception:
        return []

def load_dk_unquoted(target, path):
    p = strip_quotes(path)
    target.DK.FILE_NAME = p
    target.DK.Read(target.PARAMS)
    target.CTRL.INPUT_BRA.DECK_FILE = quote_path(p)


# ---------------------------------------------------------------------------
# Directory and target setup
# ---------------------------------------------------------------------------

def make_ctrl_paths_relative(ctrl_in, nqmcc_dir, ctrl_out):
    nqmcc_dir = os.path.normpath(strip_quotes(nqmcc_dir))
    ctrl_in   = strip_quotes(ctrl_in)
    ctrl_out  = strip_quotes(ctrl_out)

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
        util.WORKING_DIR = f"{util.WORKING_DIR}{util.NAME}-dynamic-{datetime.now().strftime('%Y-%m-%d_%H-%M')}/"
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
    """Build a fresh wavefunction_t for the given potential pair."""
    base_ctrl_in = strip_quotes(util.CTRL_FILE)
    patched_ctrl = f"{pot_dir}ctrl/base.ctrl"
    make_ctrl_paths_relative(base_ctrl_in, util.NQMCC_DIR, patched_ctrl)

    target = wavefunction_t(patched_ctrl, util.NQMCC_DIR, util.BIN_DIR, util.RUN_CMD)
    target.CTRL.FILE_NAME          = f"{pot_dir}ctrl/target.ctrl"
    target.CTRL.NUM_BLOCKS         = util.NUM_BLOCKS
    target.CTRL.BLOCK_SIZE         = util.BLOCK_SIZE
    target.CTRL.WALKERS_PER_NODE   = util.WALKERS_PER_NODE
    target.CTRL.NUM_OPT_EVALUATIONS = util.NUM_OPT_EVALUATIONS

    target.CTRL.CONST_FILE = quote_path(f"{util.NQMCC_DIR}constants/{util.CONSTANTS_FILES[pot_index]}")
    target.CTRL.L2BP_FILE  = quote_path(f"{util.NQMCC_DIR}pots/{util.TWO_BODY_FILES[pot_index]}")
    target.CTRL.L3BP_FILE  = quote_path(f"{util.NQMCC_DIR}pots/{util.THREE_BODY_FILES[pot_index]}")
    return target


# ---------------------------------------------------------------------------
# Alpha seeding
# ---------------------------------------------------------------------------

def find_alpha_core(alpha_cores_dir, pot_pair_name):
    """Return the path to the optimized He-4 alpha deck, or None if absent."""
    deck_name = f"{pot_pair_name}_optimized.dk"
    path = os.path.join(alpha_cores_dir, pot_pair_name, "dk", deck_name)
    return path if os.path.exists(path) else None

def generate_alpha_seeded_deck(target, pot_pair_name, output_path, alpha_cores_dir, he4_params_path, copy_esep=True):
    """
    Build a starting deck by copying optimized He-4 alpha-core parameters
    into the target deck.
    """
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


# ---------------------------------------------------------------------------
# Core optimization call
# ---------------------------------------------------------------------------

def try_optimize(target, opt_file, dk_out, log_prefix, instructions, e_ref, require_improvement=True):
    """
    Run one optimization pass and return (e, v, ok, msg).

    Writes stdout from the optimizer to log_prefix.txt.
    Returns ok=False if the energy is non-finite, the output deck is missing,
    or (when require_improvement=True) the energy did not decrease.
    """
    opt_file   = strip_quotes(opt_file)
    dk_out     = strip_quotes(dk_out)
    log_prefix = strip_quotes(log_prefix)

    ensure_dir(os.path.dirname(opt_file))
    ensure_dir(os.path.dirname(dk_out))
    ensure_dir(os.path.dirname(log_prefix))

    opt = GenerateOptFile(target.PARAMS, target.DK, quote_path(opt_file), instructions)
    opt.UpdateFloats(target.PARAMS, 6)

    log_txt = f"{log_prefix}.txt"
    buf = StringIO()
    try:
        with redirect_stdout(buf):
            e, v = target.Optimize(opt, quote_path(dk_out), True, log_prefix)
    except FileNotFoundError as ex:
        with open(log_txt, "w") as f:
            f.write(buf.getvalue())
            f.write(f"\nEXCEPTION: {ex}\n")
        if dk_out in str(ex):
            return None, None, False, f"No DK output: {dk_out}"
        return None, None, False, f"FileNotFoundError: {ex}"
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


# ---------------------------------------------------------------------------
# Instruction builders
# ---------------------------------------------------------------------------

# Fractional scales applied to opt_scale/flat_scale per phase.
# Stage 1 uses p0_main; Stage 2 uses the rest.
PHASE_SCALES = {
    "p0_main"    : 0.50,
    "p1_ss0"     : 1.00,
    "p2_ss_rest" : 0.50,
    "p3_full"    : 1.00,
}

SS_KEYS = [
    "SPU", "SPV", "SPR", "SPA", "SPB", "SPC", "SPK", "SPL",
    "PPU", "PPV", "PPR", "PPA", "PPB", "PPC", "PPK", "PPL",
    "WSE", "WSV", "WSR", "WSA", "WBRHO", "WBALPH",
]

def build_instructions_base(esep_scale, opt_scale, opt_3b, eta_flat=0.0, target_dk=None):
    """
    Build the main-block instruction list (ESEP + global correlations).
    """
    inst = [
        {"ss": False, "key": "ESEP", "idx": 0, "scale": esep_scale, "flat": 0.0},
        {"ss": False, "key": "ESEP", "idx": 1, "scale": esep_scale, "flat": 0.0},
        {"ss": False, "key": "ESEP", "idx": 2, "scale": esep_scale, "flat": 0.0},
        {"ss": False, "key": "ESEP", "idx": 3, "scale": esep_scale, "flat": 0.0},

        {"ss": False, "key": "ETA",   "all": True, "scale": 0.0,       "flat": float(eta_flat)},
        {"ss": False, "key": "ZETA",  "all": True, "scale": opt_scale,  "flat": 0.0},
        {"ss": False, "key": "FSCAL", "all": True, "scale": opt_scale,  "flat": 0.0},
        {"ss": False, "key": "ALPHA", "all": True, "scale": opt_scale,  "flat": 0.0},
        {"ss": False, "key": "BETA",  "all": True, "scale": opt_scale,  "flat": 0.0},

        {"ss": False, "key": "AC",    "all": True, "scale": opt_scale,  "flat": 0.0},
        {"ss": False, "key": "AA",    "all": True, "scale": opt_scale,  "flat": 0.0},
        {"ss": False, "key": "AR",    "all": True, "scale": opt_scale,  "flat": 0.0},

        {"ss": False, "key": "DELTA",   "scale": 0, "flat": opt_3b},
        {"ss": False, "key": "EPSILON", "scale": 0, "flat": opt_3b},
        {"ss": False, "key": "THETA",   "scale": 0, "flat": opt_3b},
        {"ss": False, "key": "UPSILON", "scale": 0, "flat": opt_3b},
        {"ss": False, "key": "RSCAL",   "scale": 0, "flat": opt_3b},
        {"ss": False, "key": "USCAL",   "scale": 0, "flat": opt_3b},
    ]

    q_params = [
        "QPS1",  "QPS2",  "QSSS1", "QSSS2",
        "QSSP1", "QSSP2", "QSPP1", "QSPP2", "QPPP1", "QPPP2",
        "QSPD1", "QSPD2", "QSSD1", "QSSD2", "QPPD1", "QPPD2",
        "QSDD1", "QSDD2", "QPDD1", "QPDD2", "QDDD1", "QDDD2",
    ]
    for q_key in q_params:
        if target_dk is None or hasattr(target_dk, q_key):
            inst.append({"ss": False, "key": q_key, "scale": opt_scale, "flat": 0.0})

    return inst

def build_instructions_ss(ss_indices, flat_value):
    """Build SS-block instruction list for the given SS indices."""
    inst = []
    for ssi in ss_indices:
        inst += [
            {"ss": True, "ss_idx": ssi, "key": "SPU",    "scale": 0.0, "flat": flat_value},
            {"ss": True, "ss_idx": ssi, "key": "SPV",    "scale": 0.0, "flat": flat_value},
            {"ss": True, "ss_idx": ssi, "key": "SPR",    "scale": 0.0, "flat": flat_value},
            {"ss": True, "ss_idx": ssi, "key": "SPA",    "scale": 0.0, "flat": flat_value},
            {"ss": True, "ss_idx": ssi, "key": "SPB",    "scale": 0.0, "flat": flat_value},
            {"ss": True, "ss_idx": ssi, "key": "SPC",    "scale": 0.0, "flat": flat_value},
            {"ss": True, "ss_idx": ssi, "key": "SPK",    "scale": 0.0, "flat": flat_value},
            {"ss": True, "ss_idx": ssi, "key": "SPL",    "scale": 0.0, "flat": flat_value},

            {"ss": True, "ss_idx": ssi, "key": "PPU",    "scale": 0.0, "flat": flat_value},
            {"ss": True, "ss_idx": ssi, "key": "PPV",    "scale": 0.0, "flat": flat_value},
            {"ss": True, "ss_idx": ssi, "key": "PPR",    "scale": 0.0, "flat": flat_value},
            {"ss": True, "ss_idx": ssi, "key": "PPA",    "scale": 0.0, "flat": flat_value},
            {"ss": True, "ss_idx": ssi, "key": "PPB",    "scale": 0.0, "flat": flat_value},
            {"ss": True, "ss_idx": ssi, "key": "PPC",    "scale": 0.0, "flat": flat_value},
            {"ss": True, "ss_idx": ssi, "key": "PPK",    "scale": 0.0, "flat": flat_value},
            {"ss": True, "ss_idx": ssi, "key": "PPL",    "scale": 0.0, "flat": flat_value},

            {"ss": True, "ss_idx": ssi, "key": "WSE",    "scale": 0.0, "flat": flat_value},
            {"ss": True, "ss_idx": ssi, "key": "WSV",    "scale": 0.0, "flat": flat_value},
            {"ss": True, "ss_idx": ssi, "key": "WSR",    "scale": 0.0, "flat": flat_value},
            {"ss": True, "ss_idx": ssi, "key": "WSA",    "scale": 0.0, "flat": flat_value},
            {"ss": True, "ss_idx": ssi, "key": "WBRHO",  "scale": 0.0, "flat": flat_value},
            {"ss": True, "ss_idx": ssi, "key": "WBALPH", "scale": 0.0, "flat": flat_value},
        ]
    return inst


# ---------------------------------------------------------------------------
# Phase execution helpers
# ---------------------------------------------------------------------------

def run_adaptive_phase(target_obj, phase_name, phase_tag, start_e, start_dk,
                       max_retries, inst_builder_fn, pot_dir, pair_name, attempts):
    """
    Run a phase with exponential backoff on opt/flat scales.

    On each failure, both multipliers are halved and the safe start_dk is
    reloaded before retrying — this prevents a blown-up deck from poisoning
    the next attempt. Returns on the first success or after max_retries.
    """
    opt_m = 1.0
    flat_m = 1.0

    for retry in range(max_retries):
        load_dk_unquoted(target_obj, start_dk)

        inst     = inst_builder_fn(opt_m, flat_m)
        tag      = f"{phase_tag}_r{retry}" if retry > 0 else phase_tag
        opt_file = f"{pot_dir}opt/{pair_name}_opt_{tag}_{phase_name}.opt"
        dk_out   = f"{pot_dir}dk/{pair_name}_opt_{tag}_{phase_name}.dk"
        log_pref = f"{pot_dir}logs/{pair_name}.opt_{tag}_{phase_name}"

        e, v, ok, msg = try_optimize(
            target_obj, opt_file, dk_out, log_pref, inst, start_e, require_improvement=True
        )

        attempts.append({
            "phase": phase_name, "tag": tag, "retry": retry,
            "opt_mult": opt_m, "flat_mult": flat_m,
            "ok": ok, "energy": e if ok else None, "variance": v if ok else None,
            "dk": dk_out, "opt": opt_file, "msg": msg if not ok else "",
        })

        if ok:
            print(f"  ✅ {phase_name} (retry={retry}, scale={opt_m:.3f}x): E={e:.6f}")
            return True, float(e), float(v) if is_finite(v) else float("nan"), dk_out, opt_file

        print(f"  ⚠️  {phase_name} (retry={retry}) failed: {msg} — halving scales")
        opt_m  *= 0.5
        flat_m *= 0.5

    return False, start_e, float("nan"), None, None


def run_grid_phase(target_obj, phase_name, base_tag, e_ref, dk_ref,
                   pot_dir, pair_name, attempts, multipliers, inst_builder_fn):
    """
    Grid search over scale multipliers for a single optimization phase.

    Each multiplier is tested independently from a clean dk_ref starting point.
    The multiplier producing the lowest energy wins. This is used for Stage 2
    phases where the optimal scale is not known in advance.
    """
    best_e    = e_ref
    best_v    = float("nan")
    best_dk   = None
    best_opt  = None
    best_mult = None

    for idx, mult in enumerate(multipliers):
        load_dk_unquoted(target_obj, dk_ref)

        inst     = inst_builder_fn(mult)
        tag      = f"{base_tag}_m{str(mult).replace('.', '_')}"
        opt_file = f"{pot_dir}opt/{pair_name}_opt_{tag}_{phase_name}.opt"
        dk_out   = f"{pot_dir}dk/{pair_name}_opt_{tag}_{phase_name}.dk"
        log_pref = f"{pot_dir}logs/{pair_name}.opt_{tag}_{phase_name}"

        e, v, ok, msg = try_optimize(
            target_obj, opt_file, dk_out, log_pref, inst, e_ref, require_improvement=True
        )

        attempts.append({
            "phase": phase_name, "tag": tag, "retry": idx,
            "opt_mult": mult, "flat_mult": mult,
            "ok": ok, "energy": e if ok else None, "variance": v if ok else None,
            "dk": dk_out, "opt": opt_file, "msg": msg if not ok else "",
        })

        if ok:
            print(f"  ✅ {phase_name} (mult={mult}x): E={e:.6f}")
            if e < best_e:
                best_e, best_v = float(e), float(v) if is_finite(v) else float("nan")
                best_dk, best_opt, best_mult = dk_out, opt_file, mult
        else:
            print(f"  ⚠️  {phase_name} (mult={mult}x) failed: {msg}")

    ok = best_dk is not None
    if ok:
        print(f"  Best {phase_name}: E={best_e:.6f} (mult={best_mult}x)")
    return ok, best_e, best_v, best_dk, best_opt, best_mult


# ---------------------------------------------------------------------------
# Main optimizer
# ---------------------------------------------------------------------------

def BoundStateOptimizeDynamic(util, pair_name, pot_dir, pot_index, step=0.2, start_scale=0.0):
    print("=" * 72)
    print(f"{pair_name}: Initializing Dynamic Bound State Optimizer")
    print("=" * 72)

    # --- Build initial target and optionally seed from a He-4 alpha core ---
    target0 = build_target(util, pot_dir, pot_index)

    starting_deck_label = "ctrl-default"
    if getattr(util, "ALPHA_DIR", None):
        he4_params_path  = f"{util.NQMCC_DIR}nuclei/params/he4.params"
        seeded_deck_path = f"{pot_dir}dk/{pair_name}_alpha_seed.dk"
        try:
            generate_alpha_seeded_deck(
                target0, pair_name, seeded_deck_path, util.ALPHA_DIR, he4_params_path, copy_esep=True
            )
            target0.DK = deck_t(target0.PARAMS, seeded_deck_path)
            target0.CTRL.INPUT_BRA.DECK_FILE = quote_path(seeded_deck_path)
            starting_deck_label = "alpha-seeded"
            print(f"Alpha seeding successful: {seeded_deck_path}")
        except Exception as ex:
            print(f"Alpha seeding failed ({ex}), continuing with ctrl-default deck")

    # --- Initial evaluation to establish the baseline energy ---
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

    # --- Read optimization parameters ---
    esep_scale_max  = float(util.ESEP_SCALE)
    opt_scale_base  = float(util.OPT_SCALE)
    eta_flat        = float(util.ETA_FLAT)
    three_body_flat = float(util.THREE_BODY_FLAT)
    ss_flat         = float(util.SS_FLAT)

    n_ss = len(get_ss_list(target0.DK))
    if n_ss > 0:
        print(f"Spatial symmetries detected: {n_ss} SS blocks")
    else:
        print("No SS blocks: Stage 2 will run main-only Phase 3")

    # Build the ESEP scale grid
    scales = []
    x = float(start_scale)
    while x <= esep_scale_max + 1e-12:
        scales.append(round(x, 10))
        x += float(step)
    if scales and scales[-1] < esep_scale_max - 1e-12:
        scales.append(esep_scale_max)
    print(f"ESEP scan scales: {scales}")

    attempts = []

    # -----------------------------------------------------------------------
    # STAGE 1: ESEP basin search using Phase 0 (main correlations only)
    # -----------------------------------------------------------------------
    print(f"\n--- STAGE 1: ESEP basin search (Phase 0 only) ---")
    best_stage1 = None

    for scale in scales:
        tag_scale = scale_tag(scale)
        print(f"\nScanning ESEP scale={scale} (tag={tag_scale})")

        target_s1 = build_target(util, pot_dir, pot_index)
        load_dk_unquoted(target_s1, dk_backup)

        def builder_p0(o_m, f_m):
            o_scale = opt_scale_base  * PHASE_SCALES["p0_main"] * o_m
            t_flat  = three_body_flat * PHASE_SCALES["p0_main"] * f_m
            inst    = build_instructions_base(esep_scale_max, o_scale, t_flat, eta_flat, target_s1.DK)
            return set_esep_scale(inst, scale)

        # Allow up to 3 halving retries so an aggressive opt_scale doesn't skip ESEP scale
        ok, e_res, v_res, dk_res, opt_res = run_adaptive_phase(
            target_s1, "p0_main", tag_scale, e0, dk_backup,
            max_retries=3, inst_builder_fn=builder_p0,
            pot_dir=pot_dir, pair_name=pair_name, attempts=attempts,
        )

        if ok and (best_stage1 is None or e_res < best_stage1["energy"]):
            best_stage1 = {
                "scale": scale, "tag": tag_scale,
                "energy": e_res, "variance": v_res,
                "dk": dk_res, "opt": opt_res,
            }

    final_dk  = f"{pot_dir}dk/{pair_name}_optimized.dk"
    final_opt = f"{pot_dir}opt/{pair_name}_opt.opt"

    if best_stage1 is None:
        print("\nStage 1 found no improvement across any ESEP scale. Keeping original deck.")
        shutil.copy2(dk_backup, final_dk)
        best_e, best_v, best_scale = e0, v0, 0.0

    else:
        # -------------------------------------------------------------------
        # STAGE 2: Structural fine-tuning at the winning ESEP scale
        #
        # Each phase runs a grid search over multipliers applied to both
        # opt_scale and flat_scale. The winning multiplier from Phase 1 is
        # folded into effective_ss_flat so Phases 2/3 start from a
        # calibrated baseline rather than the raw util.SS_FLAT.
        # -------------------------------------------------------------------
        best_scale = best_stage1["scale"]
        best_tag   = best_stage1["tag"]
        e_curr     = best_stage1["energy"]
        v_curr     = best_stage1["variance"]
        dk_curr    = best_stage1["dk"]
        opt_curr   = best_stage1["opt"]

        print(f"\n--- STAGE 2: fine-tuning at ESEP scale={best_scale} (E={e_curr:.6f}) ---")

        target_s2 = build_target(util, pot_dir, pot_index)
        load_dk_unquoted(target_s2, dk_curr)

        # effective_ss_flat evolves 
        effective_ss_flat = ss_flat

        GRID = [0.5, 1.0, 1.5, 2.0]

        if n_ss > 0:
            # --- Phase 1: first SS block only ---
            print(f"\nPhase 1: first SS block")
            def builder_p1(mult):
                s_flat = effective_ss_flat * PHASE_SCALES["p1_ss0"] * mult
                return build_instructions_ss([0], flat_value=s_flat)

            ok1, e_curr, v_curr, dk_curr, opt_curr, mult1 = run_grid_phase(
                target_s2, "p1_ss0", best_tag, e_curr, dk_curr,
                pot_dir, pair_name, attempts, GRID, builder_p1,
            )
            if ok1:
                effective_ss_flat *= mult1
                print(f"  Carrying effective_ss_flat={effective_ss_flat:.6f} forward")
                load_dk_unquoted(target_s2, dk_curr)

            # --- Phase 2: remaining SS blocks + main (if >1 SS block) ---
            rest = list(range(1, n_ss))
            if rest:
                print(f"\nPhase 2: remaining SS blocks ({rest}) + main")
                def builder_p2(mult):
                    o_scale = opt_scale_base  * PHASE_SCALES["p2_ss_rest"] * mult
                    t_flat  = three_body_flat * PHASE_SCALES["p2_ss_rest"] * mult
                    s_flat  = effective_ss_flat * PHASE_SCALES["p2_ss_rest"] * mult
                    inst = set_esep_scale(
                        build_instructions_base(esep_scale_max, o_scale, t_flat, eta_flat, target_s2.DK),
                        best_scale,
                    )
                    return inst + build_instructions_ss(rest, flat_value=s_flat)

                ok2, e_curr, v_curr, dk_curr, opt_curr, mult2 = run_grid_phase(
                    target_s2, "p2_ss_rest", best_tag, e_curr, dk_curr,
                    pot_dir, pair_name, attempts, GRID, builder_p2,
                )
                if ok2:
                    effective_ss_flat *= mult2
                    print(f"  Carrying effective_ss_flat={effective_ss_flat:.6f} forward")
                    load_dk_unquoted(target_s2, dk_curr)

        # --- Phase 3: full optimization (main + all SS) ---
        print(f"\nPhase 3: full optimization (main + all {n_ss} SS blocks)")
        def builder_p3(mult):
            o_scale = opt_scale_base  * PHASE_SCALES["p3_full"] * mult
            t_flat  = three_body_flat * PHASE_SCALES["p3_full"] * mult
            s_flat  = effective_ss_flat * PHASE_SCALES["p3_full"] * mult
            inst = set_esep_scale(
                build_instructions_base(esep_scale_max, o_scale, t_flat, eta_flat, target_s2.DK),
                best_scale,
            )
            return inst + build_instructions_ss(list(range(n_ss)), flat_value=s_flat)

        ok3, e_curr, v_curr, dk_curr, opt_curr, _ = run_grid_phase(
            target_s2, "p3_full", best_tag, e_curr, dk_curr,
            pot_dir, pair_name, attempts, GRID, builder_p3,
        )

        shutil.copy2(dk_curr, final_dk)
        if opt_curr and os.path.exists(opt_curr):
            shutil.copy2(opt_curr, final_opt)

        best_e, best_v = e_curr, v_curr

    result = {
        "POTENTIAL_PAIR":    pair_name,
        "STARTING_DECK_MODE": starting_deck_label,
        "ALPHA_CORES_DIR":   getattr(util, "ALPHA_DIR", None),
        "E_INITIAL":         e0,
        "V_INITIAL":         v0,
        "E_OPTIMIZED":       best_e,
        "V_OPTIMIZED":       best_v,
        "DELTA_E":           best_e - e0,
        "BEST_ESEP_SCALE":   best_scale,
        "FINAL_DK":          final_dk,
        "FINAL_OPT":         final_opt,
        "DK_BACKUP":         dk_backup,
        "ATTEMPTS":          attempts,
    }
    with open(f"{pot_dir}{pair_name}_results.json", "w") as f:
        json.dump(result, f, indent=2)

    return result


# ---------------------------------------------------------------------------
# Run loop
# ---------------------------------------------------------------------------

def run(util, step=0.2, start_scale=0.0):
    setup_main_directory(util)

    # Read radius constraint settings from the ctrl file via a temporary target.
    pair_name_temp = create_pair_name(util.TWO_BODY_FILES[0], util.THREE_BODY_FILES[0])
    pot_dir_temp   = setup_potential_pair_directory(util, pair_name_temp)
    target_temp    = build_target(util, pot_dir_temp, 0)

    limit_radii = str(getattr(target_temp.CTRL, "LIMIT_CHARGE_RADII", ".false.")).strip().lower() == ".true."
    if limit_radii:
        radius_settings = {
            "limit_charge_radii":  True,
            "neutron_radius_limit": float(getattr(target_temp.CTRL, "NEUTRON_RADIUS_LIMIT", None)),
            "proton_radius_limit":  float(getattr(target_temp.CTRL, "PROTON_RADIUS_LIMIT", None)),
        }
        print(f"LIMIT_CHARGE_RADII: True")
        print(f"NEUTRON_RADIUS_LIMIT: {radius_settings['neutron_radius_limit']}")
        print(f"PROTON_RADIUS_LIMIT:  {radius_settings['proton_radius_limit']}")
    else:
        radius_settings = {"limit_charge_radii": False}
        print(f"LIMIT_CHARGE_RADII: False")

    metadata = {
        "timestamp":          datetime.now().isoformat(),
        "optimization_mode":  "two_stage_dynamic",
        "optimization_settings": {
            "opt_scale":          util.OPT_SCALE,
            "num_opt_evaluations": util.NUM_OPT_EVALUATIONS,
            "esep_scale":         float(util.ESEP_SCALE),
            "step_size":          step,
            "start_scale":        start_scale,
            "phase_scales":       PHASE_SCALES,
            "eta_flat":           util.ETA_FLAT,
            "three_body_flat":    util.THREE_BODY_FLAT,
            "ss_flat":            util.SS_FLAT,
        },
        "radius_settings": radius_settings,
    }

    # Capture print output to a log file
    opt_output_lines = []

    class _CapturingPrint:
        def __init__(self, real_print):
            self._real = real_print
        def __call__(self, *args, **kwargs):
            opt_output_lines.append(" ".join(str(a) for a in args))
            self._real(*args, **kwargs)

    import builtins
    original_print = builtins.print
    builtins.print = _CapturingPrint(original_print)
    globals()["print"] = builtins.print

    try:
        results = []
        for i in range(len(util.TWO_BODY_FILES)):
            pair_name = create_pair_name(util.TWO_BODY_FILES[i], util.THREE_BODY_FILES[i])
            pot_dir   = setup_potential_pair_directory(util, pair_name)

            print("\n" + "=" * 72)
            print(f"PROCESSING {i+1}/{len(util.TWO_BODY_FILES)}: {pair_name}")
            print("=" * 72)

            try:
                res = BoundStateOptimizeDynamic(
                    util, pair_name, pot_dir, i, step=step, start_scale=start_scale
                )
                results.append(res)
            except Exception as ex:
                print(f"ERROR optimizing {pair_name}: {ex}")
                results.append({"POTENTIAL_PAIR": pair_name, "ERROR": str(ex)})
    finally:
        builtins.print      = original_print
        globals()["print"]  = original_print

    # Write combined results 
    processed_results = [{k: v for k, v in r.items() if k != "ATTEMPTS"} for r in results]
    combined_data = {"run_metadata": metadata, "results": processed_results}

    out = f"{util.WORKING_DIR}combined_results.json"
    with open(out, "w") as f:
        json.dump(combined_data, f, indent=2)
    print("\nSaved:", out)

    out_file = f"{util.WORKING_DIR}optimization_log.out"
    with open(out_file, "w") as f:
        f.write("# Optimization Log (Two-Stage Dynamic)\n")
        f.write(f"# Generated:           {datetime.now().isoformat()}\n")
        f.write(f"# OPT_SCALE:           {util.OPT_SCALE}\n")
        f.write(f"# NUM_OPT_EVALUATIONS: {util.NUM_OPT_EVALUATIONS}\n")
        f.write(f"# ESEP_SCALE:          {util.ESEP_SCALE}\n")
        f.write(f"# STEP_SIZE:           {step}\n")
        f.write(f"# START_SCALE:         {start_scale}\n")
        f.write(f"# PHASE_SCALES:        {PHASE_SCALES}\n")
        f.write(f"# ETA_FLAT:            {util.ETA_FLAT}\n")
        f.write(f"# THREE_BODY_FLAT:     {util.THREE_BODY_FLAT}\n")
        f.write(f"# SS_FLAT:             {util.SS_FLAT}\n")
        f.write(f"# LIMIT_CHARGE_RADII:  {radius_settings['limit_charge_radii']}\n")
        if radius_settings["limit_charge_radii"]:
            f.write(f"# NEUTRON_RADIUS_LIMIT: {radius_settings['neutron_radius_limit']}\n")
            f.write(f"# PROTON_RADIUS_LIMIT:  {radius_settings['proton_radius_limit']}\n")
        f.write("#" + "=" * 70 + "\n\n")
        for line in opt_output_lines:
            f.write(line + "\n")
    print("Saved optimization log:", out_file)

    return results


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--utility", required=True)
    p.add_argument("--step",    type=float, default=0.2)
    p.add_argument("--start",   type=float, default=0.0)
    args = p.parse_args()

    util = utility_t(args.utility)
    run(util, step=args.step, start_scale=args.start)
