"""
bopt.py
Bound-state optimizer with ESEP scan and alpha seeding
"""
#-----------------------------------------------------------------------
import os
import sys
import json
import shutil
from datetime import datetime
from pathlib import Path
import re
import numpy as np
from contextlib import redirect_stdout
from io import StringIO
#-----------------------------------------------------------------------
sys.path.append(os.path.join(os.path.dirname(__file__), "src"))
#-----------------------------------------------------------------------
from utility import utility_t
from deck import deck_t, GenerateOptFile
from wavefunction import wavefunction_t, InitNShellBoundWF
from parameters import parameters_t
#-----------------------------------------------------------------------
"""
---
Bound State Optimization Logic
---
1. Parse utility file for ESEP_SCALE, OPT_SCALE, optional ALPHA_DIR
2. Setup working directories for each potential pair
3. Build target wavefunction from control file
4. Optional: Apply alpha-core seeding if ALPHA_DIR provided
5. Perform initial evaluation to get baseline energy
6. ESEP scan from 0 to ESEP_SCALE with configurable step size (default 0.2)
7. For each ESEP scale: optimize all correlations simultaneously
8. Track best optimization result across all scales
9. Save optimized deck and results in JSON format
---
Alpha Seeding Logic (Optional)
---
1. Check if ALPHA_DIR exists and is accessible
2. Find pre-optimized He-4 alpha-core deck for current potential pair
3. Copy optimized alpha parameters to target deck
    - If seeding fails: continue with standard unseeded optimization
    - If successful: use seeded deck as starting point for better convergence
---

"""

#-----------------------------------------------------------------------

def is_finite(x):
    try:
        return x is not None and np.isfinite(float(x))
    except Exception:
        return False
#-----------------------------------------------------------------------
def strip_quotes(s):
    s = str(s).strip()
    if len(s) >= 2 and s[0] == s[-1] and s[0] in ("'", '"'):
        return s[1:-1]
    return s
#-----------------------------------------------------------------------
def quote_path(path):
    return f"'{path}'"
#-----------------------------------------------------------------------
def ensure_dir(p):
    if p:
        Path(p).mkdir(parents=True, exist_ok=True)
#-----------------------------------------------------------------------
def create_pair_name(pot2b, pot3b):
    return f"{pot2b.split('.')[0]}.{pot3b.split('.')[0]}"
#-----------------------------------------------------------------------
def set_esep_scale(instructions, esep_scale):
    out = []
    for inst in instructions:
        d = dict(inst)
        if d.get("key") == "ESEP":
            d["scale"] = float(esep_scale)
        out.append(d)
    return out
#-----------------------------------------------------------------------
def scale_tag(scale):
    return f"s{int(round(scale * 100)):03d}"
#-----------------------------------------------------------------------
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
            raw = tok[1:-1]
        else:
            raw = tok

        norm = os.path.normpath(raw)
        if norm.startswith(nqmcc_dir):
            rel = os.path.relpath(norm, nqmcc_dir)
            if quote:
                return f"{quote}{rel}{quote}"
            return rel
        return token

    text = re.sub(r"'[^']*'", lambda m: rewrite_token(m.group(0)), text)
    text = re.sub(r'"[^"]*"', lambda m: rewrite_token(m.group(0)), text)

    ensure_dir(os.path.dirname(ctrl_out))
    with open(ctrl_out, "w") as f:
        f.write(text)
    return ctrl_out
#-----------------------------------------------------------------------
def setup_main_directory(util):
    try:
        os.mkdir(util.WORKING_DIR)
    except FileExistsError:
        util.WORKING_DIR = f"{util.WORKING_DIR}{util.NAME}-{datetime.now().strftime('%Y-%m-%d_%H-%M')}/"
        os.mkdir(util.WORKING_DIR)
#-----------------------------------------------------------------------
def setup_potential_pair_directory(util, pair_name):
    pot_dir = f"{util.WORKING_DIR}{pair_name}/"
    ensure_dir(pot_dir)
    ensure_dir(f"{pot_dir}ctrl")
    ensure_dir(f"{pot_dir}logs")
    ensure_dir(f"{pot_dir}dk")
    ensure_dir(f"{pot_dir}opt")
    ensure_dir(f"{util.NQMCC_DIR}walks")
    return pot_dir
#-----------------------------------------------------------------------
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
#-----------------------------------------------------------------------
def find_alpha_core(alpha_cores_dir, pot_pair_name):
    deck_name = f"{pot_pair_name}_optimized.dk"
    path = os.path.join(alpha_cores_dir, pot_pair_name, "dk", deck_name)
    if not os.path.exists(path):
        raise FileNotFoundError(f"Could not find optimized alpha deck: {path}")
    return path
#-----------------------------------------------------------------------
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
#-----------------------------------------------------------------------
def try_optimize(target, opt_file, dk_out, log_prefix, instructions, e_ref, require_improvement=True):
    opt_file = strip_quotes(opt_file)
    dk_out = strip_quotes(dk_out)
    log_prefix = strip_quotes(log_prefix)

    ensure_dir(os.path.dirname(opt_file))
    ensure_dir(os.path.dirname(dk_out))
    ensure_dir(os.path.dirname(log_prefix))

    try:
        opt = GenerateOptFile(target.PARAMS, target.DK, quote_path(opt_file), instructions)
        opt.UpdateFloats(target.PARAMS, 6)
        e, v = target.Optimize(opt, quote_path(dk_out), True, log_prefix)
    except FileNotFoundError as ex:
        if dk_out in str(ex):
            return None, None, False, f"No DK output: {dk_out}"
        return None, None, False, f"FileNotFoundError: {ex}"
    except Exception as ex:
        return None, None, False, f"{ex}"

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

#-----------------------------------------------------------------------
def BoundStateOptimize(util, pair_name, pot_dir, pot_index, step=0.1):
    print("=" * 72)
    print(f"{pair_name}: building target + optional alpha-seeding + initial Evaluate")
    print("=" * 72)

    target0 = build_target(util, pot_dir, pot_index)

    starting_deck_label = "ctrl-default"
    if util.ALPHA_DIR:
        he4_params_path = f"{util.NQMCC_DIR}nuclei/params/he4.params"
        seeded_deck_path = f"{pot_dir}dk/{pair_name}_alpha_seed.dk"
        try:
            generate_alpha_seeded_deck(target0, pair_name, seeded_deck_path, util.ALPHA_DIR, he4_params_path)
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

    esep_scale = float(util.ESEP_SCALE)
    opt_scale = float(util.OPT_SCALE)
    opt_3b = opt_scale + 0.2
    
    instructions_base = [
        {"ss": False, "key": "ESEP", "idx": 0, "scale": esep_scale, "flat": 0.0},
        {"ss": False, "key": "ESEP", "idx": 1, "scale": esep_scale, "flat": 0.0},
        {"ss": False, "key": "ESEP", "idx": 2, "scale": esep_scale, "flat": 0.0},
        {"ss": False, "key": "ESEP", "idx": 3, "scale": esep_scale, "flat": 0.0},
        {"ss": False, "key": "FSCAL", "all": True, "scale": opt_scale, "flat": 0.0},
        {"ss": False, "key": "ALPHA", "all": True, "scale": opt_scale, "flat": 0.0},
        {"ss": False, "key": "BETA",  "all": True, "scale": opt_scale, "flat": 0.0},
        {"ss": False, "key": "DELTA",   "scale": 0, "flat": opt_3b},
        {"ss": False, "key": "EPSILON", "scale": 0, "flat": opt_3b},
        {"ss": False, "key": "THETA",   "scale": 0, "flat": opt_3b},
        {"ss": False, "key": "UPSILON", "scale": 0, "flat": opt_3b},
        {"ss": False, "key": "RSCAL",   "scale": 0, "flat": opt_3b},
        {"ss": False, "key": "USCAL",   "scale": 0, "flat": opt_3b},
        {"ss": False, "key": "QPS1",  "scale": opt_scale, "flat": 0.0},
        {"ss": False, "key": "QPS2",  "scale": opt_scale, "flat": 0.0},
        {"ss": False, "key": "QSSS1", "scale": opt_scale, "flat": 0.0},
    ]

    # Check if spatial symmetries exist and optimize all of them
    instructions = instructions_base
    if hasattr(target0.DK, 'SS') and len(target0.DK.SS) > 0:
        instructions_ss = []
        for ssi in range(len(target0.DK.SS)):
            instructions_ss += [
                {"ss":True,"ss_idx":ssi,"key":"SPU"  ,"scale":0.0, "flat":opt_scale},
                {"ss":True,"ss_idx":ssi,"key":"SPV"  ,"scale":0.0, "flat":opt_scale},
                {"ss":True,"ss_idx":ssi,"key":"SPA"  ,"scale":0.0, "flat":opt_scale},
                {"ss":True,"ss_idx":ssi,"key":"SPB"  ,"scale":0.0, "flat":opt_scale},
                {"ss":True,"ss_idx":ssi,"key":"SPC"  ,"scale":0.0, "flat":opt_scale},
                {"ss":True,"ss_idx":ssi,"key":"SPK"  ,"scale":0.0, "flat":opt_scale},
                {"ss":True,"ss_idx":ssi,"key":"SPL"  ,"scale":0.0, "flat":opt_scale},
                {"ss":True,"ss_idx":ssi,"key":"WSR"  ,"scale":0.0, "flat":opt_scale},
                {"ss":True,"ss_idx":ssi,"key":"WSA"  ,"scale":0.0, "flat":opt_scale},
                {"ss": True, "ss_idx": ssi, "key": "BETALSN", "scale": 0.0, "flat": opt_scale}
            ]
        instructions = instructions_base + instructions_ss
        print(f"Spatial symmetries detected: optimizing {len(target0.DK.SS)} SS parameters")


    max_scale = esep_scale
    scales = []
    x = 0.2
    while x <= max_scale + 1e-12:
        scales.append(round(x, 10))
        x += float(step)
    if scales[-1] < max_scale - 1e-12:
        scales.append(max_scale)

    print(f"ESEP scan scales: {scales}")

    attempts = []
    best = None

    for scale in scales:
        tag = scale_tag(scale)
        print(f"\n{pair_name}: TRY scale={scale} (tag={tag})")

        target = build_target(util, pot_dir, pot_index)
        target.DK.FILE_NAME = quote_path(dk_backup)
        target.DK.Read(target.PARAMS)

        opt_file = f"{pot_dir}opt/{pair_name}_opt_{tag}.opt"
        dk_out   = f"{pot_dir}dk/{pair_name}_opt_{tag}.dk"
        log_pref = f"{pot_dir}logs/{pair_name}.opt_{tag}"

        e, v, ok, msg = try_optimize(target, opt_file, dk_out, log_pref, set_esep_scale(instructions, scale), e0, True)

        attempt = {
            "scale": scale, "tag": tag, "ok": ok,
            "energy": e if ok else None, "variance": v if ok else None,
            "dk": dk_out, "opt": opt_file, "msg": msg if not ok else "",
        }
        attempts.append(attempt)

        if ok:
            dE = e - e0
            print(f"  ✅ success: E={e:.6f}  ΔE={dE:.6f}")
            if best is None or e < best["energy"]:
                best = {"scale": scale, "tag": tag, "energy": e, "variance": v, "dk": dk_out, "opt": opt_file}
        else:
            print(f"  ⚠️ failed: {msg}")

    if best is None:
        print(f"\n{pair_name}: No successful run. Using ORIGINAL as optimized.")
        final_dk = f"{pot_dir}dk/{pair_name}_optimized.dk"
        shutil.copy2(dk_backup, final_dk)
        best_e, best_v, best_scale = e0, v0, None
        final_opt = None
    else:
        print(f"\n{pair_name}: BEST scale={best['scale']}  E={best['energy']:.6f}  ΔE={best['energy']-e0:.6f}")
        final_dk = f"{pot_dir}dk/{pair_name}_optimized.dk"
        final_opt = f"{pot_dir}opt/{pair_name}_opt.opt"
        shutil.copy2(best["dk"], final_dk)
        if os.path.exists(best["opt"]):
            shutil.copy2(best["opt"], final_opt)
        best_e, best_v, best_scale = best["energy"], best["variance"], best["scale"]

    result = {
        "POTENTIAL_PAIR": pair_name, "STARTING_DECK_MODE": starting_deck_label,
        "ALPHA_CORES_DIR": util.ALPHA_DIR,
        "E_INITIAL": e0, "V_INITIAL": v0, "E_OPTIMIZED": best_e, "V_OPTIMIZED": best_v,
        "DELTA_E": best_e - e0, "BEST_ESEP_SCALE": best_scale,
        "FINAL_DK": final_dk, "FINAL_OPT": final_opt, "DK_BACKUP": dk_backup, "ATTEMPTS": attempts,
    }
    with open(f"{pot_dir}{pair_name}_results.json", "w") as f:
        json.dump(result, f, indent=2)
    return result
#-----------------------------------------------------------------------
def run(util, step=0.2):
    setup_main_directory(util)

    # Create metadata dictionary
    run_metadata = {
        "timestamp": datetime.now().isoformat(),
        "optimization_settings": {
            "opt_scale": util.OPT_SCALE,
            "num_opt_evaluations": util.NUM_OPT_EVALUATIONS,
            "esep_scale": float(util.ESEP_SCALE),
            "step_size": step
        }
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
        for i in range(util.NUM_POTS):
            pair_name = create_pair_name(util.TWO_BODY_FILES[i], util.THREE_BODY_FILES[i])
            pot_dir = setup_potential_pair_directory(util, pair_name)

            print("\n" + "=" * 72)
            print(f"PROCESSING {i+1}/{util.NUM_POTS}: {pair_name}")
            print("=" * 72)

            cwd = os.getcwd()
            os.chdir(pot_dir)
            try:
                res = BoundStateOptimize(util, pair_name, pot_dir, i, step)
                results.append(res)
            finally:
                os.chdir(cwd)

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
        "run_metadata": run_metadata,
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
        f.write("#" + "="*70 + "\n\n")
        for line in opt_output_lines:
            f.write(line + "\n")
    print("Saved optimization log:", out_file)
#-----------------------------------------------------------------------
if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser()
    p.add_argument("--utility", required=True)
    p.add_argument("--step", type=float, default=0.2)
    args = p.parse_args()

    util = utility_t(args.utility)
    run(util, args.step)
