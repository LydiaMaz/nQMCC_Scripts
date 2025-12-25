# Multi-potential bound-state optimizer (ESEP scan + cleanup)
#
# Per potential pair:
#   - Evaluate initial
#   - Save ORIGINAL DK to disk
#   - Scan ESEP_SCALE from 0 to util.ESEP_SCALE in steps of 0.2
#       * For EVERY scale: build new wavefunction_t, load ORIGINAL DK, Optimize
#       * Keep the run with the MIN energy
#   - Promote best to dk/{pair}_optimized.dk and opt/{pair}_opt.opt
#   - Cleanup all non-best attempt files (dk/opt/logs)
#   - If all runs fail or never improve, keep ORIGINAL as optimized
# -----------------------------------------------------------------------

import os
import sys
import json
import shutil
from datetime import datetime
from pathlib import Path
import re

import numpy as np

# -----------------------------------------------------------------------
sys.path.append(os.path.join(os.path.dirname(__file__), "src"))
# -----------------------------------------------------------------------
from utility import utility_t
from deck import GenerateOptFile
from wavefunction import wavefunction_t


# =======================================================================
# Small utilities
# =======================================================================

def _is_finite(x) -> bool:
    try:
        return x is not None and np.isfinite(float(x))
    except Exception:
        return False

def _strip_quotes(s: str) -> str:
    s = str(s).strip()
    if len(s) >= 2 and s[0] == s[-1] and s[0] in ("'", '"'):
        return s[1:-1]
    return s

def _q(path: str) -> str:
    return f"'{path}'"

def _ensure_dir(p: str):
    Path(p).mkdir(parents=True, exist_ok=True)

def create_potential_pair_name(pot2b: str, pot3b: str) -> str:
    return f"{pot2b.split('.')[0]}.{pot3b.split('.')[0]}"

def _with_esep_scale(instructions: list, esep_scale: float) -> list:
    out = []
    for inst in instructions:
        d = dict(inst)
        if d.get("key") == "ESEP":
            d["scale"] = float(esep_scale)
        out.append(d)
    return out

def _fmt_scale_tag(scale: float) -> str:
    # 0.2 -> "s020", 1.0 -> "s100"
    return f"s{int(round(scale * 100)):03d}"


# =======================================================================
# CTRL patching (fixes doubled-path crash)
# =======================================================================

def _make_ctrl_paths_relative(ctrl_in: str, nqmcc_dir: str, ctrl_out: str) -> str:
    """
    Rewrite any quoted path that starts with NQMCC_DIR to be relative to NQMCC_DIR.
    Prevents internal prefixing from duplicating absolute paths.
    """
    nqmcc_dir = os.path.normpath(_strip_quotes(nqmcc_dir))
    ctrl_in = _strip_quotes(ctrl_in)
    ctrl_out = _strip_quotes(ctrl_out)

    with open(ctrl_in, "r") as f:
        text = f.read()

    def rewrite_token(token: str) -> str:
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

    def sub_single(m):
        return rewrite_token(m.group(0))
    text = re.sub(r"'[^']*'", sub_single, text)

    def sub_double(m):
        return rewrite_token(m.group(0))
    text = re.sub(r'"[^"]*"', sub_double, text)

    _ensure_dir(os.path.dirname(ctrl_out))
    with open(ctrl_out, "w") as f:
        f.write(text)

    return ctrl_out


# =======================================================================
# Setup helpers
# =======================================================================

def setup_main_directory(util: utility_t):
    try:
        os.mkdir(util.WORKING_DIR)
    except FileExistsError:
        util.WORKING_DIR = f"{util.WORKING_DIR}{util.NAME}-{datetime.now().strftime('%Y-%m-%d_%H-%M')}/"
        os.mkdir(util.WORKING_DIR)

def setup_potential_pair_directory(util: utility_t, pair_name: str) -> str:
    pot_dir = f"{util.WORKING_DIR}{pair_name}/"
    _ensure_dir(pot_dir)
    _ensure_dir(f"{pot_dir}ctrl")
    _ensure_dir(f"{pot_dir}logs")
    _ensure_dir(f"{pot_dir}dk")
    _ensure_dir(f"{pot_dir}opt")
    _ensure_dir(f"{util.NQMCC_DIR}walks")
    return pot_dir

def build_target(util: utility_t, pot_dir: str, pot_index: int) -> wavefunction_t:
    """
    Build a wavefunction_t from a patched ctrl and set per-pair potentials.
    """
    base_ctrl_in = _strip_quotes(util.CTRL_FILE)
    patched_ctrl = f"{pot_dir}ctrl/base.ctrl"
    _make_ctrl_paths_relative(base_ctrl_in, util.NQMCC_DIR, patched_ctrl)

    target = wavefunction_t(patched_ctrl, util.NQMCC_DIR, util.BIN_DIR, util.RUN_CMD)

    target.CTRL.FILE_NAME = f"{pot_dir}ctrl/target.ctrl"
    target.CTRL.NUM_BLOCKS = util.NUM_BLOCKS
    target.CTRL.BLOCK_SIZE = util.BLOCK_SIZE
    target.CTRL.WALKERS_PER_NODE = util.WALKERS_PER_NODE
    target.CTRL.NUM_OPT_EVALUATIONS = util.NUM_OPT_EVALUATIONS

    target.CTRL.CONST_FILE = _q(f"{util.NQMCC_DIR}constants/{util.CONSTANTS_FILES[pot_index]}")
    target.CTRL.L2BP_FILE  = _q(f"{util.NQMCC_DIR}pots/{util.TWO_BODY_FILES[pot_index]}")
    target.CTRL.L3BP_FILE  = _q(f"{util.NQMCC_DIR}pots/{util.THREE_BODY_FILES[pot_index]}")

    return target


# =======================================================================
# Optimize logic
# =======================================================================

def _try_optimize(target: wavefunction_t,
                  opt_file: str,
                  dk_out: str,
                  log_prefix: str,
                  instructions: list,
                  e_ref: float,
                  require_improvement: bool = True):
    """
    Returns (e, v, ok, msg).
    If require_improvement=True, only accepts runs with e < e_ref.
    """
    opt_file = _strip_quotes(opt_file)
    dk_out = _strip_quotes(dk_out)
    log_prefix = _strip_quotes(log_prefix)

    _ensure_dir(os.path.dirname(opt_file))
    _ensure_dir(os.path.dirname(dk_out))
    _ensure_dir(os.path.dirname(log_prefix))

    try:
        opt = GenerateOptFile(target.PARAMS, target.DK, _q(opt_file), instructions)
        opt.UpdateFloats(target.PARAMS, 6)
        e, v = target.Optimize(opt, _q(dk_out), True, log_prefix)
    except FileNotFoundError as ex:
        if dk_out in str(ex):
            return None, None, False, f"Optimize produced no DK output: {dk_out}"
        return None, None, False, f"FileNotFoundError: {ex}"
    except Exception as ex:
        return None, None, False, f"{ex}"

    if not _is_finite(e):
        return None, None, False, "Optimize returned non-finite energy"
    if not os.path.exists(dk_out):
        return None, None, False, f"Optimize produced no DK output: {dk_out}"

    e = float(e)
    v = float(v) if _is_finite(v) else float("nan")

    if require_improvement:
        dE = e - float(e_ref)
        if dE >= 0.0:
            return None, None, False, f"Rejected (no improvement): ΔE={dE:.6f}"

    return e, v, True, ""


# =======================================================================
# Cleanup helpers
# =======================================================================

def _cleanup_attempt_artifacts(pot_dir: str, pair_name: str, keep_tags: set):
    """
    Delete attempt files for tags not in keep_tags.
    We consider:
      - dk/{pair}_opt_{tag}.dk
      - opt/{pair}_opt_{tag}.opt
      - logs/{pair}.opt_{tag}*  (prefix cleanup)
    We also delete the corresponding dk/opt attempt names used below.
    """
    dk_dir = f"{pot_dir}dk"
    opt_dir = f"{pot_dir}opt"
    log_dir = f"{pot_dir}logs"

    # DK cleanup
    if os.path.isdir(dk_dir):
        for fn in os.listdir(dk_dir):
            path = os.path.join(dk_dir, fn)
            if not os.path.isfile(path):
                continue
            # match attempt dk
            m = re.match(rf"^{re.escape(pair_name)}_opt_(s\d{{3}})\.dk$", fn)
            if m and m.group(1) not in keep_tags:
                try:
                    os.remove(path)
                except Exception:
                    pass

    # OPT cleanup
    if os.path.isdir(opt_dir):
        for fn in os.listdir(opt_dir):
            path = os.path.join(opt_dir, fn)
            if not os.path.isfile(path):
                continue
            m = re.match(rf"^{re.escape(pair_name)}_opt_(s\d{{3}})\.opt$", fn)
            if m and m.group(1) not in keep_tags:
                try:
                    os.remove(path)
                except Exception:
                    pass

    # LOG cleanup (prefix-based)
    if os.path.isdir(log_dir):
        for fn in os.listdir(log_dir):
            path = os.path.join(log_dir, fn)
            if not os.path.isfile(path):
                continue
            # anything starting with "{pair}.opt_sXYZ"
            m = re.match(rf"^{re.escape(pair_name)}\.opt_(s\d{{3}})", fn)
            if m and m.group(1) not in keep_tags:
                try:
                    os.remove(path)
                except Exception:
                    pass


# =======================================================================
# Main optimization per pair (ESEP scan)
# =======================================================================

def BoundStateOptimize(util: utility_t, pair_name: str, pot_dir: str, pot_index: int,
                      step: float = 0.2):
    print("=" * 72)
    print(f"{pair_name}: building target + initial Evaluate")
    print("=" * 72)

    # Build once for initial evaluation
    target0 = build_target(util, pot_dir, pot_index)

    log_init = f"{pot_dir}logs/{pair_name}.initial"
    e0, v0 = target0.Evaluate(True, log_init)
    if not _is_finite(e0):
        raise RuntimeError("Initial evaluation failed (Evaluate returned non-finite energy)")
    e0 = float(e0)
    v0 = float(v0) if _is_finite(v0) else float("nan")
    print(f"INITIAL: E={e0:.6f}  V={v0:.6f}")

    # Save ORIGINAL DK to disk (all attempts load this)
    dk_backup = f"{pot_dir}dk/{pair_name}_original.dk"
    target0.DK.Write(target0.PARAMS, _q(dk_backup))
    print(f"Saved ORIGINAL DK: {dk_backup}")

    # Instruction template (always keep ESEP entries; only scale changes)
    instructions_base = [
        {"ss": False, "key": "ESEP", "idx": 0, "scale": float(util.ESEP_SCALE), "flat": 0.0},
        {"ss": False, "key": "ESEP", "idx": 1, "scale": float(util.ESEP_SCALE), "flat": 0.0},
        {"ss": False, "key": "ESEP", "idx": 2, "scale": float(util.ESEP_SCALE), "flat": 0.0},
        {"ss": False, "key": "ESEP", "idx": 3, "scale": float(util.ESEP_SCALE), "flat": 0.0},

        {"ss": False, "key": "FSCAL", "all": True, "scale": float(util.OPT_SCALE), "flat": 0.0},
        {"ss": False, "key": "ALPHA", "all": True, "scale": float(util.OPT_SCALE), "flat": 0.0},
        {"ss": False, "key": "BETA",  "all": True, "scale": float(util.OPT_SCALE), "flat": 0.0},

        {"ss": False, "key": "DELTA",   "scale": 0, "flat": float(util.OPT_SCALE)},
        {"ss": False, "key": "EPSILON", "scale": 0, "flat": float(util.OPT_SCALE)},
        {"ss": False, "key": "THETA",   "scale": 0, "flat": float(util.OPT_SCALE)},
        {"ss": False, "key": "UPSILON", "scale": 0, "flat": float(util.OPT_SCALE)},
        {"ss": False, "key": "RSCAL",   "scale": 0, "flat": float(util.OPT_SCALE)},
        {"ss": False, "key": "USCAL",   "scale": 0, "flat": float(util.OPT_SCALE)},

        {"ss": False, "key": "QPS1",  "scale": float(util.OPT_SCALE), "flat": 0.0},
        {"ss": False, "key": "QPS2",  "scale": float(util.OPT_SCALE), "flat": 0.0},
        {"ss": False, "key": "QSSS1", "scale": float(util.OPT_SCALE), "flat": 0.0},
    ]

    # Build scale list: start at 0, increase by step, include max
    max_scale = float(util.ESEP_SCALE)
    scales = []
    x = 0.0
    while x <= max_scale + 1e-12:
        scales.append(round(x, 10))
        x += float(step)
    if scales[-1] < max_scale - 1e-12:
        scales.append(max_scale)

    print(f"ESEP scan scales: {scales}")

    attempts = []
    best = None  # dict with keys: scale, tag, e, v, dk, opt

    # Scan
    for scale in scales:
        tag = _fmt_scale_tag(scale)
        print(f"\n{pair_name}: TRY scale={scale} (tag={tag})")

        # Fresh target for every attempt
        target = build_target(util, pot_dir, pot_index)
        target.DK.FILE_NAME = _q(dk_backup)
        target.DK.Read(target.PARAMS)

        opt_file = f"{pot_dir}opt/{pair_name}_opt_{tag}.opt"
        dk_out   = f"{pot_dir}dk/{pair_name}_opt_{tag}.dk"
        log_pref = f"{pot_dir}logs/{pair_name}.opt_{tag}"

        e, v, ok, msg = _try_optimize(
            target=target,
            opt_file=opt_file,
            dk_out=dk_out,
            log_prefix=log_pref,
            instructions=_with_esep_scale(instructions_base, scale),
            e_ref=e0,
            require_improvement=True,   # keep your “must improve” safeguard
        )

        attempt = {
            "scale": scale,
            "tag": tag,
            "ok": ok,
            "energy": e if ok else None,
            "variance": v if ok else None,
            "dk": dk_out,
            "opt": opt_file,
            "msg": msg if not ok else "",
        }
        attempts.append(attempt)

        if ok:
            dE = e - e0
            print(f"  ✅ success: E={e:.6f}  ΔE={dE:.6f}")
            if best is None or e < best["energy"]:
                best = {
                    "scale": scale,
                    "tag": tag,
                    "energy": e,
                    "variance": v,
                    "dk": dk_out,
                    "opt": opt_file,
                }
        else:
            print(f"  ⚠️ failed: {msg}")

    # Decide final outcome
    if best is None:
        print(f"\n{pair_name}: No successful improving run. Using ORIGINAL (unoptimized).")
        final_dk = f"{pot_dir}dk/{pair_name}_optimized.dk"
        final_opt = None

        # promote original to optimized
        shutil.copy2(dk_backup, final_dk)

        keep_tags = set()  # no successful tags to keep
        _cleanup_attempt_artifacts(pot_dir, pair_name, keep_tags)

        best_e, best_v, best_scale = e0, v0, None

    else:
        print(f"\n{pair_name}: BEST scale={best['scale']}  E={best['energy']:.6f}  ΔE={best['energy']-e0:.6f}")

        # Promote best to canonical names
        final_dk = f"{pot_dir}dk/{pair_name}_optimized.dk"
        final_opt = f"{pot_dir}opt/{pair_name}_opt.opt"
        shutil.copy2(best["dk"], final_dk)
        if os.path.exists(best["opt"]):
            shutil.copy2(best["opt"], final_opt)

        # Cleanup all non-best attempt artifacts
        keep_tags = {best["tag"]}
        _cleanup_attempt_artifacts(pot_dir, pair_name, keep_tags)

        best_e, best_v, best_scale = best["energy"], best["variance"], best["scale"]

    # Save per-pair results
    result = {
        "POTENTIAL_PAIR": pair_name,
        "E_INITIAL": e0,
        "V_INITIAL": v0,
        "E_OPTIMIZED": best_e,
        "V_OPTIMIZED": best_v,
        "DELTA_E": best_e - e0,
        "BEST_ESEP_SCALE": best_scale,
        "FINAL_DK": final_dk,
        "FINAL_OPT": (f"{pot_dir}opt/{pair_name}_opt.opt" if best is not None else None),
        "DK_BACKUP": dk_backup,
        "ATTEMPTS": attempts,
    }
    with open(f"{pot_dir}{pair_name}_results.json", "w") as f:
        json.dump(result, f, indent=2)

    return result


# =======================================================================
# Driver
# =======================================================================

def run(util: utility_t):
    setup_main_directory(util)

    results = []
    for i in range(util.NUM_POTS):
        pair_name = create_potential_pair_name(util.TWO_BODY_FILES[i], util.THREE_BODY_FILES[i])
        pot_dir = setup_potential_pair_directory(util, pair_name)

        print("\n" + "=" * 72)
        print(f"PROCESSING {i+1}/{util.NUM_POTS}: {pair_name}")
        print("=" * 72)

        cwd = os.getcwd()
        os.chdir(pot_dir)
        try:
            res = BoundStateOptimize(util, pair_name, pot_dir, pot_index=i, step=0.2)
            results.append(res)
        finally:
            os.chdir(cwd)

    out = f"{util.WORKING_DIR}combined_results.json"
    with open(out, "w") as f:
        json.dump(results, f, indent=2)
    print("\nSaved:", out)


# =======================================================================
# CLI
# =======================================================================

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--utility", required=True)
    args = p.parse_args()

    util = utility_t(args.utility)
    run(util)