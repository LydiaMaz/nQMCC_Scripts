#-----------------------------------------------------------------------
import numpy as np
import json
import sys
import os
from datetime import datetime
#-----------------------------------------------------------------------
sys.path.append(os.path.join(os.path.dirname(__file__), "src"))
#-----------------------------------------------------------------------
from utility import utility_t
from deck import deck_t,GenerateOptFile
from wavefunction import wavefunction_t
#-----------------------------------------------------------------------
def setup(util: utility_t):
    BREAK="="*72
#-----------------------------------------------------------------------
    print("SETTING UP WORKING ENVIORMENT")
    try:
        os.mkdir(util.WORKING_DIR)
        os.chdir(util.WORKING_DIR)
    except FileExistsError:
        print("***WORKING DIRECTORY EXISTS***")
        util.WORKING_DIR=f"{util.WORKING_DIR}{util.NAME}-{datetime.now().strftime('%Y-%m-%d_%H-%M')}/"
        print(f"CURRENT WORKING DIRECTORY: {util.WORKING_DIR}")
        os.mkdir(util.WORKING_DIR)
        os.chdir(util.WORKING_DIR)
#-----------------------------------------------------------------------
    os.mkdir("ctrl")
    os.mkdir("logs")
    os.mkdir("dk")
    os.mkdir("opt")
#-----------------------------------------------------------------------
    print("... DONE")
    print(BREAK)
#-----------------------------------------------------------------------
    print("SETTING UP TARGET WAVEFUNCTION")
    # Ensure required nQMCC directories exist
    os.makedirs(f"{util.NQMCC_DIR}walks", exist_ok=True)
    target = wavefunction_t(util.CTRL_FILE,util.NQMCC_DIR,util.BIN_DIR,util.RUN_CMD)
    target_label=target.DK.NAME.strip("\'")
    target.CTRL.FILE_NAME=f"{util.WORKING_DIR}target.ctrl"
    target.CTRL.NUM_BLOCKS=util.NUM_BLOCKS
    target.CTRL.BLOCK_SIZE=util.BLOCK_SIZE
    target.CTRL.WALKERS_PER_NODE=util.WALKERS_PER_NODE
    target.CTRL.NUM_OPT_EVALUATIONS=util.NUM_OPT_EVALUATIONS
    print("... DONE")
#-----------------------------------------------------------------------
    for const,pot2b,pot3b in zip(util.CONSTANTS_FILES,util.TWO_BODY_FILES,util.THREE_BODY_FILES):
        pot2b_label=pot2b.split(".")[0]
        pot3b_label=pot3b.split(".")[0]
        tname=f"{target_label}.{pot2b_label}.{pot3b_label}"
        target.CTRL.CONST_FILE=f"'{util.NQMCC_DIR}constants/{const}'"
        target.CTRL.L2BP_FILE=f"'{util.NQMCC_DIR}pots/{pot2b}'"
        target.CTRL.L3BP_FILE=f"'{util.NQMCC_DIR}pots/{pot3b}'"
        print(f"CONSTANTS: {const}")
        print(f"V2B: {pot2b}")
        print(f"V3B: {pot3b}")
#-----------------------------------------------------------------------
    return target

def BoundStateOptimize(target: wavefunction_t, util: utility_t, label: str):
    """
    Optimize bound state wavefunction by varying ESEP separation energies
    
    Steps:
    1. Initial evaluation  
    2. Set up optimization file for ESEP values
    3. Optimize ESEP parameters
    4. Final evaluation and logging
    5. Return results
    """
    BREAK="="*72
    work_dir = util.WORKING_DIR
    opt_scale = util.OPT_SCALE
    esep_scale = util.ESEP_SCALE
    
    print(f" 🔥 BOUND STATE OPTIMIZATION: {label} 🔥")
    print(BREAK)
    
    # Step 1: Initial evaluation
    log_initial = f"{work_dir}logs/{label}.esep.initial"
    print(f"Step 1: INITIAL EVALUATION: {log_initial}.energy")
    e_initial, v_initial = target.Evaluate(True, log_initial)
    
    if e_initial is None:
        raise RuntimeError("Initial evaluation failed - could not compute energy")
    
    print(f" ⚛ INITIAL E = {e_initial:.4f} +- {v_initial:.4f}")
    esep_initial = target.DK.ESEP
    print(f"INITIAL ESEP VALUES: {esep_initial}")
    print(BREAK)
    
    # Step 2: Set up optimization file for ESEP values  
    print("Step 2: SETTING UP ESEP OPTIMIZATION FILE")
    opt_esep_file_name = f"'{work_dir}opt/{label}_esep.opt'"
    
    # Instructions for optimizing ESEP only
    instructions_esep = [
        {"ss": False, "key": "ESEP", "idx": 0, "scale": esep_scale, "flat": 0.0},  # ESEP[0]
        {"ss": False, "key": "ESEP", "idx": 1, "scale": esep_scale, "flat": 0.0},  # ESEP[1]
        {"ss": False, "key": "ESEP", "idx": 2, "scale": esep_scale, "flat": 0.0},  # ESEP[2]
        {"ss": False, "key": "ESEP", "idx": 3, "scale": esep_scale, "flat": 0.0},  # ESEP[3]
    ]

    # Instructions for correlations
    instructions_corr = [
        {"ss": False, "key": "ETA", "all": True, "scale": opt_scale, "flat": 0.0},
        {"ss": False, "key": "ZETA", "all": True, "scale": opt_scale, "flat": 0.0},
        {"ss": False, "key": "FSCAL", "all": True, "scale": opt_scale, "flat": 0.0},
        {"ss": False, "key": "AC", "all": True, "scale": opt_scale, "flat": 0.0},
        {"ss": False, "key": "AA", "all": True, "scale": opt_scale, "flat": 0.0},
        {"ss": False, "key": "AR", "all": True, "scale": opt_scale, "flat": 0.0},
        {"ss": False, "key": "ALPHA", "all": True, "scale": opt_scale, "flat": 0.0},
        {"ss": False, "key": "BETA", "all": True, "scale": opt_scale, "flat": 0.0},
        {"ss": False, "key": "GAMMA", "all": True, "scale": opt_scale, "flat": 0.0},
        {"ss": False, "key": "UUR", "scale": opt_scale, "flat": 0.0},
        {"ss": False, "key": "UUA", "scale": opt_scale, "flat": 0.0},
        {"ss": False, "key": "UUW", "scale": opt_scale, "flat": 0.0},
        {"ss": False, "key": "QPS1", "scale": opt_scale, "flat": 0.0},
        {"ss": False, "key": "QPS2", "scale": opt_scale, "flat": 0.0},
        {"ss": False, "key": "QSSS1", "scale": opt_scale, "flat": 0.0},
        {"ss": False, "key": "QSSS2", "scale": opt_scale, "flat": 0.0}
    ]

    # Generate SS instructions if SS objects exist
    instructions_ss = []
    if len(target.DK.SS) > 0:
        
        # SS float attributes that can be optimized (based on ss class)
        ss_float_attrs = [
            "SPU", "SPV", "SPR", "SPA", "SPB", "SPC", "SPK", "SPL",  # SP correlations
            "PPU", "PPV", "PPR", "PPA", "PPB", "PPC", "PPK", "PPL",  # PP correlations (if NPPART >= 2)
            "WSE", "WSV", "WSR", "WSA", "WBRHO", "WBALPH",            # Woods-Saxon parameters
            "SDU", "SDV", "SDR", "SDA", "SDB", "SDC", "SDK", "SDL",  # SD correlations (if NDPART != 0)
            "PDU", "PDV", "PDR", "PDA", "PDB", "PDC", "PDK", "PDL",  # PD correlations (if NDPART != 0)
            "DDU", "DDV", "DDR", "DDA", "DDB", "DDC", "DDK", "DDL",  # DD correlations (if NDPART >= 2)
            "D_WSE", "D_WSV", "D_WSR", "D_WSA", "D_WBRHO", "D_WBALPH"  # D Woods-Saxon (if NDPART != 0)
        ]
        
        for ss_idx, ss_obj in enumerate(target.DK.SS):
            for attr in ss_float_attrs:
                if hasattr(ss_obj, attr):
                    try:
                        value = getattr(ss_obj, attr)
                        # Check if it's a float
                        if isinstance(value, str) and "." in value:
                            instruction = {
                                "ss": True,
                                "ss_idx": ss_idx,
                                "key": attr,
                                "scale": opt_scale,
                                "flat": 0.0
                            }
                            instructions_ss.append(instruction)
                            # print(f"  {attr}: {value} -> Added to optimization")
                        else:
                            print(f"  {attr}: {value} (skipped - not float)")
                    except Exception as e:
                        print(f"  {attr}: <could not access: {e}>")

    # Combined instructions for final optimization
    instructions_all = instructions_esep + instructions_corr + instructions_ss
    
    # Generate ESEP-only optimization file
    opt_esep = GenerateOptFile(target.PARAMS, target.DK, opt_esep_file_name, instructions_esep)
    opt_esep.UpdateFloats(target.PARAMS, 6)
    print("... DONE")
    print(BREAK)
    
    # Step 3a: Optimize ESEP parameters only
    print("Step 3a: ESEP-ONLY PARAMETER OPTIMIZATION")
    dk_name_esep = f"'{work_dir}dk/{label}_esep_only.dk'"
    log_optimize_esep = f"{work_dir}logs/{label}.esep_only.optimize"
    print(f"BEGIN ESEP-ONLY OPTIMIZATION: {log_optimize_esep}.optimize")
    
    e_esep_only, v_esep_only = target.Optimize(opt_esep, dk_name_esep, True, log_optimize_esep)
    
    print(f" ⚛ ESEP-ONLY E = {e_esep_only:.4f} +- {v_esep_only:.4f}")
    print(f"ESEP-ONLY IMPROVEMENT: {e_esep_only - e_initial:.4f} MeV")
    print(f"ESEP VALUES AFTER ESEP-ONLY: {target.DK.ESEP}")
    print(BREAK)
    
    # Step 3b: Optimize ESEP + correlations
    print("Step 3b: ESEP + CORRELATIONS OPTIMIZATION")
    
    # Load the ESEP-optimized deck as the base for correlation optimization
    esep_optimized_deck = deck_t(target.PARAMS, dk_name_esep)
    
    opt_all_file_name = f"'{work_dir}opt/{label}_all.opt'"
    opt_all = GenerateOptFile(target.PARAMS, esep_optimized_deck, opt_all_file_name, instructions_all)
    opt_all.UpdateFloats(target.PARAMS, 6)
    
    dk_name_all = f"'{work_dir}dk/{label}_esep_plus_corr.dk'"
    log_optimize_all = f"{work_dir}logs/{label}.esep_plus_corr.optimize"
    print(f"BEGIN ESEP + CORRELATIONS OPTIMIZATION: {log_optimize_all}.optimize")
    
    e_optimized, v_optimized = target.Optimize(opt_all, dk_name_all, True, log_optimize_all)
    
    print(f" ⚛ FINAL OPTIMIZED E = {e_optimized:.4f} +- {v_optimized:.4f}")
    print(f"CORRELATIONS IMPROVEMENT: {e_optimized - e_esep_only:.4f} MeV")
    print(f"FINAL ESEP VALUES: {target.DK.ESEP}")
    print(BREAK)
    
    # Step 4: Final evaluation and logging
    print("Step 4: FINAL EVALUATION")
    total_improvement = e_optimized - e_initial
    esep_contribution = e_esep_only - e_initial
    corr_contribution = e_optimized - e_esep_only
    
    print(f"TOTAL ENERGY IMPROVEMENT: {total_improvement:.4f} MeV")
    print(f"  - ESEP contribution: {esep_contribution:.4f} MeV")
    print(f"  - Correlations contribution: {corr_contribution:.4f} MeV")
    print(f"OPTIMIZATION SCALES:")
    print(f"  - ESEP scale: {esep_scale}")
    print(f"  - OPT scale: {opt_scale}")
    
    # Step 5: Return results
    result = {
        "TYPE": "BOUND_STATE",
        "LABEL": label,
        "ESEP_INITIAL": [float(val) for val in esep_initial],
        "ESEP_FINAL": [float(val) for val in target.DK.ESEP],
        "E_INITIAL": e_initial,
        "V_INITIAL": v_initial,
        "E_ESEP_ONLY": e_esep_only,
        "V_ESEP_ONLY": v_esep_only,
        "E_OPTIMIZED": e_optimized,
        "V_OPTIMIZED": v_optimized,
        "TOTAL_IMPROVEMENT": total_improvement,
        "ESEP_CONTRIBUTION": esep_contribution,
        "CORR_CONTRIBUTION": corr_contribution,
        "ESEP_SCALE": esep_scale,
        "OPT_SCALE": opt_scale,
        "NUM_SS_INSTRUCTIONS": len(instructions_ss),
        "DECK_PATH_ESEP_ONLY": dk_name_esep,
        "DECK_PATH_FINAL": dk_name_all
    }
    
    print("OPTIMIZATION SUMMARY:")
    print(json.dumps(result, indent=2))
    print(BREAK)
    return result

def testBoundStateOpt(util: utility_t):
    """
    Test function for bound state optimization
    """
    # Set up working environment and get configured wavefunction
    target = setup(util)
    
    # Run bound state optimization
    label = "test_bound_state"
    result = BoundStateOptimize(target, util, label)
    
    # Save results to JSON file
    with open(f"{util.WORKING_DIR}{label}_results.json", "w") as f:
        json.dump(result, f, indent=2)
    
    return result
#-----------------------------------------------------------------------

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Test bound state optimization using ESEP parameters",
        epilog="Example: python3 optFileTest.py --utility examples/he4.util"
    )
    parser.add_argument('--utility', required=True, help="Path to utility file")
    args = parser.parse_args()
    
    print("="*72)
    print("BOUND STATE OPTIMIZATION TEST")
    print("="*72)
    
    util = utility_t(args.utility)
    result = testBoundStateOpt(util)
    
    print("="*72)
    print("TEST COMPLETED SUCCESSFULLY!")
    print("="*72)
#-----------------------------------------------------------------------
