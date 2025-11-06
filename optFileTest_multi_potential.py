#-----------------------------------------------------------------------
import numpy as np
import json
import sys
import os
from datetime import datetime
from pathlib import Path
#-----------------------------------------------------------------------
sys.path.append(os.path.join(os.path.dirname(__file__), "src"))
#-----------------------------------------------------------------------
from utility import utility_t
from deck import deck_t,GenerateOptFile
from wavefunction import wavefunction_t
#-----------------------------------------------------------------------

def create_potential_pair_name(pot2b, pot3b):
    """Create a clean directory name from potential pair"""
    pot2b_clean = pot2b.split(".")[0]
    pot3b_clean = pot3b.split(".")[0]
    return f"{pot2b_clean}.{pot3b_clean}"

def setup_main_directory(util: utility_t):
    """Set up the main working directory structure"""
    BREAK="="*72
    print("SETTING UP MAIN WORKING ENVIRONMENT")
    
    # Create main working directory
    try:
        os.mkdir(util.WORKING_DIR)
        print(f"Created main directory: {util.WORKING_DIR}")
    except FileExistsError:
        print("***MAIN WORKING DIRECTORY EXISTS***")
        util.WORKING_DIR = f"{util.WORKING_DIR}{util.NAME}-{datetime.now().strftime('%Y-%m-%d_%H-%M')}/"
        print(f"Using timestamped directory: {util.WORKING_DIR}")
        os.mkdir(util.WORKING_DIR)
    
    # Save run metadata
    run_metadata = {
        "timestamp": datetime.now().isoformat(),
        "nucleus": util.NAME,
        "system_type": util.SYSTEM_TYPE,
        "num_potential_pairs": util.NUM_POTS,
        "potential_pairs": [
            {
                "index": i,
                "2b_potential": util.TWO_BODY_FILES[i],
                "3b_potential": util.THREE_BODY_FILES[i],
                "constants": util.CONSTANTS_FILES[i],
                "pair_name": create_potential_pair_name(util.TWO_BODY_FILES[i], util.THREE_BODY_FILES[i])
            }
            for i in range(util.NUM_POTS)
        ],
        "optimization_settings": {
            "num_blocks": util.NUM_BLOCKS,
            "block_size": util.BLOCK_SIZE,
            "walkers_per_node": util.WALKERS_PER_NODE,
            "opt_scale": util.OPT_SCALE,
            "num_opt_evaluations": util.NUM_OPT_EVALUATIONS,
            "esep_scale": util.ESEP_SCALE
        }
    }
    
    with open(f"{util.WORKING_DIR}run_metadata.json", "w") as f:
        json.dump(run_metadata, f, indent=2)
    
    print("... DONE")
    print(BREAK)
    return run_metadata

def setup_potential_pair_directory(util: utility_t, pot_pair_name: str, pot_index: int):
    """Set up directory structure for a specific potential pair"""
    pot_dir = f"{util.WORKING_DIR}{pot_pair_name}/"
    
    print(f"Setting up potential pair directory: {pot_pair_name}")
    os.makedirs(pot_dir, exist_ok=True)
    os.makedirs(f"{pot_dir}ctrl", exist_ok=True)
    os.makedirs(f"{pot_dir}logs", exist_ok=True) 
    os.makedirs(f"{pot_dir}dk", exist_ok=True)
    os.makedirs(f"{pot_dir}opt", exist_ok=True)
    
    # Create and configure wavefunction for this potential pair
    os.makedirs(f"{util.NQMCC_DIR}walks", exist_ok=True)
    target = wavefunction_t(util.CTRL_FILE, util.NQMCC_DIR, util.BIN_DIR, util.RUN_CMD)
    target_label = target.DK.NAME.strip("\'")
    
    # Set up control file for this potential pair (write into the pair's ctrl/ directory)
    target.CTRL.FILE_NAME = f"{pot_dir}ctrl/target.ctrl"
    target.CTRL.NUM_BLOCKS = util.NUM_BLOCKS
    target.CTRL.BLOCK_SIZE = util.BLOCK_SIZE
    target.CTRL.WALKERS_PER_NODE = util.WALKERS_PER_NODE
    target.CTRL.NUM_OPT_EVALUATIONS = util.NUM_OPT_EVALUATIONS
    
    # Configure potentials for this pair
    const_file = util.CONSTANTS_FILES[pot_index]
    pot2b_file = util.TWO_BODY_FILES[pot_index]
    pot3b_file = util.THREE_BODY_FILES[pot_index]
    
    target.CTRL.CONST_FILE = f"'{util.NQMCC_DIR}constants/{const_file}'"
    target.CTRL.L2BP_FILE = f"'{util.NQMCC_DIR}pots/{pot2b_file}'"
    target.CTRL.L3BP_FILE = f"'{util.NQMCC_DIR}pots/{pot3b_file}'"
    
    print(f"  CONSTANTS: {const_file}")
    print(f"  V2B: {pot2b_file}")
    print(f"  V3B: {pot3b_file}")
    
    return target, pot_dir

def BoundStateOptimize(target: wavefunction_t, util: utility_t, pot_pair_name: str, pot_dir: str):
    """
    Optimize bound state wavefunction for a specific potential pair
    """
    BREAK="="*72
    opt_scale = util.OPT_SCALE
    esep_scale = util.ESEP_SCALE
    
    print(f" 🔥 BOUND STATE OPTIMIZATION: {pot_pair_name} 🔥")
    print(BREAK)
    
    # Step 1: Initial evaluation
    log_initial = f"{pot_dir}logs/{pot_pair_name}.esep.initial"
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
    opt_esep_file_name = f"'{pot_dir}opt/{pot_pair_name}_esep.opt'"
    
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
        ss_float_attrs = [
            "SPU", "SPV", "SPR", "SPA", "SPB", "SPC", "SPK", "SPL",  # SP correlations
            "PPU", "PPV", "PPR", "PPA", "PPB", "PPC", "PPK", "PPL",  # PP correlations
            "WSE", "WSV", "WSR", "WSA", "WBRHO", "WBALPH",            # Woods-Saxon parameters
            "SDU", "SDV", "SDR", "SDA", "SDB", "SDC", "SDK", "SDL",  # SD correlations
            "PDU", "PDV", "PDR", "PDA", "PDB", "PDC", "PDK", "PDL",  # PD correlations
            "DDU", "DDV", "DDR", "DDA", "DDB", "DDC", "DDK", "DDL",  # DD correlations
            "D_WSE", "D_WSV", "D_WSR", "D_WSA", "D_WBRHO", "D_WBALPH"  # D Woods-Saxon
        ]
        
        for ss_idx, ss_obj in enumerate(target.DK.SS):
            for attr in ss_float_attrs:
                if hasattr(ss_obj, attr):
                    try:
                        value = getattr(ss_obj, attr)
                        if isinstance(value, str) and "." in value:
                            instruction = {
                                "ss": True,
                                "ss_idx": ss_idx,
                                "key": attr,
                                "scale": opt_scale,
                                "flat": 0.0
                            }
                            instructions_ss.append(instruction)
                    except Exception as e:
                        pass

    # Combined instructions for final optimization
    instructions_all = instructions_esep + instructions_corr + instructions_ss
    
    # Generate ESEP-only optimization file
    opt_esep = GenerateOptFile(target.PARAMS, target.DK, opt_esep_file_name, instructions_esep)
    opt_esep.UpdateFloats(target.PARAMS, 6)
    print("... DONE")
    print(BREAK)
    
    # Step 3a: Optimize ESEP parameters only
    print("Step 3a: ESEP-ONLY PARAMETER OPTIMIZATION")
    dk_name_esep = f"'{pot_dir}dk/{pot_pair_name}_esep_only.dk'"
    log_optimize_esep = f"{pot_dir}logs/{pot_pair_name}.esep_only.optimize"
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
    
    opt_all_file_name = f"'{pot_dir}opt/{pot_pair_name}_all.opt'"
    opt_all = GenerateOptFile(target.PARAMS, esep_optimized_deck, opt_all_file_name, instructions_all)
    opt_all.UpdateFloats(target.PARAMS, 6)
    
    dk_name_all = f"'{pot_dir}dk/{pot_pair_name}_esep_plus_corr.dk'"
    log_optimize_all = f"{pot_dir}logs/{pot_pair_name}.esep_plus_corr.optimize"
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
    
    # Step 5: Return results
    result = {
        "TYPE": "BOUND_STATE",
        "POTENTIAL_PAIR": pot_pair_name,
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
    
    # Save individual result
    with open(f"{pot_dir}{pot_pair_name}_results.json", "w") as f:
        json.dump(result, f, indent=2)
    
    print("OPTIMIZATION SUMMARY:")
    print(json.dumps(result, indent=2))
    print(BREAK)
    return result

def run_multi_potential_optimization(util: utility_t):
    """
    Run optimization for all potential pairs
    """
    # Set up main directory and get metadata
    run_metadata = setup_main_directory(util)
    
    all_results = []
    
    # Process each potential pair
    for i in range(util.NUM_POTS):
        pot_pair_name = create_potential_pair_name(util.TWO_BODY_FILES[i], util.THREE_BODY_FILES[i])
        
        print(f"\n{'='*72}")
        print(f"PROCESSING POTENTIAL PAIR {i+1}/{util.NUM_POTS}: {pot_pair_name}")
        print(f"{'='*72}")
        
        # Set up directory and wavefunction for this potential pair
        target, pot_dir = setup_potential_pair_directory(util, pot_pair_name, i)

        # If an old target.ctrl exists at the pair root (from previous runs), move it into the
        # newly created ctrl/ folder so the control files are consistently stored there.
        old_target = f"{pot_dir}target.ctrl"
        new_target = f"{pot_dir}ctrl/target.ctrl"
        try:
            if os.path.exists(old_target) and not os.path.exists(new_target):
                os.replace(old_target, new_target)
                print(f"Moved existing target.ctrl into {pot_pair_name}/ctrl/")
        except Exception as e:
            print(f"Warning: could not move existing target.ctrl for {pot_pair_name}: {e}")
        
        # Change to potential pair directory
        original_cwd = os.getcwd()
        os.chdir(pot_dir)
        
        try:
            # Run optimization for this potential pair
            result = BoundStateOptimize(target, util, pot_pair_name, pot_dir)
            all_results.append(result)
            
            print(f"✅ Successfully completed optimization for {pot_pair_name}")
            
        except Exception as e:
            error_result = {
                "TYPE": "BOUND_STATE",
                "POTENTIAL_PAIR": pot_pair_name,
                "ERROR": str(e),
                "STATUS": "FAILED"
            }
            all_results.append(error_result)
            print(f"❌ Failed optimization for {pot_pair_name}: {e}")
            
        finally:
            # Return to original directory
            os.chdir(original_cwd)
    
    # Create combined results summary
    combined_results = {
        "run_metadata": run_metadata,
        "timestamp_completed": datetime.now().isoformat(),
        "total_potential_pairs": util.NUM_POTS,
        "successful_optimizations": len([r for r in all_results if "ERROR" not in r]),
        "failed_optimizations": len([r for r in all_results if "ERROR" in r]),
        "individual_results": all_results,
        "energy_comparison": [
            {
                "potential_pair": r["POTENTIAL_PAIR"],
                "initial_energy": r.get("E_INITIAL", None),
                "final_energy": r.get("E_OPTIMIZED", None),
                "total_improvement": r.get("TOTAL_IMPROVEMENT", None)
            }
            for r in all_results if "ERROR" not in r
        ]
    }
    
    # Save combined results
    with open(f"{util.WORKING_DIR}combined_results.json", "w") as f:
        json.dump(combined_results, f, indent=2)
    
    # Print summary
    print(f"\n{'='*72}")
    print("MULTI-POTENTIAL OPTIMIZATION SUMMARY")
    print(f"{'='*72}")
    print(f"Total potential pairs processed: {util.NUM_POTS}")
    print(f"Successful optimizations: {combined_results['successful_optimizations']}")
    print(f"Failed optimizations: {combined_results['failed_optimizations']}")
    
    if combined_results['energy_comparison']:
        print(f"\nEnergy improvements:")
        for comp in combined_results['energy_comparison']:
            if comp['total_improvement'] is not None:
                print(f"  {comp['potential_pair']}: {comp['total_improvement']:.4f} MeV")
    
    print(f"\nResults saved in: {util.WORKING_DIR}")
    print(f"{'='*72}")
    
    return combined_results

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Test bound state optimization for multiple potential pairs",
        epilog="Example: python3 optFileTest_multi_potential.py --utility examples/he4.util"
    )
    parser.add_argument('--utility', required=True, help="Path to utility file")
    args = parser.parse_args()
    
    print("="*72)
    print("MULTI-POTENTIAL BOUND STATE OPTIMIZATION")
    print("="*72)
    
    util = utility_t(args.utility)
    results = run_multi_potential_optimization(util)
    
    print("="*72)
    print("ALL OPTIMIZATIONS COMPLETED!")
    print("="*72)