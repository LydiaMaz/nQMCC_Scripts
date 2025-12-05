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
            "base_opt_scale": util.OPT_SCALE,           # Base values from util file
            "base_num_opt_evaluations": util.NUM_OPT_EVALUATIONS,
            "base_esep_scale": util.ESEP_SCALE
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
    
    # Get potential-specific parameters
    scales = get_potential_specific_scales(pot_pair_name, util)
    
    print(f"Setting up potential pair directory: {pot_pair_name}")
    print(f"  Using potential-specific parameters:")
    print(f"    {scales['multiplier_info']}")
    print(f"    ESEP scale: {scales['esep_scale']} (util: {util.ESEP_SCALE})")
    print(f"    OPT scale: {scales['opt_scale']} (util: {util.OPT_SCALE})")
    print(f"    Evaluations: {scales['num_evaluations']} (util: {util.NUM_OPT_EVALUATIONS})")
    
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
    target.CTRL.NUM_OPT_EVALUATIONS = scales['num_evaluations']  # Use potential-specific value
    
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

def get_potential_specific_scales(pot_pair_name, util):
    """
    Get optimization scales specific to potential family, using util file values as base.
    This allows control from util file while applying potential-specific multipliers.
    Uses more conservative scaling to avoid NaN generation.
    """
    base_esep_scale = util.ESEP_SCALE
    base_opt_scale = util.OPT_SCALE
    base_num_evaluations = util.NUM_OPT_EVALUATIONS
    
    if 'av18' in pot_pair_name.lower():
        # AV18 potentials: use util file values directly (already optimized)
        return {
            'esep_scale': base_esep_scale,
            'opt_scale': base_opt_scale,
            'num_evaluations': base_num_evaluations,
            'multi_stage': False,  # Single-stage optimization is sufficient
            'multiplier_info': 'AV18 (1.0x util values)'
        }
    elif 'nv2' in pot_pair_name.lower():
        # NV2 potentials: use ultra-conservative scaling to prevent gamma NaN issues
        # The gamma parameter is very sensitive with NV2 - use minimal scaling
        return {
            'esep_scale': base_esep_scale * 1.01,    # Ultra-conservative: only 1% increase
            'opt_scale': base_opt_scale * 1.005,     # Ultra-conservative: only 0.5% increase  
            'num_evaluations': base_num_evaluations,  # Keep same as util file
            'multi_stage': True,   # Use multi-stage optimization for gradual improvement
            'multiplier_info': 'NV2 (1.01x esep, 1.005x opt, 1.0x evals from util - ultra-conservative to prevent gamma NaN)'
        }
    else:
        # Default for unknown potentials: moderate scaling from util values
        return {
            'esep_scale': base_esep_scale * 1.2,
            'opt_scale': base_opt_scale * 1.1,
            'num_evaluations': int(base_num_evaluations * 1.2),
            'multi_stage': True,
            'multiplier_info': 'Other (1.2x esep, 1.1x opt, 1.2x evals from util)'
        }

def check_for_optimization_failure(energy, variance, deck_obj, step_name):
    """
    Check if optimization produced NaN values or failed in other ways.
    Returns True if optimization failed, False if successful.
    """
    # Check for NaN in energy or variance
    if energy is None or variance is None:
        print(f"❌ {step_name} FAILED: Energy or variance is None")
        return True
        
    if np.isnan(energy) or np.isnan(variance):
        print(f"❌ {step_name} FAILED: Energy or variance is NaN")
        return True
    
    # Check for NaN in ESEP parameters
    if hasattr(deck_obj, 'ESEP') and deck_obj.ESEP is not None:
        try:
            esep_values = deck_obj.ESEP
            esep_array = np.array(esep_values, dtype=float)
            if np.any(np.isnan(esep_array)):
                print(f"❌ {step_name} FAILED: ESEP contains NaN values: {esep_values}")
                return True
        except (TypeError, ValueError):
            print(f"⚠️ {step_name}: Could not validate ESEP parameters: {deck_obj.ESEP}")
    
    # Check for NaN in gamma parameters (critical for wavefunction stability)
    try:
        if hasattr(deck_obj, 'OPT_OBJECTS') or (hasattr(deck_obj, 'target') and hasattr(deck_obj.target, 'OPT_OBJECTS')):
            opt_objects = getattr(deck_obj, 'OPT_OBJECTS', None) or getattr(deck_obj.target, 'OPT_OBJECTS', None)
            if opt_objects:
                for i, obj in enumerate(opt_objects):
                    if hasattr(obj, 'gamma'):
                        try:
                            gamma_val = float(obj.gamma)
                            if np.isnan(gamma_val) or np.isinf(gamma_val):
                                print(f"❌ {step_name} FAILED: gamma[{i}] is NaN/Inf: {gamma_val}")
                                return True
                            if abs(gamma_val) > 100.0:  # Unreasonably large gamma
                                print(f"❌ {step_name} FAILED: gamma[{i}] too large: {gamma_val}")
                                return True
                        except (TypeError, ValueError):
                            print(f"❌ {step_name} FAILED: gamma[{i}] not numeric: {obj.gamma}")
                            return True
                    if hasattr(obj, 'beta'):
                        try:
                            beta_val = float(obj.beta)
                            if np.isnan(beta_val) or np.isinf(beta_val):
                                print(f"❌ {step_name} FAILED: beta[{i}] is NaN/Inf: {beta_val}")
                                return True
                        except (TypeError, ValueError):
                            print(f"❌ {step_name} FAILED: beta[{i}] not numeric: {obj.beta}")
                            return True
    except Exception as e:
        print(f"⚠️ {step_name}: Could not validate correlation parameters: {e}")
    
    # Check for unreasonably large variance (sign of instability)
    if variance > 1.0:  # Adjust threshold as needed
        print(f"⚠️  {step_name} WARNING: Large variance detected ({variance:.4f})")
        print("This may indicate optimization instability")
    
    return False


def validate_gamma_parameters(target, context=""):
    """
    Validate gamma parameters before optimization to prevent NaN generation.
    Returns True if valid, False if problematic.
    """
    try:
        if hasattr(target, 'OPT_OBJECTS'):
            for i, obj in enumerate(target.OPT_OBJECTS):
                if hasattr(obj, 'gamma'):
                    try:
                        gamma_val = float(obj.gamma)
                        if np.isnan(gamma_val) or np.isinf(gamma_val):
                            print(f"❌ {context}: Pre-optimization gamma[{i}] is NaN/Inf: {gamma_val}")
                            return False
                        if abs(gamma_val) > 50.0:  # Flag potentially problematic gamma
                            print(f"⚠️ {context}: Pre-optimization gamma[{i}] is large: {gamma_val}")
                    except (TypeError, ValueError):
                        print(f"❌ {context}: Pre-optimization gamma[{i}] not numeric: {obj.gamma}")
                        return False
        return True
    except Exception as e:
        print(f"⚠️ {context}: Could not validate gamma parameters: {e}")
        return True  # Assume valid if we can't check


def safe_optimize(target, opt_obj, deck_name, write_deck, log_name, step_name):
    """
    Safely perform optimization with parameter validation before and after
    """
    print(f"  🔧 {step_name}: Starting optimization...")
    
    # Pre-optimization gamma validation
    if not validate_gamma_parameters(target, step_name):
        print(f"❌ {step_name}: Pre-optimization gamma validation failed - aborting")
        return None, None
    
    # Validate parameters before optimization
    try:
        if hasattr(target, 'OPT_OBJECTS'):
            for i, obj in enumerate(target.OPT_OBJECTS):
                try:
                    if hasattr(obj, 'gamma'):
                        gamma_val = float(obj.gamma)
                        if np.isnan(gamma_val):
                            raise ValueError(f"Pre-optimization: gamma parameter {i} is NaN")
                    if hasattr(obj, 'beta'):
                        beta_val = float(obj.beta)
                        if np.isnan(beta_val):
                            raise ValueError(f"Pre-optimization: beta parameter {i} is NaN")
                except (TypeError, ValueError, AttributeError):
                    print(f"⚠️ Could not validate optimization object {i} parameters before {step_name}")
    except Exception as e:
        print(f"⚠️ Pre-optimization parameter check failed for {step_name}: {e}")
    
    # Perform optimization
    try:
        energy, variance = target.Optimize(opt_obj, deck_name, write_deck, log_name)
    except FileNotFoundError as e:
        print(f"❌ {step_name}: File not found during optimization: {e}")
        print(f"   Deck path: {deck_name}")
        print(f"   Log path: {log_name}")
        return None, None
    except PermissionError as e:
        print(f"❌ {step_name}: Permission denied during optimization: {e}")
        print(f"   Deck path: {deck_name}")
        print(f"   Log path: {log_name}")
        return None, None
    except Exception as e:
        print(f"❌ {step_name}: Unexpected error during optimization: {e}")
        print(f"   Error type: {type(e).__name__}")
        print(f"   Deck path: {deck_name}")
        print(f"   Log path: {log_name}")
        return None, None
    
    # Validate results - handle None values first
    if energy is None or variance is None:
        print(f"❌ {step_name}: Optimization returned None values: E={energy}, V={variance}")
        return None, None
    
    # Then check for NaN
    if np.isnan(energy) or np.isnan(variance):
        print(f"❌ {step_name}: Optimization produced NaN values: E={energy}, V={variance}")
        return None, None
    
    # Validate deck parameters after optimization
    if hasattr(target, 'DK') and hasattr(target.DK, 'ESEP'):
        try:
            esep_values = target.DK.ESEP
            esep_array = np.array(esep_values, dtype=float)
            if np.any(np.isnan(esep_array)):
                print(f"❌ {step_name}: Post-optimization ESEP contains NaN values: {esep_values}")
                return None, None
        except (TypeError, ValueError):
            print(f"⚠️ Could not validate ESEP parameters after {step_name}")
    
    # Post-optimization gamma validation
    if not validate_gamma_parameters(target, f"{step_name} POST"):
        print(f"❌ {step_name}: Post-optimization gamma validation failed")
        return None, None
    
    print(f"  ✓ {step_name}: E = {energy:.4f} +- {variance:.4f}")
    return energy, variance


def BoundStateOptimize(target: wavefunction_t, util: utility_t, pot_pair_name: str, pot_dir: str):
    """
    Optimize bound state wavefunction for a specific potential pair
    """
    BREAK="="*72
    
    # Get potential-specific optimization parameters
    scales = get_potential_specific_scales(pot_pair_name, util)
    opt_scale = scales['opt_scale']
    esep_scale = scales['esep_scale']
    num_evaluations = scales['num_evaluations']
    
    print(f" 🔥 BOUND STATE OPTIMIZATION: {pot_pair_name} 🔥")
    print(f"Using potential-specific scales:")
    print(f"  {scales['multiplier_info']}")
    print(f"  ESEP scale: {esep_scale}")
    print(f"  OPT scale: {opt_scale}")
    print(f"  Evaluations: {num_evaluations}")
    print(BREAK)
    
    # Validate initial parameters before optimization
    print("🔍 Checking initial wavefunction parameters...")
    if hasattr(target, 'DK') and hasattr(target.DK, 'ESEP'):
        initial_esep = target.DK.ESEP
        try:
            # Convert to numpy array and check for NaN/Inf
            esep_array = np.array(initial_esep, dtype=float)
            if np.any(np.isnan(esep_array)) or np.any(np.isinf(esep_array)):
                raise ValueError(f"Initial ESEP parameters contain NaN/Inf values: {initial_esep}")
            print(f"✓ Initial ESEP parameters are valid: {initial_esep}")
        except (TypeError, ValueError) as e:
            print(f"⚠️ Could not validate ESEP parameters as numeric: {initial_esep}")
            print(f"   Error: {e}")
            print("   Proceeding with caution...")
    
    # Check for basic wavefunction parameter validity
    try:
        # Basic parameter validation - attempt to access key attributes
        if hasattr(target, 'OPT_OBJECTS'):
            for i, obj in enumerate(target.OPT_OBJECTS):
                try:
                    if hasattr(obj, 'gamma'):
                        gamma_val = float(obj.gamma)
                        if np.isnan(gamma_val):
                            raise ValueError(f"Initial gamma parameter {i} is NaN")
                    if hasattr(obj, 'beta'):
                        beta_val = float(obj.beta)
                        if np.isnan(beta_val):
                            raise ValueError(f"Initial beta parameter {i} is NaN")
                except (TypeError, ValueError, AttributeError):
                    print(f"⚠️ Could not validate optimization object {i} parameters")
        print("✓ Initial optimization parameters appear valid")
    except Exception as e:
        print(f"⚠️ Warning during parameter validation: {e}")
        # Continue but with extra caution
    
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
    
    # Adaptive optimization: try standard scale first, then increase if needed
    e_esep_only, v_esep_only = safe_optimize(target, opt_esep, dk_name_esep, True, log_optimize_esep, "ESEP-ONLY")
    
    # Check for optimization failure
    if e_esep_only is None or v_esep_only is None or check_for_optimization_failure(e_esep_only, v_esep_only, target.DK, "ESEP-ONLY"):
        print("❌ ESEP optimization failed - attempting ultra-conservative recovery")
        
        # Try with much more conservative scale as fallback
        conservative_esep_scale = util.ESEP_SCALE * 0.25  # Ultra-ultra-conservative: quarter the util file value
        instructions_esep_conservative = [
            {"ss": False, "key": "ESEP", "idx": 0, "scale": conservative_esep_scale, "flat": 0.0},
            {"ss": False, "key": "ESEP", "idx": 1, "scale": conservative_esep_scale, "flat": 0.0},
            {"ss": False, "key": "ESEP", "idx": 2, "scale": conservative_esep_scale, "flat": 0.0},
            {"ss": False, "key": "ESEP", "idx": 3, "scale": conservative_esep_scale, "flat": 0.0},
        ]
        
        opt_esep_recovery_file = f"'{pot_dir}opt/{pot_pair_name}_esep_recovery.opt'"
        opt_esep_recovery = GenerateOptFile(target.PARAMS, target.DK, opt_esep_recovery_file, instructions_esep_conservative)
        opt_esep_recovery.UpdateFloats(target.PARAMS, 6)
        
        dk_name_esep_recovery = f"'{pot_dir}dk/{pot_pair_name}_esep_recovery.dk'"
        log_optimize_esep_recovery = f"{pot_dir}logs/{pot_pair_name}.esep_recovery.optimize"
        print(f"🔧 RECOVERY ESEP OPTIMIZATION (scale={conservative_esep_scale:.1f}): {log_optimize_esep_recovery}")
        
        e_esep_recovery, v_esep_recovery = safe_optimize(target, opt_esep_recovery, dk_name_esep_recovery, True, log_optimize_esep_recovery, "ESEP-RECOVERY")
        
        if e_esep_recovery is not None and v_esep_recovery is not None and not check_for_optimization_failure(e_esep_recovery, v_esep_recovery, target.DK, "ESEP-RECOVERY"):
            print("✅ Recovery optimization successful")
            e_esep_only, v_esep_only = e_esep_recovery, v_esep_recovery
            dk_name_esep = dk_name_esep_recovery
            esep_scale = conservative_esep_scale  # Update for metadata
        else:
            raise RuntimeError(f"ESEP optimization failed completely for {pot_pair_name} - unable to recover")
    
    # Ensure we have valid values before proceeding
    if e_esep_only is None or v_esep_only is None:
        raise RuntimeError(f"ESEP optimization returned None for {pot_pair_name}")
    
    initial_esep_improvement = abs(e_esep_only - e_initial)
    
    # If improvement is small for NV2 potentials, try more aggressive optimization
    if 'nv2' in pot_pair_name.lower() and initial_esep_improvement < 0.5:
        print(f"⚠️  Small ESEP improvement ({initial_esep_improvement:.3f} MeV) detected for NV2 potential")
        print("🔄 Trying more aggressive ESEP optimization...")
        
        # Create more aggressive optimization file (but not too aggressive to avoid NaNs)
        aggressive_esep_scale = esep_scale * 1.5  # Modest increase from already conservative scale
        instructions_esep_aggressive = [
            {"ss": False, "key": "ESEP", "idx": 0, "scale": aggressive_esep_scale, "flat": 0.0},
            {"ss": False, "key": "ESEP", "idx": 1, "scale": aggressive_esep_scale, "flat": 0.0},
            {"ss": False, "key": "ESEP", "idx": 2, "scale": aggressive_esep_scale, "flat": 0.0},
            {"ss": False, "key": "ESEP", "idx": 3, "scale": aggressive_esep_scale, "flat": 0.0},
        ]
        
        opt_esep_aggressive_file = f"'{pot_dir}opt/{pot_pair_name}_esep_aggressive.opt'"
        opt_esep_aggressive = GenerateOptFile(target.PARAMS, target.DK, opt_esep_aggressive_file, instructions_esep_aggressive)
        opt_esep_aggressive.UpdateFloats(target.PARAMS, 6)
        
        dk_name_esep_aggressive = f"'{pot_dir}dk/{pot_pair_name}_esep_aggressive.dk'"
        log_optimize_esep_aggressive = f"{pot_dir}logs/{pot_pair_name}.esep_aggressive.optimize"
        print(f"AGGRESSIVE ESEP OPTIMIZATION (scale={aggressive_esep_scale:.1f}): {log_optimize_esep_aggressive}")
        
        e_esep_aggressive, v_esep_aggressive = safe_optimize(target, opt_esep_aggressive, dk_name_esep_aggressive, True, log_optimize_esep_aggressive, "ESEP-AGGRESSIVE")
        
        # Check for failure in aggressive optimization
        if e_esep_aggressive is None or v_esep_aggressive is None or check_for_optimization_failure(e_esep_aggressive, v_esep_aggressive, target.DK, "ESEP-AGGRESSIVE"):
            print("❌ Aggressive optimization failed (NaN detected), keeping standard result")
            print(f"✓ Using stable standard result: E = {e_esep_only:.4f} +- {v_esep_only:.4f}")
            # Keep the standard result (e_esep_only, v_esep_only, dk_name_esep already set correctly)
        else:
            aggressive_improvement = abs(e_esep_aggressive - e_initial)
            print(f"Standard improvement: {initial_esep_improvement:.3f} MeV")
            print(f"Aggressive improvement: {aggressive_improvement:.3f} MeV")
            
            if e_esep_aggressive < e_esep_only:  # Use aggressive result if better
                print("✅ Aggressive optimization successful, using aggressive result")
                e_esep_only, v_esep_only = e_esep_aggressive, v_esep_aggressive
                dk_name_esep = dk_name_esep_aggressive
                esep_scale = aggressive_esep_scale  # Update for metadata
            else:
                print("⚠️  Aggressive optimization did not improve, keeping standard result")
                print(f"Standard: {e_esep_only:.4f} MeV vs Aggressive: {e_esep_aggressive:.4f} MeV")
    
    print(f" ⚛ ESEP-ONLY E = {e_esep_only:.4f} +- {v_esep_only:.4f}")
    print(f"ESEP-ONLY IMPROVEMENT: {e_esep_only - e_initial:.4f} MeV")
    print(f"ESEP VALUES AFTER ESEP-ONLY: {target.DK.ESEP}")
    
    # Final safety check: ensure we have valid ESEP results before proceeding
    if e_esep_only is None or v_esep_only is None:
        raise RuntimeError(f"ESEP optimization returned invalid results for {pot_pair_name}")
    
    # Multi-stage ESEP optimization for difficult potentials
    if scales['multi_stage'] and abs(e_esep_only - e_initial) > 0.5:
        print("🔄 APPLYING MULTI-STAGE ESEP OPTIMIZATION (large improvement detected)")
        
        # Stage 2: More aggressive ESEP optimization with larger scale (but not too large)
        stage2_esep_scale = esep_scale * 1.2  # Conservative increase
        
        # For NV2 potentials, be even more careful with stage 2 scaling
        if 'nv2' in pot_pair_name.lower():
            stage2_esep_scale = esep_scale * 1.05  # Very small increase for NV2
            print(f"🔧 Using extra-conservative stage 2 scaling for NV2: {stage2_esep_scale:.3f}")
        
        instructions_esep_stage2 = [
            {"ss": False, "key": "ESEP", "idx": 0, "scale": stage2_esep_scale, "flat": 0.0},
            {"ss": False, "key": "ESEP", "idx": 1, "scale": stage2_esep_scale, "flat": 0.0},
            {"ss": False, "key": "ESEP", "idx": 2, "scale": stage2_esep_scale, "flat": 0.0},
            {"ss": False, "key": "ESEP", "idx": 3, "scale": stage2_esep_scale, "flat": 0.0},
        ]
        
        # CRITICAL: Load the Stage 1 optimized deck as starting point for Stage 2
        stage1_deck_path = dk_name_esep.strip("'\"")  # Remove quotes
        if os.path.exists(stage1_deck_path):
            print(f"✓ Loading Stage 1 deck for Stage 2: {stage1_deck_path}")
            stage1_optimized_deck = deck_t(target.PARAMS, dk_name_esep)
            
            # CRUCIAL: Update the target wavefunction with Stage 1 optimized parameters
            # This ensures gamma/beta parameters start from Stage 1 values, not original values
            target.DK = stage1_optimized_deck
            print(f"✓ Updated target wavefunction with Stage 1 parameters")
            print(f"   Stage 1 ESEP: {target.DK.ESEP}")
            
            opt_esep_stage2_file = f"'{pot_dir}opt/{pot_pair_name}_esep_stage2.opt'"
            opt_esep_stage2 = GenerateOptFile(target.PARAMS, stage1_optimized_deck, opt_esep_stage2_file, instructions_esep_stage2)
            opt_esep_stage2.UpdateFloats(target.PARAMS, 6)
            
            dk_name_esep_stage2 = f"'{pot_dir}dk/{pot_pair_name}_esep_stage2.dk'"
            log_optimize_esep_stage2 = f"{pot_dir}logs/{pot_pair_name}.esep_stage2.optimize"
            print(f"STAGE 2 ESEP OPTIMIZATION (scale={stage2_esep_scale:.1f}): {log_optimize_esep_stage2}")
            
            e_esep_stage2, v_esep_stage2 = safe_optimize(target, opt_esep_stage2, dk_name_esep_stage2, True, log_optimize_esep_stage2, "ESEP-STAGE2")
        else:
            print(f"❌ Stage 1 deck not found: {stage1_deck_path}")
            print("⚠️ Skipping Stage 2 optimization")
            # Set stage 2 results to None to indicate failure
            e_esep_stage2, v_esep_stage2 = None, None
        
        
        # Check for failure in stage 2 optimization
        if e_esep_stage2 is None or v_esep_stage2 is None or check_for_optimization_failure(e_esep_stage2, v_esep_stage2, target.DK, "ESEP-STAGE2"):
            print("❌ Stage 2 optimization failed (NaN detected), keeping stage 1 result")
            print(f"✓ Falling back to stable Stage 1: E = {e_esep_only:.4f} +- {v_esep_only:.4f}")
            # Keep dk_name_esep as the stage 1 deck (already set correctly)
        elif e_esep_stage2 < e_esep_only:  # Only keep if improvement
            print(f" ⚛ STAGE 2 E = {e_esep_stage2:.4f} +- {v_esep_stage2:.4f}")
            print(f"STAGE 2 ADDITIONAL IMPROVEMENT: {e_esep_stage2 - e_esep_only:.4f} MeV")
            print("✅ Stage 2 successful, using Stage 2 result")
            e_esep_only, v_esep_only = e_esep_stage2, v_esep_stage2
            dk_name_esep = dk_name_esep_stage2  # Use stage 2 deck for correlation opt
        else:
            print("⚠️  Stage 2 did not improve energy, keeping stage 1 result")
            print(f"Stage 1: {e_esep_only:.4f} MeV vs Stage 2: {e_esep_stage2:.4f} MeV")
            # Keep dk_name_esep as the stage 1 deck (already set correctly)
    
    print(BREAK)
    
    # Step 3b: Optimize ESEP + correlations
    print("Step 3b: ESEP + CORRELATIONS OPTIMIZATION")
    
    # Verify the ESEP deck exists before trying to load it
    esep_deck_path = dk_name_esep.strip("'\"")  # Remove quotes if present
    print(f"🔍 Checking for ESEP deck: {esep_deck_path}")
    
    if not os.path.exists(esep_deck_path):
        print(f"❌ ESEP deck not found: {esep_deck_path}")
        print("🔧 Using current target wavefunction state for correlations")
        # Use the current target.DK which should have the ESEP-optimized parameters
        esep_optimized_deck = target.DK
        
        # Write the current state to a fallback deck file for consistency
        fallback_deck_path = f"'{pot_dir}dk/{pot_pair_name}_current_state.dk'"
        print(f"💾 Saving current wavefunction state: {fallback_deck_path}")
        target.DK.WriteDeck(fallback_deck_path.strip("'\""))
        dk_name_for_correlations = fallback_deck_path
    else:
        # Load the ESEP-optimized deck as the base for correlation optimization
        try:
            print(f"📖 Loading ESEP-optimized deck: {esep_deck_path}")
            esep_optimized_deck = deck_t(target.PARAMS, dk_name_esep)
            dk_name_for_correlations = dk_name_esep
            print(f"✓ Successfully loaded ESEP-optimized deck")
        except Exception as e:
            print(f"❌ Failed to load ESEP deck: {e}")
            print("🔧 Using current target wavefunction state as fallback")
            esep_optimized_deck = target.DK
            
            # Write the current state to a fallback deck file
            fallback_deck_path = f"'{pot_dir}dk/{pot_pair_name}_fallback_state.dk'"
            print(f"💾 Saving fallback state: {fallback_deck_path}")
            target.DK.WriteDeck(fallback_deck_path.strip("'\""))
            dk_name_for_correlations = fallback_deck_path
    
    opt_all_file_name = f"'{pot_dir}opt/{pot_pair_name}_all.opt'"
    opt_all = GenerateOptFile(target.PARAMS, esep_optimized_deck, opt_all_file_name, instructions_all)
    opt_all.UpdateFloats(target.PARAMS, 6)
    
    dk_name_all = f"'{pot_dir}dk/{pot_pair_name}_esep_plus_corr.dk'"
    log_optimize_all = f"{pot_dir}logs/{pot_pair_name}.esep_plus_corr.optimize"
    print(f"BEGIN ESEP + CORRELATIONS OPTIMIZATION: {log_optimize_all}.optimize")
    print(f"🎯 Final deck will be saved to: {dk_name_all}")
    print(f"🔍 Current working directory: {os.getcwd()}")
    
    # Debug: Check if the dk directory exists and is writable
    dk_dir = f"{pot_dir}dk"
    print(f"🔍 Checking dk directory: {dk_dir}")
    print(f"   Exists: {os.path.exists(dk_dir)}")
    print(f"   Is directory: {os.path.isdir(dk_dir)}")
    print(f"   Writable: {os.access(dk_dir, os.W_OK)}")
    
    # Also check the absolute path
    abs_dk_dir = os.path.abspath(dk_dir)
    print(f"🔍 Absolute dk directory: {abs_dk_dir}")
    print(f"   Abs exists: {os.path.exists(abs_dk_dir)}")
    print(f"   Abs writable: {os.access(abs_dk_dir, os.W_OK)}")
    
    # Ensure the dk directory exists
    os.makedirs(dk_dir, exist_ok=True)
    
    # Test write a dummy file to check permissions
    test_file = f"{dk_dir}/test_write.tmp"
    try:
        with open(test_file, 'w') as f:
            f.write("test")
        os.remove(test_file)
        print(f"✓ Directory write test passed")
    except Exception as e:
        print(f"❌ Directory write test failed: {e}")
        raise RuntimeError(f"Cannot write to dk directory {dk_dir}: {e}")
    
    # Try using relative path for the deck name (QMC might expect relative paths)
    dk_name_all_relative = f"dk/{pot_pair_name}_esep_plus_corr.dk"
    print(f"🔄 Using relative deck path: {dk_name_all_relative}")
    
    e_optimized, v_optimized = safe_optimize(target, opt_all, dk_name_all_relative, True, log_optimize_all, "ESEP+CORR")
    
    # Check for failure in final optimization
    if e_optimized is None or v_optimized is None or check_for_optimization_failure(e_optimized, v_optimized, target.DK, "FINAL"):
        print("❌ Final optimization failed (NaN detected), using ESEP-only result")
        e_optimized, v_optimized = e_esep_only, v_esep_only
        print(f"⚠️  Falling back to ESEP-only result: E = {e_optimized:.4f} +- {v_optimized:.4f}")
    else:
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
        "NUM_OPT_EVALUATIONS": num_evaluations,
        "BASE_ESEP_SCALE": util.ESEP_SCALE,      # Original util file values
        "BASE_OPT_SCALE": util.OPT_SCALE,
        "BASE_NUM_OPT_EVALUATIONS": util.NUM_OPT_EVALUATIONS,
        "SCALE_MULTIPLIER_INFO": scales['multiplier_info'],
        "OPT_SCALE": opt_scale,
        "NUM_SS_INSTRUCTIONS": len(instructions_ss),
        "DECK_PATH_ESEP_ONLY": dk_name_esep,
        "DECK_PATH_FINAL": dk_name_all_relative,  # Use the relative path that was actually used
        "DECK_PATH_FOR_CORRELATIONS": dk_name_for_correlations  # Track which deck was actually used
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


base_esep_scale = util.ESEP_SCALE   # Start exactly from util file

trial_multipliers = [1.0, 1.5, 2.0, 3.0]
best = None

print("="*72)
print(f"NV ADAPTIVE ESEP: base util scale = {base_esep_scale}")
print(f"Trial multipliers: {trial_multipliers}")
print("="*72)

for m in trial_multipliers:
    current_scale = base_esep_scale * m
    label = f"nv_m{m:.1f}"

    instructions_esep = [
        {"ss": False, "key": "ESEP", "idx": 0, "scale": current_scale, "flat": 0.0},
        {"ss": False, "key": "ESEP", "idx": 1, "scale": current_scale, "flat": 0.0},
        {"ss": False, "key": "ESEP", "idx": 2, "scale": current_scale, "flat": 0.0},
        {"ss": False, "key": "ESEP", "idx": 3, "scale": current_scale, "flat": 0.0},
    ]

    opt_esep_file_name = f"'{pot_dir}opt/{pot_pair_name}_esep_{label}.opt'"
    dk_name_esep = f"'{pot_dir}dk/{pot_pair_name}_esep_{label}.dk'"
    log_optimize_esep = f"{pot_dir}logs/{pot_pair_name}.esep_{label}.optimize"

    print(f" → Trying scale = {current_scale:.4f} ({label})")

    opt_esep = GenerateOptFile(target.PARAMS, target.DK, opt_esep_file_name, instructions_esep)
    opt_esep.UpdateFloats(target.PARAMS, 6)

    e_try, v_try = safe_optimize(
        target,
        opt_esep,
        dk_name_esep,
        True,
        log_optimize_esep,
        f"ESEP-NV-{label}",
    )

    if e_try is None:
        print("   ✖ Invalid result, skipping.")
        continue

    print(f"   ✓ E = {e_try:.4f}")

    # stop when energy increases → minimum found
    if best is not None and e_try > best["E"]:
        print(f"   ↑ Energy increased at m={m}. Minimum reached — stopping NV scan.")
        break

    best = {
        "E": e_try,
        "V": v_try,
        "scale": current_scale,
        "multiplier": m,
        "dk_name": dk_name_esep,
    }
