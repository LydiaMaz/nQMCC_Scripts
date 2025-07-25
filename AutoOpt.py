#!/usr/bin/env python
import argparse
import copy
import subprocess
import re
import numpy as np
import os
import sys
import shutil
from pathlib import Path
from scipy.optimize import minimize_scalar
#-----------------------------------------------------------------------
src_path = os.path.join(os.path.dirname(__file__), "src")
sys.path.append(src_path)
#-----------------------------------------------------------------------
from control import Control
from deck import read_params_and_deck
from zero_params import zero_var_params, save_opt_file
from nuc_system import NuclearSystem
from utility import Utility
#-----------------------------------------------------------------------
"""
todo
"""
def change_working_directory(working_dir):
    os.chdir(working_dir)
    print(f"Changed working directory to: {os.getcwd()}")
#-----------------------------------------------------------------------
"""
Create necessary directories (`min`, `opt`, `scratch`) in the current working directory.
"""
def setup_environment(working_dir):
    change_working_directory(working_dir)
    cwd = Path.cwd()
    (cwd / "logs").mkdir(exist_ok=True) # save output from commands
    (cwd / "min").mkdir(exist_ok=True) # directory to store optimized decks
    (cwd / "opt").mkdir(exist_ok=True) # directory to store opt files
    (cwd / "scratch").mkdir(exist_ok=True) # directory to store temp files
#-----------------------------------------------------------------------
"""
todo
"""
def load_system(util_file, nuclear_system_file):
    util = Utility(util_file)
    setup_environment(util.working_dir)
    util.copy_control_files()  # Copy control files to working directory
#-----------------------------------------------------------------------
    nuclear_system = NuclearSystem(nuclear_system_file)
#-----------------------------------------------------------------------
    target_control = util.target_control
    scattering_control = util.scattering_control
#-----------------------------------------------------------------------
# logic for looping through 2b_potentials
# logic for changing 3b_potentials
#-----------------------------------------------------------------------
    existing_path = target_control.parameters['3b_file'].strip("'\"")
    new_filename = nuclear_system.parameters['3b_potentials'][0]
    print(f"Updating 3b_file from {existing_path} to {new_filename}")
    updated_3b_file = os.path.join(os.path.dirname(existing_path), new_filename)
    target_control.parameters['3b_file'] = f"'{updated_3b_file}'"
    scattering_control.parameters['3b_file'] = f"'{updated_3b_file}'"
    target_control.write_control()
    scattering_control.write_control()
#-----------------------------------------------------------------------
# logic for looping through channels
#-----------------------------------------------------------------------
    return util, nuclear_system
#-----------------------------------------------------------------------
"""
AutoOpt Program
"""
def main(util: Utility, system: NuclearSystem):
#-----------------------------------------------------------------------
    WSE_SHIFT=0.5 #these needs to be in system file
    BSCAT_INIT_MAX=15
    MAX_BSCAT_SLOPE=2
#-----------------------------------------------------------------------
# Setup nQMCC Inputs
#-----------------------------------------------------------------------
    BIN_PATH = util.binary_dir
    cmd = util.run_cmd
#-----------------------------------------------------------------------
    target_control = util.target_control
    target_ctrl = os.path.join(util.working_dir, target_control.Name + ".ctrl")
    scattering_control = util.scattering_control
    scattering_ctrl = os.path.join(util.working_dir, scattering_control.Name + ".ctrl")
#-----------------------------------------------------------------------
    target_param_file = target_control.parameters['wf'].param_file.strip("'\"")
#-----------------------------------------------------------------------
# STEP 1: Optimize target deck (he4) to get optimized parameters
#-----------------------------------------------------------------------
    target_energy, target_energy_var = optimize_deck(target_ctrl, cmd, BIN_PATH) # saves opt deck to min/
    print(f"Optimized target energy: {target_energy:.4f} MeV, variance: {target_energy_var:.4f}")
    opt_target_deck, _ = read_params_and_deck(target_param_file,\
      target_control.parameters['optimized_deck'].strip("'\""))
#-----------------------------------------------------------------------
# make temp file for working deck (he5)
#-----------------------------------------------------------------------
    wf = scattering_control.parameters.get('wf', None)
    if wf is None:
        raise RuntimeError("No wavefunction block found in control file.")
    shutil.copy(wf.deck_file.strip("'\""), "temp.dk")  # use temp.dk to update params instead of original deck
    setattr(wf, 'deck_file', "'temp.dk'")              # update deck file in control
    scattering_control.write_control()
    working_deck_file = wf.deck_file.strip("'\"")
    print(f"Working deck file: {working_deck_file}")
    param_file = wf.param_file.strip("'\"")
#-----------------------------------------------------------------------
    scattering_deck, scattering_deck_name = read_params_and_deck(param_file, working_deck_file)
#-----------------------------------------------------------------------
# copy optimized params from target deck to scattering deck (lines 4-19)
#-----------------------------------------------------------------------
    params_to_skip = ['Name', 'PI', 'J', 'MJ', 'T', 'MT', 'lwf', 'lsc', 'l3bc', 'lcut']
    for target_param, scattering_param in zip(opt_target_deck.parameters.keys()\
        ,scattering_deck.parameters.keys()):
        if target_param not in params_to_skip \
          and target_param == scattering_param \
          or target_param in ['cutR', 'cutA', 'cutW', 'delta', 'eps', 'theta', 'ups',]: # 'diff names for [] in opt_deck
            scattering_deck.parameters[scattering_param] = opt_target_deck.parameters[target_param]
#-----------------------------------------------------------------------
    scattering_deck.write_deck(scattering_deck_name) 
#-----------------------------------------------------------------------
# STEP 2: Find starting bscat value (bscat yielding E_rel ~ 3 MeV)
#-----------------------------------------------------------------------
    ss = str(system.parameters['spatial_symmetry'])
    scattering_deck.parameters[ss].wse = system.parameters['e_start'] - WSE_SHIFT
    scattering_deck.write_deck(scattering_deck_name) # Set wse based on closest E_rel
#-----------------------------------------------------------------------
    b_0 = minimize_scalar(find_starting_bscat, \
      args=(ss, target_energy, scattering_deck, \
      scattering_deck_name, system.parameters['e_start'],\
      scattering_ctrl, cmd, BIN_PATH), bounds=(-0.15, 0.15),\
      method='bounded', options={'maxiter': BSCAT_INIT_MAX}).x  # Find bscat that gives E_rel close to E_start
#-----------------------------------------------------------------------
    print(f"Initial bscat found: {b_0:.4f}")
    scattering_deck.parameters[ss].bscat = b_0
    scattering_deck.write_deck(scattering_deck_name) # Update bscat in deck
#-----------------------------------------------------------------------
# optimize deck with initial bscat
#-----------------------------------------------------------------------
    opt_E(b_0, ss, scattering_ctrl, target_energy, "scratch", cmd, BIN_PATH)
#-----------------------------------------------------------------------
    control = Control()
    control.read_control(scattering_ctrl)
    opt_deck_file = control.parameters['wf'].deck_file.strip("'\"")
    optimized_scattering_deck, optimized_scattering_deck_name = read_params_and_deck(param_file, opt_deck_file)
#-----------------------------------------------------------------------
# do bscat search again
#-----------------------------------------------------------------------
    optimized_scattering_deck.parameters[ss].wse = system.parameters['e_start'] - WSE_SHIFT
    optimized_scattering_deck.write_deck(optimized_scattering_deck_name) # Set wse based on closest E_rel
#-----------------------------------------------------------------------
    b_0 = minimize_scalar(find_starting_bscat, \
      args=(ss, target_energy, optimized_scattering_deck,\
      opt_deck_file, system.parameters['e_start'], scattering_ctrl, \
      cmd, BIN_PATH), bounds=(-0.15, 0.15) \
      ,method='bounded', options={'maxiter': BSCAT_INIT_MAX}).x  # Find bscat that gives E_rel close to E_start
#-----------------------------------------------------------------------
    print(f"Bscat after 2nd search: {b_0:.4f}")
    optimized_scattering_deck.parameters[ss].bscat = b_0
    optimized_scattering_deck.write_deck(opt_deck_file.strip(".dk")) # Update bscat in deck
    optimized_scattering_deck.parameters[ss].wse = \
      calculate_energy_com(\
      extract_final_energy(\
      run_energy(scattering_ctrl, cmd, BIN_PATH)), target_energy)-WSE_SHIFT
    optimized_scattering_deck.write_deck(opt_deck_file.strip(".dk")) # Set wse based on closest E_rel
    print(f"Setting wse to {optimized_scattering_deck.parameters[ss].wse:.4f} in {opt_deck_file}")
#-----------------------------------------------------------------------
# STEP 3: Run bscat-step-algorithm with initial bscat value
#-----------------------------------------------------------------------
    from bscat_optimizer import run_bscat_scan
    run_bscat_scan(
        control_file = scattering_ctrl,
        b_0 = b_0,
        ss = ss,
        input_db = system.parameters['input_db'] * b_0,     # x% of b_0
        input_de = system.parameters['input_de'], 
        E_min = system.parameters['e_min'], 
        E_max = system.parameters['e_max'], 
        E_start = system.parameters['e_start'], 
        target_energy = target_energy, 
        slope_max = MAX_BSCAT_SLOPE,
        cmd = cmd,
        BIN_PATH = BIN_PATH
    )
#-----------------------------------------------------------------------
    clean_dir(scattering_ctrl, target_ctrl)  # Clean up working directory after processing
#-----------------------------------------------------------------------
    print("Bscat scan completed.")
    print("optimized_decks.json file created at", Path("optimized_decks.json").resolve())
    print("Optimized decks saved to", Path("min").resolve())
    print("Logs saved to", Path("logs").resolve())
#-----------------------------------------------------------------------
"""
RunCommand
"""
def run_command(cmd, input_file, bscat=None):
    with open(input_file, "r") as ctrl:
        result = subprocess.run(cmd, stdin=ctrl, capture_output=True, text=True)
#-----------------------------------------------------------------------
    if result.returncode != 0:
        raise RuntimeError(f"Command failed: {' '.join(cmd)}\n{result.stderr}")
#-----------------------------------------------------------------------
# write to per-bscat log file
#-----------------------------------------------------------------------
    if bscat is not None:
        Path("logs").mkdir(exist_ok=True)
        log_file = Path("logs") / f"{bscat:.4f}.txt"
        with open(log_file, "a") as f:
            f.write(f"\n===== Command: {' '.join(cmd)} =====\n")
            f.write(f"Input file: {input_file}\n")
            f.write(result.stdout)
            f.write("\n")
#-----------------------------------------------------------------------
    return result.stdout
#-----------------------------------------------------------------------
"""
"""
def optimize_deck(control_file, cmd, BIN_PATH, bscat=None):
    OPTIMIZE_BIN = os.path.join(BIN_PATH, "optimize")
    cmd = cmd + [OPTIMIZE_BIN]
    output_text = run_command(cmd, control_file, bscat=bscat)
    energy, variance = extract_opt_energy(output_text)
    return energy, variance
#-----------------------------------------------------------------------
"""
"""
def run_energy(ctrl_file, cmd, BIN_PATH, bscat=None):
    ENERGY_BIN = os.path.join(BIN_PATH, "energy")
    cmd = cmd + [ENERGY_BIN]
    output_text = run_command(cmd, ctrl_file, bscat=bscat)
    return output_text
#-----------------------------------------------------------------------
"""
"""
def extract_opt_energy(output_text):
    match = re.search(r'OPTIMIZED ENERGY:\s*([-+]?[0-9]*\.?[0-9]+)\s*\(([-+]?[0-9]*\.?[0-9]+)\)', output_text)
    if match:
        energy = float(match.group(1))
        variance = float(match.group(2))
        return energy, variance
    else:
        print("Could not find optimized energy and variance.")
        return None, None
#-----------------------------------------------------------------------
"""
"""
def extract_final_energy(output):
    matches = re.findall(r'H\s*=\s*([-+]?\d+\.\d+)', output)
    return float(matches[-1]) if matches else None
"""
"""
def calculate_energy_com(energy_system, target_energy):
    return energy_system - target_energy
#-----------------------------------------------------------------------
"""
"""
def find_starting_bscat(bscat, ss, target_energy, scattering_deck \
     ,scattering_deck_file, e_start, scattering_ctrl_file, cmd, BIN_PATH):
    scattering_deck.parameters[ss].bscat = bscat                   # Update bscat in deck
    scattering_deck.write_deck(scattering_deck_file.strip(".dk"))  # Write updated deck with new bscat
    E_system = extract_final_energy(run_energy(scattering_ctrl_file, cmd, BIN_PATH))  # Run energy calculation
    E_rel = calculate_energy_com(E_system, target_energy)
    print(bscat, E_rel)
    return abs(E_rel - e_start)
#-----------------------------------------------------------------------
"""
"""
def opt_E(bscat, ss, control_file, target_energy, scratch_dir, cmd, BIN_PATH):
#-----------------------------------------------------------------------
    WSE_SHIFT=0.5
    WSE_OPT_VALUE=0.5
    WSE_NUM_OPTS=5
    WSE_NUM_EVALS=5
    OPT_SCALE=0.2
#-----------------------------------------------------------------------
    print(f"\n🔁 Bscat = {bscat:.4f}")
    control = Control()
    control.read_control(control_file)
    param_file = control.parameters['wf'].param_file.strip("'\"")
    deck_file = control.parameters['wf'].deck_file.strip("'\"")
    deck, deck_name = read_params_and_deck(param_file, deck_file)
    num_evals=control.parameters['num_opt_evaluations']
    num_walks=control.parameters['num_opt_walks']
    print(f"Using deck file: {deck_file}")
#-----------------------------------------------------------------------
# 1. Load and update deck
#-----------------------------------------------------------------------
    deck.parameters[ss].bscat = bscat
    deck.write_deck(deck_file.strip(".dk"))  # Write updated deck with new bscat
    print(f"Initial wse: {deck.parameters[ss].wse:.4f}")
    E_scattering = extract_final_energy(run_energy(control.Name + ".ctrl", cmd, BIN_PATH, bscat=bscat))
    E_rel_start = calculate_energy_com(E_scattering, target_energy)
    wse_val = E_rel_start-WSE_SHIFT
    print(f"E_rel = {E_rel_start:.4f}, setting wse = {wse_val:.4f}")
#-----------------------------------------------------------------------
# 2. Set wse in deck
#-----------------------------------------------------------------------
    deck.parameters[ss].wse = wse_val
    deck.write_deck(deck_file.strip(".dk"))  # Write updated deck with new wse
#-----------------------------------------------------------------------
# 3. Generate opt input for wse only
#-----------------------------------------------------------------------
    wse_corr, _ = zero_var_params(
        param_file,
        deck_file,
        ss,
        correlation_groups=[
            {'params': ['wse'], 'mode': 'set', 'value': WSE_OPT_VALUE}
        ]
    )
#-----------------------------------------------------------------------
    wse_opt_path = save_opt_file(bscat, wse_corr, "./opt", suffix="_wse")
    print(f"Generated wse opt file: {wse_opt_path}")
#-----------------------------------------------------------------------
# 4. Update control file with new optimization input
#-----------------------------------------------------------------------
    deck_basename = f"opt_he4n_av18_{bscat:.4f}.dk"
    optimized_deck_path = Path(scratch_dir) / deck_basename
    control.parameters['optimized_deck'] = f"'{optimized_deck_path}'"
    control.parameters['optimization_input'] = f"'{wse_opt_path}'"
    control.parameters['num_opt_walks']=WSE_NUM_OPTS
    control.parameters['num_opt_evaluations']=WSE_NUM_EVALS
    control.write_control()
#-----------------------------------------------------------------------
# 5. Estimate WSE
#-----------------------------------------------------------------------
    energy, variance = optimize_deck(control.Name + ".ctrl", cmd, BIN_PATH, bscat=bscat)
    E_rel_WSE = calculate_energy_com(energy, target_energy)
    WSE_EDIFF = E_rel_WSE-E_rel_start
    print(f"OPT WSE: E_rel = {E_rel_WSE:.4f}, var = {variance:.4f}, wse = {deck.parameters[ss].wse:.4f}")
    print(f"ENERGY DIFFERENCE: {WSE_EDIFF:.4f} MeV")
#-----------------------------------------------------------------------
# 6. Generate opt file for scattered nucleon correlations
#-----------------------------------------------------------------------
    spu_corrs, _ = zero_var_params(
        param_file,
        deck_file,
        ss,
        correlation_groups=[
            {'params': ['qssp1', 'qssp2'], 'mode': 'scale', 'value': OPT_SCALE},
            {'params': ['spu', 'spv', 'spr', 'spa', 'spb', 'spc', 'spk', 'spl'], 'mode': 'scale', 'value': OPT_SCALE},
            {'params': ['wsr', 'wsa'], 'mode': 'scale', 'value': OPT_SCALE}
        ]
    )
    opt_path_he5 = save_opt_file(bscat, spu_corrs, "./opt")
    print(f"Generated opt file for He5 correlations: {opt_path_he5}")
#-----------------------------------------------------------------------
# 7. Save opt input for He5 and update control
#-----------------------------------------------------------------------
    control.parameters['wf'].deck_file = f"'{optimized_deck_path}'" # use optimized deck
    control.parameters['optimization_input'] = f"'{opt_path_he5}'"
    control.parameters['num_opt_walks']=num_walks
    control.parameters['num_opt_evaluations']=num_evals
    control.write_control()
    print(f"Updated control file with optimized deck: {optimized_deck_path}")
#-----------------------------------------------------------------------
# 8. Run optimization with He5 correlations
#-----------------------------------------------------------------------
    energy, variance = optimize_deck(control.Name + ".ctrl", cmd, BIN_PATH, bscat=bscat)
    if energy is None:
        print("Optimization failed, no energy returned.")
        return None, None, None, None
    E_rel_CORR = calculate_energy_com(energy, target_energy)
    CORR_EDIFF_T = E_rel_CORR-E_rel_start
    CORR_EDIFF = E_rel_CORR-E_rel_WSE
    print(f"OPT CORR: E_rel = {E_rel_CORR:.4f}, var = {variance:.4f}")
    print(f"ENERGY DIFFERENCE: TOTAL: {CORR_EDIFF_T:.4f}, FROM WSE: {CORR_EDIFF:.4f} MeV")
#-----------------------------------------------------------------------
# 9. Update WSE to erel-wse_shift then optimize again
#-----------------------------------------------------------------------
    all_corrs, _ = zero_var_params(
        param_file,
        deck_file,
        ss,
        correlation_groups=[
            {'params': ['wse'], 'mode': 'set', 'value': WSE_OPT_VALUE},
            {'params': ['qssp1', 'qssp2'], 'mode': 'scale', 'value': OPT_SCALE},
            {'params': ['spu', 'spv', 'spr', 'spa', 'spb', 'spc', 'spk', 'spl'], 'mode': 'scale', 'value': OPT_SCALE},
            {'params': ['wsr', 'wsa'], 'mode': 'scale', 'value': OPT_SCALE}
        ]
    )
    opt_path_all = save_opt_file(bscat, all_corrs, "./opt", suffix="_all")
    print(f"Generated opt file for He5 correlations: {opt_path_all}")
#-----------------------------------------------------------------------
    control.parameters['optimization_input'] = f"'{opt_path_all}'"
    control.write_control()
    print(f"Updated control file with optimized input: {opt_path_all}")
#-----------------------------------------------------------------------
    energy, variance = optimize_deck(control.Name + ".ctrl", cmd, BIN_PATH, bscat=bscat)
    if energy is None:
        print("Optimization failed, no energy returned.")
        return None, None, None, None
#-----------------------------------------------------------------------
# 10. Final optimization; updates otimized deck
#-----------------------------------------------------------------------
    final_E_rel = calculate_energy_com(energy, target_energy)
    final_deck, _ = read_params_and_deck(param_file, optimized_deck_path)
    FINAL_EDIFF_T = final_E_rel-E_rel_start
    FINAL_EDIFF_C = final_E_rel-E_rel_CORR
    FINAL_EDIFF_W = final_E_rel-E_rel_WSE
    print(f"FINAL OPT: final E_rel = {final_E_rel:.4f}, var = {variance:.4f}, wse = {final_deck.parameters[ss].wse:.4f}")
    print(f"ENERGY DIFFERENCE: TOTAL: {FINAL_EDIFF_T:.4f}, FROM CORR: {FINAL_EDIFF_C:.4f}, FROM WSE: {FINAL_EDIFF_W:.4f} MeV")
#-----------------------------------------------------------------------
# 11. Save final deck 
#-----------------------------------------------------------------------
    folder_path = Path("min")
    dest_path = folder_path / deck_basename
    shutil.copy(optimized_deck_path, dest_path)
#-----------------------------------------------------------------------
    return final_E_rel, variance, str(dest_path), final_deck.parameters[ss].wse
#-----------------------------------------------------------------------
"""
"""
def clean_dir(scattering_ctrl_file, target_ctrl_file):
    os.remove("temp.dk")
    os.remove(scattering_ctrl_file)
    os.remove(target_ctrl_file)
#-----------------------------------------------------------------------
"""
Command Line Interface (CLI)
"""
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Run AutoOpt using system parameters from the command line.",
        epilog="""\
            Example:
            python AutoOpt.py --utility /path/to/he4n_.util \
                --system /path/to/he4n.sys 
            """
    )
#-----------------------------------------------------------------------
    parser.add_argument('--utility', required=True, help="Path to the utility file.")
    parser.add_argument('--system', required=True, help="Path to the nuclear system file.")
    args = parser.parse_args()
#-----------------------------------------------------------------------
# Create the System object from CLI parameters
#-----------------------------------------------------------------------
    utility_obj, system_obj = load_system(args.utility, args.system)
    print(f"System '{system_obj.parameters['name']}' successfully initialized.")
#-----------------------------------------------------------------------
# Execute the main auto-optimization routine
#-----------------------------------------------------------------------
    main(utility_obj, system_obj)
#-----------------------------------------------------------------------