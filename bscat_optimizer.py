from pathlib import Path
import shutil
import sys
import os
import copy
# import matplotlib.pyplot as plt
import json

src_path = os.path.join(os.path.dirname(__file__), "src")
sys.path.append(src_path)

from control import Control 
from deck import read_params_and_deck
from AutoOpt import opt_E

# --- Helper Functions ---
def estimate_initial_db_de(b_0, ss, input_db, control_file, og_deck, og_deck_file, target_energy, cmd, BIN_PATH):
    print("input_db =", input_db)
    b_plus = b_0 + input_db
    b_minus = b_0 - input_db

    control = Control()
    control.read_control(control_file)

    print("\nBACKWARD (b_minus)")
    E_minus, *_ = opt_E(b_minus, ss, control_file, target_energy, "scratch", cmd, BIN_PATH) # <- deck file changes 

    # reset control deck 
    control.parameters['wf'].deck_file = f"'{og_deck_file}'"
    control.write_control()
    og_deck.write_deck(og_deck_file.strip(".dk"))
    
    print("\nFORWARD (b_plus)")
    E_plus, *_ = opt_E(b_plus, ss, control_file, target_energy, "scratch", cmd, BIN_PATH)

    # reset control deck 
    control.parameters['wf'].deck_file = f"'{og_deck_file}'"
    control.write_control()
    og_deck.write_deck(og_deck_file.strip(".dk"))

    db_de = (b_plus - b_minus) / (E_plus - E_minus)

    print(f"Initial db/de = {db_de:.5f} (b_plus = {b_plus:.5f}, E_plus = {E_plus:.5f}, b_minus = {b_minus:.5f}, E_minus = {E_minus:.5f})")

    return db_de

def update_db_de(b_new, b_prev, E_new, E_prev):
    return (b_new - b_prev) / (E_new - E_prev)

def log_step_info(step, direction, bscat, db_de, E_next):
    print(f"[{direction} Step {step}] bscat = {bscat:.5f}, db/de = {db_de:.5f}, E_rel = {E_next:.5f}")


def run_bscat_scan(control_file, b_0, ss, input_db, input_de, E_min, E_max, E_start, target_energy, slope_max, cmd, BIN_PATH):
    
    E_0, var_0, path_0, wse_0 = opt_E(b_0, ss, control_file, target_energy, "scratch", cmd, BIN_PATH)

    results = [{
        "bscat": b_0,
        "E_rel": E_0,
        "variance": var_0,
        "deck_path": path_0,
        "wse": wse_0
    }]

    # Create copy of the control file for backward walk

    control = Control()
    control.read_control(control_file)

    param_file = control.parameters['wf'].param_file.strip("'\"")
    og_deck_file = control.parameters['wf'].deck_file.strip("'\"")
    og_deck, _ = read_params_and_deck(param_file, og_deck_file)

    db_de_initial = estimate_initial_db_de(b_0, ss, input_db, control_file, og_deck, og_deck_file, target_energy, cmd, BIN_PATH)

    # --- Forward Walk ---
    bscat = b_0
    b_prev = b_0
    E_prev = E_0
    num_samples_up = int((E_max - E_start) / input_de)
    db_de = db_de_initial

    print("Starting forward walk")

    for step in range(num_samples_up):
        b_next = bscat + (db_de * input_de)
        E_next, var_next, path_next, wse_next = opt_E(b_next, ss, control_file, target_energy, "scratch", cmd, BIN_PATH)

        if E_next is None:
            break
        db_de = update_db_de(b_next, b_prev, E_next, E_prev)
        log_step_info(step, "Forward", b_next, db_de, E_next)

        results.append({
            "bscat": b_next,
            "E_rel": E_next,
            "variance": var_next,
            "deck_path": path_next,
            "wse": wse_next
        })

        if E_next < E_min or E_next > E_max or abs(db_de) > slope_max:
            break

        bscat = b_next
        b_prev = b_next
        E_prev = E_next

    # --- Backward Walk ---

    backward_control = copy.deepcopy(control)
    backward_control.parameters['wf'].deck_file = f"'{og_deck_file}'"
    og_deck.write_deck(og_deck_file.strip(".dk")) # Reset deck 
    backward_control.write_control() # overwrite control file
    
    bscat = b_0
    b_prev = b_0
    E_prev = E_0
    num_samples_down = int((E_start - E_min) / input_de)
    db_de = db_de_initial

    print("Starting backward walk")

    for step in range(num_samples_down):
        b_next = bscat - (db_de * input_de)
        E_next, var_next, path_next, wse_next = opt_E(b_next, ss, control_file, target_energy, "scratch", cmd, BIN_PATH)
        if E_next is None:
            break
        db_de = update_db_de(b_next, b_prev, E_next, E_prev)

        log_step_info(step, "Backward", b_next, db_de, E_next)

        results.append({
            "bscat": b_next,
            "E_rel": E_next,
            "variance": var_next,
            "deck_path": path_next,
            "wse": wse_next
        })

        if E_next < E_min or E_next > E_max or abs(db_de) > slope_max:
            break

        bscat = b_next
        b_prev = b_next
        E_prev = E_next

    # --- Save Results ---
    results_sorted = sorted(results, key=lambda r: r["E_rel"])
    with open("optimized_decks.json", "w") as f:
        json.dump(results_sorted, f, indent=2)

    # --- Plotting Results ---
    # plt.figure(figsize=(10, 6))
    # plt.plot([r['E_rel'] for r in results_sorted], [r['bscat'] for r in results_sorted], marker='o')
    # plt.xlabel("E_rel (MeV)")
    # plt.ylabel("bscat")
    # plt.title("bscat vs E_rel (Forward + Backward Walks)")
    # plt.grid(True)
    # plt.show()

    # --- Clean up working dir ---
    clean_up_dir()

    return results_sorted

def clean_up_dir():
    shutil.rmtree(Path("scratch"), ignore_errors=True)

# control_file = "/Users/lydiamazeeva/QMC/nQMCC/nQMCC/external/nQMCC_Scripts/he4n_copy.ctrl"
# run_bscat_scan(control_file, 0.0768, '1S[0]', 0.01536, 0.25, 0.5, 5.5, 3, -23.8819, 2, '/Users/lydiamazeeva/QMC/nQMCC/nQMCC/build/bin')
