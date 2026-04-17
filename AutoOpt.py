"""
AutoOpt.py
Automatic Optimazation Interface for nQMCC
Quantum Monte Carlo Group @ Washington University in St. Louis
08/01/2025
"""
import sys
import os
import json
import argparse
import builtins
from datetime import datetime
#-----------------------------------------------------------------------
sys.path.append(os.path.join(os.path.dirname(__file__), "src"))
#-----------------------------------------------------------------------
from utility import utility_t
from wavefunction import wavefunction_t,InitPShellScattWF
from bscat import SingleChannelScan
from bound import BoundStateOptimizePhases
#-----------------------------------------------------------------------
def BoundStates(util: utility_t):
#-----------------------------------------------------------------------
    BREAK="="*72
#-----------------------------------------------------------------------
    print("SETTING UP WORKING ENVIRONMENT")
    try:
        os.mkdir(util.WORKING_DIR)
    except FileExistsError:
        print("***WORKING DIRECTORY EXISTS***")
        util.WORKING_DIR=f"{util.WORKING_DIR}{util.NAME}-{datetime.now().strftime('%Y-%m-%d_%H-%M')}/"
        print(f"CURRENT WORKING DIRECTORY: {util.WORKING_DIR}")
        os.mkdir(util.WORKING_DIR)
    os.makedirs(f"{util.NQMCC_DIR}walks", exist_ok=True)
    print("... DONE")
    print(BREAK)
#-----------------------------------------------------------------------
    # Capture all print output to be written to a log file
    opt_output_lines=[]
    class _CapturingPrint:
        def __init__(self,real_print):
            self._real=real_print
        def __call__(self,*args,**kwargs):
            opt_output_lines.append(" ".join(str(a) for a in args))
            self._real(*args,**kwargs)
    original_print=builtins.print
    builtins.print=_CapturingPrint(original_print)
#-----------------------------------------------------------------------
    try:
        results=[]
        for i,(const,pot2b,pot3b) in enumerate(zip(util.CONSTANTS_FILES,util.TWO_BODY_FILES,util.THREE_BODY_FILES)):
            pair_name=f"{pot2b.split('.')[0]}.{pot3b.split('.')[0]}"
            pot_dir=f"{util.WORKING_DIR}{pair_name}/"
            os.makedirs(f"{pot_dir}ctrl", exist_ok=True)
            os.makedirs(f"{pot_dir}logs", exist_ok=True)
            os.makedirs(f"{pot_dir}dk",   exist_ok=True)
            os.makedirs(f"{pot_dir}opt",  exist_ok=True)
            print(f"\n{BREAK}")
            print(f"PROCESSING {i+1}/{len(util.TWO_BODY_FILES)}: {pair_name}")
            print(f"CONSTANTS: {const}")
            print(f"V2B: {pot2b}")
            print(f"V3B: {pot3b}")
            print(BREAK)
#-----------------------------------------------------------------------
            try:
                res=BoundStateOptimizePhases(util,pair_name,pot_dir,i)
                results.append(res)
            except Exception as ex:
                print(f"ERROR optimizing {pair_name}: {ex}")
                results.append({"POTENTIAL_PAIR":pair_name,"ERROR":str(ex)})
#-----------------------------------------------------------------------
    finally:
        builtins.print=original_print
#-----------------------------------------------------------------------
    # Write combined results across all potential pairs
    processed_results=[{k:v for k,v in r.items() if k!="ATTEMPTS"} for r in results]
    combined_data={
        "run_metadata":{
            "timestamp":datetime.now().isoformat(),
            "optimization_mode":"phased",
            "optimization_settings":{
                "opt_scale":util.OPT_SCALE,
                "num_opt_evaluations":util.NUM_OPT_EVALUATIONS,
                "esep_scale":f"[{util.ESEP_START}, {util.ESEP_STOP}] n={util.ESEP_NUM}",
                "eta_flat":util.ETA_FLAT,
                "three_body_flat":util.THREE_BODY_FLAT,
                "ss_flat":util.SS_FLAT,
            },
        },
        "results":processed_results,
    }
    out=f"{util.WORKING_DIR}combined_results.json"
    with open(out,"w") as f:
        json.dump(combined_data,f,indent=2)
    print("\nSaved:",out)
#-----------------------------------------------------------------------
    # Write full log with settings header
    out_file=f"{util.WORKING_DIR}optimization_log.out"
    with open(out_file,"w") as f:
        f.write("# Optimization Log (Phased)\n")
        f.write(f"# Generated:           {datetime.now().isoformat()}\n")
        f.write(f"# OPT_SCALE:           {util.OPT_SCALE}\n")
        f.write(f"# NUM_OPT_EVALUATIONS: {util.NUM_OPT_EVALUATIONS}\n")
        f.write(f"# ESEP_SCALE:          [{util.ESEP_START}, {util.ESEP_STOP}] n={util.ESEP_NUM}\n")
        f.write(f"# ETA_FLAT:            {util.ETA_FLAT}\n")
        f.write(f"# THREE_BODY_FLAT:     {util.THREE_BODY_FLAT}\n")
        f.write(f"# SS_FLAT:             {util.SS_FLAT}\n")
        f.write("#" + "="*70 + "\n\n")
        for line in opt_output_lines:
            f.write(line+"\n")
    print("Saved optimization log:",out_file)


def SingleChannelScattering(util: utility_t):
#-----------------------------------------------------------------------
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
    target = wavefunction_t(util.CTRL_FILE,util.NQMCC_DIR,util.BIN_DIR,util.RUN_CMD)
    target_label=target.DK.NAME.strip("\'")
    target.CTRL.FILE_NAME=f"{util.WORKING_DIR}target.ctrl"
    target.CTRL.NUM_BLOCKS=util.NUM_BLOCKS
    target.CTRL.BLOCK_SIZE=util.BLOCK_SIZE
    target.CTRL.WALKERS_PER_NODE=util.WALKERS_PER_NODE
    target.CTRL.CALC_TYPE=util.CALC_TYPE
    target.CTRL.CALC_FILE=f"'{util.WORKING_DIR}ctrl/target.calc'"
    target.SetCalc(util.CALC_TYPE,f"{util.NQMCC_DIR}ctrl/calc/{util.CALC_FILE}")
    print("... DONE")
    print(BREAK)
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
        print(BREAK)
#-----------------------------------------------------------------------
        if util.OPTIMIZE_TARGET:
            print(f"OPTIMIZING TARGET: {tname}")
            ecore,vcore=[0,1]
        else:
            print(f"EVALUATING TARGET: {tname}")
            ecore,vcore=target.Evaluate(True,tname)
            print(f"E = {ecore:.4f} +- {vcore:.4f}")
#-----------------------------------------------------------------------
        print(BREAK)
#-----------------------------------------------------------------------
        for scat_ctrl,ssi in zip(util.SCATTERING_CTRL_FILES,util.SS_INDEXS):
#-----------------------------------------------------------------------
            print("SETTING UP SCATTERING WAVEFUNCTION")
            scatter = wavefunction_t(scat_ctrl,util.NQMCC_DIR,util.BIN_DIR,util.RUN_CMD)
            scatter_label=scatter.DK.NAME.strip("\'")
            scatter.CTRL.FILE_NAME=f"{util.WORKING_DIR}scatter.ctrl"
            scatter.CTRL.NUM_BLOCKS=util.NUM_BLOCKS
            scatter.CTRL.BLOCK_SIZE=util.BLOCK_SIZE
            scatter.CTRL.WALKERS_PER_NODE=util.WALKERS_PER_NODE
            scatter.CTRL.CALC_TYPE=util.CALC_TYPE
            scatter.CTRL.CALC_FILE=f"'{util.WORKING_DIR}ctrl/scatter.calc'"
            scatter.SetCalc(util.CALC_TYPE,f"{util.NQMCC_DIR}ctrl/calc/{util.CALC_FILE}")
            scatter.CTRL.CONST_FILE=f"'{util.NQMCC_DIR}constants/{const}'"
            scatter.CTRL.L2BP_FILE=f"'{util.NQMCC_DIR}pots/{pot2b}'"
            scatter.CTRL.L3BP_FILE=f"'{util.NQMCC_DIR}pots/{pot3b}'"
            sname=f"{scatter_label}.{pot2b_label}.{pot3b_label}"
            print(f"Scattering wave function: {sname}")
            print("... DONE")
            print(BREAK)
#-----------------------------------------------------------------------
            #Copy Core deck to scatter deck, write deck to file
            InitPShellScattWF(scatter,target,f"\'{util.WORKING_DIR}temp.dk\'")
            print(f"SCANNING BSCAT")
            print(BREAK)
#-----------------------------------------------------------------------
            SingleChannelScan(util,sname,scatter,ssi,ecore,vcore)
#-----------------------------------------------------------------------
def AutoOptAPI(util: utility_t):
#-----------------------------------------------------------------------
    match util.SYSTEM_TYPE.lower():
        case "bound":
            BoundStates(util)
        case "sc_scattering":
            SingleChannelScattering(util)
#-----------------------------------------------------------------------
if __name__ == '__main__':
#-----------------------------------------------------------------------
    BREAK="="*72
#-----------------------------------------------------------------------
    parser = argparse.ArgumentParser(
        description="Run AutoOpt using a utility file",
        epilog="""
            Example:
            python3 AutoOpt.py --utility /test/test.util
            """
    )
#-----------------------------------------------------------------------
    parser.add_argument('--utility', required=True, help="/path/to/utility/file")
    args = parser.parse_args()
#-----------------------------------------------------------------------
    print(BREAK)
    util = utility_t(args.utility)
    print(BREAK)
#-----------------------------------------------------------------------
    AutoOptAPI(util)
#-----------------------------------------------------------------------