
"""
wavefunction.py
AutoOpt paths + command line interface
"""
#-----------------------------------------------------------------------
from pathlib import Path
from subprocess import run
import os
import re
#-----------------------------------------------------------------------
from control import control_t
#-----------------------------------------------------------------------
class utility_t:
    FILE_NAME=""
    SYSTEM_TYPE=""
    RUN_CMD=[]
    NQMCC_DIR=""
    BIN_DIR=""
    WORKING_DIR=""
    CTRL_FILE=""
    NUM_POTS=1
    CONSTANTS_FILES=[]
    TWO_BODY_FILES=[]
    THREE_BODY_FILES=[]
    NUM_BLOCKS=""
    BLOCK_SIZE=""
    WALKERS_PER_NODE=""
    NUM_OPT_SAMPLES=""
#-----------------------------------------------------------------------
#   single channel scattering
#-----------------------------------------------------------------------
    SCATTERING_CTRL_FILE=""
    OPTIMIZE_TARGET=False
    SS_IDX=0
    ENERGY_LOWER_BOUND=0.5
    ENERGY_UPPER_BOUND=6.0
    DELTA_ENERGY=0.25
    INITIAL_DB=0.25
#-----------------------------------------------------------------------
    def __init__(self, file_name_):
        self.FILE_NAME = file_name_
        self.Read(self.FILE_NAME)
#-----------------------------------------------------------------------
    def Read(self, filename):
        file = open(self.FILE_NAME.strip("\'"), 'r')
        data = [(l.strip().split()) for l in file.readlines()]
        file.close()
        self.NAME=data[0][0]
        self.SYSTEM_TYPE=data[1][0]
        self.RUN_CMD=data[2]
        self.NQMCC_DIR=data[3][0]
        self.BIN_DIR=data[4][0]
        self.WORKING_DIR=data[5][0]
        self.CTRL_FILE=data[6][0]
        self.NUM_POTS=int(data[7][0])
        self.CONSTANTS_FILES=data[8][:self.NUM_POTS]
        self.TWO_BODY_FILES=data[9][:self.NUM_POTS]
        self.THREE_BODY_FILES=data[10][:self.NUM_POTS]
        self.NUM_BLOCKS,self.BLOCK_SIZE,self.WALKERS_PER_NODE=data[11][:3]
        self.NUM_OPT_SAMPLES=data[12][0]
#-----------------------------------------------------------------------
        if self.SYSTEM_TYPE.lower() == "sc_scattering":
            data=data[13:]
            self.SCATTERING_CTRL_FILE=data[0][0]
            self.OPTIMIZE_TARGET = 1 == int(data[1][0])
            self.SS_IDX=int(data[2][0])
            self.ENERGY_LOWER_BOUND, self.ENERGY_UPPER_BOUND, self.DELTA_ENERGY=[float(d) for d in data[3][:3]]
            self.INITIAL_DELTA_BSCAT=float(data[4][0])

#-----------------------------------------------------------------------
def nQMCC(binary: str, ctrl: control_t, bin_dir: str, runner: list, write_log=False, log_name=""):
        cmd = f"{" ".join(runner)} {bin_dir}{binary}".split()
        ctrl.Write(ctrl.FILE_NAME)
        with open(ctrl.FILE_NAME, "r") as ctrl:
            result = run(cmd, stdin=ctrl, capture_output=True, text=True)
        if write_log:
            Path("logs").mkdir(exist_ok=True)
            log_file = Path("logs") / f"{log_name}.{binary}"
            with open(log_file, "a") as f:
                f.write(f"\n===== Command: {' '.join(cmd)} =====\n")
                f.write(result.stdout)
        return result.stdout