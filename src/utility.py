
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
#-----------------------------------------------------------------------
    def __init__(self, file_name_):
        self.FILE_NAME = file_name_
        self.Read(self.FILE_NAME)
#-----------------------------------------------------------------------
        print("UTILITY INPUTS")
        print("="*72)
        for key,val in self.__dict__.items():
            print(f"{key} :: {val}")
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
        self.OPT_SCALE,self.NUM_OPT_EVALUATIONS=[float(data[12][0]),data[12][1]]
#-----------------------------------------------------------------------
        if self.SYSTEM_TYPE.lower() == "sc_scattering":
            data=data[13:]
            self.NUM_CHANNELS=int(data[0][0])
            self.SCATTERING_CTRL_FILES=data[1][:self.NUM_CHANNELS]
            self.SS_INDEXS=[int(d) for d in data[2][:self.NUM_CHANNELS]]
            self.OPTIMIZE_TARGET = 1 == int(data[3][0])
            self.ENERGY_LOWER_BOUND, self.ENERGY_UPPER_BOUND, self.DELTA_ENERGY=[float(d) for d in data[4][:3]]
            self.INITIAL_BSCAT,self.INITIAL_DELTA_BSCAT,self.MAX_BSCAT_SLOPE=[float(d) for d in data[5][:3]]
            self.MAX_SCAN_COUNT=int(data[6][0])
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
                f.write(f"===== Command: {' '.join(cmd)} =====\n")
                f.write(result.stdout)
        return result.stdout