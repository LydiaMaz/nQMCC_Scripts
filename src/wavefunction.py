"""
wavefunction.py
nQMCC variational w.f. interface
"""
from re import findall
#-----------------------------------------------------------------------
from control import control_t
from parameters import parameters_t
from deck import deck_t
from utility import nQMCC
#-----------------------------------------------------------------------
class wavefunction_t:
#-----------------------------------------------------------------------
    def __init__(self,ctrl_file_name,bin_dir_,run_cmd_):
        self.CTRL=control_t(ctrl_file_name)
        self.PARAMS=parameters_t(self.CTRL.INPUT_BRA.PARAM_FILE)
        self.DK=deck_t(self.PARAMS,self.CTRL.INPUT_BRA.DECK_FILE)
        self.BIN_DIR=bin_dir_
        self.RUN_CMD=run_cmd_
#-----------------------------------------------------------------------
    def Evaluate(self,write_log,log_name):
        log = nQMCC("energy", self.CTRL, self.BIN_DIR, self.RUN_CMD, write_log, log_name)
        rx=r'H\s*=\s*(-?\d+\.\d+)\s*\((\d+\.\d+)\)'
        energy,var = findall(rx, log)[-1]
        return float(energy),float(var)
#-----------------------------------------------------------------------
    def Optimize(self,opt: deck_t,odk_file_name,write_log,log_name):
        opt.Write(self.PARAMS,opt.FILE_NAME)
        self.CTRL.OPTIMIZATION_INPUT_FILE=opt.FILE_NAME
        self.CTRL.OPTIMIZED_DECK_FILE=odk_file_name
        log = nQMCC("optimize", self.CTRL, self.BIN_DIR, self.RUN_CMD, write_log, log_name)
        self.DK=deck_t(self.PARAMS,self.CTRL.OPTIMIZED_DECK_FILE)
        self.CTRL.INPUT_BRA.DECK_FILE=self.CTRL.OPTIMIZED_DECK_FILE
        rx=r' OPTIMIZED ENERGY: (-?\d+\.\d+) \((\d+\.\d+)\)'
        opt_e,opt_v = findall(rx, log)[-1]
        return opt_e,opt_v