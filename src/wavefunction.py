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
    def __init__(self,ctrl_file_name,nqmcc_dir_,bin_dir_,run_cmd_):
        self.CTRL=control_t(ctrl_file_name,nqmcc_dir_)
        self.PARAMS=parameters_t(self.CTRL.INPUT_BRA.PARAM_FILE)
        self.DK=deck_t(self.PARAMS,self.CTRL.INPUT_BRA.DECK_FILE)
        self.NQMCC_DIR=nqmcc_dir_
        self.BIN_DIR=bin_dir_
        self.RUN_CMD=run_cmd_
#-----------------------------------------------------------------------
    def Evaluate(self,write_log,log_name):
        log = nQMCC("energy", self.CTRL, self.BIN_DIR, self.RUN_CMD, write_log, log_name)
        rx=r'H\s*=\s*(-?\d+\.\d+)\s*\((\d+\.\d+)\)'
        try:
            energy_str,var_str = findall(rx, log)[-1]
            return float(energy_str),float(var_str)
        except IndexError:
            print("***Warning Evaluate did not find a valid energy***")
            return None, None
#-----------------------------------------------------------------------
    def Optimize(self,opt: deck_t,odk_file_name,write_log,log_name):
        opt.Write(self.PARAMS,opt.FILE_NAME)
        self.CTRL.OPTIMIZATION_INPUT_FILE=opt.FILE_NAME
        self.CTRL.OPTIMIZED_DECK_FILE=odk_file_name
        log = nQMCC("optimize", self.CTRL, self.BIN_DIR, self.RUN_CMD, write_log, log_name)
        self.DK=deck_t(self.PARAMS,self.CTRL.OPTIMIZED_DECK_FILE)
        self.CTRL.INPUT_BRA.DECK_FILE=self.CTRL.OPTIMIZED_DECK_FILE
        rx=r' OPTIMIZED ENERGY: (-?\d+\.\d+) \((\d+\.\d+)\)'
        try:
            opt_e_str,opt_v_str = findall(rx, log)[-1]
            return float(opt_e_str),float(opt_v_str)
        except IndexError:
            print("***Warning Optimize did not find a valid energy***")
            return None, None
#-----------------------------------------------------------------------
def InitPShellScattWF(scatter: wavefunction_t,target: wavefunction_t,out_file):
    for key,val in target.DK.__dict__.items():
        if (key not in ["FILE_NAME","NAME","SS"]):
            if isinstance(val,list):
                scatter.DK.__dict__[key]=val
            else:
                if "." in val: scatter.DK.__dict__[key]=val
    scatter.DK.FILE_NAME=out_file
    scatter.DK.Write(scatter.PARAMS,out_file)
    scatter.CTRL.INPUT_BRA.DECK_FILE=out_file
#-----------------------------------------------------------------------
def InitNShellBoundWF(nucleus: wavefunction_t, alpha: wavefunction_t,
                      out_file: str, copy_esep=True):
    """
    Initialize a bound-state wavefunction (e.g. Li-6, Li-7) by
    copying the optimized alpha-core correlation parameters
    from 'alpha' into 'nucleus'.

    Only the correlation block corresponding to deck lines 5–18
    (ZETA through QSSS2) is replaced, plus ESEP if desired.

    This is the bound-state analogue of InitPShellScattWF.
    """

    dst = nucleus.DK       # deck for Li-6, Li-7, etc.
    src = alpha.DK         # optimized He-4 deck (for same potential)

    # Optionally copy ESEP (potential-dependent separation energy)
    if copy_esep:
        dst.ESEP = src.ESEP[:]

    # Correlation block: deck lines 5–18
    corr_keys = [
        "ETA",
        "ZETA",
        "FSCAL",
        "AC",
        "AA",
        "AR",
        "ALPHA",
        "BETA",
        "GAMMA",
        "UUR", "UUA", "UUW",
        "SSH_LFSP", "SSH_E_OR_V", "SSH_LSCAT",
        "SSH_LMU1", "SSH_LMU2", "SSH_LQNUM", "SSH_LNODES",
        "SSH_WSE", "SSH_WSV", "SSH_WSR", "SSH_WSA",
        "SSH_WBRHO", "SSH_WBALPH",
        "DELTA", "EPSILON", "THETA", "UPSILON", "RSCAL", "USCAL",
        "QPS1", "QPS2",
        "QSSS1", "QSSS2",
    ]

    # Copy these from alpha → nucleus
    for key in corr_keys:
        dst.__dict__[key] = src.__dict__[key][:] if isinstance(src.__dict__[key], list) \
                            else src.__dict__[key]

    # Do NOT touch nucleus-specific parts:
    # NAME, PARITY, J/T, LCUT/LOPC, orbitals, SS block, etc.

    # Write new deck file
    dst.FILE_NAME = out_file
    dst.Write(nucleus.PARAMS, out_file)

    # Update control file's active deck
    nucleus.CTRL.INPUT_BRA.DECK_FILE = out_file

