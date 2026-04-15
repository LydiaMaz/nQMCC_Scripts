"""
control.py
Control Files for nQMCC
"""
#-----------------------------------------------------------------------
class wavefunction_input_t:
#-----------------------------------------------------------------------
    def __init__(self,wf_type,list_wf_input,ridx=0):
        match wf_type.strip("\'"):
#-----------------------------------------------------------------------
            case "variational":
#-----------------------------------------------------------------------
                self.PARAM_FILE  = list_wf_input[0][0]
                self.DECK_FILE   = list_wf_input[1][0]
                self.RW_SPIN,  self.SPIN_FILE, self.ISOSPIN_FILE = list_wf_input[2][:3]
                self.RW_YLM_PHI,self.YLM_FILE, self.PHI_FILE     = list_wf_input[3][:3]
                self.RW_CONFIG, self.CONFIG_FILE                  = list_wf_input[4][:2]
                self.IDX=5
#-----------------------------------------------------------------------
            case "product":
#-----------------------------------------------------------------------
                self.PRODUCT_PARAM_FILE = list_wf_input[0][0]
                self.PRODUCT_DECK_FILE  = list_wf_input[1][0]
                self.PARAM_FILE         = list_wf_input[2][0]
                self.DECK_FILE          = list_wf_input[3][0]
                self.RW_SPIN,  self.SPIN_FILE, self.ISOSPIN_FILE = list_wf_input[4][:3]
                self.RW_YLM_PHI,self.YLM_FILE, self.PHI_FILE     = list_wf_input[5][:3]
                self.RW_CONFIG, self.CONFIG_FILE                  = list_wf_input[6][:2]
                self.IDX=7
#-----------------------------------------------------------------------
            case _:
#-----------------------------------------------------------------------
                self.IDX=0
#-----------------------------------------------------------------------
    def Write(self,wf_type,file):
        match wf_type.strip("\'"):
#-----------------------------------------------------------------------
            case "variational":
#-----------------------------------------------------------------------
                file.write(self.PARAM_FILE+"\n")
                file.write(self.DECK_FILE+"\n")
                file.write(" ".join([self.RW_SPIN,self.SPIN_FILE,self.ISOSPIN_FILE])+"\n")
                file.write(" ".join([self.RW_YLM_PHI,self.YLM_FILE,self.PHI_FILE])+"\n")
                file.write(" ".join([self.RW_CONFIG,self.CONFIG_FILE])+"\n")
#-----------------------------------------------------------------------
            case "product":
#-----------------------------------------------------------------------
                file.write(self.PRODUCT_PARAM_FILE+"\n")
                file.write(self.PRODUCT_DECK_FILE+"\n")
                file.write(self.PARAM_FILE+"\n")
                file.write(self.DECK_FILE+"\n")
                file.write(" ".join([self.RW_SPIN,self.SPIN_FILE,self.ISOSPIN_FILE])+"\n")
                file.write(" ".join([self.RW_YLM_PHI,self.YLM_FILE,self.PHI_FILE])+"\n")
                file.write(" ".join([self.RW_CONFIG,self.CONFIG_FILE])+"\n")
#-----------------------------------------------------------------------
            case _:
#-----------------------------------------------------------------------
                pass
#-----------------------------------------------------------------------
    def AddPrefix(self,wf_type,prefix):
#-----------------------------------------------------------------------
        self.PARAM_FILE   =JoinPath(prefix,self.PARAM_FILE)
        self.DECK_FILE    =JoinPath(prefix,self.DECK_FILE)
        self.SPIN_FILE    =JoinPath(prefix,self.SPIN_FILE)
        self.ISOSPIN_FILE =JoinPath(prefix,self.ISOSPIN_FILE)
        self.YLM_FILE     =JoinPath(prefix,self.YLM_FILE)
        self.PHI_FILE     =JoinPath(prefix,self.PHI_FILE)
        self.CONFIG_FILE  =JoinPath(prefix,self.CONFIG_FILE)
#-----------------------------------------------------------------------
        if wf_type.strip("\'") == "product":
            self.PRODUCT_PARAM_FILE=JoinPath(prefix,self.PRODUCT_PARAM_FILE)
            self.PRODUCT_DECK_FILE =JoinPath(prefix,self.PRODUCT_DECK_FILE)
#-----------------------------------------------------------------------
class control_t:
#-----------------------------------------------------------------------
    def __init__(self,file_name_,prefix=""):
#-----------------------------------------------------------------------
        self.FILE_NAME = file_name_
        self.Read()
        self.AddPrefix(prefix)
#-----------------------------------------------------------------------
    def Read(self):
#-----------------------------------------------------------------------
        file = open(self.FILE_NAME.strip("\'"), 'r')
        data = [(l.strip().split()) for l in file.readlines()]
        file.close()
#----------------------------------------------------------------------
        self.BASIS                      = data[0][0]
        self.CONST_FILE                 = data[1][0]
        self.L2BP_FILE, self.L3BP_FILE  = data[2][:2]
        self.LKE, self.LEMP             = data[3][:2]
        self.CALC_TYPE, self.CALC_FILE  = data[4][:2]
        self.BRA_EQ_KET                 = data[5][0]
        self.BRA_TYPE                   = data[6][0]
#----------------------------------------------------------------------
        self.INPUT_BRA = wavefunction_input_t(self.BRA_TYPE,data[7:])
        data_index = 7 + self.INPUT_BRA.IDX
#----------------------------------------------------------------------
        if self.BRA_EQ_KET != ".true.":
#----------------------------------------------------------------------
            self.KET_TYPE  = data[data_index][0]
            self.INPUT_KET = wavefunction_input_t(self.KET_TYPE,data[data_index+1:])
            data_index     = data_index + 1 + self.INPUT_KET.IDX
        else:
            self.KET_TYPE  = self.BRA_TYPE
            self.INPUT_KET = self.INPUT_BRA
#----------------------------------------------------------------------
        self.RW_WALK,  self.WALK_FILE                           = data[data_index][:2]
        self.RNG_SEED                                           = data[data_index+1][0]
        self.NUM_BLOCKS,self.BLOCK_SIZE,self.NUM_WALKERS_PER_NODE = data[data_index+2][:3]
        self.BURN_IN_COUNT, self.NUM_MOVES_BETWEEN              = data[data_index+3][:2]
        self.PARTICLE_MAX_DX                                    = data[data_index+4][0]
        self.SAMPLE_L2, self.RSAM, self.PSAM                    = data[data_index+5][:3]
        self.NPTS, self.NPTS_IN_ONE_FERMI, self.FD_FACTOR       = data[data_index+6][:3]
        self.SS_LIMIT, self.SP_LIMIT, self.SD_LIMIT             = data[data_index+7][:3]
        self.BOX_SIZE                                           = data[data_index+8][0]
        self.BIN_WIDTH, self.BIN_MAXR                           = data[data_index+9][:2]
        self.SCRATCH_DIR                                        = data[data_index+10][0]
#-----------------------------------------------------------------------
    def Write(self, out_file):
        file = open(out_file.strip("\'"), 'w')
        file.write(self.BASIS+"\n")
        file.write(self.CONST_FILE+"\n")
        file.write(" ".join([self.L2BP_FILE,self.L3BP_FILE])+"\n")
        file.write(" ".join([self.LKE,self.LEMP])+"\n")
        file.write(" ".join([self.CALC_TYPE,self.CALC_FILE])+"\n")
        file.write(self.BRA_EQ_KET+"\n")
        file.write(self.BRA_TYPE+"\n")
        self.INPUT_BRA.Write(self.BRA_TYPE,file)
#----------------------------------------------------------------------
        if self.BRA_EQ_KET != ".true.":
#----------------------------------------------------------------------
            file.write(self.KET_TYPE+"\n")
            self.INPUT_KET.Write(self.KET_TYPE,file)
#----------------------------------------------------------------------
        file.write(" ".join([self.RW_WALK,self.WALK_FILE])+"\n")
        file.write(self.RNG_SEED+"\n")
        file.write(" ".join([self.NUM_BLOCKS,self.BLOCK_SIZE,self.NUM_WALKERS_PER_NODE])+"\n")
        file.write(" ".join([self.BURN_IN_COUNT,self.NUM_MOVES_BETWEEN])+"\n")
        file.write(self.PARTICLE_MAX_DX+"\n")
        file.write(" ".join([self.SAMPLE_L2,self.RSAM,self.PSAM])+"\n")
        file.write(" ".join([self.NPTS,self.NPTS_IN_ONE_FERMI,self.FD_FACTOR])+"\n")
        file.write(" ".join([self.SS_LIMIT,self.SP_LIMIT,self.SD_LIMIT])+"\n")
        file.write(self.BOX_SIZE+"\n")
        file.write(" ".join([self.BIN_WIDTH,self.BIN_MAXR])+"\n")
        file.write(self.SCRATCH_DIR)
        file.close()
#----------------------------------------------------------------------
    def AddPrefix(self,prefix):
        self.WALK_FILE=JoinPath(prefix,self.WALK_FILE)
        self.INPUT_BRA.AddPrefix(self.BRA_TYPE,prefix)
        if self.BRA_EQ_KET != ".true.":
            self.INPUT_KET.AddPrefix(self.KET_TYPE,prefix)
#-----------------------------------------------------------------------
# CALC_FILE_TYPE integer codes (matches Fortran SELECT CASE)
# 0=none  1=optimize  2=gfmc  3=phi  4=spec  5=anc  6=momd  7=emff
#-----------------------------------------------------------------------
class calc_t:
#-----------------------------------------------------------------------
    def __init__(self,calc_type,file_name_):
        self.CALC_TYPE = int(calc_type)
        self.FILE_NAME = file_name_
        self.Read()
#-----------------------------------------------------------------------
    def Read(self):
        if self.CALC_TYPE == 0: return
        file = open(self.FILE_NAME.strip("\'"), 'r')
        data = [(l.strip().split()) for l in file.readlines()]
        file.close()
        match self.CALC_TYPE:
#-----------------------------------------------------------------------
            case 1: # optimize
#-----------------------------------------------------------------------
                self.NLOPT_METHOD,self.NUM_OPT_WALKS,self.NUM_OPT_EVALUATIONS        = data[0][:3]
                self.LIMIT_CHARGE_RADII,self.NEUTRON_RADIUS_LIMIT,self.PROTON_RADIUS_LIMIT = data[1][:3]
                self.OPTIMIZATION_INPUT_FILE                                           = data[2][0]
                self.OPTIMIZED_DECK_FILE                                               = data[3][0]
#-----------------------------------------------------------------------
            case _:
#-----------------------------------------------------------------------
                pass
#-----------------------------------------------------------------------
    def Write(self,out_file):
        if self.CALC_TYPE == 0: return
        file = open(out_file.strip("\'"), 'w')
        match self.CALC_TYPE:
#-----------------------------------------------------------------------
            case 1: # optimize
#-----------------------------------------------------------------------
                file.write(" ".join([self.NLOPT_METHOD,self.NUM_OPT_WALKS,self.NUM_OPT_EVALUATIONS])+"\n")
                file.write(" ".join([self.LIMIT_CHARGE_RADII,self.NEUTRON_RADIUS_LIMIT,self.PROTON_RADIUS_LIMIT])+"\n")
                file.write(self.OPTIMIZATION_INPUT_FILE+"\n")
                file.write(self.OPTIMIZED_DECK_FILE)
#-----------------------------------------------------------------------
            case _:
#-----------------------------------------------------------------------
                pass
        file.close()
#----------------------------------------------------------------------
def JoinPath(p1,p2):
    return f"\'{p1.strip("\'")}{p2.strip("\'")}\'"
#----------------------------------------------------------------------
