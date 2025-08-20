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
                self.PARAM_FILE=list_wf_input[0][0]
                self.DECK_FILE =list_wf_input[1][0]
                self.RW_SPIN,self.SPIN_FILE=list_wf_input[2][:2]
                self.RW_YLM_PHI,self.YLM_FILE,self.PHI_FILE=list_wf_input[3][:3]
                self.RW_CONFIG,self.CONFIG_FILE=list_wf_input[4][:2]
                self.IDX=5
#-----------------------------------------------------------------------
            case "product":
#-----------------------------------------------------------------------
                self.PRODUCT_PARAM_FILE=list_wf_input[0][0]
                self.PRODUCT_DEC_FILE=list_wf_input[1][0]
                self.PARAM_FILE =list_wf_input[2][0]
                self.DECK_FILE  =list_wf_input[3][0]
                self.RW_SPIN,self.SPIN_FILE=list_wf_input[4][:2]
                self.RW_YLM_PHI,self.YLM_FILE,self.PHI_FILE=list_wf_input[5][:3]
                self.RW_CONFIG,self.CONFIG_FILE=list_wf_input[6][:2]
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
                file.write(" ".join([self.RW_SPIN,self.SPIN_FILE])+"\n")
                file.write(" ".join([self.RW_YLM_PHI,self.YLM_FILE,self.PHI_FILE])+"\n")
                file.write(" ".join([self.RW_CONFIG,self.CONFIG_FILE])+"\n")
#-----------------------------------------------------------------------
            case "product":
#-----------------------------------------------------------------------
                file.write(self.PRODUCT_PARAM_FILE+"\n")
                file.write(self.PRODUCT_DEC_FILE+"\n")
                file.write(self.PARAM_FILE+"\n")
                file.write(self.DECK_FILE+"\n")
                file.write(" ".join([self.RW_SPIN,self.SPIN_FILE])+"\n")
                file.write(" ".join([self.RW_YLM_PHI,self.YLM_FILE,self.PHI_FILE])+"\n")
                file.write(" ".join([self.RW_CONFIG,self.CONFIG_FILE])+"\n")
#-----------------------------------------------------------------------
            case _:
#-----------------------------------------------------------------------
                pass
#-----------------------------------------------------------------------
    def AddPrefix(self,wf_type,prefix):
#-----------------------------------------------------------------------
        self.PARAM_FILE =JoinPath(prefix,self.PARAM_FILE)
        self.DECK_FILE  =JoinPath(prefix,self.DECK_FILE)
        self.SPIN_FILE  =JoinPath(prefix,self.SPIN_FILE)
        self.YLM_FILE   =JoinPath(prefix,self.YLM_FILE)
        self.PHI_FILE   =JoinPath(prefix,self.PHI_FILE)
        self.CONFIG_FILE=JoinPath(prefix,self.CONFIG_FILE)
#-----------------------------------------------------------------------
        if wf_type.strip("\'") == "product":
            self.PRODUCT_PARAM_FILE=JoinPath(prefix,self.PRODUCT_PARAM_FILE)
            self.PRODUCT_DECK_FILE =JoinPath(prefix,self.PRODUCT_DECK_FILE )
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
        self.BASIS=data[0][0]
        self.BRA_EQ_KET=data[1][0]
        self.BRA_TYPE=data[2][0]
#----------------------------------------------------------------------
        self.INPUT_BRA=wavefunction_input_t(self.BRA_TYPE,data[3:10])
        data=data[3+self.INPUT_BRA.IDX:]
#----------------------------------------------------------------------
        if (self.BRA_EQ_KET != ".true."):
#----------------------------------------------------------------------
            self.KET_TYPE=data[0][0]
            self.INPUT_KET=wavefunction_input_t(self.KET_TYPE,data[1:8])
            data=data[1+self.INPUT_KET.IDX:]
        else:
            self.KET_TYPE=self.BRA_TYPE
            self.INPUT_KET=self.INPUT_BRA
#----------------------------------------------------------------------
        self.RW_WALK,self.WALK_FILE=data[0][:2]
        self.CONST_FILE=data[1][0]
        self.L2BP_FILE,self.L3BP_FILE=data[2][:2]
        self.LKE,self.LEMP=data[3][:2]
        self.RNG_SEED=data[4][0]
        self.NUM_BLOCKS,self.BLOCK_SIZE,self.NUM_WALKERS_PER_NODE=data[5][:3]
        self.BURN_IN_COUNT,self.NUM_MOVES_BETWEEN=data[6][:2]
        self.PARTICLE_MAX_DX=data[7][0]
        self.NPTS,self.NPTS_IN_ONE_FERMI,self.FD_FACTOR=data[8][:3]
        self.SS_LIMIT,self.SP_LIMIT,self.SD_LIMIT=data[9][:3]
        self.BOX_SIZE=data[10][0]
        self.SAMPLE_L2,self.RSAM,self.PSAM=data[11][:3]
        self.BIN_WIDTH,self.BIN_MAXR=data[12][:2]
        self.LASTP_SAMPLE_TYPE,self.LPS_A,self.LPS_B,self.EXTRA_NORM,self.EXTRA_NORM_ERROR=data[13][:5]
        self.ANC_ENERGY,self.REG_GAMMA,self.ECLSTR1,self.DCLSTR1,self.ECLSTR2,self.DCLSTR2=data[14][:6]
        self.DO_GROUP,self.GROUP_FILE,self.NORTAB_FILE=data[15][:3]
        self.NLOPT_METHOD,self.NUM_OPT_WALKS,self.NUM_OPT_EVALUATIONS=data[16][:3]
        self.OPTIMIZATION_INPUT_FILE=data[17][0]
        self.OPTIMIZED_DECK_FILE=data[18][0]
        self.SCRATCH_DIR=data[19][0]
#-----------------------------------------------------------------------
    def Write(self, out_file):
        file = open(out_file.strip("\'"), 'w')
        file.write(self.BASIS+"\n")
        file.write(self.BRA_EQ_KET+"\n")
        file.write(self.BRA_TYPE+"\n")
        self.INPUT_BRA.Write(self.BRA_TYPE,file)
#----------------------------------------------------------------------
        if (not self.BRA_EQ_KET):
#----------------------------------------------------------------------
            file.write(self.KET_TYPE+"\n")
            self.INPUT_KET.Write(self.KET_TYPE,file)
#----------------------------------------------------------------------
        file.write(" ".join([self.RW_WALK,self.WALK_FILE])+"\n")
        file.write(self.CONST_FILE+"\n")
        file.write(" ".join([self.L2BP_FILE,self.L3BP_FILE])+"\n")
        file.write(" ".join([self.LKE,self.LEMP])+"\n")
        file.write(self.RNG_SEED+"\n")
        file.write(" ".join([self.NUM_BLOCKS,self.BLOCK_SIZE,self.NUM_WALKERS_PER_NODE])+"\n")
        file.write(" ".join([self.BURN_IN_COUNT,self.NUM_MOVES_BETWEEN])+"\n")
        file.write(self.PARTICLE_MAX_DX+"\n")
        file.write(" ".join([self.NPTS,self.NPTS_IN_ONE_FERMI,self.FD_FACTOR])+"\n")
        file.write(" ".join([self.SS_LIMIT,self.SP_LIMIT,self.SD_LIMIT])+"\n")
        file.write(self.BOX_SIZE+"\n")
        file.write(" ".join([self.SAMPLE_L2,self.RSAM,self.PSAM])+"\n")
        file.write(" ".join([self.BIN_WIDTH,self.BIN_MAXR])+"\n")
        file.write(" ".join([self.LASTP_SAMPLE_TYPE,self.LPS_A,self.LPS_B,self.EXTRA_NORM,self.EXTRA_NORM_ERROR])+"\n")
        file.write(" ".join([self.ANC_ENERGY,self.REG_GAMMA,self.ECLSTR1,self.DCLSTR1,self.ECLSTR2,self.DCLSTR2])+"\n")
        file.write(" ".join([self.DO_GROUP,self.GROUP_FILE,self.NORTAB_FILE])+"\n")
        file.write(" ".join([self.NLOPT_METHOD,self.NUM_OPT_WALKS,self.NUM_OPT_EVALUATIONS])+"\n")
        file.write(self.OPTIMIZATION_INPUT_FILE+"\n")
        file.write(self.OPTIMIZED_DECK_FILE+"\n")
        file.write(self.SCRATCH_DIR)
        file.close()
#----------------------------------------------------------------------
    def AddPrefix(self,prefix):
        self.WALK_FILE=JoinPath(prefix,self.WALK_FILE)
        self.INPUT_BRA.AddPrefix(self.BRA_TYPE,prefix)
        if not self.BRA_EQ_KET: self.INPUT_KET.AddPrefix(self.KET_TYPE,prefix)
#----------------------------------------------------------------------
def JoinPath(p1,p2):
    return f"\'{p1.strip("\'")}{p2.strip("\'")}\'"
#----------------------------------------------------------------------