"""
control.py
Control Files for nQMCC
"""
#-----------------------------------------------------------------------
class wavefunction_input_t:
#-----------------------------------------------------------------------
    IDX=0 #NUMBER OF LINES TO READ
    PARAM_FILE ="" # PARAM FILES
    DECK_FILE  ="" # VARIATIONAL PARAMS
    RW_SPIN    ="" # READ / WRITE SPIN
    RW_YLM_PHI ="" # READ / WRITE PHI
    RW_CONFIG  ="" # READ / WRITE CONFIG (SPIN - ISOSPIN MATRIX)
    SPIN_FILE  ="" # SPIN TABLE FILE
    YLM_FILE   ="" # KALMAN-YLM FILE
    PHI_FILE   ="" # KALMAN-PHI FILE
    CONFIG_FILE="" # CONFIG FILE (PLACE HOLDER -- USED TO READ GFMC CONFIGS)
    PRODUCT_PARAM_FILE=""  # PARAMETERS FOR A-BODY CLUSTER SYSTEM
    PRODUCT_DECK_FILE =""  # DECK FOR CLUSTERS -- COUPLING INFORMATION
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
class control_t:
#-----------------------------------------------------------------------
    FILE_NAME=""
    BASIS =""                                    #
    BRA_EQ_KET=""
    BRA_TYPE=""                                 # WAVEFUNCTION TYPE FOR BRA
    KET_TYPE=""                                 # WAVEFUNCTION TYPE FOR KET
    INPUT_BRA=wavefunction_input_t("init",[])
    INPUT_KET=wavefunction_input_t("init",[])
    RW_WALK=""                                  # READ / WRITE COORDINATES + PAIR ORDER + RNG + PSISQ
    WALK_FILE=""                                # NAME OF STORED WALK
    LKE=""                                      # SWITCH FOR KINETIC ENERGY
    LEMP=""                                     # SWITCH FOR ELECTROMAGNETIC POTENTIAL
    CONST_FILE=""                               # NAME OF CONSTANTS FILE
    L2BP_FILE=""                                # NAME OF TWO-BODY POT FILE
    L3BP_FILE=""                                # NAME OF THREE-BODY POT FILE
    RNG_SEED=""                                 # RANDOM SEED
    NUM_BLOCKS=""                               # HOW MANY BLOCKS
    BLOCK_SIZE=""                               # TOTAL EVAVLUATIONS IN A BLOCK
    NUM_WALKERS_PER_NODE=""                     # NUMBER OF WALKERS PER MPI_SHARED_MEMORY_REGION
    BURN_IN_COUNT=""                            # HOW MANY CONFIGS TO BURN IN
    NUM_MOVES_BETWEEN=""                        # HOW MANY STEPS BETWEEN EVALUATIONS
    PARTICLE_MAX_DX=""                          # MAX DISTANCE A PARTICLE CAN MOVE
    NPTS=""                                     # GRID SIZE
    NPTS_IN_ONE_FERMI=""                        # NUMBER OF GRID POINTS IN ONE FERMI
    FD_FACTOR=""                                # FINITE DIFFERENCE GRID FACTOR
    BOX_SIZE=""                                 # MAX DISTANCE BEWTEEN TWO CLUSTERS
    SS_LIMIT=""                                 # MAX DISTANCE BEWTEEN s to s shell
    SP_LIMIT=""                                 # MAX DISTANCE BEWTEEN s to p shell
    SD_LIMIT=""                                 # MAX DISTANCE BEWTEEN s+p to sd shell
    SAMPLE_L2=""                                # SAMPLE L^2 POTENTIAL -- massive computational savings for A>=10
    RSAM=""                                     # distance between nucleons -- determines rejection L2 sampling
    PSAM=""                                     # always accept min with this probability
    BIN_WIDTH=""
    BIN_MAXR=""                                 # DENSITY BINS
    LASTP_SAMPLE_TYPE=""                        # SAMPLING FUNCTION SWITCH FOR P(R_sp)
    LPS_A=""
    LPS_B=""                                    # PARAMS FOR LAST PARTICLE SAMPLING FUNCTIONS
    EXTRA_NORM=""                               # NORMS FOR GFMC CALCS (MIXED ESTIMATES)
    EXTRA_NORM_ERROR=""                         # 
    ANC_ENERGY=""                               # RELATIVE ENERGY BETWEEN CLUSTERS
    REG_GAMMA=""                                # GAMMA FOR REGULARIZER
    ECLSTR1=""
    ECLSTR2=""                                  # ENERGY FOR A CLUSTER
    DCLSTR1=""
    DCLSTR2=""                                  # ST.DEV. FOR THE ENERGY OF THE CLUSTER
    DO_GROUP=""                                 # DO GROUP THEORY OR READ IT
    GROUP_FILE=""                               # FILE FOR PHIGEN
    NORTAB_FILE=""                              # FILE FOR PHIGEN
    NLOPT_METHOD=""                             # NLOPT ALGO
    NUM_OPT_WALKS=""                            # HOW MANY TIMES TO LOOP OPTIMIZATION SEARCH
    NUM_OPT_EVALUATIONS=""                      # HOW MANY ENERGY EVALUATIONS TO DO IN OPTIMIZATION
    OPTIMIZATION_INPUT_FILE=""                  # HOW MUCH TO CHANGE EACH VARIATIONAL PARAMETER
    OPTIMIZED_DECK_FILE=""                      # DECK WITH OPTIMAL PARAMETERS
    SCRATCH_DIR=""                              # TEMPORARY STORAGE DIRECTORY
#-----------------------------------------------------------------------
    def __init__(self,file_name_):
#-----------------------------------------------------------------------
        self.FILE_NAME = file_name_
        self.Read()
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