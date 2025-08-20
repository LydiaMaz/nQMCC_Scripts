"""
deck.py
nQMCC variational w.f. parameters_t
"""
from parameters import parameters_t
from copy import deepcopy
#-----------------------------------------------------------------------
class spatial_symmetry_t:
#-----------------------------------------------------------------------
    def __init__(self,params: parameters_t,data_):
#-----------------------------------------------------------------------
        data=deepcopy(data_)
#-----------------------------------------------------------------------
        self.IDX=5
        self.BLSNAME=data[0][0]
#-----------------------------------------------------------------------
        if params.PHI_TYPE != 0:
            self.IDX=self.IDX+2
            self.L_CORE,self.S_CORE,self.J_CORE,self.T_CORE,self.TZ_CORE=data[1][:5]
            self.L_SDSH,self.S_SDSH,self.J_SDSH,self.T_SDSH,self.TZ_SDSH=data[2][:5]
            data=data[3:]
        else:
            data=data[1:]
#-----------------------------------------------------------------------
        self.BETALSN=data[0][0]
        self.LFCSP,self.SPU,self.SPV,self.SPR,self.SPA,self.SPB,self.SPC,self.SPK,self.SPL=data[1][:9]
#-----------------------------------------------------------------------
        if params.NPPART >= 2:
            self.IDX=self.IDX+1
            self.LFCPP,self.PPU,self.PPV,self.PPR,self.PPA,self.PPB,self.PPC,self.PPK,self.PPL=data[2][:9]
            data=data[3:]
        else:
            data=data[2:]
#-----------------------------------------------------------------------
        self.LFSP,self.E_OR_V,self.LSCAT,self.LMU1,self.LMU2,self.LQNUM,self.LNODES=data[0][:7]
#-----------------------------------------------------------------------
        if int(self.LSCAT) == 1:
            self.IDX=self.IDX+1
            self.BSCAT=data[1][0]
            data=data[2:]
        else:
            data=data[1:]
#-----------------------------------------------------------------------
        self.WSE,self.WSV,self.WSR,self.WSA,self.WBRHO,self.WBALPH=data[0][:6]
#-----------------------------------------------------------------------
        if params.NDPART != 0:
            self.IDX=self.IDX+4
            data=data[1:]
#-----------------------------------------------------------------------
            self.LFCSD,self.SDU,self.SDV,self.SDR,self.SDA,self.SDB,self.SDC,self.SDK,self.SDL=data[0][:9]
            self.LFCPD,self.PDU,self.PDV,self.PDR,self.PDA,self.PDB,self.PDC,self.PDK,self.PDL=data[1][:9]
#-----------------------------------------------------------------------
            if params.NDPART >= 2: 
                self.IDX=self.IDX+1
                self.LFCDD,self.DDU,self.DDV,self.DDR,self.DDA,self.DDB,self.DDC,self.DDK,self.DDL=data[2][:9]
                data=data[3:]
            else:
                data=data[2:]
#-----------------------------------------------------------------------
            self.D_LFSP,self.D_E_OR_V,self.D_LSCAT,self.D_LMU1,self.D_LMU2,self.D_LQNUM,self.D_LNODES=data[0][:7]
#-----------------------------------------------------------------------
            if int(self.LSCAT) == 1:
                self.IDX=self.IDX+1
                self.D_BSCAT=data[1][0]
                data=data[2:]
            else:
                data=data[1:]
#-----------------------------------------------------------------------
            self.D_WSE,self.D_WSV,self.D_WSR,self.D_WSA,self.D_WBRHO,self.D_WBALPH=data[0][:6]
#-----------------------------------------------------------------------
    def Write(self,params: parameters_t, file, last_char="\n"):
#-----------------------------------------------------------------------
        file.write(self.BLSNAME+"\n")
#-----------------------------------------------------------------------
        if params.PHI_TYPE != 0:
            file.write(" ".join([self.L_CORE,self.S_CORE,self.J_CORE,self.T_CORE,self.TZ_CORE])+"\n")
            file.write(" ".join([self.L_SDSH,self.S_SDSH,self.J_SDSH,self.T_SDSH,self.TZ_SDSH])+"\n")
#-----------------------------------------------------------------------
        file.write(self.BETALSN+"\n")
        file.write(" ".join([self.LFCSP,self.SPU,self.SPV,self.SPR,self.SPA,self.SPB,self.SPC,self.SPK,self.SPL])+"\n")
#-----------------------------------------------------------------------
        if params.NPPART >= 2:
            file.write(" ".join([self.LFCPP,self.PPU,self.PPV,self.PPR,self.PPA,self.PPB,self.PPC,self.PPK,self.PPL])+"\n")
#-----------------------------------------------------------------------
        file.write(" ".join([self.LFSP,self.E_OR_V,self.LSCAT,self.LMU1,self.LMU2,self.LQNUM,self.LNODES])+"\n")
#-----------------------------------------------------------------------
        if int(self.LSCAT) == 1:
            file.write(self.BSCAT+"\n")
#-----------------------------------------------------------------------
        file.write(" ".join([self.WSE,self.WSV,self.WSR,self.WSA,self.WBRHO,self.WBALPH])+"\n")
#-----------------------------------------------------------------------
        if params.NDPART != 0:
#-----------------------------------------------------------------------
            file.write(" ".join([self.LFCSD,self.SDU,self.SDV,self.SDR,self.SDA,self.SDB,self.SDC,self.SDK,self.SDL])+"\n")
            file.write(" ".join([self.LFCPD,self.PDU,self.PDV,self.PDR,self.PDA,self.PDB,self.PDC,self.PDK,self.PDL])+"\n")
#-----------------------------------------------------------------------
            if params.NDPART >= 2: 
                file.write(" ".join([self.LFCDD,self.DDU,self.DDV,self.DDR,self.DDA,self.DDB,self.DDC,self.DDK,self.DDL])+"\n")
#-----------------------------------------------------------------------
            file.write(" ".join([self.D_LFSP,self.D_E_OR_V,self.D_LSCAT,self.D_LMU1,self.D_LMU2,self.D_LQNUM,self.D_LNODES])+"\n")
#-----------------------------------------------------------------------
            if int(self.D_LSCAT) == 1:
                file.write(self.D_BSCAT+"\n")
#-----------------------------------------------------------------------
            file.write(" ".join([self.D_WSE,self.D_WSV,self.D_WSR,self.D_WSA,self.D_WBRHO,self.D_WBALPH])+last_char)
#-----------------------------------------------------------------------
class deck_t:
#-----------------------------------------------------------------------
    def __init__(self, params: parameters_t, file_name_: str, read=True):
#-----------------------------------------------------------------------
        self.FILE_NAME = file_name_
        if read: self.Read(params)
#-----------------------------------------------------------------------
    def Read(self, params: parameters_t):
#-----------------------------------------------------------------------
        file = open(self.FILE_NAME.strip("\'"), 'r')
        data = [(l.strip().split()) for l in file.readlines()]
        file.close()
        self.NAME=data[0][0]
        if self.NAME[-1] != "\'": self.NAME=self.NAME+"\'"
        self.PARITY,self.TOTAL_J,self.TOTAL_JZ,self.TOTAL_T,self.TOTAL_TZ=data[1][:5]
        self.LWF,self.LSC,self.LOPC,self.LCUT=data[2][:4]
        self.ESEP=data[3][:4]
        self.ETA=data[4][:2]
        self.ZETA=data[5][:2]
        self.FSCAL=data[6][:8]
        self.AC=data[7][:8]
        self.AA=data[8][:8]
        self.AR=data[9][:8]
        self.ALPHA=data[10][:8]
        self.BETA=data[11][:8]
        self.GAMMA=data[12][:8]
        self.UUR,self.UUA,self.UUW=data[13][:3]
        self.SSH_LFSP,self.SSH_E_OR_V,self.SSH_LSCAT,self.SSH_LMU1,self.SSH_LMU2,self.SSH_LQNUM,self.SSH_LNODES=data[14][:7]
        self.SSH_WSE,self.SSH_WSV,self.SSH_WSR,self.SSH_WSA,self.SSH_WBRHO,self.SSH_WBALPH=data[15][:6]
        self.DELTA,self.EPSILON,self.THETA,self.UPSILON,self.RSCAL,self.USCAL=data[16][:6]
        self.QPS1,self.QPS2=data[17][:2]
        self.QSSS1,self.QSSS2=data[18][:2]
# ----------------------------------------------------------------------
        if (params.NPPART >= 1):
            data=data[19:]
            self.QSSP1,self.QSSP2=data[0][:2]
            data=data[1:]
        if (params.NPPART >= 2):
            self.QSPP1,self.QSPP2=data[0][:2]
            data=data[1:]
        if (params.NPPART >= 3):
            self.QPPP1,self.QPPP2=data[0][:2]
            data=data[1:]
# ----------------------------------------------------------------------
        if (params.NDPART >= 1):
            self.QSPD1,self.QSPD2=data[0][:2]
            self.QSSD1,self.QSSD2=data[1][:2]
            self.QPPD1,self.QPPD2=data[2][:2]
            data=data[3:]
        if (params.NDPART >= 2):
            self.QSDD1,self.QSDD2=data[0][:2]
            self.QPDD1,self.QPDD2=data[1][:2]
            data=data[2:]
        if (params.NDPART >= 3):
            self.QDDD1,self.QDDD2=data[0][:2]
            data=data[1:]
# ----------------------------------------------------------------------
        self.SS=[]
# ----------------------------------------------------------------------
        if params.NPPART >= 1:
            for b in range(params.NBETA):
                self.SS.append(spatial_symmetry_t(params,data))
                if b != params.NBETA: data=data[self.SS[b].IDX:]
#-----------------------------------------------------------------------
    def Write(self,params: parameters_t,out_file):
# ----------------------------------------------------------------------
        file = open(out_file.strip("\'"), 'w')
        file.write(self.NAME+"\n")
        file.write(" ".join([self.PARITY,self.TOTAL_J,self.TOTAL_JZ,self.TOTAL_T,self.TOTAL_TZ])+"\n")
        file.write(" ".join([self.LWF,self.LSC,self.LOPC,self.LCUT])+"\n")
        file.write(" ".join(self.ESEP)+"\n")
        file.write(" ".join(self.ETA)+"\n")
        file.write(" ".join(self.ZETA)+"\n")
        file.write(" ".join(self.FSCAL)+"\n")
        file.write(" ".join(self.AC)+"\n")
        file.write(" ".join(self.AA)+"\n")
        file.write(" ".join(self.AR)+"\n")
        file.write(" ".join(self.ALPHA)+"\n")
        file.write(" ".join(self.BETA)+"\n")
        file.write(" ".join(self.GAMMA)+"\n")
        file.write(" ".join([self.UUR,self.UUA,self.UUW])+"\n")
        file.write(" ".join([self.SSH_LFSP,self.SSH_E_OR_V,self.SSH_LSCAT,self.SSH_LMU1,self.SSH_LMU2,self.SSH_LQNUM,self.SSH_LNODES])+"\n")
        file.write(" ".join([self.SSH_WSE,self.SSH_WSV,self.SSH_WSR,self.SSH_WSA,self.SSH_WBRHO,self.SSH_WBALPH])+"\n")
        file.write(" ".join([self.DELTA,self.EPSILON,self.THETA,self.UPSILON,self.RSCAL,self.USCAL])+"\n")
        file.write(" ".join([self.QPS1,self.QPS2])+"\n")
        file.write(" ".join([self.QSSS1,self.QSSS2])+"\n")
# ----------------------------------------------------------------------
        if (params.NPPART >= 1):
            file.write(" ".join([self.QSSP1,self.QSSP2])+"\n")
        if (params.NPPART >= 2):
            file.write(" ".join([self.QSPP1,self.QSPP2])+"\n")
        if (params.NPPART >= 3):
            file.write(" ".join([self.QPPP1,self.QPPP2])+"\n")
# ----------------------------------------------------------------------
        if (params.NDPART >= 1):
            file.write(" ".join([self.QSPD1,self.QSPD2])+"\n")
            file.write(" ".join([self.QSSD1,self.QSSD2])+"\n")
            file.write(" ".join([self.QPPD1,self.QPPD2])+"\n")
        if (params.NDPART >= 2):
            file.write(" ".join([self.QSDD1,self.QSDD2])+"\n")
            file.write(" ".join([self.QPDD1,self.QPDD2])+"\n")
        if (params.NDPART >= 3):
            file.write(" ".join([self.QDDD1,self.QDDD2])+"\n")
# ----------------------------------------------------------------------
        if params.NPPART >= 1:
            for b in range(params.NBETA):
                if b < params.NBETA: self.SS[b].Write(params,file)
        file.close()
# ----------------------------------------------------------------------
    def UpdateFloats(self, params: parameters_t, precision: int):
# ----------------------------------------------------------------------
        for key,val in self.__dict__.items():
            if (key not in ["FILE_NAME","NAME","SS"]):
                if isinstance(val,list):
                    self.__dict__[key]=[f"{float(i):.{precision}f}" for i in self.__dict__[key]]
                else:
                    if "." in val: self.__dict__[key]=f"{float(val):.{precision}f}"
# ----------------------------------------------------------------------
        if params.NBETA >=1:
            for ss in self.SS:
                for key,val in ss.__dict__.items():
                    if isinstance(val,list):
                        ss.__dict__[key]=[f"{float(i):.{precision}f}" for i in self.__dict__[key]]
                    elif isinstance(val,int):
                        pass
                    else:
                        if "." in val: ss.__dict__[key]=f"{float(val):.{precision}f}"
# ----------------------------------------------------------------------
def GenerateZeroDeck(params: parameters_t, file_name):
# ----------------------------------------------------------------------
    dk=deck_t(params,file_name)
# ----------------------------------------------------------------------
    for key,val in dk.__dict__.items():
        if (key not in ["FILE_NAME","NAME","SS"]):
            if isinstance(val,list):
                dk.__dict__[key]=["0." for i in dk.__dict__[key]]
            else:
                if "." in val: dk.__dict__[key]="0."
# ----------------------------------------------------------------------
    if params.NBETA >=1:
        for ss in dk.SS:
            for key,val in ss.__dict__.items():
                if (key not in ["BLSNAME"]):
                    if isinstance(val,list):
                        ss.__dict__[key]=["0." for i in dk.__dict__[key]]
                    elif isinstance(val,int):
                        pass
                    else:
                        if "." in val: ss.__dict__[key]="0."
# ----------------------------------------------------------------------
    return dk
# ----------------------------------------------------------------------
def GenerateOptFile(params: parameters_t, dk_: deck_t, file_name, instructions: list):
# ----------------------------------------------------------------------
    opt=GenerateZeroDeck(params,dk_.FILE_NAME)
    opt.FILE_NAME=file_name
    try:
        for i in instructions:
            if i["ss"]:
                opt.SS[i["ss_idx"]].__dict__[i["key"]]=str(float(dk_.SS[i["ss_idx"]].__dict__[i["key"]])*i["scale"]+i["flat"])
            else:
                if isinstance(opt.__dict__[i["key"]],list):
                    if i["all"]:
                        opt.__dict__[i["key"]]=[str(float(v)*i["scale"]+i["flat"]) for v in dk_.__dict__[i["key"]]]
                    else:
                        opt.__dict__[i["key"]][i["idx"]]=str(float(dk_.__dict__[i["key"]][i["idx"]])*i["scale"]+i["flat"])
                else:
                    opt.__dict__[i["key"]]=str(float(dk_.__dict__[i["key"]])*i["scale"]+i["flat"])
    except KeyError:
        print(f"GenerateOptFile *** BAD KEY: {i["key"]} ***")
    opt.Write(params,file_name)
    return opt
# ----------------------------------------------------------------------