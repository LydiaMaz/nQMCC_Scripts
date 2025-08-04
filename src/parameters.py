"""
parameters.py
Nuclei Param Files for nQMCC
"""
#-----------------------------------------------------------------------
class parameters_t:
#-----------------------------------------------------------------------
    def __init__(self,file_name_):
#-----------------------------------------------------------------------
        self.FILE_NAME = file_name_
        self.Read()
#-----------------------------------------------------------------------
    def Read(self):
#-----------------------------------------------------------------------
        file = open(self.FILE_NAME.strip("\'"), 'r')
        data = [int((l.strip().split())[0]) for l in file.readlines()]
        file.close()
        self.NPART=data[0]
        self.NPROT=data[1]
        self.ISOSPIN_NT=data[2]
        self.NSPART=data[3]
        self.NPPART=data[4]
        self.NDPART=data[5]
        self.PHI_TYPE=data[6]
        self.NBETA=data[7]
        self.NPHIM=data[8]
        self.NS0=data[9]
        self.MXL=data[10]
        self.MNYOUNG=data[11]
        self.M_SPIN_LIM=data[12]
        self.IMA_LIM=data[13]
        self.NY_LIM=data[14]
        self.KMACO_LIM=data[15]
        self.M_L_LIM=data[16]
        self.NTERM=data[17]
        self.NSTA0=data[18]
        self.NSC0=data[19]
        self.NORTAB_DIM=data[20]
        self.MAX_ICY1=data[21]
#-----------------------------------------------------------------------
    def Write(self, out_file):
        file = open(out_file.strip("\'"), 'w')
        file.write(str(self.NPART)+"\n")
        file.write(str(self.NPROT)+"\n")
        file.write(str(self.ISOSPIN_NT)+"\n")
        file.write(str(self.NSPART)+"\n")
        file.write(str(self.NPPART)+"\n")
        file.write(str(self.NDPART)+"\n")
        file.write(str(self.PHI_TYPE)+"\n")
        file.write(str(self.NBETA)+"\n")
        file.write(str(self.NPHIM)+"\n")
        file.write(str(self.NS0)+"\n")
        file.write(str(self.MXL)+"\n")
        file.write(str(self.MNYOUNG)+"\n")
        file.write(str(self.M_SPIN_LIM)+"\n")
        file.write(str(self.IMA_LIM)+"\n")
        file.write(str(self.NY_LIM)+"\n")
        file.write(str(self.KMACO_LIM)+"\n")
        file.write(str(self.M_L_LIM)+"\n")
        file.write(str(self.NTERM)+"\n")
        file.write(str(self.NSTA0)+"\n")
        file.write(str(self.NSC0)+"\n")
        file.write(str(self.NORTAB_DIM)+"\n")
        file.write(str(self.MAX_ICY1))
        file.close()
#-----------------------------------------------------------------------