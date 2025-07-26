"""
parameters.py
Nuclei Param Files for nQMCC
"""
#-----------------------------------------------------------------------
class PARAMETERS:
#-----------------------------------------------------------------------
    def __init__(self,file_name_):
#-----------------------------------------------------------------------
        self.FILE_NAME = file_name_ #os.path.splitext(os.path.basename(file_name_))[0]
        self.Read()
#-----------------------------------------------------------------------
    def Read(self):
#-----------------------------------------------------------------------
        file = open(self.FILE_NAME, 'r')
        param_data = [int((l.strip().split())[0]) for l in file.readlines()]
        file.close()
        self.NPART=param_data[0]
        self.NPROT=param_data[1]
        self.ISOSPIN_NT=param_data[2]
        self.NSPART=param_data[3]
        self.NPPART=param_data[4]
        self.NSDPART=param_data[5]
        self.PHI_TYPE=param_data[6]
        self.NBETA=param_data[7]
        self.NPHIM=param_data[8]
        self.NS0=param_data[9]
        self.MXL=param_data[10]
        self.MNYOUNG=param_data[11]
        self.M_SPIN_LIM=param_data[12]
        self.IMA_LIM=param_data[13]
        self.NY_LIM=param_data[14]
        self.KMACO_LIM=param_data[15]
        self.M_L_LIM=param_data[16]
        self.NTERM=param_data[17]
        self.NSTA0=param_data[18]
        self.NSC0=param_data[19]
        self.NORTAB_DIM=param_data[20]
        self.MAX_ICY1=param_data[21]
#-----------------------------------------------------------------------
    def Write(self, out_file):
        file = open(out_file, 'w')
        file.write(str(self.NPART)+"\n")
        file.write(str(self.NPROT)+"\n")
        file.write(str(self.ISOSPIN_NT)+"\n")
        file.write(str(self.NSPART)+"\n")
        file.write(str(self.NPPART)+"\n")
        file.write(str(self.NSDPART)+"\n")
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