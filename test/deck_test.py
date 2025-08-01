import sys
import os
#-----------------------------------------------------------------------
sys.path.append(os.path.join(os.path.dirname(__file__), "../src"))
#-----------------------------------------------------------------------
from deck import deck_t,GenerateOptFile
from parameters import parameters_t
#-----------------------------------------------------------------------
if __name__=="__main__":
  params=parameters_t("test.params")
  dk=deck_t(params,"test.dk")
  dk.Write(params,"dummy.out")
  instr=[{"ss":True,"ss_idx":3,"key":"WSE","scale":1, "flat":-171231.5}\
  ,{"ss":False,"key":"SSH_WSE","scale":1, "flat":-2342.27777777777777342}\
  ,{"ss":False,"key":"ALPHA","all":True,"scale":1, "flat":-69.420}\
  ,{"ss":False,"key":"BETA","all":False,"idx":2,"scale":1, "flat":-14}]
  opt=GenerateOptFile(params, dk, "dummy.out",instr)
  opt.UpdateFloats(params,4)
  opt.Write(params,"dummy.out")
#-----------------------------------------------------------------------