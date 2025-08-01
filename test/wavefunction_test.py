import sys
import os
#-----------------------------------------------------------------------
sys.path.append(os.path.join(os.path.dirname(__file__), "../src"))
#-----------------------------------------------------------------------
from deck import deck_t
from wavefunction import wavefunction_t
#-----------------------------------------------------------------------
if __name__=="__main__":
  wf=wavefunction_t("he4.ctrl","./build/bin/",["mpirun", "-np", "4"])
  #e,v=wf.Evaluate(True,"dummy.out")
  opt=deck_t(wf.PARAMS,wf.CTRL.OPTIMIZATION_INPUT_FILE)
  e,v=wf.Optimize(opt,"dummy.dk",True,"dummy.log")
  print(e,v,sep=" : ")
#-----------------------------------------------------------------------