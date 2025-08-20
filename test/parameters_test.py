import sys
import os
#-----------------------------------------------------------------------
sys.path.append(os.path.join(os.path.dirname(__file__), "../src"))
#-----------------------------------------------------------------------
from parameters import parameters_t
#-----------------------------------------------------------------------
if __name__=="__main__":
  params=parameters_t("test.params")
  params.Write("dummy.out")
#-----------------------------------------------------------------------