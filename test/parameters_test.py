import sys
import os
#-----------------------------------------------------------------------
sys.path.append(os.path.join(os.path.dirname(__file__), "../src"))
#-----------------------------------------------------------------------
from parameters import PARAMETERS
#-----------------------------------------------------------------------
if __name__=="__main__":
  params=PARAMETERS("test.params")
  params.Read()
  params.Write("dummy.out")
#-----------------------------------------------------------------------