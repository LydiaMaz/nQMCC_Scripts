import sys
import os
#-----------------------------------------------------------------------
sys.path.append(os.path.join(os.path.dirname(__file__), "../src"))
#-----------------------------------------------------------------------
from utility import utility_t
#-----------------------------------------------------------------------
if __name__=="__main__":
  util=utility_t("test.util")
  for key,val in util.__dict__.items():
    print(key,val)
#-----------------------------------------------------------------------