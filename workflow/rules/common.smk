from glob import glob
import pandas as pd

def get_zz_folders():
    #return glob("phastest_results/ZZ_*")
    return [folder for folder in os.listdir("phastest_results") if os.path.isdir(os.path.join("phastest_results",folder))]
regions = {}

