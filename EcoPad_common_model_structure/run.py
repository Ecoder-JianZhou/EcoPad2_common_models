# all model must have the same form:
#     ./modelName pars.txt forcing.txt outputDir other_setting

from ecopadLibs import ecopadObj
# from .teco_spruce_v2 import myflib
# from teco_spruce_v2 import myflib
import teco 

def run(modelName, mode):
    '''
        arg: modelName means which model to run
             mode      means which mode to run
                       modes: simulation
                              spinup
                              data_assimilation
                              forecasting
    '''
    runCase = ecopadObj("teco_spruce_v2")
    print(runCase.modelName)
    runCase.run_data_assimilation()

run("a", "b")
print("call teco ...")
teco.TECO_simu("a","b")