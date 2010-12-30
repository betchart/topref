import configuration
#####################################
for module in configuration.stepsFiles() :
    exec("from steps%s import *"%module)
#####################################
def adjustSteps(inSteps, dataOrMc = None) :
    outSteps = []
    blackList = getattr(configuration, "stepsToDisableFor%s"%dataOrMc)()
    histoBlackList = getattr(configuration, "histogramsToDisableFor%s"%dataOrMc)()
    for step in inSteps :
        disable = False
        name = step.__class__.__name__
        if name in blackList : disable = True
        if name == "histogrammer" :
            for item in histoBlackList :
                if item in step.var : disable = True
        outSteps.append(copy.deepcopy(step))
        if disable : outSteps[-1].disable()
    return outSteps
#####################################
def adjustStepsForData(inSteps) : return adjustSteps(inSteps, "Data")
def adjustStepsForMc(inSteps)   : return adjustSteps(inSteps, "Mc")
#####################################
def insertPtHatFilter(inSteps,value) :
    inSteps.insert(0,ptHatFilter(value))
    inSteps[0].ignore()
#####################################
