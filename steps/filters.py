from supy import analysisStep
#####################################
class runLsEvent(analysisStep) :
    def __init__(self, fileName) :
        self.moreName = "run:ls:event in %s"%fileName
        file = open(fileName)
        self.tuples = [ eval("(%s,%s,%s)"%tuple(line.replace("\n","").split(":"))) for line in file]
        file.close()

    def select (self,eventVars) :
        return (eventVars["run"], eventVars["lumiSection"], eventVars["event"]) in self.tuples
#####################################
class run(analysisStep) :
    def __init__(self,runs,acceptRatherThanReject) :
        self.runs = runs
        self.accept = acceptRatherThanReject

        self.moreName = "run%s in list %s" % ( ("" if self.accept else " not"),str(runs) )
        
    def select (self,eventVars) :
        return not ((eventVars["run"] in self.runs) ^ self.accept)
#####################################
