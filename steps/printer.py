import time, ROOT as r
from supy import analysisStep,utils
#####################################
class eventPrinter(analysisStep) :

    def __init__(self) :
        self.nHyphens=56

    def uponAcceptance(self,eventVars) :
        print
        print "".ljust(self.nHyphens,"-")
        outString ="run %7d"%eventVars["run"]
        outString+="  event %10d"%eventVars["event"]
        outString+="  ls %#5d"%eventVars["lumiSection"]
        outString+="  bx %4d"%eventVars["bunch"]
        print outString
#####################################
class particleP4Printer(analysisStep) :

    def __init__(self,cs) :
        self.cs = cs
        self.moreName="%s %s" %self.cs
        self.nHyphens=56
        self.p4Name = "%sP4%s"%self.cs
    def select (self,eventVars) :
        p4Vector=eventVars[self.p4Name]

        nParticles=len(p4Vector)
        for iParticle in range(nParticles) :
            particle=p4Vector[iParticle]

            outString = "%s%s %2d" %(self.cs[0],self.cs[1],iParticle)
            outString+="; pT %#6.1f GeV"      %particle.pt()
            outString+="; eta %#6.3f"         %particle.eta()
            outString+="; phi %#6.3f"         %particle.phi()
            print outString

        if nParticles>0 : print
        else :            print "no %s%s found"%self.cs

        return True
#####################################
class metPrinter(analysisStep) :

    def __init__(self,collections) :
        self.collections=collections
        self.moreName = str(self.collections)
        self.nHyphens=56

    def select (self,eventVars) :
        print
        for met in self.collections :
            metVector=eventVars[met]
            outString=met.ljust(15)
            outString+=" pT %#6.1f GeV"%metVector.pt()
            outString+="; phi %#4.1f"  %metVector.phi()
            print outString
        print
        return True
#####################################
