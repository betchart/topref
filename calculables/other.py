from supy import wrappedChain,utils
import configuration
import math,ROOT as r
try: import numpy as np
except: np = None

#####################################
class pthatLess(wrappedChain.calculable) :
    def __init__(self, maxPtHat = None) :
        self.fixes = ("","%d"%maxPtHat)
        self.maxPtHat = maxPtHat
    def update(self,ignored) : self.value = None if self.maxPtHat < self.source["genpthat"] else 1
##############################
class Mt(wrappedChain.calculable) :
    @property
    def name(self) :
        return "%sMt%s%s"%(self.fixes[0], self.fixes[1], self.met)
    
    def __init__(self, collection = None, met = None, byHand = True , allowNonIso = False, isSumP4=False) :
        self.met = met
        self.fixes = collection
        self.stash(["Indices","IndicesNonIso","P4"])
        self.byHand = byHand
        self.isSumP4 = isSumP4
        self.allowNonIso = allowNonIso
        self.moreName = "%s%s, %s, byHand=%d"%(collection[0], collection[1], met, byHand)

    def update(self, ignored) :
        if (not self.source[self.Indices]) and not (self.allowNonIso or self.source[self.indicesNonIso]) :
            self.value= -1.0
            return
        index = self.source[self.Indices][0] if self.source[self.Indices] else \
                self.source[self.IndicesNonIso][0]
        lep = self.source[self.P4][index]
        met = self.source[self.met] * (-1 if self.isSumP4 else 1)

        if self.byHand :
            self.value = math.sqrt( 2.0*lep.pt()*met.pt()*(1.0 - math.cos(r.Math.VectorUtil.DeltaPhi(lep, met))) )
        else :
            self.value = (lep+met).Mt()
#####################################
class jsonWeight(wrappedChain.calculable) :
    def __init__(self, fileName = "", acceptFutureRuns = False) :
        self.moreName = "run:ls in %s"%fileName
        self.acceptFutureRuns = acceptFutureRuns
        if self.acceptFutureRuns : self.moreName += " OR future runs"

        self.json = {}
        self.runs = []
        self.maxRunInJson = -1
        if fileName :
            file = open(fileName)
            self.makeIntJson(eval(file.readlines()[0].replace("\n","")))
            file.close()

    def makeIntJson(self, json) :
        for key,value in json.iteritems() :
            self.json[int(key)] = value
        self.maxRunInJson = max(self.json.keys())
        self.runs = self.json.keys()

    def inJson(self) :
        run = self.source["run"]
        if self.acceptFutureRuns and run>self.maxRunInJson : return True
        if not (run in self.runs) : return False
        lumiRanges = self.json[run]
        ls = self.source["lumiSection"]
        for lumiRange in lumiRanges :
            if (ls>=lumiRange[0] and ls<=lumiRange[1]) : return True
        return False

    def update(self, ignored) :
        self.value = 1.0 if self.inJson() else None
#####################################
class PtSorted(wrappedChain.calculable) :
    def __init__(self, collection) :
        self.fixes = collection
        self.stash(["Pt"])
    def update(self,_) :
        pt = self.source[self.Pt]
        self.value = all(i>j for i,j in zip(pt[:-1],pt[1:]))
#####################################
class Covariance(wrappedChain.calculable) :
    def __init__(self, collection) :
        self.fixes = collection
        self.stash(['SigmaXX','SigmaXY','SigmaYY'])
    def update(self,_) :
        self.value = np.array([[self.source[self.SigmaXX],self.source[self.SigmaXY]],
                               [self.source[self.SigmaXY],self.source[self.SigmaYY]]])
#####################################
class KarlsruheDiscriminant(wrappedChain.calculable) :
    def __init__(self, jet, met) :
        self.stash(['M3'],jet)
        self.met = met
        self.moreName = "-8.met if met<40 else m3; %s%s; %s"%(jet+(met,))
    def update(self,_) :
        met = self.source[self.met].pt()
        self.value = -8*met if met<40 else self.source[self.M3]
