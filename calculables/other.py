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
class jw(wrappedChain.calculable) :
    def __init__(self, fileName = "") :
        self.moreName = "run:ls in %s"%fileName
        with open(fileName) as jsonFile :
            json = eval(''.join(l.strip() for l in jsonFile.readlines()))
        self.json = dict([(int(k),sorted(v)) for k,v in json.items()])

    def inJson(self) :
        run,ls = self.source["run"],self.source['lumiSection']
        return run in self.json and next( ( lo<=ls for lo,hi in self.json[run] if ls<=hi) , False)

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
class KarlsruheDiscriminant(wrappedChain.calculable) :
    def __init__(self, jet, met) :
        self.stash(['M3'],jet)
        self.met = met
        self.moreName = "-8.met if met<40 else m3; %s%s; %s"%(jet+(met,))
    def update(self,_) :
        met = self.source[self.met].pt()
        self.value = -8*met if met<40 else self.source[self.M3]
