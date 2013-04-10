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
#####################################
class pileUpRatios(wrappedChain.calculable) :
    def __init__(self, var, fileNames) :
        self.var = var
        hists = []
        for i,name in enumerate(fileNames) :
            tfile = r.TFile.Open(name)
            hists.append(tfile.Get('pileup').Clone())
            hists[-1].Scale(1.0/hists[-1].Integral())
            hists[-1].SetDirectory(0)
            tfile.Close()
        for h in reversed(hists) : h.Divide(hists[0])
        self.hists = hists

    def update(self,_) :
        val = self.source[self.var]
        self.value = [h.GetBinContent(h.FindFixBin(val)) for h in self.hists[1:]]


class ScaleFactors(wrappedChain.calculable):
    def __init__(self,nontrivial): raise Exception("Virtual class")
    def update(self,_):
        i = self.source[self.Indices]
        if not i: row, col = None, None
        else:
            p4 = self.source[self.P4][i[0]]
            col = next((i for i,(lo,hi) in enumerate(self.columns) if lo <= p4.pt() <= hi), None)
            row = next((i for i,(lo,hi) in enumerate(self.rows) if lo<= abs(p4.eta()) <=hi), None)

        if None in [row,col]: self.value = 3*[1.0]
        else:
            cent = self.central[row][col]
            self.value = [cent, cent + self.deltaUp[row][col], cent + self.deltaDn[row][col]]
