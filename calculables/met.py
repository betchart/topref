import math,ROOT as r
from supy import wrappedChain,utils
try: import numpy as np
except: pass

##############################
class AdjustedP4(wrappedChain.calculable) :
    def __init__(self, met = None, jet = None, smear="", djec = 0) :
        self.fixes = met
        self.stash(["P4"])
        self.djecfactor = djec
        self.dJEC = "DeltaMETJEC".join(jet)
        self.dSmear = ("DeltaMET"+smear).join(jet) if smear else None
        self.moreName = " + ".join(filter(None,[self.P4,self.dSmear,'dJEC(%d)'%djec]))

    def update(self,_) :
        djec = utils.LorentzV() if not self.djecfactor else self.source[self.dJEC]*self.djecfactor
        djer = utils.LorentzV() if self.source['isRealData'] or not self.dSmear else self.source[self.dSmear]
        self.value = self.source[self.P4] + djec + djer
#####################################
class Covariance(wrappedChain.calculable) :
    def __init__(self, collection) :
        self.fixes = collection
        self.stash(['SigmaXX','SigmaXY','SigmaYY'])
    def update(self,_) :
        XY = self.source[self.SigmaXY]
        self.value = (np.array([[1000,0],[0,1000]]) if math.isnan(XY) else
                      np.array([[self.source[self.SigmaXX],XY ],
                                [XY ,self.source[self.SigmaYY]]]) )
#####################################
class MetMt(wrappedChain.calculable) :
    def __init__(self, lepton = None, met = None, byHand = True) :
        self.fixes = lepton
        self.stash(["Indices","P4"])
        self.met = met
        self.byHand = byHand
        self.moreName = "%s, byHand=%d"%(met, byHand)

    def update(self,_) :
        index = next(iter(self.source[self.Indices]), None)
        if index==None :
            self.value = -1.0
            return
        lep = self.source[self.P4][index]
        met = self.source[self.met]

        if self.byHand :
            self.value = math.sqrt( 2.0*lep.pt()*met.pt()*(1.0 - math.cos(r.Math.VectorUtil.DeltaPhi(lep, met))) )
        else :
            self.value = (lep+met).Mt()
