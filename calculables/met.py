import math,ROOT as r
from supy import wrappedChain,utils
##############################
class AdjustedP4(wrappedChain.calculable) :
    def __init__(self, met = None, jet = None, smear="", z = None) :
        self.fixes = met
        self.stash(["P4"])
        self.dSmear = ("DeltaMET"+smear).join(jet) if smear else None
        self.z = z
        self.moreName = " + ".join(filter(None,[self.P4,self.dSmear]))

    def update(self,_) :
        met = self.source[self.P4]
        self.value = utils.LorentzV(met.pt(),met.eta(),met.phi(),met.mass())
        if not self.source["isRealData"] :
            smear = self.source[self.dSmear]
            self.value += utils.LorentzV(smear.pt(),smear.eta(),smear.phi(),smear.mass())
#####################################
class Covariance(wrappedChain.calculable) :
    def __init__(self, collection) :
        self.fixes = collection
        self.stash(['SigmaXX','SigmaXY','SigmaYY'])
    def update(self,_) :
        self.value = np.array([[self.source[self.SigmaXX],self.source[self.SigmaXY]],
                               [self.source[self.SigmaXY],self.source[self.SigmaYY]]])
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
