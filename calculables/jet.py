import math,bisect,itertools,ROOT as r
from supy import wrappedChain,calculables,utils
try: import numpy as np
except: np = None

##############################
class AdjustedP4(wrappedChain.calculable) :
    def __init__(self, collection = None, smear = "", jesAbs = 1, jesRel = 0) :
        self.fixes = collection
        self.jesAbs,self.jesRel = jesAbs,jesRel
        self.Smear = smear.join(collection)
        self.stash(['P4'])
        self.moreName = "p4 * %s * %.2f * (1 + %.2f|eta|)" % (smear, jesAbs, jesRel)

    @staticmethod
    def smear(p4,smear) : return p4*smear
    def jes(self,p4) : return p4 * (self.jesAbs * (1+self.jesRel*abs(p4.eta())))

    def update(self,_) :
        self.value = ( utils.hackMap(self.jes, self.source[self.P4]) if self.source['isRealData'] else
                       utils.hackMap(self.smear, self.source[self.P4], self.source[self.Smear]) )
##############################
class Pt(wrappedChain.calculable) :
    def __init__(self,collection=None) :
        self.fixes = collection
        self.stash(["AdjustedP4"])
    def update(self,ignored) :
        self.value = utils.hackMap( lambda x: x.pt(), self.source[self.AdjustedP4] )
##############################
class Indices(wrappedChain.calculable) :
    def __init__(self, collection = None, ptMin = None):
        self.fixes = collection
        self.stash(["Pt"])
        self.ptMin = ptMin
        self.moreName = "pT>=%.1f GeV" % ptMin

    def update(self,_) :
        pts = self.source[self.Pt]
        self.value = sorted( [ i for i,pt in enumerate(pts) if pt>= self.ptMin ],
                             key = pts.__getitem__, reverse = True)
###################################
class IndicesBtagged(wrappedChain.calculable) :
    def __init__(self,collection,tag) :
        self.fixes = collection
        self.stash(["Indices"])
        self.Tag = tag.join(collection)
        self.moreName = "Ordered by %s"%self.Tag
    def update(self,_) :
        self.value = sorted(self.source[self.Indices],
                            key = self.source[self.Tag].__getitem__, reverse = True )
###################################
class BIndices(wrappedChain.calculable) :
    def __init__(self, collection, nMax) :
        self.fixes = collection
        self.stash(["IndicesBtagged"])
        self.nMax = nMax
        self.moreName = "First %d %s"%(nMax,self.IndicesBtagged)
    def update(self,_) : self.value = self.source[self.IndicesBtagged][:self.nMax]
###################################
class IndicesGenB(wrappedChain.calculable) :
    def __init__(self,collection) :
        self.fixes = collection
        self.stash(["Indices","GenFlavor"])
    def update(self,ignored) :
        flav = self.source[self.GenFlavor]
        self.value = [i for i in self.source[self.Indices] if abs(flav[i])==5]
###################################
class IndicesGenWqq(wrappedChain.calculable) :
    def __init__(self,collection) :
        self.fixes = collection
        self.stash(["Indices","AdjustedP4",'GenFlavor'])
    
    def matchesGenWqq(self,index) :
        genP4s = self.source["genP4"]
        p4s = self.source[self.AdjustedP4]
        for iGenQ in self.source["genIndicesWqq"] :
            qP4 = genP4s[iGenQ]
            p4 = p4s[index]
            if r.Math.VectorUtil.DeltaR(p4,qP4) < 0.5 and abs(p4.pt()-qP4.pt()) / qP4.pt() < 0.4 : return True
        return False
    def update(self,ignored) :
        flav = self.source[self.GenFlavor]
        self.value = [i for i in self.source[self.Indices] if abs(flav[i]) not in [5,21] and self.matchesGenWqq(i)]
##############################
class M3(wrappedChain.calculable) :
    def __init__(self, collection) :
        self.fixes = collection
        self.stash(["AdjustedP4","Indices"])
    def update(self,ignored) :
        p4 = self.source[self.AdjustedP4]
        sumP4s = sorted([p4[i]+p4[j]+p4[k] for i,j,k in itertools.combinations(self.source[self.Indices], 3)], key = lambda sumP4 : -sumP4.pt() )
        self.value = sumP4s[0].M() if sumP4s else None
####################################
class CovariantResolution2(wrappedChain.calculable) :
    '''[[xx,xy],[xy,yy]] in the transverse plane.'''

    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(["AdjustedP4","Resolution"])

    @staticmethod
    def matrix(p4,res) :
        phi = p4.phi()
        return p4.Perp2() * res**2 * np.outer(*(2*[[math.cos(phi),math.sin(phi)]]))

    def update(self,ignored) : 
        self.value = utils.hackMap(self.matrix , self.source[self.AdjustedP4] , self.source[self.Resolution] )
#####################################
class ProbabilityGivenBQN(calculables.secondary) :
    def __init__(self, collection = None, bvar = None, binning = (0,0,0), samples = ('',[]), tag = None,) :
        self.fixes = collection
        self.__name = (bvar+self.__class__.__name__).join(self.fixes)
        self.bvar = bvar.join(collection)
        for item in ['binning','samples','tag'] : setattr(self,item,eval(item))
        self.stash(['Indices','IndicesGenB','IndicesGenWqq'])
        self.moreName = (tag if tag!=None else '') + '; ' + ','.join(samples[1] if samples[1] else [samples[0]])
    @property
    def name(self) : return self.__name

    def onlySamples(self) : return [self.samples[0]]

    def setup(self,*_) :
        hists = self.fromCache([self.samples[0]],['B','Q','N'], tag = self.tag)
        self.histsBQN = [hists[self.samples[0]][jetType] for jetType in ['B','Q','N']]
        for hist in filter(None,self.histsBQN) : hist.Scale(1./hist.Integral(0,hist.GetNbinsX()+1),"width")
        
    def uponAcceptance(self,ev) :
        if ev['isRealData'] : return
        indices = ev[self.Indices]
        iB = ev[self.IndicesGenB]
        iQ = ev[self.IndicesGenWqq]
        bvar = ev[self.bvar]
        for i in indices :
            jetType = "B" if i in iB else "Q" if i in iQ else "N"
            self.book.fill(bvar.at(i), jetType, *self.binning, title = ";%s (%s);events / bin"%(self.bvar,jetType))
    
    def update(self,_) :
        self.value = [tuple(hist.GetBinContent(hist.FindFixBin(bvar)) if hist else 0 for hist in self.histsBQN)
                      for bvar in self.source[self.bvar]]
        
    def organize(self,org) :
        if org.tag != self.tag : return
        if self.samples[1] :
            missing = [s for s in self.samples[1] if s not in [ss['name'] for ss in org.samples]]
            if missing: print self.name, "-- no such samples :\n", missing
            org.mergeSamples( targetSpec = {'name':self.samples[0]}, sources = self.samples[1] )
        else: org.mergeSamples( targetSpec = {'name':self.samples[0]}, allWithPrefix = self.samples[0] )
#####################################
class Kt(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(["IndicesMinDeltaR","AdjustedP4"])
        self.moreName = "min(pt_i,pt_j) * dR(i,j); ij with minDR; %s%s"%collection
    def update(self,ignored) :
        p4 = self.source[self.AdjustedP4]
        i,j = self.source[self.IndicesMinDeltaR]
        self.value = min(p4[i].pt(),p4[j].pt()) * r.Math.VectorUtil.DeltaR(p4[i],p4[j]) if j else -1
######################################
class HTopCandidateIndices(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(["Indices","BIndices"])
    def update(self,_) :
        bIndices = self.source[self.BIndices]
        self.value = [(p,q,h) for p,q,h in itertools.permutations(self.source[self.Indices],3) if h in bIndices and q>p]
#####################################
class ComboPQBRawMassWTop(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(['HTopCandidateIndices','AdjustedP4'])
    def update(self,_) :
        p4 = self.source[self.AdjustedP4]
        self.value = {}
        for iPQB in self.source[self.HTopCandidateIndices] :
            _,W,t = np.cumsum([p4[i] for i in iPQB])
            self.value[iPQB] = (W.M(),t.M())
######################################
class ComboPQBDeltaRawMassWTop(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(['ComboPQBRawMassWTop'])
    def update(self,_) : self.value = dict([(key,(val[0]-80.4,val[1]-172.5)) for key,val in self.source[self.ComboPQBRawMassWTop].iteritems()])
######################################

######################################
class __value__(wrappedChain.calculable) :
    def __init__(self, jets = None, index = 0, Btagged = True ) :
        self.fixes = ("%s%s%d"%(jets[0],'B' if Btagged else '', index), jets[1])
        self.stash(["AdjustedP4","Indices","IndicesBtagged"],jets)
        self.index = index
        self.Btagged = Btagged
    def update(self,_) :
        p4 = self.source[self.AdjustedP4]
        indices = self.source[self.IndicesBtagged if self.Btagged else self.Indices]
        self.value = self.function(p4[indices[self.index]]) if len(indices)>self.index else 0
######################################
class pt(__value__) :
    def function(self, x) : return x.pt()
######################################
class absEta(__value__) :
    def function(self,x) : return abs(x.eta())
######################################
class eta(__value__) :
    def function(self,x) : return x.eta()
######################################

######################################
class FourJetPtThreshold(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(['Pt', 'Indices'])
    def update(self,_):
        pt = self.source[self.Pt]
        indices = self.source[self.Indices]
        idPt = [pt[i] for i in indices]
        self.value = 0 if len(idPt)<4 else idPt[3]
######################################
class FourJetAbsEtaThreshold(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(['AdjustedP4', 'Indices'])
    def update(self,_):
        p4 = self.source[self.AdjustedP4]
        indices = self.source[self.Indices]
        idAbsEta = sorted([abs(p4.at(i).eta()) for i in indices])
        self.value = 0 if len(idAbsEta)<4 else idAbsEta[3]
######################################
class MaxAbsEta(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(['AdjustedP4', 'Indices'])
    def update(self,_):
        p4 = self.source[self.AdjustedP4]
        indices = self.source[self.Indices]
        self.value = max(abs(p4.at(i).eta()) for i in indices)
######################################

