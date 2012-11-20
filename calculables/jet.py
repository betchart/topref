import math,bisect,itertools,ROOT as r
from supy import wrappedChain,calculables,utils
try: import numpy as np
except: np = None

##############################
class Indices(wrappedChain.calculable) :
    def __init__(self, collection = None, ptMin = None, etaMax = None):
        self.fixes = collection
        self.stash(["P4"],collection)
        self.pt2Min = ptMin*ptMin
        self.moreName = "pT>=%.1f GeV; |eta|<%.1f"% (ptMin, etaMax)

    def update(self,_) :
        self.value = []
        pt2s    = []

        for i,p4 in enumerate(self.source[self.P4]) :
            pt2 = p4s.at(i).Perp2()
            pt2s.append(pt2)
            if pt2 < self.pt2Min : continue
            elif abs(p4s.at(i).eta()) < self.etaMax :
                self.value.append(i)
            else: other.append(i)
        self.value.sort( key = pt2s.__getitem__, reverse = True)
###################################
class IndicesBtagged(wrappedChain.calculable) :
    '''
    CMS PAS BTV-09-001
    CMS PAS BTV-10-001
    '''
    def __init__(self,collection,tag) :
        self.fixes = collection
        self.stash(["Indices"])
        self.tag = ("%s"+tag+"%s") % xcStrip(collection)
        self.moreName = "Ordered by %s; %s%s"%((tag,)+collection)
    def update(self,ignored) :
        self.value = sorted(self.source[self.Indices],
                            key = self.source[self.tag].__getitem__, reverse = True )
###################################
class IndicesBtagged2(wrappedChain.calculable) :
    '''
    CMS PAS BTV-09-001
    CMS PAS BTV-10-001
    '''
    def __init__(self, collection, tag, threshold = None) :
        self.fixes = collection
        self.stash(["Indices"])
        self.tag = ("%s"+tag+"%s") % xcStrip(collection)
        self.moreName = " (>%g)"%threshold
        self.threshold = threshold

    def update(self,ignored) :
        self.value = []
        for i in self.source[self.Indices] :
            v = self.source[self.tag].at(i)
            if v > self.threshold :
                self.value.append(i)
#####################################
class IndicesGenB(wrappedChain.calculable) :
    def __init__(self,collection) :
        self.fixes = collection
        self.stash(["Indices","CorrectedP4"])
    
    def matchesGenB(self,index) :
        genP4s = self.source["genP4"]
        p4s = self.source[self.CorrectedP4]
        for iGenB in self.source["genIndicesB"] :
            bP4 = genP4s[iGenB]
            p4 = p4s[index]
            if r.Math.VectorUtil.DeltaR(p4,bP4) < 0.6 : return True #and abs(p4.pt()-bP4.pt()) / bP4.pt() < 0.4 : return True
        return False
    def update(self,ignored) : self.value = filter(self.matchesGenB, self.source[self.Indices])
###################################
class IndicesGenWqq(wrappedChain.calculable) :
    def __init__(self,collection) :
        self.fixes = collection
        self.stash(["Indices","CorrectedP4"])
    
    def matchesGenWqq(self,index) :
        genP4s = self.source["genP4"]
        p4s = self.source[self.CorrectedP4]
        for iGenQ in self.source["genIndicesWqq"] :
            qP4 = genP4s[iGenQ]
            p4 = p4s[index]
            if r.Math.VectorUtil.DeltaR(p4,qP4) < 0.5 and abs(p4.pt()-qP4.pt()) / qP4.pt() < 0.4 : return True
        return False
    def update(self,ignored) : self.value = filter(self.matchesGenWqq, self.source[self.Indices])
###################################
class LeadingPt(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(["CorrectedP4","Indices"])
    def update(self,ignored) :
        p4s = self.source[self.CorrectedP4]
        indices = self.source[self.Indices]
        self.value = p4s.at(indices[0]).pt() if len(indices) else None
##############################
class Pt(wrappedChain.calculable) :
    def __init__(self,collection=None) :
        self.fixes = collection
        self.stash(["CorrectedP4"])
    def update(self,ignored) :
        p4s = self.source[self.CorrectedP4]
        self.value = [p4s.at(i).pt() for i in range(len(p4s))]
##############################
class M3(wrappedChain.calculable) :
    def __init__(self, collection) :
        self.fixes = collection
        self.stash(["CorrectedP4","Indices"])
    def update(self,ignored) :
        p4 = self.source[self.CorrectedP4]
        sumP4s = sorted([p4[i]+p4[j]+p4[k] for i,j,k in itertools.combinations(self.source[self.Indices], 3)], key = lambda sumP4 : -sumP4.pt() )
        self.value = sumP4s[0].M() if sumP4s else None
####################################
class CovariantResolution2(wrappedChain.calculable) :
    '''[[xx,xy],[xy,yy]] in the transverse plane.'''

    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(["CorrectedP4","Resolution"])

    @staticmethod
    def matrix(p4,res) :
        phi = p4.phi()
        return p4.Perp2() * res**2 * np.outer(*(2*[[math.cos(phi),math.sin(phi)]]))

    def update(self,ignored) : 
        self.value = utils.hackMap(self.matrix , self.source[self.CorrectedP4] , self.source[self.Resolution] )
#####################################
class ProbabilityGivenBQN(calculables.secondary) :
    def __init__(self, collection = None, bvar = None, binning = (0,0,0), samples = ('',[]), tag = None,) :
        self.fixes = collection
        self.__name = ('%s'+bvar+self.__class__.__name__+'%s')%self.fixes
        self.bvar = ("%s"+bvar+"%s")%xcStrip(collection)
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
        self.stash(["IndicesMinDeltaR","CorrectedP4"])
        self.moreName = "min(pt_i,pt_j) * dR(i,j); ij with minDR; %s%s"%collection
    def update(self,ignored) :
        p4 = self.source[self.CorrectedP4]
        i,j = self.source[self.IndicesMinDeltaR]
        self.value = min(p4[i].pt(),p4[j].pt()) * r.Math.VectorUtil.DeltaR(p4[i],p4[j]) if j else -1
######################################
class ComboPQBRawMassWTop(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(['Indices','CorrectedP4'])
    def update(self,_) :
        p4 = self.source[self.CorrectedP4]
        self.value = {}
        for iPQB in itertools.permutations(self.source[self.Indices],3) :
            if iPQB[0]>iPQB[1] : continue
            _,W,t = np.cumsum([p4[i] for i in iPQB])
            self.value[iPQB] = (W.M(),t.M())
######################################
class ComboPQBDeltaRawMassWTop(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(['ComboPQBRawMassWTop'])
    def update(self,_) : self.value = dict([(key,(val[0]-80.4,val[1]-172)) for key,val in self.source[self.ComboPQBRawMassWTop].iteritems()])
######################################

######################################
class __value__(wrappedChain.calculable) :
    def __init__(self, jets = None, index = 0, Btagged = True ) :
        self.fixes = ("%s%s%d"%(jets[0],'B' if Btagged else '', index), jets[1])
        self.stash(["CorrectedP4","Indices","IndicesBtagged"],jets)
        self.index = index
        self.Btagged = Btagged
    def update(self,_) :
        p4 = self.source[self.CorrectedP4]
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
        self.stash(['CorrectedP4', 'Indices'])
    def update(self,_):
        p4 = self.source[self.CorrectedP4]
        indices = self.source[self.Indices]
        idAbsEta = sorted([abs(p4.at(i).eta()) for i in indices])
        self.value = 0 if len(idAbsEta)<4 else idAbsEta[3]
######################################
class MaxAbsEta(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(['CorrectedP4', 'Indices'])
    def update(self,_):
        p4 = self.source[self.CorrectedP4]
        indices = self.source[self.Indices]
        self.value = max(abs(p4.at(i).eta()) for i in indices)
######################################

