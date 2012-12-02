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
    def update(self,_) :
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
class IndicesGenB(wrappedChain.calculable) :
    def __init__(self,collection) :
        self.fixes = collection
        self.stash(["Indices","GenFlavor"])
    def update(self,_) :
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
    def update(self,_) :
        flav = self.source[self.GenFlavor]
        self.value = [i for i in self.source[self.Indices] if abs(flav[i]) not in [5,21] and self.matchesGenWqq(i)]
##############################
class M3(wrappedChain.calculable) :
    def __init__(self, collection) :
        self.fixes = collection
        self.stash(["AdjustedP4","Indices"])
    def update(self,_) :
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

    def update(self,_) : 
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
    def update(self,_) :
        p4 = self.source[self.AdjustedP4]
        i,j = self.source[self.IndicesMinDeltaR]
        self.value = min(p4[i].pt(),p4[j].pt()) * r.Math.VectorUtil.DeltaR(p4[i],p4[j]) if j else -1
######################################
class BScaling(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(["AdjustedP4"])
        self.factor = 1.1
        self.moreName = "Uniform factor %.1f"%self.factor
    def update(self,_) :
        self.value = [self.factor] * len(self.source[self.AdjustedP4])
######################################
class TopCandidateIndices(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(["Indices"])
    def update(self,_) :
        self.value = [ PQ+HL
                       for PQ in itertools.combinations(sorted(self.source[self.Indices]),2)
                       for HL in itertools.permutations(self.source[self.Indices],2)
                       if len(set(PQ+HL))==4 ]
######################################
class HTopCandidateIndices(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(["TopCandidateIndices"])
    def update(self,_) :
        self.value = sorted(set(iPQHL[:3] for iPQHL in self.source[self.TopCandidateIndices]))
#####################################
class RawMassWTopCorrectPQB(wrappedChain.calculable) :
    def __init__(self, collection) :
        self.fixes = collection
        self.stash(['IndicesGenTopPQHL','RawMassWTopPQB'])
    def update(self,_) :
        iPQB = self.source[self.IndicesGenTopPQHL][:3]
        self.value = self.source[self.RawMassWTopPQB][iPQB] if len(set(iPQB))==3 and None not in iPQB else ()
#####################################
class RawMassWTopPQB(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(['HTopCandidateIndices','AdjustedP4','BScaling'])
    def update(self,_) :
        p4 = self.source[self.AdjustedP4]
        bscale = self.source[self.BScaling]
        self.value = {}
        for iPQB in self.source[self.HTopCandidateIndices] :
            _,W,t = np.cumsum([p4[i]*(bscale[i] if i==2 else 1) for i in iPQB])
            self.value[iPQB] = (W.M(),t.M())
######################################
class HTopSigmasPQB(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(['RawMassWTopPQB','RawMassWTopCorrectPQB'])
        self.matrix = None
    def update(self,_) :
        if self.matrix==None:
            self.matrix = dict.__getitem__(self.source.someDict if type(self.source).__name__=='keyTracer' else self.source,
                                           self.RawMassWTopCorrectPQB+'TwoDChiSquared').matrix
        rawM = self.source[self.RawMassWTopPQB]
        self.value = dict( ( iPQH, math.sqrt(np.dot( v, np.dot(self.matrix, v) )) ) for iPQH,v in
                           [( i3, ms+(1,) ) for i3,ms in rawM.items()])
######################################
class HTopSigmasLikelihoodRatioPQB(calculables.secondary) :
    def __init__(self, collection, samples, tag) :
        self.fixes = collection
        for item in ['samples','tag'] : setattr(self,item,eval(item))
        self.stash(['HTopSigmasPQB','IndicesGenTopPQHL'])
        self.max = 5 # sigma

    def update(self,_) :
        self.value = dict( (iPQB,
                            self.likeR.Interpolate(level) if self.likeR else 1)
                           for iPQB,level in self.source[self.HTopSigmasPQB].items() )
    
    def uponAcceptance(self,ev) :
        if not ev['genTopTTbar'] : return
        iCorrect = ev[self.IndicesGenTopPQHL][:3]
        for iPQB,level in ev[self.HTopSigmasPQB].items() :
            self.book.fill( min(level,self.max-1e-6), '%scorrect'%('' if iCorrect==iPQB else 'in'), 250,0,self.max, title = ';HTopSigmasPQB;')

    def setup(self,*_) :
        hists = self.fromCache(['merged'],['%scorrect'%s for s in ['','in']], tag = self.tag)['merged']
        if None in hists.values() :
            self.likeR = None
            print self.name, ': Histograms not found'
        else :
            for h in hists.values() : h.Scale(1./h.Integral(0,h.GetNbinsX()+1))
            self.likeR = hists['correct'].Clone('likelihoodR')
            self.likeR.Divide(hists['incorrect'])
            return hists

    def onlySamples(self) : return ['merged']
    def baseSamples(self) : return self.samples
    def organize(self,org) :
        if org.tag ==self.tag :
            org.mergeSamples(targetSpec = {"name":'merged'}, sources = self.samples )
    def reportCache(self) :
        hists = self.setup()
        if not hists :
            print '%s.setup() failed'%self.name
            return
        fileName = '/'.join(self.outputFileName.split('/')[:-1]+[self.name]) + '.pdf'
        c = r.TCanvas()
        c.Print(fileName +'[')
        hists['correct'].SetLineColor(r.kRed)
        hists['correct'].Draw('hist')
        hists['incorrect'].Draw('hist same')
        c.Print(fileName)
        self.likeR.SetTitle("LikelihoodRatio, Correct:Incorrect; HTopSigmasPQB")
        self.likeR.Draw('hist')
        c.Print(fileName)
        c.Print(fileName +']')
        print 'Wrote : %s'%fileName
######################################
class ProbabilityHTopMasses(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.prior = 0.5
        self.pfactor = (1-self.prior) / self.prior
        self.stash(['HTopSigmasLikelihoodRatioPQB'])
    def update(self,_):
        LRs = self.source[self.HTopSigmasLikelihoodRatioPQB].values()
        sumLRs = sum(LRs)
        self.value = sumLRs / (sumLRs + len(LRs) * self.pfactor )
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
