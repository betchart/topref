import math,bisect,itertools,ROOT as r
from supy import wrappedChain,calculables,utils
try: import numpy as np
except: np = None

##############################
class AdjustedP4(wrappedChain.calculable) :
    def __init__(self, collection = None, smear = "", jec = 0 ) :
        self.fixes = collection
        self.Smear = smear.join(collection)
        self.jec = jec
        self.stash(['P4','JecUnc','JecFactor'])
        self.moreName = "p4 * %s * (1 + %d * jecUnc/jecFac)" % (smear, jec)

    @staticmethod
    def smear(p4,smear) : return p4*smear
    def jes(self,p4,fac,unc) : return p4 * ( 1 + self.jec * unc / fac )

    def update(self,_) :
        self.value = ( utils.hackMap(self.jes, self.source[self.P4], self.source[self.JecFactor], self.source[self.JecUnc]) if self.source['isRealData'] else
                       utils.hackMap(self.smear, self.source[self.P4], self.source[self.Smear]) )
##############################
class DeltaMETJEC(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stashed = ['P4','JecFactor','JecUnc']
        self.stash(self.stashed)

    @staticmethod
    def delta(p4,fac,unc) : return p4 * (-unc/fac)

    def update(self,_) :
        self.value = sum(utils.hackMap(self.delta, *[self.source[getattr(self,item)] for item in self.stashed]), utils.LorentzV())
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
class GenIndex(wrappedChain.calculable) :
    def __init__(self,collection) :
        self.fixes = collection
        self.stash(["AdjustedP4","GenFlavor"])

    def genIndex(self,jet,jflav) :
        id = self.source['genPdgId']
        p4 = self.source['genP4']
        return next((i for i in range(len(id)) if id[i]==jflav and r.Math.VectorUtil.DeltaR(p4[i],jet)<0.5), None)

    def update(self,_) :
        self.value = utils.hackMap(self.genIndex, self.source[self.AdjustedP4], self.source[self.GenFlavor])
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
##############################
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

    def reportCache(self) :
        optStat = r.gStyle.GetOptStat()
        r.gStyle.SetOptStat(0)
        self.setup(None)
        if None in self.histsBQN :
            print '%s.setup() failed'%self.name
            r.gStyle.SetOptStat(optStat)
            return
        fileName = '/'.join(self.outputFileName.split('/')[:-1]+[self.name]) + '.pdf'
        c = r.TCanvas()
        c.Print(fileName +'[')
        leg = r.TLegend(0.4,0.55,0.7,0.85)
        leg.SetHeader("jet flavor")
        for h in self.histsBQN :
            h.Fill(h.GetBinCenter(1), h.GetBinContent(0))
            h.SetBinContent(0,0)
        height = 1.1 * max(h.GetMaximum() for h in self.histsBQN)
        for i,(f,color) in enumerate(zip('BQN',[r.kRed,r.kBlue,r.kGreen])) :
            h = self.histsBQN[i]
            h.SetTitle(";%s;"%h.GetXaxis().GetTitle().split()[0])
            h.SetLineColor(color)
            h.SetLineWidth(2)
            h.SetMarkerColor(color)
            h.SetMaximum(height)
            h.SetMinimum(0)
            h.Draw("hist" + ("same" if i else ""))
            leg.AddEntry(h,{"B":"B jets","Q":"jets from W","N":"other (gluon)"}[f],'l')
        leg.Draw()
        c.Print(fileName)
        c.Print(fileName +']')
        print 'Wrote : %s'%fileName
        r.gStyle.SetOptStat(optStat)
######################################
class ScalingBQN(calculables.secondary) :
    def __init__(self, collection = None , samples = [], tag = None) :
        self.fixes = collection
        self.stash(['Indices','IndicesGenB','IndicesGenWqq','AdjustedP4','GenIndex'])
        self.samples = samples
        self.tag = tag
        self.etasMax = sorted([ 1.131, 1.653, 2.51 ])

    def etaBin(self,absEta) : return next( j for j,eMax in enumerate(self.etasMax) if absEta<eMax)

    def update(self,_) :
        jets = self.source[self.AdjustedP4]
        def scales(p4) :
            etaBin = self.etaBin( abs(p4.eta()) )
            pt = p4.pt()
            lRs = (self.hists[f+str(etaBin)].Interpolate(pt) for f in 'BQN')
            return [math.exp(lR) for lR in lRs] if self.hists else [1,1,1]

        self.value = dict(zip('BQN', zip(*[scales(jets[i]) for i in range(len(jets))]) ) )

    def uponAcceptance(self,ev) :
        if ev['isRealData'] : return
        indices = ev[self.Indices]
        iB = ev[self.IndicesGenB]
        iQ = ev[self.IndicesGenWqq]
        p4 = ev[self.AdjustedP4]
        iGen = ev[self.GenIndex]
        genP4 = ev['genP4']
        for i in indices :
            if iGen[i] == None : continue
            jetType = "B" if i in iB else "Q" if i in iQ else "N"
            etaBin = self.etaBin( abs(p4[i].eta()) )
            ptRatio = math.log( genP4[iGen[i]].E() / p4[i].E() )
            binPt = min( 300-1e-6, p4[i].pt())
            self.book.fill( (binPt,ptRatio), jetType+str(etaBin)+'_2', (50,80), (0,-0.4), (300,0.4), title = "%s, %.3f<|#eta|<%.3f;pt;<log(E.parton/E.jet)>"%(jetType,
                                                                                                                                                              self.etasMax[etaBin-1] if etaBin else 0,
                                                                                                                                                              self.etasMax[etaBin]))
    def setup(self,*_) :
        names = [f+str(i) for f in "BQN" for i in range(len(self.etasMax))]
        names2 = [n+"_2" for n in names]
        hists = self.fromCache(['merged'],names2, tag = self.tag)['merged']
        if None in hists.values() :
            print self.name, ": Histograms not found."
            self.hists = {}
            return

        def median(h) :
            if not h.Integral() : return 0
            q = np.array([0.5],'d')
            v = np.array([-100.0],'d')
            h.GetQuantiles(1,v,q)
            return v[0]

        self.hists = {}
        for n,n2 in zip(names,names2) :
            h = hists[n2].ProjectionX()
            self.hists[n] = r.TH1D(n,';'.join([h.GetTitle(),h.GetXaxis().GetTitle(),h.GetYaxis().GetTitle()]), h.GetNbinsX(), h.GetBinLowEdge(1), h.GetBinLowEdge(1+h.GetNbinsX()))
            for i in range(1,1+h.GetNbinsX()) :
                bins = ("",i,i) if i < 15 else ("",i-1,i+1)
                self.hists[n].SetBinContent(i, median(hists[n2].ProjectionY(*bins)))
                self.hists[n].SetBinError(i, h.GetBinError(i))

    def reportCache(self) :
        optStat = r.gStyle.GetOptStat()
        r.gStyle.SetOptStat(0)
        self.setup(None)
        if not self.hists :
            print '%s.setup() failed'%self.name
            r.gStyle.SetOptStat(optStat)
            return
        fileName = '/'.join(self.outputFileName.split('/')[:-1]+[self.name]) + '.pdf'
        c = r.TCanvas()
        c.Print(fileName +'[')
        for f in 'BQN' :
            leg = r.TLegend(0.6,0.6,0.9,0.9)
            leg.SetHeader("#eta range")
            for i,color in enumerate([r.kRed,r.kBlue,r.kGreen]) :
                h = self.hists[f+str(i)]
                label = h.GetTitle().split(',')[1]
                h.SetTitle(h.GetTitle().split(',')[0])
                h.SetLineColor(color)
                h.SetMarkerColor(color)
                h.SetMaximum(0.4)
                h.SetMinimum(-0.4)
                h.Draw("same" if i else "")
                leg.AddEntry(h,label,'l')
            leg.Draw()
            c.Print(fileName)
        c.Print(fileName +']')
        print 'Wrote : %s'%fileName
        r.gStyle.SetOptStat(optStat)

    def onlySamples(self) : return ['merged']
    def baseSamples(self) : return self.samples
    def organize(self,org) :
        if org.tag == self.tag :
            org.mergeSamples(targetSpec = {"name":'merged'}, sources = self.samples )
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
