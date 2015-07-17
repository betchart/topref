import math,bisect,itertools,operator,ROOT as r
from supy import wrappedChain,calculables,utils,whereami
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
class Moments2Sum(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(['Eta2Moment','Phi2Moment'])
        self.moreName = '%s%s; eta2moment+phi2moment'%collection
    def update(self,_) :
        self.value = utils.hackMap(operator.add,self.source[self.Eta2Moment],self.source[self.Phi2Moment])
##############################
class Unsmeared(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(['Eta2Moment'])
    def update(self,_) :
        self.value = [1.0]*len(self.source[self.Eta2Moment])
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
    def __init__(self, collection = None, bvar = None, binning = (0,0,0), samples = [], tag = None,) :
        self.fixes = collection
        self.__name = (bvar+self.__class__.__name__).join(self.fixes)
        self.bvar = bvar.join(collection)
        for item in ['binning','samples','tag'] : setattr(self,item,eval(item))
        self.stash(['Indices','IndicesGenB','IndicesGenWqq'])
        self.moreName = (tag if tag!=None else '') + '; ' + ','.join(samples)
    @property
    def name(self) : return self.__name

    def baseSamples(self) : return self.samples

    def setup(self,*_) :
        hists = self.fromCache(['merged'],['B','Q','N'], tag = self.tag)
        self.histsBQN = [hists['merged'][jetType] for jetType in ['B','Q','N']]
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
        missing = [s for s in self.samples if s not in [ss['name'] for ss in org.samples]]
        if missing: print self.name, "-- no such samples :\n", missing
        org.mergeSamples( targetSpec = {'name':'merged'}, sources = self.samples )

    def reportCache(self) :
        optStat = r.gStyle.GetOptStat()
        r.gStyle.SetOptStat(0)
        r.gROOT.ProcessLine(".L %s/cpp/tdrstyle.C"%whereami())
        r.setTDRStyle()
        self.setup(None)
        for hist in filter(None,self.histsBQN) :
            hist.Scale(1./hist.Integral(0,hist.GetNbinsX()+1))
        if None in self.histsBQN :
            print '%s.setup() failed'%self.name
            r.gStyle.SetOptStat(optStat)
            return
        fileName = '/'.join(self.outputFileName.split('/')[:-1]+[self.name]) + '.pdf'
        c = r.TCanvas()
        c.Print(fileName +'[')
        leg = r.TLegend(0.18,0.45,0.95,0.75)
        #leg.SetHeader("jet flavor")
        leg.SetFillColor(r.kWhite)
        leg.SetBorderSize(0)
        leg.SetTextFont(42)
        for h in self.histsBQN :
            h.Fill(h.GetBinCenter(1), h.GetBinContent(0))
            h.SetBinContent(0,0)
        height = 1.1 * max(h.GetMaximum() for h in self.histsBQN)
        for i,(f,color,style) in enumerate(zip('BQN',[r.kBlack,r.kRed,r.kBlue],[1,7,8])) :
            h = self.histsBQN[i]
            h.UseCurrentStyle()
            h.SetTitle(";%s;Probability / %.2f"%(h.GetXaxis().GetTitle().split()[0].replace('jet',''),(self.binning[2]-self.binning[1]) / self.binning[0]))
            h.SetLineColor(color)
            h.SetLineWidth(2 if style==1 else 3)
            h.SetLineStyle(style)
            h.SetMaximum(height)
            h.SetMinimum(0)
            h.Draw("hist" + ("same" if i else ""))
            leg.AddEntry(h,{"B":"Jets from b quark hadronization","Q":"Jets from W boson decay","N":"Other jets"}[f],'l')
        leg.Draw()
        r.gPad.RedrawAxis()

        stamp = r.TText()
        ssize = stamp.GetTextSize()
        stamp.SetTextFont(62)
        stamp.SetTextSize(ssize)
        stamp.DrawTextNDC(0.20 ,0.88,"CMS")
        stamp.SetTextFont(52)
        stamp.SetTextSize(0.8 * ssize)
        stamp.DrawTextNDC(0.20, 0.83, "Simulation")
        stamp.SetTextFont(42)
        stamp.SetTextSize(ssize)
        stamp.DrawTextNDC(0.84, 0.96, "(8 TeV)")


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
        r.gROOT.ProcessLine(".L %s/cpp/tdrstyle.C"%whereami())
        r.setTDRStyle()
        r.tdrStyle.SetPadRightMargin(0.06)
        self.setup(None)
        if not self.hists :
            print '%s.setup() failed'%self.name
            r.gStyle.SetOptStat(optStat)
            return
        fileName = '/'.join(self.outputFileName.split('/')[:-1]+[self.name]) + '.pdf'
        c = r.TCanvas()
        c.Print(fileName +'[')

        stamp = r.TText()
        ssize = stamp.GetTextSize()

        for f in 'BQN' :
            leg = r.TLegend(0.31,0.65,0.9,0.95)
            #leg.SetHeader("#eta range")
            leg.SetFillColor(r.kWhite)
            leg.SetBorderSize(0)
            leg.SetTextFont(42)
            for i,(color,style) in enumerate(zip([r.kBlack,r.kRed,r.kBlue],[1,7,8])) :
                h = self.hists[f+str(i)]
                label = h.GetTitle().split(',')[1].replace("0.000<","").replace("<"," < ")
                letter = h.GetTitle().split(',')[0]
                h.SetTitle(';#lower[0.2]{p_{#lower[-0.25]{T}}^{meas} (GeV)};#lower[-0.15]{Median log(E^{#lower[0.4]{gen}}/E^{#lower[0.4]{meas}})}')
                h.SetLineColor(color)
                h.SetLineStyle(style)
                h.SetLineWidth(3 if style!=1 else 2)
                h.SetMaximum(0.2)
                h.SetMinimum(-0.2)
                h.Draw("histsame" if i else "hist")
                stamp.SetTextFont(42)
                stamp.SetTextSize(ssize)
                stamp.DrawTextNDC(0.2,0.2, {'B':'Jets from b quark hadronization', 'Q':'Jets from W boson decay', 'N':'Other jets'}[letter])
                leg.AddEntry(h,label,'l')
            leg.Draw()
            r.gPad.RedrawAxis()

            stamp.SetTextFont(62)
            stamp.SetTextSize(ssize)
            stamp.DrawTextNDC(0.16 ,0.96,"CMS")
            stamp.SetTextFont(52)
            stamp.SetTextSize(0.8 * ssize)
            stamp.DrawTextNDC(0.27, 0.96, "Simulation")
            stamp.SetTextFont(42)
            stamp.SetTextSize(ssize)
            stamp.DrawTextNDC(0.8, 0.96, "(8 TeV)")

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
class CSV(calculables.secondary) :
    def __init__(self, collection = None, binning = (0,0,0), samples = [], tag = None, activated = False) :
        self.fixes = collection
        self.leaf = 'combinedSecondaryVertex'.join(collection)
        for item in ['binning','samples','tag','activated'] : setattr(self,item,eval(item))
        self.stash(['Indices','IndicesGenB'])
        self.moreName = ('CSV sf-adjusted') + '; ' + ','.join(samples)
        self.wpNames = ['zero','CSVL','CSVM','CSVT','one']
        self.workingPoints = [0.0, 0.244, 0.679, 0.898, 1.0]
        # SF reference: http://cds.cern.ch/record/1494669/files/1748-0221_8_04_P04013.pdf
        self.sfPointsBN = ([1.0, 1.01, 0.97, 0.96, 1.0], # Table 10, b-jet sf fot tt-like events
                           [1.0, 1.10, 1.11, 1.17, 1.0]) # Table 12, misidentification sf
        self.sfBN = [r.TGraph(5, np.array(self.workingPoints), np.array(points)) for points in self.sfPointsBN]

    def baseSamples(self) : return self.samples

    def setup(self,*_) :
        hists = self.fromCache(['merged'],['B','N'], tag = self.tag)
        self.histsBN = [hists['merged'][jetType] for jetType in ['B','N']]
        if not any(self.histsBN):
            self.activated = False
        if not self.activated: return
        for hist in filter(None,self.histsBN) : hist.Scale(1./hist.Integral(),"width")
        def makeCDF(h):
            cdf = h.Clone(h.GetName()+"_cdf")
            cdf.Reset()
            for i in range(1,2+h.GetNbinsX()):
                cdf.SetBinContent(i, h.Integral(1,i,"width"))
            return cdf
        self.cdfBN = map(makeCDF, self.histsBN)
        self.cdfBN_d = [h.Clone(h.GetName()+"_d") for h in self.cdfBN]
        for h,sf in zip(self.cdfBN_d, self.sfBN):
            for i in range(1,1+h.GetNbinsX()):
                P_MC = h.GetBinContent(i)
                SF = sf.Eval(h.GetBinLowEdge(i)+h.GetBinWidth(i))
                h.SetBinContent(i, 1-SF*(1-P_MC))

        def interpolate(x,h):
            i = h.FindFixBin(x)
            P = h.GetBinContent(i)
            P_ = h.GetBinContent(i-1)
            X = h.GetBinLowEdge(i+1)
            X_ = h.GetBinLowEdge(i)
            return P_ + (x-X_) * (P-P_) / (X-X_)
        points = np.arange(0,1+1e-8,0.005)
        cdfs = [np.array([interpolate(p,cdf) for p in points]) for cdf in self.cdfBN]
        self.testsBN = [r.TGraph(len(points),points, Ps) for Ps in cdfs]

        def inverse(p,h):
            i = h.FindFirstBinAbove(p)
            if i==-1: return 1
            P = h.GetBinContent(i)
            P_ = h.GetBinContent(i-1)
            X = h.GetBinLowEdge(i+1)
            X_ = h.GetBinLowEdge(i)
            return X_ + (X-X_) * (p-P_) / (P-P_)
        funcs = [np.array([inverse(interpolate(p,cdf),cdf_d) for p in points ]) for cdf,cdf_d in zip(self.cdfBN, self.cdfBN_d)]
        self.funcsBN = [r.TGraph(len(points),points,func) for func in funcs]
        
        
    def uponAcceptance(self,ev) :
        if ev['isRealData'] : return
        indices = ev[self.Indices]
        iB = ev[self.IndicesGenB]
        bvar = ev[self.leaf]
        for i in indices :
            jetType = "B" if i in iB else "N"
            self.book.fill(bvar.at(i), jetType, *self.binning, title = ";CSV (%s);events / bin"%(jetType))
    
    def update(self,_) :
        csv = self.source[self.leaf]
        if self.source['isRealData'] or not self.activated:
            self.value = csv
            return
        iB = self.source[self.IndicesGenB]
        self.value = utils.vector( [ self.funcsBN[0 if i in iB else 1].Eval(c) for i,c in enumerate(csv) ] )
        
    def organize(self,org) :
        if org.tag != self.tag : return
        missing = [s for s in self.samples if s not in [ss['name'] for ss in org.samples]]
        if missing: print self.name, "-- no such samples :\n", missing
        org.mergeSamples( targetSpec = {'name':'merged'}, sources = self.samples )

    def reportCache(self) :
        if not self.activated:
            print '%s: Not activated' % self.name
            return
        optStat = r.gStyle.GetOptStat()
        r.gStyle.SetOptStat(0)
        r.gROOT.ProcessLine(".L %s/cpp/tdrstyle.C"%whereami())
        r.setTDRStyle()
        self.setup(None)
        if None in self.histsBN :
            print '%s.setup() failed'%self.name
            r.gStyle.SetOptStat(optStat)
            return
        fileName = '/'.join(self.outputFileName.split('/')[:-1]+[self.name]) + '.pdf'
        c = r.TCanvas()
        c.Print(fileName +'[')
        histlist = [self.histsBN, self.cdfBN, self.cdfBN_d, self.sfBN, self.funcsBN]
        titles = [";CSV';pdf",";CSV';cdf",";CSV;cdf",";CSV Working Point;Scale Factor",";CSV';CSV\""]
        for hists,title in zip(histlist, titles):
            leg = r.TLegend(0.36,0.55,0.6,0.7)
            leg.SetHeader("jet flavor")
            height = 1.1 * max(h.GetMaximum() for h in hists)
            if height<0: height = 1.1 * max(sum(self.sfPointsBN,[]))
            for i,(f,color) in enumerate(zip('BN',[r.kRed,r.kBlue,r.kGreen])) :
                h = hists[i]
                h.SetTitle(title)
                h.UseCurrentStyle()
                h.SetLineColor(color)
                h.SetLineWidth(2)
                h.SetMarkerColor(color)
                h.SetMaximum(height)
                h.SetMinimum(0)
                h.Draw("hist" + ("same" if i else "")) if type(h)!=r.TGraph else h.Draw("APL"[1 if i else None:])
                leg.AddEntry(h,{"B":"B jets","N":"other jets"}[f],'l')
            leg.Draw()
            c.Print(fileName)
        c.Print(fileName +']')
        print 'Wrote : %s'%fileName
        r.gStyle.SetOptStat(optStat)
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
