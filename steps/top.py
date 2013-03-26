import math,itertools,ROOT as r
from supy import analysisStep,steps,calculables
try:
    import numpy as np
    import scipy.stats
except:  pass

#####################################

class channelClassification(analysisStep) :
    labels = ['','ee','mm','tt','em','et','mt','ej','mj','tj','jj']
    nbins = len(labels)

    def uponAcceptance(self,ev) :
        iBin = self.labels.index(ev['ttDecayMode'])
        if iBin :
            self.book.fill( iBin, 'ttDecayMode', self.nbins, 0, self.nbins, xAxisLabels = self.labels )

#####################################
class kinFitLook(analysisStep) :
    def __init__(self,indexName) : self.moreName = indexName
    def uponAcceptance(self,ev) :
        index = ev[self.moreName]
        if index<0 : return
        topReco = ev["TopReconstruction"][index]
        residuals = topReco["residuals"]
        lepTopM = topReco['lepTopP4'].M()
        hadTopM = topReco['hadTopP4'].M()

        for name,val in residuals.iteritems() :
            self.book.fill(val, "residual_%s"%name+self.moreName, 50,-7,7, title = ";residual %s;events / bin"%name)
            self.book.fill(scipy.stats.norm.cdf(val), "residualCDF_%s"%name+self.moreName, 100,0,1, title = ";cdf(residual) %s;events / bin"%name)

        self.book.fill( (lepTopM, hadTopM), "topMassesFit"+self.moreName, (100,100),(0,0),(300,300),
                        title = ";fit mass_{top} (leptonic); fit mass_{top} (hadronic);events / bin",)
        
        self.book.fill( topReco['chi2'], "topReco_Chi2"+self.moreName, 50, 0 , 100, title = ';ttbar kin. fit #chi^{2};events / bin')
#####################################

class combinatorialFrequency(analysisStep) :
    def uponAcceptance(self,ev) :
        if not ev['genTopTTbar'] : return
        recos = ev['TopReconstruction']
        iP,iQ,iH,iL = ev['IndicesGenTopPQHL']

        igen=ev['genTopRecoIndex']

        iGenHadB = next((i for i,R in enumerate(recos) if R['iPQHL'][2]==iH) , -1)
        iGenLepB = next((i for i,R in enumerate(recos) if R['iPQHL'][3]==iL) , -1)
        iGenHadPQB  = next((i for i,R in enumerate(recos) if R['iPQHL'][:3]==(iP,iQ,iH)), -1 )
        iGenJets = next((i for i,R in enumerate(recos) if R['iPQHL']==(iP,iQ,iH,iL)), -1)

        for i,item in enumerate(['iGenHadB','iGenLepB','iGenHadPQB','iGenJets','igen']) :
            self.book.fill(eval(item), str(i)+item, 10,-1.5,8.5, title=";%s;events / bin"%item )
#####################################
class kinematics(analysisStep) :
    def __init__(self,indexName) : self.moreName = indexName
    def uponAcceptance(self,ev) :
        index = ev["%sRecoIndex"%self.moreName]
        if index < 0 : return
        topReco = ev["TopReconstruction"]
        self.book.fill(ev[self.moreName+"TtxMass"], "TTX.mass", 50,300,1300, title = ";ttx invariant mass;events / bin")
        self.book.fill(ev[self.moreName+"PtSum"],   "TT.pt",   100,  0, 200, title = ";ttbar.pt;events / bin")
        self.book.fill(ev[self.moreName+"RapiditySum"], "TT.y", 50, 0, 3, title = ";|t#bar{t}.y|;events / bin" )
        self.book.fill( topReco[index]['hadTraw'].mass(), "rawHadTopMass", 100, 100,300, title = ";%s raw hadronic top mass;events / bin"%self.moreName)
        self.book.fill( topReco[index]['hadWraw'].mass(), "rawHadWMass", 100, 0,200, title = ";%s raw hadronic W mass;events / bin"%self.moreName)

        for low,var in  zip([0,0,-1],["fitTopPtOverSumPt", "fitTopTanhRapiditySum", "TridiscriminantQQggQg"]) :
            lo = low
            up = 1
            v = max(lo,min(ev[var],up-1e-6))
            self.book.fill( v, var, 100, lo, up, title = ';%s;events / bin'%var )

        i5 = ev['fitTopFifthJetIndex']
        triD = ev['TridiscriminantWTopQCD']
        if i5!=None:
            jets = ev['TopJets']
            vars = [item.join(jets) for item in ['Pt','Moments2Sum']]
            los = (0,0)
            ups = (400,0.11)
            for var,lo,up in zip(vars,los,ups) :
                v = max(lo, min(up-1e-8, ev[var][i5]))
                self.book.fill( (v,triD), '%s_triD'%var, (100,100), (lo,-1), (up,1))
                self.book.fill( v, var, 100, lo, up, title = ';%s;events / bin'%var)
#####################################
class kinematics3D(analysisStep) :
    def __init__(self,indexName) : self.moreName = indexName
    def uponAcceptance(self,ev) :
        index = ev["%sRecoIndex"%self.moreName]
        if index < 0 : return
        topReco = ev["TopReconstruction"]

        for low,var in  zip([0,0,-1],["fitTopPtOverSumPt", "fitTopTanhRapiditySum", "TridiscriminantQQggQg"]) :
            lo = (low,-1)
            up = (1, 1)
            v = tuple( max(L,min(val,U-1e-6)) for val,L,U  in zip((ev[var],ev['TridiscriminantWTopQCD']),lo, up) )
            self.book.fill( v, '%s_triD'%var, (50,5), lo, up)
#####################################
class resolutions(analysisStep) :
    def __init__(self,indexName) : self.moreName = indexName
    def uponAcceptance(self,ev) :
        genTTbar = ev["genTopTTbar"]
        if not genTTbar : return

        topReco = ev["TopReconstruction"]

        index = ev[self.moreName]
        self.book.fill(index, '0_'+self.moreName, 21, -1.5, 19.5, title = ';%s;events / bin'%self.moreName)
        if index<0 : return
        
        if ev['genTTbarIndices']['semi'] :
            for s in ['lep','nu','bLep','bHad','q'] :
                self.book.fill(ev['%sDeltaRTopRecoGen'%s][index], '1_'+s+'DeltaRTopRecoGen', 50,0,2, title = ';%s DeltaR reco gen;events / bin'%s)
        
        gsP4 = ev['genSumP4']
        self.book.fill( topReco[index]['ttx'].pt() - gsP4.pt(), "resolution_pt", 100, -100, 100, title = ";%s #Delta_{reco-gen} ttx.pt;events / bin"%self.moreName )
        self.book.fill( topReco[index]['ttx'].mass() - gsP4.mass(), "resolution_m", 100, -200, 200, title = ";%s #Delta_{reco-gen} ttx.m;events / bin"%self.moreName )

        iLep = min(0,topReco[index]["lepCharge"])
        gen =  (ev["genP4"][genTTbar[0]], ev["genP4"][genTTbar[1]])
        reco = (topReco[index]['top'],topReco[index]['tbar'])
        unfit = (topReco[index]['lepTraw'], topReco[index]['hadTraw'])[::topReco[index]["lepCharge"]]

        for func in ['Rapidity'] :
            genFunc = (getattr(gen[0],func)(), getattr(gen[1],func)())
            recoFunc = (getattr(reco[0],func)(),getattr(reco[1],func)())
            unfitFunc = (getattr(unfit[0],func)(), getattr(unfit[1],func)())
            for f,fit in [('fit',recoFunc),('unfit',unfitFunc)] :
                self.book.fill( fit[iLep]     - genFunc[iLep],       "d%sLepTop_%s"%(func,f), 100,-1,1, title=";lep top #Delta %s_{%s reco gen};events / bin"%(func,f))
                self.book.fill( fit[not iLep] - genFunc[not iLep],   "d%sHadTop_%s"%(func,f), 100,-1,1, title=";had top #Delta %s_{%s reco gen};events / bin"%(func,f))
                self.book.fill( fit[0]-fit[1] - (genFunc[0]-genFunc[1]), "dd%sTTbar_%s"%(func,f), 100,-1,1, title = ";#Delta %s_{t#bar{t} %s reco}-#Delta %s_{t#bar{t} gen};events / bin"%(func,f,func))
