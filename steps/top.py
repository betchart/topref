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
class leptonSigned(analysisStep) :
    def __init__(self, var, binning = (100,0,1)):
        self.moreName = "+/- lepton of %s"%var
        self.var = var
        self.binning = binning

    def uponAcceptance(self,ev) :
        q = 'Pos' if ev['fitTopLeptonCharge']>0 else 'Neg'
        self.book.fill(ev[self.var], self.var+q, *self.binning, title = ';%s (%s lepton);events / bin'%(self.var,q))
#####################################
class Asymmetry(analysisStep) :
    def __init__(self, collection, bins = 18 ) :
        self.collection = collection
        for item in ["LeptonCharge","SignedLeptonRapidity","RelativeLeptonRapidity",
                     "DeltaAbsYttbar","DirectedDeltaYttbar","Beta","DirectedDeltaYLHadt",
                     "DirectedLTopRapidity","DirectedHTopRapidity"] :
            setattr(self,item,("%s"+item+"%s")%collection)
        self.bins = bins
        self.moreName = "with %d bins."%bins

    def uponAcceptance(self,ev) :
        for charge in ["",["Negative","Positive"][max(0,ev[self.LeptonCharge])]][:1] :
            self.book.fill(ev[self.SignedLeptonRapidity], "leptonSignedY"+charge, self.bins,-3,3, title = "%s;leptonSignedY;events / bin"%charge)
            self.book.fill(ev[self.RelativeLeptonRapidity], "leptonRelativeY"+charge, self.bins,-3,3, title = "%s;#Delta y;events / bin"%charge)
            #self.book.fill(ev[self.DirectedLTopRapidity], "dirLtopY"+charge, self.bins,-3,3, title = "%s;y_{ltop};events / bin"%charge)
            #self.book.fill(ev[self.DirectedHTopRapidity], "dirHtopY"+charge, self.bins,-3,3, title = "%s;y_{htop};events / bin"%charge)

        self.book.fill( ev[self.DeltaAbsYttbar],      'ttbarDeltaAbsY',    self.bins, -3, 3, title = ';#Delta|Y|_{ttbar};events / bin' )
        self.book.fill( ev[self.DirectedDeltaYttbar], 'ttbarSignedDeltaY', self.bins, -4, 4, title = ';sumP4dir * #Delta Y_{ttbar};events / bin' )
        self.book.fill( ev[self.DirectedDeltaYLHadt], 'lHadtDeltaY',       self.bins, -4, 4, title = ';sumP4dir * #Delta Y_{lhadt};events / bin')
        #self.book.fill( ev[self.Beta],                'ttbarBeta',         self.bins, -math.sqrt(2), math.sqrt(2), title = ';#beta_{ttbar};events / bin')
        self.book.fill( ev["TTbarSignExpectation"], 'ttbarSignExpectation', self.bins, -1, 1, title = ";<sign #Delta y>_{t#bar{t}};events / bin" )
#####################################
class signCheck(analysisStep) :
    def uponAcceptance(self,ev) :
        if ev["isRealData"] : return
        if not ev["genQQbar"] : return
        if not ev["genTopTTbar"] : return
        genTT_DDY = ev["genTopDirectedDeltaYttbar"]
        signGenTT_DDY = 1 if genTT_DDY > 0 else -1 if genTT_DDY<0 else 0
        self.book.fill( ev["TTbarSignExpectation"] * signGenTT_DDY, "TTbarSignExpectation_trueSign", 41, -1, 1, title = ";<sign #Delta y>_{t#bar{t}} . trueSign;events / bin")
        self.book.fill( ev["TTbarSignExpectation"] * signGenTT_DDY, "TTbarSignExpectation_trueSign_2bin", 2, -1, 1, title = ";<sign #Delta y>_{t#bar{t}} . trueSign;events / bin")

#####################################
class Spin(analysisStep) :
    def __init__(self, collection) :
        self.collection = collection
        for item in ['CosHelicityThetaL', 'CosHelicityThetaQ'] :
            setattr(self,item,('%s'+item+'%s')%collection)
        self.bins = 18

    def uponAcceptance(self,ev) :
        cosTL = ev[self.CosHelicityThetaL]
        cosTQ = ev[self.CosHelicityThetaQ]
        self.book.fill( cosTL, 'CosHelicityThetaL', self.bins, -1, 1, title = ';CosHelicityThetaL;events / bin' )
        self.book.fill( cosTQ, 'CosHelicityThetaQ', self.bins, -1, 1, title = ';CosHelicityThetaQ;events / bin' )
        self.book.fill( cosTL*cosTQ, 'helicityCos2', self.bins, -1, 1, title = ';helicityCos2;events / bin' )
        self.book.fill( (cosTL, cosTQ), 'vs_CosHelicityThetaL_CosHelicityThetaQ', (self.bins,self.bins), (-1,-1), (1,1),
                        title = ';CosHelicityThetaL;CosHelicityThetaQ;events / bin')
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
        #lepWM = topReco['lepW'].M()
        #hadWM = topReco['hadW'].M()
        #rawHadWM = topReco['hadWraw'].M()

        for name,val in residuals.iteritems() :
            self.book.fill(val, "residual_%s"%name+self.moreName, 50,-7,7, title = ";residual %s;events / bin"%name)
            self.book.fill(scipy.stats.norm.cdf(val), "residualCDF_%s"%name+self.moreName, 100,0,1, title = ";cdf(residual) %s;events / bin"%name)

        #self.book.fill( (topReco["dS"],topReco["dL"]), "topKinFit_DSoverDL"+self.moreName, (100,100), (0,0), (30,30), title = ";ds;dL;events / bin")
        #self.book.fill( (topReco["sigmaS"],topReco["sigmaL"]), "topKinFit_SigmaSoverSigmaL"+self.moreName, (100,100), (0,0), (30,30), title = ";#sigma_{s};#sigma_{L};events / bin")
        #self.book.fill((residuals["hadP"],residuals["hadQ"]), "topKinFit_residual_had_PQ"+self.moreName, (100,100),(-5,-5),(5,5), title = ';residual hadP;residual hadQ;events / bin')
        #self.book.fill((residuals["lepS"],residuals["lepL"]), "topKinFit_residual_lep_SL"+self.moreName, (100,100),(-5,-5),(5,5), title = ';residual lepS;residual lepL;events / bin')

        #self.book.fill( lepWM, "wMassLepFit"+self.moreName, 60, 0, 180, title = ';fit mass_{W} (leptonic);events / bin')
        #self.book.fill( hadWM, "wMassHadFit"+self.moreName, 60, 0, 180, title = ';fit mass_{W} (hadronic);events / bin')
        #self.book.fill( rawHadWM, "wMassHadRaw"+self.moreName, 60, 0, 180, title = ';raw mass_{W} (hadronic);events / bin')
        #self.book.fill( lepTopM, "topMassLepFit"+self.moreName, 100,0,300, title = ";fit mass_{top} (leptonic);events / bin" )
        #self.book.fill( hadTopM, "topMassHadFit"+self.moreName, 100,0,300, title = ";fit mass_{top} (hadronic);events / bin" )
        self.book.fill( (lepTopM, hadTopM), "topMassesFit"+self.moreName, (100,100),(0,0),(300,300),
                        title = ";fit mass_{top} (leptonic); fit mass_{top} (hadronic);events / bin",)
        
        self.book.fill( topReco['chi2'], "topReco_Chi2"+self.moreName, 50, 0 , 100, title = ';ttbar kin. fit #chi^{2};events / bin')
        #self.book.fill( math.exp(-0.5*topReco['chi2']), "topReco_L"+self.moreName, 50, 0, 1, title = ";ttbar kin. fit exp(-0.5#chi^{2});events / bin" )
        #self.book.fill( math.log(1+topReco['chi2']), "topRecoLogOnePlusChi2"+self.moreName, 50, 0 , 7, title = ';ttbar kin. fit log(1+#chi^{2});events / bin')
        #self.book.fill( math.log(1+topReco['key']), "topRecoLogOnePlusKey"+self.moreName, 50, 0 , 7, title = ';ttbar kin. fit log(1+#chi^{2}-2logP);events / bin')

        #hadX2 = math.log(1+topReco['hadChi2'])
        #lepX2 = math.log(1+topReco['lepChi2'])
        #bound = ("_bound" if topReco['lepBound'] else "_unbound")

        #self.book.fill( hadX2, "topRecoLHadX2"+self.moreName, 50, 0 , 10, title = ';ttbar kin. fit log(1+#chi^{2}_{had});events / bin')
        #self.book.fill( lepX2, "topRecoLLepX2"+self.moreName, 50, 0 , 10, title = ';ttbar kin. fit log(1+#chi^{2}_{lep});events / bin')
        #self.book.fill( lepX2, "topRecoLLepX2"+bound+self.moreName, 50, 0 , 10, title = ';ttbar kin. fit log(1+#chi^{2}_{lep});events / bin')
        #self.book.fill( (lepX2,hadX2), "topRecoVsLX2"+self.moreName, (50,50),(0,0),(10,10), title = ";log(1+#chi^{2}_{lep});log(1+#chi^{2}_{had});events / bin" )
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
        self.book.fill(ev[self.moreName+"RapiditySum"], "TT.y", 50, 0, 3, title = ";ttbar.rapidity;events / bin" )
        self.book.fill( topReco[index]['hadTraw'].mass(), "rawHadTopMass", 100, 100,300, title = ";%s raw hadronic top mass;events / bin"%self.moreName)
        self.book.fill( topReco[index]['hadWraw'].mass(), "rawHadWMass", 100, 0,200, title = ";%s raw hadronic W mass;events / bin"%self.moreName)

        lo = (-1,0)
        up = (1,1)
        v = tuple( max(L,min(val,U-1e-6)) for val,L,U  in zip((ev['TridiscriminantWTopQCD'],
                                                               ev[self.moreName+'SqrtPtOverSumPt']), lo, up) )
        self.book.fill( v, 'triD_v_sqtsumptopt', (100,100), lo, up, title = ""  )
#####################################
class resolutions(analysisStep) :
    def __init__(self,indexName) : self.moreName = indexName
    def uponAcceptance(self,ev) :
        genTTbar = ev["genTopTTbar"]
        if not genTTbar : return

        topReco = ev["TopReconstruction"]

        index = ev[self.moreName]
        self.book.fill(index, self.moreName, 21, -1.5, 19.5, title = ';%s;events / bin'%self.moreName)
        if index<0 : return
        
        if ev['genTTbarIndices']['semi'] :
            for s in ['lep','nu','bLep','bHad','q'] :
                self.book.fill(ev['%sDeltaRTopRecoGen'%s][index], s+'DeltaRTopRecoGen', 50,0,2, title = ';%s DeltaR reco gen;events / bin'%s)
        
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

        #iHad = max(0,topReco[index]["lepCharge"])
        #genLepY = ev['genP4'][max(ev['genTTbarIndices'][item] for item in ['lplus','lminus'])].Rapidity()
        #self.book.fill( recoY[iHad] - topReco[index]['lep'].Rapidity() - (genY[iHad]-genLepY), "ddRapidityLHadTop", 100,-1,1, title = ";#Delta y_{l-htop reco}-#Delta y_{l-htop gen};events / bin")
######################################
class mcTruthQDir(analysisStep) :
    def __init__(self,withLepton = False, withNu = False) :
        self.withNu = withNu and withLepton
        self.withLepton = withLepton
        
    def uponAcceptance(self,ev) :
        if ev['isRealData'] : return
        genSumPz = ev['genSumP4'].pz()
        #for sumP4 in ['genTopNuP4','genTopTTbarSumP4','mixedSumP4','mixedSumP4Nu'][:4 if self.withNu else 3 if self.withLepton else 2] :
        #    self.book.fill( (genSumPz, ev[sumP4].pz()), "genSumP4_%s_pz"%sumP4, (100,100),(-3000,-3000),(3000,3000),
        #                    title = ";genSumP4 pz;%s pz;events/bin"%sumP4)

        qqbar = ev['genQQbar']
        if qqbar :
            qdir = 1 if ev['genP4'][qqbar[0]].pz()>0 else -1
            for sumP4 in ['genSumP4','genTopSumP4','mixedSumP4','mixedSumP4Nu','fitTopSumP4'][:5 if self.withNu else 3 if self.withLepton else 2] :
                self.book.fill( qdir * ev[sumP4].pz(), "qdir_%s_pz"%sumP4, 100,-3000,3000, title = ';qdir * %s.pz;events/bin'%sumP4)
                self.book.fill( qdir * ev[sumP4].Eta(), "qdir_%s_eta"%sumP4, 100,-10,10, title = ';qdir * %s.eta;events/bin'%sumP4)
        
#####################################
class mcTruthAcceptance(analysisStep) :
    def uponAcceptance(self,ev) :
        if not ev['genTopTTbar'] : return

        indices = ev['genTTbarIndices']
        if not bool(indices['lplus'])^bool(indices['lminus']) : return
        lep = ev['genP4'][max(indices['lplus'],indices['lminus'])]
        iJets = [indices['b'],indices['bbar']] + indices['wplusChild' if indices['lminus'] else 'wminusChild']
        jets = [ev['genP4'][i] for i in iJets]
        iBlep = indices['b'] if indices['lplus'] else indices['bbar']
        
        self.book.fill(lep.eta(),"lepEta",31,-5,5, title=';#eta_{lep};events / bin')
        self.book.fill(max([abs(p4.eta()) for p4 in jets]), 'jetEtaMax', 30,0,5, title=';jet |#eta|_{max};events / bin')
        self.book.fill(max([abs(p4.eta()) for p4 in jets[:2]]), 'jetEtaMaxB', 30,0,5, title=';b jet |#eta|_{max};events / bin')
        self.book.fill(max([abs(p4.eta()) for p4 in jets[2:]]), 'jetEtaMaxLite', 30,0,5, title=';lite jet |#eta|_{max};events / bin')

        pts = [p4.pt() for p4 in jets]
        self.book.fill(min(pts), 'jetMinPt', 50,0,100, title=';jet pT_{min};events / bin')
        self.book.fill(min(pts[:2]), 'jetMinPtB', 50,0,100, title=';b jet pT_{min};events / bin')
        self.book.fill(min(pts[2:]), 'jetMinPtLite', 50,0,100, title=';lite jet pT_{min};events / bin')

        self.book.fill( max(pts[:2]) - min(pts[2:]), "diffBigBLittleQ", 50,-50,100,title=';pT_{maxb}-pT_{minq};events / bin' )
        self.book.fill( min(pts[:2]) - max(pts[2:]), "diffLittleBBigQ", 50,-50,100,title=';pT_{minb}-pT_{maxq};events / bin' )
        self.book.fill( sum(pts[:2]) - sum(pts[2:]), "diffSumBBSumQQ", 50,-50,100,title=';sumpT_{b}-sumpT_{q};events / bin' )
        
        self.book.fill(sum(pts), 'jetSumPt', 50, 0, 800, title=';#sum pT_{top jets};events / bin')
        self.book.fill(sum(pts)-ev['genP4'][iBlep].pt(), 'jetSumPtHad', 50, 0, 500, title=';hadronic #sum pT_{top jets};events / bin')

        self.book.fill( int(max(pts)==max(pts[:2])), "maxPtJetIsBjet", 2, 0 , 1, title = ';maxPt is bjet;events / bin')
        self.book.fill( int(max(pts[:2])>min(pts[2:])), "maxPtOrNextJetIsBjet", 2, 0 , 1, title = ';maxPt or next is bjet;events / bin')
        self.book.fill( int(sum(pts[:2])>sum(pts[2:])), "sumPtBB_gt_sumPtPQ", 2, 0 , 1, title = ';sumPt of bs > sumPt of pq;events / bin')
#####################################
class mcTruthTemplates(analysisStep) :
    def uponAcceptance(self,ev) :
        if not ev['genTopTTbar'] : return

        self.book.fill(ev['genTopAlpha'],'alpha',10,0,1,title=';genTopAlpha;events / bin')
        self.book.fill(math.sqrt(ev['genTopAlpha']),'alpha_sqrt',10,0,1,title=';#sqrt{#alpha};events / bin')
        alpha = '_alpha%02d'%int(10*ev['genTopAlpha'])
        self.book.fill(ev['genTopAlpha'], "alpha%s"%alpha, 100,0,1, title = ';#alpha;events / bin')

        cts,ctsb = ev['genttCosThetaStar']
        self.book.fill(cts, 'genCosT', 20, -1, 1, title = ';gen cosThetaStar;events / bin')
        self.book.fill(cts, 'genCosT%s'%alpha, 20, -1, 1, title = ';gen cosThetaStar;events / bin')
        self.book.fill(ctsb, 'genCosTbar', 20, -1, 1, title = ';gen cosThetaStarBar;events / bin')
        self.book.fill(ctsb, 'genCosTbar%s'%alpha, 20, -1, 1, title = ';gen cosThetaStarBar;events / bin')
        self.book.fill(0.5*(cts+ctsb), 'genCosTavg', 20, -1, 1, title = ';gen cosThetaStarAvg;events / bin')
        self.book.fill(0.5*(cts+ctsb), 'genCosTavg%s'%alpha, 20, -1, 1, title = ';gen cosThetaStarAvg;events / bin')

        return
        
        #self.book.fill(ev['genTopTTbarSumP4'].M(), "genttbarinvmass", 40,0,1000, title = ';ttbar invariant mass;events / bin' )
        #for i in [0,1]: self.book.fill(ev['genP4'][ev['genTopTTbar'][i]].M(), "topmass", 50, 120, 220, title = ';top mass;events / bin')

        qqbar = ev['genQQbar']
        genP4 = ev['genP4']
        qdir = 1 if qqbar and genP4[qqbar[0]].pz()>0 else -1
        genP4dir = 1 if ev['genSumP4'].pz() > 0 else -1
        
        self.book.fill(    qdir * ev['genTopDeltaYttbar'], 'genTopTrueDeltaYttbar', 31,-5,5, title = ';True Signed #Delta y_{ttbar};events / bin')
        self.book.fill(genP4dir * ev['genTopDeltaYttbar'], 'genTopMezDeltaYttbar', 31,-5,5, title = ';MEZ Signed #Delta y_{ttbar};events / bin')
        self.book.fill(        ev['genTopDeltaAbsYttbar'], 'genTopDeltaAbsYttbar', 31,-5,5, title = ';#Delta |y|_{ttbar};events / bin')

        indices = ev['genTTbarIndices']
        if indices['lplus'] and indices['lminus'] :
            dy = genP4[indices['lplus']].Rapidity() - genP4[indices['lminus']].Rapidity()
            self.book.fill(    qdir * dy, "genTopTrueDeltaYll", 31,-5,5, title = ';True Signed #Delta y_{ll};events / bin')
            self.book.fill(genP4dir * dy, "genTopMezDeltaYll", 31,-5,5, title = ';MEZ Signed #Delta y_{ll};events / bin')
        elif indices['lplus'] or indices['lminus'] :
            Q = 1 if indices['lplus'] else -1
            lRapidity = genP4[max(indices['lplus'],indices['lminus'])].Rapidity()
            dy = (lRapidity - ev['genSumP4'].Rapidity())
            for suf in ['','Positive' if Q>0 else 'Negative'] :
                self.book.fill(    qdir * Q * dy, "genTopTrueDeltaYlmiss"+suf, 31,-5,5, title = '%s;True Signed #Delta y_{lmiss};events / bin'%suf)
                self.book.fill(genP4dir * Q * dy, "genTopMezDeltaYlmiss"+suf, 31,-5,5, title = '%s;MEZ Signed #Delta y_{lmiss};events / bin'%suf)
                self.book.fill(    qdir * Q * lRapidity, "genTopTrueLRapidity"+suf, 31,-5,5, title = "%s;True Signed y_l;events / bin"%suf)
                self.book.fill(genP4dir * Q * lRapidity, "genTopMezLRapidity"+suf, 31,-5,5, title = "%s;MEZ Signed y_l;events / bin"%suf)
######################
class mcTruthAsymmetryBinned(analysisStep) :
    def __init__(self, binVar, bins, min, max, collection = ("genTop","")) :
        for item in ['bins', 'min', 'max'] : setattr(self,item,eval(item))
        self.asymmVar = "%sDeltaY%s"%collection
        self.binVar = ("%s"+binVar+"%s")%collection
        self.binName = "%s_%s"%(self.asymmVar, self.binVar) + "%03d"
        
    def uponAcceptance(self,ev) :
        qqbar = ev['genQQbar']
        genP4 = ev['genP4']
        qdir = 1 if qqbar and genP4[qqbar[0]].pz()>0 else -1

        binVar = ev[self.binVar]
        Dy = ev[self.asymmVar] * qdir
        self.book.fill(binVar, self.binVar, self.bins, self.min, self.max, title = ';%s;events / bin'%self.binVar )
        bin = min(self.book[self.binVar].FindFixBin(binVar),self.bins)
        self.book.fill(Dy, self.binName%bin, 2, -50, 50, title = ";%s %d;events / bin"%(self.asymmVar,bin))

    def outputSuffix(self) : return steps.master.outputSuffix()

    def varsToPickle(self) :
        return ["bins","min","max","binName","asymmVar","binVar"]

    @staticmethod
    def asymmetryFromHist(hist) :
        if not hist : return 0,0
        nMinus = hist.GetBinContent(1)
        nMinusE = hist.GetBinError(1)
        nPlus = hist.GetBinContent(2)
        nPlusE = hist.GetBinError(2)
        S = nPlus + nMinus
        asymm = float(nPlus - nMinus) / S
        err = 2./S**2 * math.sqrt((nPlus*nMinusE)**2+(nMinus*nPlusE)**2)
        return asymm,err

    def mergeFunc(self, products) :
        file = r.TFile.Open(self.outputFileName, "UPDATE")
        master = file.FindObjectAny("Master")
        asymm = [self.asymmetryFromHist(master.FindObjectAny(self.binName%(bin+1))) for bin in range(self.bins) ]
        binVarHist = master.FindObjectAny(self.binVar)
        binVarHist.GetDirectory().cd()

        asymmByBinVar = binVarHist.Clone("%s_%s"%(self.binVar,self.asymmVar))
        asymmByBinVar.SetTitle(";%s;%s"%(self.binVar,"A_{fb}"))
        asymmByBinVar.SetMinimum(-0.5)
        asymmByBinVar.SetMaximum(0.5)
        
        for i in range(self.bins) :
            print asymm[i]
            asymmByBinVar.SetBinContent(i+1,asymm[i][0])
            asymmByBinVar.SetBinError(i+1,asymm[i][1])
        asymmByBinVar.SetBinContent(self.bins+1,0)
        asymmByBinVar.SetBinError(self.bins+1,0)
        asymmByBinVar.Write()
        r.gROOT.cd()
        file.Close()
        #print "Output updated with %s."%asymmByBinVar.GetName()
######################
class collisionType(calculables.secondary) :
    def uponAcceptance(self,ev) :
        id = ev['genPdgId']
        iHard = ev['genIndicesHardPartons']
        iQs = [i for i in iHard if id[i]!=21]
        nGlu = 2 - len(iQs)
        self.book.fill(nGlu, 'nCollidingGluons', 3,-0.5,2.5, title = ';N colliding gluons;events / bin')
        glubar = nGlu==1 and filter(lambda i: id[i]<0,iHard)
        self.book.fill(ev['genCosThetaStar'],'cosThetaStar_%dglu'%nGlu + ('bar' if glubar else ''), 100,-1,1,
                       title = ";cos#theta* (%d glu)%s;events / bin"%(nGlu,'(bar)' if glubar else ''))
        self.book.fill(ev['genCosThetaStarBar'],'cosThetaStarBar_%dglu'%nGlu + ('bar' if glubar else ''), 100,-1,1,
                       title = ";cos#bar{#theta}* (%d glu)%s;events / bin"%(nGlu,'(bar)' if glubar else ''))
        p4 = ev['genP4']
        system = p4[iHard[0]]+p4[iHard[1]]
        self.book.fill( system.E(), 'sqrtshat_%dglu'%nGlu, 100,0,2000, title = ';#sqrt{#hat{s}} (%d colliding gluons);events / bin'%nGlu )

        iQ = max(iQs, key = lambda i : id[i]) if nGlu!=2 else iHard[0]
        self.book.fill( system.z() * (1 if p4[iQ].z() > 0 else -1) * (-1 if glubar else 1), 'systemZbyQdir_%dglu'%nGlu, 100,-3000,3000,
                        title = ";(p_{z} system) * sign(p_{z} quark) , %d colliding gluons;events / bin"%nGlu )

    def organize(self,org) : org.scale(5000)

    def reportCache(self) :
        fileName = '/'.join(self.outputFileName.split('/')[:-1]+[self.name])
        optstat = r.gStyle.GetOptStat()
        r.gStyle.SetOptStat(0)
        c = r.TCanvas()
        c.Print(fileName+'.pdf[')

        samples = ['ttj_mg','ttj_mn','ttj_ph']
        colors = [r.kBlack,r.kBlue,r.kRed]
        names = [(0,''),(1,''),(1,'bar'),(2,'')]
        for name in names :
            pattern = 'cosThetaStar%s'+'_%dglu%s'%name
            hists = self.fromCache( samples, [pattern%'',pattern%'Bar'])
            toDraw = []
            for s in samples :
                hdiff = hists[s][pattern%''].Clone('diff')
                hsum = hists[s][pattern%''].Clone('sum')
                hdiff.Add(hists[s][pattern%'Bar'],-1)
                hsum.Add(hists[s][pattern%'Bar'])
                hdiff.SetLineColor(colors[samples.index(s)])
                hdiff.SetTitle('N_{t}-N_{#bar{t}}')
                hperc = hdiff.Clone('perc')
                hperc.Divide(hsum)
                hperc.GetYaxis().SetTitle('\%')
                hperc.SetTitle('(N_{t}-N_{#bar{t}}) / (N_{t}+N_{#bar{t}})')
                toDraw+=[(hdiff,hperc)]
            for i,(h,_) in enumerate(sorted(toDraw, key = lambda hh : hh[0].GetMaximum(), reverse = True)) : h.Draw('hist' + ('same' if i else ''))
            c.Print(fileName+'.pdf')
            for i,(_,h) in enumerate(sorted(toDraw, key = lambda hh : hh[1].GetMaximum(), reverse = True)) : h.Draw('hist' + ('same' if i else ''))
            c.Print(fileName+'.pdf')
        c.Print(fileName+'.pdf]')
        r.gStyle.SetOptStat(optstat)
        print "Wrote file: %s.pdf"%fileName


class mcQuestions(analysisStep) :

    def uponAcceptance(self,ev) :
        p4 = ev['genP4']
        id = ev['genPdgId']
        status = ev['genStatus']
        iQ,iQbar = ev['genQQbar']
        iT,iTbar = ev['genTopTTbar']
        qqbar = p4[iQ] + p4[iQbar]
        ttbar = p4[iT] + p4[iTbar]
        qqcmzboost = r.Math.BoostZ( qqbar.BoostToCM().z() )
        beta = qqbar.BoostToCM()
        qqcmboost = r.Math.Boost( beta.x(), beta.y(), beta.z() )
        glus = [i for i in range(iT,len(id)) if id[i]==21 and p4[i].pt() and status[i]==3]
        nglu = len(glus)
        qqbarMttbar = qqbar - ttbar

        self.book.fill(nglu, 'N gluons', 5,-0.5,4.5, title = 'N gluons')
        self.book.fill(qqbar.pt() , 'qqbarpt%d'%nglu, 100,0,50, title = ';q#bar{q} pt, nGlu=%d;'%nglu)
        self.book.fill(ttbar.pt() , 'ttbarpt%d'%nglu, 100,0,50, title = ';t#bar{t} pt, nGlu=%d;'%nglu)
        self.book.fill(sum([p4[i] for i in glus],-qqbarMttbar).pt(), 'qqbarMttbarglus%d'%nglu, 10,0,20, title = ';(q#bar{q} - t#bar{t} - glus) pt, nGlu=%d'%nglu)
        for i in glus :
            cmzGlu = qqcmzboost( p4[i] )
            cmGlu  = qqcmboost( p4[i] )
            self.book.fill(cmGlu.P(), 'cmgluP', 200,0,200, title = ';q#bar{q}cm gluP')
            self.book.fill(cmzGlu.pt(), 'cmzgluPt', 200,0,200, title = ';q#bar{q}cm_{z} gluP_{T}')
            self.book.fill(abs(cmzGlu.eta()), 'cmzgluAbsEta', 100,0,5.5, title = ';|cmz glu #eta|')

class mcQuestions2(analysisStep) :
    def __init__(self,pred=None) :
        self.pred = pred
        self.maxes = {'genttCosThetaStar':1,
                      #'genTopCosEqualThetaZ':1,
                      #'genTopCosThetaBoost':1,
                      'genTopCosThetaBoostAlt':1,
                      #'genTopDeltaAbsYttbar':3,
                      #'genTopBetaProjection':1,
                      'genTopDeltaBetazRel':1,
                      'genTopCosPhiBoost':1}

    def uponAcceptance(self,ev) :
        if self.pred and not ev[self.pred] : return
        for var,m in self.maxes.items() :
            self.book.fill(ev[var], var          , 100, -m,m, title = ';%s;events / bin'%var)
            self.book.fill(ev[var], var+"_lowres",   2, -m,m, title = ';%s;events / bin'%var)

        for ((v1,m1),(v2,m2)) in itertools.combinations(self.maxes.items(),2) :
            self.book.fill(ev[v1]*ev[v2], "product_%s_%s"%(v1,v2), 2, -m1*m2, m1*m2, title = ";%s . %s;events / bin"%(v1,v2))

class fractions(calculables.secondary) :

    def uponAcceptance(self,ev) :
        abssum = ev['genTopAbsSumRapidities']
        absdelta = ev['genTopDeltaAbsYttbar']
        self.book.fill(  abssum,  'abssum', 100,  0, 8, w = 1)
        self.book.fill(absdelta,'absdelta', 100, -3, 3, w = 1)

    def reportCache(self) :
        fileName = '/'.join(self.outputFileName.split('/')[:-1]+[self.name])
        optstat = r.gStyle.GetOptStat()
        r.gStyle.SetOptStat(0)
        c = r.TCanvas()
        c.Print(fileName+'.pdf[')

        names = ['%s#rightarrow^{}t#bar{t} '%i for i in ['gg','qg','q#bar{q}','#bar{q}g']]
        samples = ['ttj_ph.w%s.tw.pu'%s for s in ['GG','QG','QQ','AG']]
        colors = [r.kBlack,r.kBlue,r.kRed,r.kGreen]
        hists = self.fromCache(samples, ['abssum','absdelta'])

        def arrange(name, rebin = 1) :
            coll = [hists[s][name] for s in samples]
            for h_ in coll : h_.Rebin(rebin)
            h = coll[0].Clone();
            h.Reset()
            for h_ in coll : h.Add(h_)
            return h,coll
        abssum,abssumColl = arrange('abssum',4)
        absdelta,absdeltaColl = arrange('absdelta',2)
        leg = r.TLegend(0.72,0.55,0.87,0.80)
        leg.SetTextFont(102)
        for i,h in enumerate(abssumColl) :
            leg.AddEntry(h,names[i],'l')
            h.Divide(abssum)
            h.SetTitle(';|y_{t}+y_{#bar{t}}|;fraction  t#bar{t}')
            h.SetMinimum(0.01)
            h.SetMaximum(1.1)
            h.SetLineWidth(2)
            h.SetLineColor(colors[i])
            h.Draw('hist same'[:None if i else -4])
        c.SetLogy(1)
        leg.SetTextSize(0.04)
        leg.Draw()
        c.Print(fileName+'.pdf')
        leg2 = r.TLegend(0.62,0.55,0.87,0.88)
        leg2.SetTextSize(0.04)
        leg2.SetHeader('(N^{+}-N^{-}) / (N^{+}+N^{-})  (%)')
        #leg2.SetTextAlign(22)
        leg2.SetTextFont(102)
        for i,h in enumerate(absdeltaColl) :
            h.SetTitle(';|y_{t}|-|y_{#bar{t}}|;Events / bin / pb^{-1}')
            h.SetMinimum(0)
            h.SetLineWidth(2)
            h.SetLineColor(colors[i])
            h.Draw('hist same'[:None if i else -4])
            asymm = 100 * (h.Integral(h.FindFixBin(0.0001), h.GetNbinsX()+2) - h.Integral(0, h.FindFixBin(-0.0001))) / h.Integral(0,h.GetNbinsX()+2)
            leg2.AddEntry(h, names[i]+'%+.1f'%asymm,'l')
        c.SetLogy(0)
        leg.Draw()
        leg2.Draw()
        c.Print(fileName+'.pdf')
        c.Print(fileName+'.pdf]')
        r.gStyle.SetOptStat(optstat)
        print "Wrote file: %s.pdf"%fileName
