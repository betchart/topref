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
class jetPrinter(analysisStep) :
    def uponAcceptance(self,ev) :
        jets = ev['TopJets']['fixes']
        iTT = ev['genTTbarIndices']
        gen = ev['genP4']
        print 'mu,gen b,bbar,q*'
        for i in [max(iTT['lminus'],iTT['lplus']),iTT['b'],iTT['bbar']]+iTT['q'] :
            if i==None :
                print "--"
                continue
            p4 = gen[i]
            print ("\t%.0f\t%+.1f\t%+.1f")%(p4.pt(),p4.eta(),p4.phi())
        print 'jet'
        fPU = ev['%sPileUpPtFraction%s'%jets]
        for i in ev['%sIndices%s'%jets] :
            p4 = ev['%sCorrectedP4%s'%jets][i]
            print i, '\t',
            print 'b' if i in ev['%sIndicesGenB%s'%jets] else '',
            print 'q' if i in ev['%sIndicesGenWqq%s'%jets] else '',
            print ("\t%.0f\t%+.1f\t%+.1f")%(p4.pt(),p4.eta(),p4.phi()),
            print "  %.2f  "%fPU[i] ,
            print 'p' if i==ev['%sIndicesGenTopPQHL%s'%jets][0] else '',
            print 'q' if i==ev['%sIndicesGenTopPQHL%s'%jets][1] else '',
            print 'h' if i==ev['%sIndicesGenTopPQHL%s'%jets][2] else '',
            print 'l' if i==ev['%sIndicesGenTopPQHL%s'%jets][3] else '',
            print 'g' if i in ev['%sIndicesGenTopExtra%s'%jets] else ''
        if ev['%sIndicesOther%s'%jets] : print '-id'
        for i in ev['%sIndicesOther%s'%jets] :
            p4 = ev['%sCorrectedP4%s'%jets][i]
            print ("\t%.0f\t%+.1f\t%+.1f")%(p4.pt(),p4.eta(),p4.phi())
        print
#####################################
class pileupJets(analysisStep) :
    def uponAcceptance(self,ev) :
        xcjets = ev['TopJets']['fixes']
        jets = ev['TopJets']['fixesStripped']
        indices = ev['Indices'.join(xcjets)]
        indicesPU = ev['IndicesGenPileup'.join(xcjets)]
        p4 = ev['CorrectedP4'.join(xcjets)]
        n = ev['CountwithPrimaryHighPurityTracks'.join(jets)]
        nPU = ev['CountwithPileUpHighPurityTracks'.join(jets)]
        puPtFrac = ev['PileUpPtFraction'.join(xcjets)]

        for i in indices :
            eta = abs(p4[i].eta())
            label = "pileup" if i in indicesPU else "primary"
            label2 = ('inner' if eta<1.9 else 'middle' if eta<2.6 else 'outer')  + '_' +  label
            if "outer" in label2 : continue

            self.book.fill( p4[i].pt(), "pT_"+label, 50, 0, 100, title = ';p_{T};jets / bin')
            self.book.fill( n[i], "ntracks_%s"%label2, 30, 0, 30, title = ';ntracks (%s);jets / bin'%label2 )
            self.book.fill( puPtFrac[i], "pileupPtFrac_%s"%label2, 50, 0, 1, title = ';pileup Pt Fraction (%s);jets / bin'%label2)

        indicesPrimary = [i for i in indices if i not in indicesPU]
        if len(indicesPrimary) :
            iMaxPrimary = max(indicesPrimary, key = puPtFrac.__getitem__)
            eta = abs(p4[iMaxPrimary].eta())
            label = ('inner' if eta<1.9 else 'middle' if eta<2.6 else 'outer')  + '_' +  'primaryMax'
            self.book.fill(puPtFrac[iMaxPrimary], "pileupPtFrac_%s"%label, 50, 0, 1, title = ';pilup Pt Fraction (%s);events / bin'%label)

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

class combinatorialFiltering(analysisStep) :
    def __init__(self,jets=None) :
        for item in ["IndicesBtagged","CorrectedP4","IndicesGenTopPQHL","ComboPQBDeltaRawMassWTop"] :
            setattr(self, item, item.join(jets) )
        theta = math.pi/6
        self.ellipseR = np.array([[math.cos(theta),-math.sin(theta)],[math.sin(theta), math.cos(theta)]])
        
    def uponAcceptance(self,ev) :
        iGenPQHL = ev[self.IndicesGenTopPQHL]

        iP,iQ,iH,iL = iGenPQHL

        comboDRawMass = ev[self.ComboPQBDeltaRawMassWTop]
        indices = ev[self.IndicesBtagged]

        iBless = -1 if iH not in indices or iL not in indices else max( indices.index(iH) ,
                                                                        indices.index(iL) )
        
        self.book.fill( iBless , "max_bjet_genIndex", 8, -1.5, 6.5, title = ";b-index of lesser TCHE bjet;events / bin")

        for iPQH in itertools.permutations(indices,3) :
            if iPQH[0]>iPQH[1] : continue
            rawMs = comboDRawMass[iPQH]
            passEllipse = np.dot(*(2*[self.ellipseR.dot(rawMs) / [35,70]])) <= 1 
            tag = "correct" if iPQH == iGenPQHL[:3] else "incorrect"
            for ptag in ['','pass'][:2 if passEllipse else 1] :
                self.book.fill( (min(rawMs[0],150),min(rawMs[1],250)), "combo_raw_m_%s_%s"%(tag,ptag), (100,150),(-80,-150),(180,300), title = ";(%s) #Delta raw hadronic W mass;(%s) #Delta raw hadronic top mass;events / bin"%(tag,tag) )

class combinatorialFrequency(analysisStep) :
    def uponAcceptance(self,ev) :
        if not ev['genTopTTbar'] : return
        recos = ev['TopReconstruction']
        iP,iQ,iH,iL = ev['IndicesGenTopPQHL'.join(ev['TopJets']['fixes'])]

        igen=ev['genTopRecoIndex']

        iGenHadB = next((i for i,R in enumerate(recos) if R['iPQHL'][2]==iH) , -1)
        iGenLepB = next((i for i,R in enumerate(recos) if R['iPQHL'][3]==iL) , -1)
        iGenHadPQB  = next((i for i,R in enumerate(recos) if R['iPQHL'][:3]==(iP,iQ,iH)), -1 )
        iGenJets = next((i for i,R in enumerate(recos) if R['iPQHL']==(iP,iQ,iH,iL)), -1)

        for i,item in enumerate(['iGenHadB','iGenLepB','iGenHadPQB','iGenJets','igen']) :
            self.book.fill(eval(item), str(i)+item, 10,-1.5,8.5, title=";%s;events / bin"%item )

class combinatorialBG(analysisStep) :
    def __init__(self,jets=None) : self.jets = jets        
    def uponAcceptance(self,ev) :
        maxP = ev["TopComboQQBBMaxProbability"]
        iTrue = ev['genTopRecoIndex']
        recos = ev['TopReconstruction']
        jetIndices = set()
        hadIndices = set()
        bIndices = set()
        for iReco in [iTrue]+list(set(range(len(recos)))-set([iTrue])) :
            reco = recos[iReco]
            tag = "correct" if iReco==iTrue else "incorrect"
            indicesB = ev['%sIndicesBtagged%s'%self.jets]

            iPQHL = reco['iPQHL']
            iPQH = iPQHL[:-1]
            iHL = iPQHL[2:]

            iB2 = max(iHL, key = lambda i: indicesB.index(i))
            iiB2 = indicesB.index(iB2)
            self.book.fill(math.log(1+reco['chi2']), "logOnePlusChi2_"+tag, 50,0,7, title = ';%s ttbar kin. fit log(1+#chi^{2});events / bin'%tag )
            self.book.fill(math.log(1+reco['key']), 'logOnePlusKey_'+tag, 50,0,7, title = ';%s log(1+key);events / bin'%tag)
            self.book.fill( (math.log(1+reco['chi2']),math.log(1-2*math.log(reco['probability']))), "logOnePlus_chi2p_"+tag, (50,50),(0,0),(7,4), title = ";%s log(1+#chi^{2});log(1-2log(P));events / bin"%tag)
            
            if iPQHL not in jetIndices :
                jetIndices.add(iPQHL)
                p = ev['TopComboQQBBProbability'][iPQH[:2]+tuple(sorted(iHL))]
                self.book.fill( math.log(1-2*math.log(reco['probability'])), "logOnePlus_m2p_"+tag, 50,0,4, title = ";%s log(1-2log(P));events / bin"%tag)
                for i in range(4) : self.book.fill( ev["%sIndices%s"%self.jets].index(sorted(iPQHL)[i]), "iJ_j%d_"%i+tag, 10,-0.5,9.5, title = ';%s index of topjet%d;events / bin'%(tag,i) )

            if iPQH not in hadIndices :
                hadIndices.add(iPQH)
                mW,mT = ev["%sComboPQBRawMassWTop%s"%self.jets][iPQH]
                dmW,dmT = ev["%sComboPQBDeltaRawMassWTop%s"%self.jets][iPQH]
                self.book.fill(mW, "rawMassHadW_"+tag,60,0,180, title = ";raw mass_{W} had (%s);events / bin"%tag)
                self.book.fill(mT, "rawMassHadT_"+tag,100,0,300, title = ";raw mass_{t} had (%s);events / bin"%tag)
                topMass = 172; wMass = 80.4
                self.book.fill((dmT,dmW), "rawDeltaMassHadTW_"+tag,(100,60),(-topMass,-wMass),(300-topMass,180-wMass), title= ";(%s) raw had #Delta mass_{T};raw had #Delta mass_{W};events / bin"%tag)

            if self.jets and iHL not in bIndices :
                bIndices.add(iHL)
                indicesPt = ev['%sIndices%s'%self.jets]
                self.book.fill(indicesB.index(iB2), "iB_b2_"+tag, 10,-0.5,9.5, title = ';b-ordered index of second b (%s);events / bin'%tag)
                self.book.fill(indicesPt.index(iB2), "iPt_b2_"+tag, 10,-0.5,9.5, title = ';pt-ordered index of second b (%s);events / bin'%tag)

#####################################
class topProbLook(analysisStep) :
    def __init__(self, jets) :
        self.indicesGenB = "%sIndicesGenB%s"%jets
        self.indicesGenWqq = "%sIndicesGenWqq%s"%jets
        self.indices = "%sIndices%s"%jets
    def uponAcceptance(self,ev) :
        maxProb = ev["TopComboMaxProbability"]
        trueCombo = tuple( sorted(ev[self.indicesGenB]) + sorted(ev[self.indicesGenWqq]) )
        multiplicity = len(ev[self.indices])
        for key,val in ev["TopComboProbability"].iteritems() :
            tag = "true" if key == trueCombo else "other"
            self.book.fill(val, "topProbability"+tag, 100,0,1, title = ";%s top probability;events / bin"%tag)
            self.book.fill(val/maxProb, "topRelProbability"+tag, 100,0,1, title = ";%s top probability/ maxTopProb;events / bin"%tag)
            self.book.fill((val/maxProb,multiplicity), "topRelProbabilityByMulti"+tag, (100,10),(0,0),(1,10), title = ";%s top probability/ maxTopProb;jet multiplicity;events / bin"%tag)
            self.book.fill((maxProb,val), "topProbability_vMax"+tag, (100,100),(0,0),(1,1), title = ";maxTopProb;%s top probability;events / bin"%tag)
        self.book.fill(maxProb, "TopComboMaxProbability", 100,0,1, title = ';TopComboMaxProbability;events / bin')
        self.book.fill((maxProb,multiplicity), "TopComboMaxProbabilityLen"+self.indices, (100,10), (0,-0.5), (1,9.5), title = ';TopComboMaxProbability;jet multiplicity;events / bin')
#####################################
class kinematics(analysisStep) :
    def __init__(self,indexName) : self.moreName = indexName
    def uponAcceptance(self,ev) :
        index = ev["%sRecoIndex"%self.moreName]
        if index < 0 : return
        topReco = ev["TopReconstruction"]
        self.book.fill(ev[self.moreName+"TtxMass"], "TTX.mass", 50,300,1300, title = ";ttx invariant mass;events / bin")
        self.book.fill(ev[self.moreName+"PtSum"],   "TT.pt",   100,  0, 200, title = ";ttbar.pt;events / bin")
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
class discriminateNonTop(analysisStep) :
    def __init__(self, pars) :
        obj = pars['objects']
        lepCollection = obj[pars['lepton']['name']]
        self.MT = "%sMt%s"%lepCollection+"mixedSumP4"
        self.sumPt = obj["sumPt"]
        self.HT = "%sSumPt%s"%obj["jet"]
        self.jetP4 = "%sCorrectedP4%s"%obj["jet"]
        self.iJet = "%sIndices%s"%obj["jet"]
        self.bJet = "%sIndicesBtagged%s"%obj["jet"]        
        self.lepP4 = "%sP4%s"%lepCollection
        self.iLep = "%sSemileptonicTopIndex%s"%lepCollection
        self.kT = "%sKt%s"%obj["jet"]

    def uponAcceptance(self, ev) :
        jetP4 = ev[self.jetP4]
        iJet = ev[self.iJet]
        bJet = ev[self.bJet]
        lepP4 = ev[self.lepP4][ev[self.iLep]]
        
        self.book.fill( ev["mixedSumP4"].pt(), "met", 60, 0, 180, title = ';met;events / bin')
        self.book.fill( lepP4.pt(), "leptonPt", 60, 0, 180, title = ';lepton Pt;events / bin')
        self.book.fill( abs(lepP4.eta()), "leptonEta", 50, 0, 3, title = ';lepton |#eta|;events / bin')
        self.book.fill(ev[self.MT],self.MT,30,0,180, title = ";M_{T};events / bin")
        dphiLnu = abs(r.Math.VectorUtil.DeltaPhi(lepP4,ev["mixedSumP4"]))
        self.book.fill( dphiLnu, "dphiLnu", 20,0,math.pi, title = ";#Delta#phi l#nu;events / bin" )
        self.book.fill( (dphiLnu,ev[self.MT]), "dphiLnu_v_mt", (20,30),(0,0),(math.pi,180), title = ";#Delta#phi l#nu;M_{T};events / bin" )
        self.book.fill( ev[self.kT], "kT", 30, 0, 150, title = ";k_{T};events / bin")

        self.book.fill(ev[self.sumPt],self.sumPt,50,0,1500, title = ';%s;events / bin'%self.sumPt)
        self.book.fill(ev[self.HT],self.HT,50,0,1500, title = ';%s;events / bin'%self.HT)
        
        self.book.fill(jetP4[iJet[0]].pt(), "jetPtI0", 40,0,400, title = ';pT jets[0] pt;events / bin')
        self.book.fill(jetP4[bJet[0]].pt(), "jetPtB0", 40,0,400, title = ';b- jets[0] pt;events / bin')
        self.book.fill( abs(r.Math.VectorUtil.DeltaPhi(jetP4[bJet[0]],jetP4[bJet[1]])), "dphiBjets", 20,0,math.pi,
                        title = ";#Delta#phi leading b-tagged jets;events / bin" )
        for i in range(min(4,len(iJet))) :
            self.book.fill(abs(jetP4[iJet[i]].eta()), "jetEtaI%d"%i, 40,0,4, title = ';pT jets[%d] |#eta|;events / bin'%i)

#####################################
class discriminateQQbar(analysisStep) :
    def __init__(self, collection) :
        for item in ['DeltaPhi','PtOverSumPt','SumP4','CosThetaDaggerTT'] :
            setattr(self,item,('%s'+item+'%s')%collection)
        
    @staticmethod
    def phiMod(phi) : return phi + 2*math.pi*int(phi<0)

    def uponAcceptance(self,ev) :
        if not ev['genTopTTbar'] : return

        dphi = self.phiMod(ev[self.DeltaPhi])

        ### dphi is highly correlated with PtAsym and/or PtOverSumPt, but they are mostly uncorrelated to alpha
        #self.book.fill( (dphi,ev['genTopPtAsymttbar']), 'corrDphiPtAsym', (51,51), (0,-1),(2*math.pi,1), title=';dphi;ptasymm;events / bin' )
        #self.book.fill( (dphi,ev['genTopAlpha']), 'corrDphiAlpha', (51,10), (0,0),(2*math.pi,1), title=';dphi;#alpha;events / bin' )
        #self.book.fill( (ev['genTopTTbarPtOverSumPt'],ev['genTopAlpha']), 'corrPtAsymAlpha', (50,10), (0,0),(1,1), title=';(t+tbar)_{pt}/(t_{pt}+tbar_{pt});#alpha;events / bin' )

        self.book.fill( dphi, self.DeltaPhi, 51,0,2*math.pi, title = ';#Delta #phi_{ttbar};events / bin')
        self.book.fill( ev[self.PtOverSumPt], self.PtOverSumPt, 50,0,1, title = ';(t+tbar)_{pt}/(t_{pt}+tbar_{pt});events / bin')
        self.book.fill( ev[self.CosThetaDaggerTT], self.CosThetaDaggerTT, 50,-1,1, title = ';cos#theta^{#dagger}_{tt}')
        sumP4 = ev[self.SumP4]
        self.book.fill( abs(sumP4.Rapidity()), self.SumP4+'AbsRapidity', 50,0,3, title = ';y_{ttbar};events / bin')
        self.book.fill( abs(sumP4.Eta()), self.SumP4+'AbsEta', 40,0,10, title = ';|#eta_{ttbar}|;events / bin')
        self.book.fill( abs(sumP4.Pz()), self.SumP4+'AbsPz', 50,0,3000, title = ';|pz|_{ttbar};events / bin')
        
#####################################
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
