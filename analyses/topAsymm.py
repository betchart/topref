import supy,steps,calculables,samples
import os,math,copy,itertools,ROOT as r, numpy as np

class topAsymm(supy.analysis) :
    ''' Analysis for measurement of production asymmetry in top/anti-top pairs
    
    This analysis contains several secondary calculables, which need
    to be primed in the following order:
    
    Pass#
    0: ^selection : prime sim.*   : pileup
    1: ^finegrain : prime top.top : RawMassWTopCorrectPQBTwoDChiSquared,jetCSVProbabilityGivenBQN,LTopUnfitSqrtChi2CombinationsLR
    2: ^finegrain : prime top.top : HTopSigmasPQBCombinationsLR
    3: ^finegrain : prime needed  : TridiscriminantWTopQCD
    4: ^genTop    : prime top.top : TridiscriminantGGqqGq
    5.            : *
    '''

    def parameters(self) :

        reweights = {
            'abbr'  :'pu',
            'func'  :'pileup',
            'var'   :'pileupTrueNumInteractionsBX0Target'
            }

        leptons = {
            'name'     : [                    'mu',                      'el'],
            'ptMin'    : [                   26.0 ,                     30.0 ],
            'etaMax'   : [                    2.1 ,                      2.5 ],
            'triggers' : [       self.mutriggers(),         self.eltriggers()],
            'isoNormal': [             {"max":0.12},             {'max':0.10}],
            'isoInvert': [ {"min":0.13, "max":0.20}, {'min':0.11, 'max':0.15}]
            }

        #btag working points: https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagPerformanceOP
        csvWP = {"CSVL" : 0.244, "CSVM" : 0.679, "CSVT" : 0.898 }
        bCut = {"normal"   : {"index":0, "min":csvWP['CSVM']},
                "inverted" : {"index":0, "min":csvWP['CSVL'], "max":csvWP['CSVM']}}

        return { "vary" : ['selection','lepton','toptype'],
                 "discriminant2DPlots": True,
                 "bVar" : "CSV", # "Combined Secondary Vertex"
                 "objects" : dict([(item,(item,'')) for item in ['jet','mu','el','met']]),
                 "lepton" : self.vary([ ( leptons['name'][index], dict((key,val[index]) for key,val in leptons.iteritems())) for index in range(2) if leptons['name'][index] in ['mu','el'][:2]]),
                 "reweights" : reweights,
                 "selection" : self.vary({"top" : {"bCut":bCut["normal"],  "iso":"isoNormal"},
                                          "QCD" : {"bCut":bCut["normal"],  "iso":"isoInvert"}
                                          }),
                 "toptype" : self.vary({"ph":"ph"}),
                 "topBsamples": ("ttj_%s",['ttj_%s.wGG.%s','ttj_%s.wQG.%s','ttj_%s.wAG.%s','ttj_%s.wQQ.%s']),
                 "smear" : "Smear",
                 }

    @staticmethod
    def scaleFactor() : return 0.9

    @staticmethod
    def mutriggers() :
        pattern = "HLT_IsoMu%d_eta2p1_TriCentralPF%sJet%s_v%d"
        defs = [( (17,'NoPU','30_30_20') , (1,2) ),
                ( (17,'NoPU',      '30') , (1,) ),
                ( (20,'NoPU',      '30') , (4,) ),
                ( (17,    '',      '30') , (2,3,4,5) ),
                ( (20,    '',      '30') , (2,3,4) ),
                ( (17,'NoPU','50_40_30') , (1,) ),
                ( (20,'NoPU','50_40_30') , (4,) ),
                ( (20,    '','50_40_30') , (2,3,4) ),
                ]
        return (['HLT_IsoMu24_CentralPFJet30_CentralPFJet25_v%d'%v for v in (1,2)] +
                ['HLT_Mu24_eta2p1_CentralPFJet30_CentralPFJet25_v3'] +
                ['HLT_IsoMu24_eta2p1_CentralPFJet30_CentralPFJet25_v%d'%v for v in (3,4,5)]+
                sum([[pattern%(d+(v,)) for v in vs] for d,vs in defs],[])
                )
    
    @staticmethod
    def eltriggers() :
        pattern = "HLT_Ele25_CaloIdVT_CaloIso%s_TrkId%s_TrkIsoT_TriCentralPF%sJet%s_v%d"
        defs = [ ( ( 'T', 'T','NoPU','30_30_20'), (1,2)    ),
                 ( ('VL','VL','NoPU',      '30'), (1,2)),
                 ( ( 'T', 'T','NoPU',      '30'), (5,)     ),
                 ( ( "T", "T",    '',      "30"), (8,9,10) ),
                 ( ('VL','VL','NoPU','50_40_30'), (1,2)),
                 ( ( 'T', 'T','NoPU','50_40_30'), (5,)     ),
                 ( ( 'T', 'T',    '','50_40_30'), (3,4,5)  ),
                 ]
        return sum([[pattern%(d+(v,)) for v in vs] for d,vs in defs],[])
    ########################################################################################

    def listOfSampleDictionaries(self) : return [getattr(samples,item) for item in ['lepton112', 'top112', 'ewk112', 'qcd112']]

    @staticmethod
    def muons() : return ['MuHad.2012A_1','MuHad.2012A_2','SingleMu.2012B','SingleMu.2012C']
    @staticmethod
    def electrons() : return ['ElHad.2012A_1','ElHad.2012A_2','SingleEl.2012B','SingleEl.2012C']
    @staticmethod
    def single_top() : return ['top_s_ph','top_t_ph','top_tW_ph','tbar_s_ph','tbar_t_ph','tbar_tW_ph']
    @staticmethod
    def qcd_mu(suf='') :
        bins = [15,20,30,50,80,120,170,300,470,600,800,1000,None]
        return ['qcd_mu_%s'%'_'.join(str(i) + suf for i in [lo,hi] if i) for lo,hi in zip(bins[:-1],bins[1:])]
    @staticmethod
    def qcd_el(suf='') :
        bins = [20,30,80,170,250,350,None]
        return ['qcd_em_%s'%'_'.join(str(i) + suf for i in [lo,hi] if i) for lo,hi in zip(bins[:-1],bins[1:])]

    def listOfSamples(self,pars) :
        rw = pars['reweights']['abbr']
        lname = pars['lepton']["name"]

        def ewk(eL = None) :
            dy = ['dy%dj_mg'%n for n in range(1,5)] if "QCD" not in pars['tag'] else []
            w =  ['w%dj_mg'%n for n in range(1,5)]
            return supy.samples.specify( names = ( w + dy ),
                                         effectiveLumi = eL, color = 28, weights = [rw] ) #if "QCD" not in pars['tag'] else []

        def single_top(eL = None) :
            return supy.samples.specify( names = self.single_top(),
                                         effectiveLumi = eL, color = r.kGray, weights = [rw]) #if "QCD" not in pars['tag'] else []
        
        def qcd(eL = None) :
            return supy.samples.specify( names = { 'mu' : self.qcd_mu(), 'el' : self.qcd_el(),}[pars['lepton']['name']],
                                         effectiveLumi = eL, color = 9, weights = [rw] ) if "QCD" not in pars['tag'] else []

        def ttbar(eL = None) :
            tt = pars['toptype']
            color = [r.kBlue,r.kCyan,r.kCyan+1,r.kOrange]
            wSub =  [  "wGG",  "wQG",    "wAG",    "wQQ"]
            return sum([supy.samples.specify(names = 'ttj_%s'%tt, effectiveLumi = eL,
                                             color = c, weights = [w,rw]) for w,c in zip(wSub,color)],[])
        
        def data() :
            return { "mu" : supy.samples.specify( names = self.muons()),
                     "el" : supy.samples.specify( names = self.electrons()),
                     }[pars['lepton']['name']]

        return  ( data() + ewk() + ttbar() + single_top() ) #+ qcd() )

    ########################################################################################
    def listOfCalculables(self, pars) :
        obj = pars["objects"]
        jet = obj['jet']
        met = obj['met']
        mu = obj['mu']
        el = obj['el']
        lepton = obj[pars["lepton"]["name"]]
        calcs = sum( [supy.calculables.zeroArgs(module) for module in [calculables, supy.calculables]]+
                     [supy.calculables.fromCollections(getattr(calculables,item), [obj[it]])
                      for item,it in  [('jet','jet'),('electron','el'),('muon','mu')]], [])
        calcs += supy.calculables.fromCollections(calculables.top,[('genTop',""),('fitTop',"")])
        calcs += [
            calculables.met.AdjustedP4(met, jet, pars['smear']),
            calculables.jet.AdjustedP4(jet, pars['smear']),
            calculables.jet.Indices(jet,ptMin = 20 ),
            calculables.jet.IndicesBtagged(jet,pars["bVar"]),
            supy.calculables.other.abbreviation('combinedSecondaryVertex','CSV',jet),
            calculables.muon.Indices(mu),
            calculables.electron.Indices(el),

            calculables.trigger.lowestUnPrescaledTrigger(pars["lepton"]["triggers"]),

            calculables.gen.genIndicesHardPartons(),
            calculables.top.TopJets( jet ),
            calculables.top.TopLeptons( lepton ),
            calculables.top.TopReconstruction(),

            calculables.comb.TopComboQQBBLikelihood( pars['bVar'] ),
            calculables.comb.OtherJetsLikelihood( pars['bVar'] ),
            calculables.comb.TopRatherThanWProbability( priorTop = 0.05 ),

            calculables.met.MetMt( lepton, "AdjustedP4".join(met)),
            calculables.met.Covariance( met ),
            calculables.other.KarlsruheDiscriminant( jet, "AdjustedP4".join(met) ),

            calculables.jet.pt( jet, index = 0, Btagged = True ),
            calculables.jet.pt( jet, index = 3, Btagged = False ),
            calculables.jet.absEta( jet, index = 3, Btagged = False),

            supy.calculables.other.pt( "AdjustedP4".join(met) ),
            supy.calculables.other.size( "Indices".join(jet) ),
            #supy.calculables.other.abbreviation( pars['reweights']['var'], pars['reweights']['abbr'] ),
            supy.calculables.other.fixedValue('one',1),
            supy.calculables.other.abbreviation( 'one', pars['reweights']['abbr'] ),
            ]
        return calcs
    ########################################################################################

    def listOfSteps(self, pars) :
        obj = pars["objects"]
        jet = obj['jet']
        met = obj['met']
        lname = pars["lepton"]["name"]
        lepton = obj[lname]
        otherLepton = obj[{'el':'mu','mu':'el'}[lname]]
        lIsoMinMax = pars["lepton"][pars['selection']['iso']]
        bVar = pars["bVar"].join(jet)
        rw = pars['reweights']['abbr']
        tt = pars['toptype']
        topSamples = (pars['topBsamples'][0]%tt,[s%(tt,rw) for s in pars['topBsamples'][1]])
        topTag = pars['tag'].replace("QCD","top")
        
        ssteps = supy.steps

        saDisable = 'ttj' not in pars['sample'] or not any(w in pars['sample'].split('.') for w in ['wQQ','wQG','wAG','wGG'])
        saWeights = []
        return (
            [ssteps.printer.progressPrinter()
             , ssteps.histos.value("genQ",200,0,1000,xtitle="#hat{Q} (GeV)").onlySim()
             #, steps.top.fractions().disable(saDisable)
             #, getattr(self,pars['reweights']['func'])(pars)
             #, supy.calculables.other.SymmAnti(pars['sample'],"genTopCosPhiBoost",1, inspect=True, nbins=160, weights = saWeights,
             #                                  funcEven = r.TF1('phiboost',"[0]*(1+[1]*x**2)/sqrt(1-x**2)",-1,1),
             #                                  funcOdd = r.TF1('phiboostodd','[0]*x/sqrt(1-x**2)',-1,1)).disable(saDisable)
             #, supy.calculables.other.SymmAnti(pars['sample'],"genTopDeltaBetazRel",1, inspect=True, weights = saWeights,
             #                                  funcEven = '++'.join(['(1-abs(x))']+['x**%d'%d for d in [0,2,4,6,8,10,12,14,16,18]]),
             #                                  funcOdd = '++'.join(['x**%d'%d for d in [1,3,5,7,9,11,13]])).disable(saDisable)
             , calculables.top.ttSymmAnti(pars['sample'], inspect=True).disable(saDisable)
             ####################################
             , ssteps.filters.label('selection').invert(),
             ssteps.filters.value("muHandleValid",min=True),
             ssteps.filters.multiplicity( max=1, min=1, var = "Charge".join(lepton)),
             ssteps.filters.multiplicity( max=0, var = "Charge".join(otherLepton)),
             ssteps.filters.multiplicity( min=1, var = 'Indices'.join(lepton)),
             ssteps.filters.value('PassConversionVeto'.join(lepton), min=True, indices='Indices'.join(lepton), index=0) if lname=='el' else supy.steps.filters.label('empty'),
             ssteps.filters.value("Pt".join(jet), min = 45, indices = "Indices".join(jet), index=0),
             ssteps.filters.value("Pt".join(jet), min = 45, indices = "Indices".join(jet), index=1),
             ssteps.filters.value("Pt".join(jet), min = 35, indices = "Indices".join(jet), index=2),
             ssteps.filters.value("Pt".join(jet), min = 20, indices = "Indices".join(jet), index=3),
             ssteps.filters.value(bVar, indices = "IndicesBtagged".join(jet), **pars["selection"]["bCut"]),
             ssteps.filters.value('RelIso'.join(lepton), indices='Indices'.join(lepton), index=0, **lIsoMinMax),
             steps.trigger.lowestUnPrescaledTriggerFilter().onlyData() #, steps.trigger.hltKeys().onlyData(),
             , steps.trigger.lowestUnPrescaledTriggerHistogrammer(drop = {'mu' : ['IsoMu','_eta2p1_','CentralPF','Jet'],
                                                                          'el' : ['Ele25_CaloIdVT_','CaloIso','TrkId','_TrkIsoT','TriCentralPF','Jet']}[lname]).onlyData()
             ####################################
             , ssteps.filters.label("selection complete")
             , steps.top.channelClassification().onlySim()
             , ssteps.histos.multiplicity('Indices'.join(jet))

             , ssteps.filters.label("secondaries")
             , calculables.jet.ProbabilityGivenBQN(jet, pars['bVar'], binning=(51,-0.02,1), samples = topSamples, tag = topTag)
             , calculables.jet.ScalingBQN(jet, samples = topSamples[1], tag = topTag)
             , supy.calculables.other.TwoDChiSquared('RawMassWTopCorrectPQB', samples = topSamples[1], tag = topTag, binningX = (300,0,600), binningY = (300,0,1200), labelsXY = ("Raw Hadronic m_{W}","Raw Hadronic m_{top}"), tailSuppression=0.01)
             , supy.calculables.other.CombinationsLR( var = 'HTopSigmasPQB', varMax = 5, trueKey = 'IndicesGenTopPQH', samples = topSamples[1], tag = topTag)
             , supy.calculables.other.CombinationsLR( var = 'LTopUnfitSqrtChi2',varMax = 10, trueKey = 'IndexGenTopL', samples = topSamples[1], tag = topTag)
             , self.tridiscriminant(pars)

             , ssteps.filters.label('finegrain')
             , ssteps.histos.value('MetMt'.join(lepton), 120, 0, 120)
             , ssteps.histos.value('ProbabilityHTopMasses', 100,0,1)
             , ssteps.histos.value("TopRatherThanWProbability", 100,0,1)

             , ssteps.filters.label('object pt')
             , ssteps.histos.pt('AdjustedP4'.join(met), 125, 0, 250 )
             , ssteps.histos.pt('P4'.join(lepton), 125, 0, 250, indices = 'Indices'.join(lepton), index=0)
             , ssteps.histos.value("Pt".join(jet), 125, 0, 250, indices = 'Indices'.join(jet), index=0)
             , ssteps.histos.value("Pt".join(jet), 125, 0, 250, indices = 'Indices'.join(jet), index=1)
             , ssteps.histos.value("Pt".join(jet), 125, 0, 250, indices = 'Indices'.join(jet), index=2)
             , ssteps.histos.value("Pt".join(jet), 125, 0, 250, indices = 'Indices'.join(jet), index=3)

             , ssteps.filters.label('object eta')
             , ssteps.histos.absEta("P4".join(lepton), 100,0,4, indices = "Indices".join(lepton), index = 0)
             , ssteps.histos.absEta("AdjustedP4".join(jet), 100,0,4, indices = "Indices".join(jet), index = 0)
             , ssteps.histos.absEta("AdjustedP4".join(jet), 100,0,4, indices = "Indices".join(jet), index = 1)
             , ssteps.histos.absEta("AdjustedP4".join(jet), 100,0,4, indices = "Indices".join(jet), index = 2)
             , ssteps.histos.absEta("AdjustedP4".join(jet), 100,0,4, indices = "Indices".join(jet), index = 3)

             ####################################
             #, steps.displayer.ttbar(jets=jet, met=obj['met'], muons = obj['mu'], electrons = obj['el'])
             , ssteps.filters.label('top reco')
             , steps.top.combinatorialFrequency().onlySim()
             , ssteps.histos.value('genTopRecoIndex', 10,-1.5,8.5)
             , ssteps.histos.value('TopGenLikelihoodIndex', 10,-1.5,8.5)
             , ssteps.histos.value('TopFitLikelihoodIndex',10,-1.5,8.5)
             , ssteps.histos.value('TopFitLikelihoodCorrectIndex',10,-1.5,8.5)

             , ssteps.filters.label('GGqqGq')
             , self.tridiscriminant2(pars)

             , ssteps.filters.label('genTop')
             , steps.top.kinFitLook('genTopRecoIndex')
             , steps.top.resolutions('genTopRecoIndex')

             , ssteps.filters.label('recoTop')
             , steps.top.kinFitLook('fitTopRecoIndex')
             , steps.top.resolutions('fitTopRecoIndex')
             , steps.top.kinematics('fitTop')
             ####################################

             , ssteps.filters.label('signal distributions')
             , ssteps.histos.symmAnti('genTopCosPhiBoost','genTopCosPhiBoost',100,-1,1).disable(saDisable)
             , ssteps.histos.symmAnti('genTopDeltaBetazRel','genTopDeltaBetazRel',100,-1,1).disable(saDisable)

             , ssteps.histos.symmAnti('genTopCosPhiBoost','fitTopCosPhiBoost',100,-1,1, other = ('TridiscriminantWTopQCD',100,-1,1))
             , ssteps.histos.symmAnti('genTopDeltaBetazRel','fitTopDeltaBetazRel',100,-1,1, other = ('TridiscriminantWTopQCD',100,-1,1))

             ])
    ########################################################################################

    @classmethod
    def pileup(cls,pars) :
        rw = pars['reweights']['abbr']
        tt = pars['toptype']
        return supy.calculables.other.Target("pileupTrueNumInteractionsBX0", thisSample = pars['baseSample'],
                                             target = ("data/pileup_true_Cert_160404-180252_7TeV_ReRecoNov08_Collisions11_JSON.root","pileup"),
                                             groups = [('qcd_mu',[]),
                                                       ('dy1j_ll_mg',[]),
                                                       ("w2j_mg",[]),("w3j_mg",[]),("w4j_mg",[]),
                                                       ('single_top', ['%s.%s'%(s,rw) for s in cls.single_top()]),
                                                       ('ttj_%s'%tt,['ttj_%s.%s.%s'%(tt,s,rw) for s in ['',
                                                                                                           'wGG','wQG','wAG','wQQ']])]).onlySim()

    def tridiscriminant(self,pars) :
        rw = pars['reweights']['abbr']
        lname = pars['lepton']['name']
        tt = pars['toptype']
        tops = ['ttj_%s.%s.%s'%(tt,s,rw) for s in ['wGG','wQG','wAG','wQQ']]
        others = ['%s.%s'%(o,rw) for o in self.single_top() + ['w%dj_mg'%d for d in [1,2,3,4]]]
        datas = {"mu" : self.muons(),
                 "el": self.electrons()}[lname]
        sf = self.scaleFactor()
        return supy.calculables.other.Tridiscriminant( fixes = ("","WTopQCD"),
                                                       zero = {"pre":"ttj_%s"%tt, "tag":"top_%s_%s"%(lname,tt), "samples": tops},
                                                       pi23 = {"pre":"Multijet", "tag":"QCD_%s_%s"%(lname,tt), "samples":['data']+tops+others, 'sf':[1] + [-sf]*len(tops+others)},
                                                       pi43 = {"pre":"wj", "tag":"top_%s_%s"%(lname,tt), "samples":['w%dj_mg.%s'%(n,rw) for n in [1,2,3,4]]},
                                                       correlations = pars['discriminant2DPlots'],
                                                       otherSamplesToKeep = datas,
                                                       dists = {"TopRatherThanWProbability" : (20,0,1),
                                                                'ProbabilityHTopMasses' : (20,0,1),
                                                                "MetMt".join(pars['objects'][lname]) : (20,0,100),
                                                                })
    @staticmethod
    def tridiscriminant2(pars) :
        rw = pars['reweights']['abbr']
        lname = pars['lepton']['name']
        tt = pars['toptype']
        jets = pars["objects"]['jet']
        return supy.calculables.other.Tridiscriminant( fixes = ("","GGqqGq"),
                                                       zero = {"pre":"gg", "tag":"top_%s_%s"%(lname,tt), "samples":['ttj_%s.wGG.%s'%(tt,rw)]},
                                                       pi23 = {"pre":"qq", "tag":"top_%s_%s"%(lname,tt), "samples":['ttj_%s.wQQ.%s'%(tt,rw)]},
                                                       pi43 = {"pre":"qg", "tag":"top_%s_%s"%(lname,tt), "samples":['ttj_%s.%s.%s'%(tt,s,rw) for s in ['wQG','wAG']]},
                                                       correlations = pars['discriminant2DPlots'],
                                                       dists = { "fitTopPtPlusSumPt" : (20,0,600),
                                                                 "fitTopPtOverSumPt" : (20,0,1),
                                                                 "fitTopAbsSumRapidities" : (20, 0, 4),
                                                                 #"fitTopSqrtPtOverSumPt" : (10,0,1),
                                                                 #"fitTopSumP4AbsEta" : (20,0,6),
                                                                 #"fitTopPartonLo" : (20,-0.2,0.2),
                                                                 #"fitTopPartonHi" : (20,0,0.4),
                                                                 #"fitTopPartonXlo" : (20,0,0.12),              # 0.036
                                                                 #"fitTopPartonXhi" : (20,0.04,0.4),           # 0.003
                                                                 })
    ########################################################################################
    def concludeAll(self) :
        self.orgMelded = {}
        self.sensitivityPoints = {}
        self.rowcolors = 2*[13] + 2*[45]
        super(topAsymm,self).concludeAll()
        for tt,rw,lname in set([(pars['toptype'],pars['reweights']['abbr'],pars['lepton']['name']) for pars in self.readyConfs]) :
            self.meldScale(rw,lname,tt)
            self.plotMeldScale(rw,lname,tt)
            #self.PEcurves(rw,lname, tt)
        for tt,rw in set([(pars['toptype'],pars['reweights']['abbr']) for pars in self.readyConfs]) :
            self.measureAmplitudes(rw,tt)
        #self.sensitivity_graphs()
        #self.grant_proposal_plots()

    def conclude(self,pars) :
        rw = pars['reweights']['abbr']
        tt = pars['toptype']
        org = self.organizer(pars, verbose = True )

        if pars['lepton']['name']=='mu' :
            org.mergeSamples(targetSpec = {"name":"Mu.2012", "color":r.kBlack, "markerStyle":20}, sources = self.muons() )
        else:
            org.mergeSamples(targetSpec = {"name":"El.2012", "color":r.kBlack, "markerStyle":20}, sources = self.electrons())
            
        #org.mergeSamples(targetSpec = {"name":"t#bar{t}", "color":r.kViolet}, sources=["ttj_%s.%s.%s"%(tt,s,rw) for s in ['wQQ','wQG','wAG','wGG']])#, keepSources = True)
        org.mergeSamples(targetSpec = {"name":"W+jets", "color":28}, allWithPrefix = 'w')
        org.mergeSamples(targetSpec = {"name":"DY+jets", "color":r.kYellow}, allWithPrefix="dy")
        org.mergeSamples(targetSpec = {"name":"Single top", "color":r.kGray}, sources = ["%s.%s"%(s,rw) for s in self.single_top()])
        org.mergeSamples(targetSpec = {"name":"St.Model", "color":r.kGreen+2}, sources = ["t#bar{t}","W+jets","DY+jets","Single top","Multijet"], keepSources = True)

        #self.skimStats(org)

        org.scale( lumiToUseInAbsenceOfData = 5814 )

        names = [ss["name"] for ss in org.samples]
        kwargs = {"detailedCalculables": False,
                  "blackList":["lumiHisto","xsHisto","nJobsHisto"],
                  "samplesForRatios" : next(iter(filter(lambda x: x[0] in names and x[1] in names, [("Mu.2012","St.Model"),
                                                                                                    ("El.2012","St.Model")])), ("","")),
                  "sampleLabelsForRatios" : ("data","s.m."),
                  "detailedCalculables" : True,
                  "rowColors" : self.rowcolors,
                  "rowCycle" : 100,
                  "omit2D" : False,
                  }
        
        #supy.plotter(org, pdfFileName = self.pdfFileName(org.tag+"_log"),  doLog = True, pegMinimum = 0.01, **kwargs ).plotAll()
        supy.plotter(org, pdfFileName = self.pdfFileName(org.tag+"_nolog"), doLog = False, **kwargs ).plotAll()

    def skimStats(self,org) :
        statsname = {'DY+jets':'dy',
                     'W+jets': 'wj',
                     'Mullijet': 'qcd',
                     't#bar{t}':'tt',
                     'ttj_ph.wGG.pu':'ttgg',
                     'ttj_ph.wQG.pu':'ttqg',
                     'ttj_ph.wAG.pu':'ttag',
                     'ttj_ph.wQQ.pu':'ttqq',
                     'Single top' : 'st',
                     'Mu.2012': 'mu',
                     'El.2012': 'el'
                     }

        fileName = '%s/stats_%s.root'%(self.globalStem,org.tag)
        tfile = r.TFile.Open(fileName,'RECREATE')
        grab = ['beamHaloCSCLooseHaloId', 'triD_v_sqtsumptopt']
        for g in grab :
            tfile.mkdir(g).cd()
            for ss,hist in zip( org.samples,
                                org.steps[next(org.indicesOfStepsWithKey(g))][g] ) :
                if not hist or ss['name']=='Standard Model': continue
                h = hist.Clone(statsname[ss['name']] if ss['name'] in statsname else ss['name'])
                h.Write()
        tfile.Close()
        print 'Wrote: ', fileName

    def plotMeldScale(self,rw,lname,tt) :
        if (lname,rw, tt) not in self.orgMelded : print "run meldScale() before plotMeldScale()"; return
        melded = copy.deepcopy(self.orgMelded[(lname,rw,tt)])
        for s in ['top.ttj_%s.%s.%s'%(tt,s,rw) for s in ['wQQ','wQG','wAG','wGG']] :
            melded.drop(s)
        for log,label in [(False,""),(True,"_log")][:1] : 
            pl = supy.plotter(melded, pdfFileName = self.pdfFileName(melded.tag + label),
                              doLog = log,
                              blackList = ["lumiHisto","xsHisto","nJobsHisto"],
                              rowColors = self.rowcolors,
                              #samplesForRatios = ("top.Data 2012","S.M."),
                              #sampleLabelsForRatios = ('data','s.m.'),
                              rowCycle = 100,
                              omit2D = True,
                              pageNumbers = False,
                              ).plotAll()

    def meldScale(self,rw,lname,tt) :
        meldSamples = {"top_%s_%s"%(lname,tt) :( { 'mu': self.muons(),
                                                    'el': self.electrons()}[lname]+
                                                 ["ttj_%s"%tt]+self.single_top()+
                                                 ["dy%dj_mg"%n for n in [1,2,3,4]]+
                                                 ["w%dj_mg"%n for n in [1,2,3,4]]),
                       "QCD_%s_%s"%(lname,tt) : ( { 'mu':self.muons(),
                                                   'el':self.electrons()}[lname] +
                                                  ["ttj_%s"%tt] ) }

        organizers = [supy.organizer(tag, [s for s in self.sampleSpecs(tag) if any(item in s['name'] for item in meldSamples[tag])])
                      for tag in [p['tag'] for p in self.readyConfs if p["tag"] in meldSamples]]

        if len(organizers) < len(meldSamples) : return
        for org in organizers :
            org.mergeSamples(targetSpec = {"name":"t#bar{t}", "color":r.kViolet}, sources=["ttj_%s.%s.%s"%(tt,s,rw) for s in ['wQQ','wQG','wAG','wGG']], keepSources = 'top' in org.tag)
            org.mergeSamples(targetSpec = {"name":"W", "color":r.kRed}, sources = ["w%dj_mg.%s"%(n,rw) for n in [1,2,3,4]] )
            org.mergeSamples(targetSpec = {"name":"DY", "color":28}, allWithPrefix = "dy")
            org.mergeSamples(targetSpec = {"name":"Single", "color":r.kGray}, sources = ["%s.%s"%(s,rw) for s in self.single_top()], keepSources = False )
            org.mergeSamples(targetSpec = {"name":"Data 2012", "color":r.kBlack, "markerStyle":20}, sources={'mu':self.muons(),'el':self.electrons()}[lname])
            org.scale()
            if "QCD_" in org.tag :
                org.mergeSamples(targetSpec = {"name":"multijet","color":r.kBlue},
                                 sources=["Data 2012",'t#bar{t}'],
                                 scaleFactors = [1,-self.scaleFactor()],
                                 force=True, keepSources = False)

        self.orgMelded[(lname,rw,tt)] = supy.organizer.meld(organizers = organizers)
        org = self.orgMelded[(lname,rw,tt)]
        #self.skimStats(org)
        templateSamples = ['top.t#bar{t}','top.W','QCD.multijet']
        baseSamples = ['top.Single','top.DY']

        mfCanvas = r.TCanvas()
        mfFileName = "%s/%s_measuredFractions"%(self.globalStem, lname )
        supy.utils.tCanvasPrintPdf( mfCanvas, mfFileName, option = '[', verbose = False)
        with open(mfFileName+'.txt','w') as file : print >> file, ""
        
        def measureFractions(dist, rebin = 1) :
            before = next(org.indicesOfStep('label','selection complete'))
            distTup = org.steps[next(iter(filter(lambda i: before<i, org.indicesOfStepsWithKey(dist))))][dist]

            templates = [None] * len(templateSamples)
            bases = []
            for ss,hist in zip(org.samples,distTup) :
                contents = supy.utils.binValues(hist)
                if rebin!=1 : contents = contents[:1] + [sum(bins) for bins in zip(*[contents[1:-1][i::rebin] for i in range(rebin)])] + contents[-1:]
                if ss['name'] == "top.Data 2012" :
                    observed = contents
                elif ss['name'] in templateSamples :
                    templates[templateSamples.index(ss['name'])] = contents
                elif ss['name'] in baseSamples :
                    bases.append(contents)
                else : pass

            from supy.utils.fractions import componentSolver,drawComponentSolver
            cs = componentSolver(observed, templates, 1e4, base = np.sum(bases, axis=0) )
            stuff = drawComponentSolver( cs, mfCanvas, distName = dist,
                                         templateNames = [t.replace("top.ttj_%s.wQQ.%s"%(tt,rw),"q#bar{q}-->t#bar{t}").replace("top.ttj_%s.wQG.%s"%(tt,rw),"qg-->t#bar{t}").replace("top.ttj_%s.wAG.%s"%(tt,rw),"#bar{q}g-->t#bar{t}").replace("top.ttj_%s.wGG.%s"%(tt,rw),"gg-->t#bar{t}").replace("QCD.Data 2012","Multijet").replace("top.W","W+jets").replace('top.',"") for t in  templateSamples])
            supy.utils.tCanvasPrintPdf( mfCanvas, mfFileName, verbose = False)
            with open(mfFileName+'.txt','a') as file : print >> file, "\n",dist+"\n", cs
            return distTup,cs

        def mf2(dist) : return measureFractions(dist,1)

        distTup,cs = map(mf2,["KarlsruheDiscriminant","TridiscriminantWTopQCD"][1:])[-1]

        iTT = next(i for i,ss in enumerate(org.samples) if ss['name']=='top.t#bar{t}')
        nTT = distTup[iTT].Integral(0,distTup[iTT].GetNbinsX()+1)
        fractions = dict(zip(templateSamples,cs.fractions))
        for iSample,ss in enumerate(org.samples) :
            if ss['name'] in fractions :
                f = fractions[ss['name']]
                n = distTup[iSample].Integral(0,distTup[iSample].GetNbinsX()+1)
            elif ss['name'] in ['top.ttj_%s.%s.%s'%(tt,s,rw) for s in ['wQQ','wQG','wAG','wGG']] :
                f = fractions['top.t#bar{t}']
                n = nTT
            else : continue
            org.scaleOneRaw(iSample, f * sum(cs.observed) / n )

        org.mergeSamples(targetSpec = {"name":"bg", "color":r.kBlack,"fillColor":r.kGray, "markerStyle":1, "goptions":"hist"}, sources = set(baseSamples + templateSamples) - set(['top.t#bar{t}']), keepSources = True, force = True)
        baseSamples = ['bg']
        templateSamples = ['top.ttj_%s.%s.%s'%(tt,s,rw) for s in ['wGG','wQQ','wQG','wAG']]
        org.mergeSamples(targetSpec = {"name":'qgag'}, sources = templateSamples[2:], keepSources = True, force = True)
        templateSamples = templateSamples[:-2] + ['qgag']
        distTup,cs = map(measureFractions,["fitTopPtOverSumPt","fitTopPtPlusSumPt","fitTopAbsSumRapidities","TridiscriminantGGqqGq"])[-1]
        supy.utils.tCanvasPrintPdf( mfCanvas, mfFileName, option = ']')
        org.drop('qgag')

        templateSamples = ['top.t#bar{t}'] # hack !!
        org.mergeSamples(targetSpec = {"name":"S.M.", "color":r.kGreen+2}, sources = templateSamples + baseSamples , keepSources = True, force = True)
        org.drop('bg')

    def measureAmplitudes(self,rw,tt) :
        lnames = ['mu','el']
        if not all((lname,rw, tt) in self.orgMelded for lname in lnames):
            print "run meldScale() before measureAmplitudes()"
            return
        orgMuEl = [self.orgMelded[(lname,rw,tt)] for lname in lnames]
        omitSamples = ["S.M.","top.Data 2012","top.t#bar{t}"]
        for org in orgMuEl :
            print ", ".join(ss["name"] for ss in org.samples if ss["name"] not in omitSamples)
        canvas = r.TCanvas()
        maFileName = "%s/measuredAmplitudes"%(self.globalStem)
        supy.utils.tCanvasPrintPdf( canvas, maFileName, option = '[', verbose = False)
        with open(maFileName+'.txt','w') as file : print >> file, ""

        def BV(h, rebin=10) :
            contents = supy.utils.binValues(h)
            if rebin!=1 : contents = contents[:1] + [sum(bins) for bins in zip(*[contents[1:-1][i::rebin] for i in range(rebin)])] + contents[-1:]
            return contents

        SA = supy.utils.symmAnti
        def measure(fitvar,genvar, antisamples, sumTemplates = False) :
            steps = [org.steps[ next(org.indicesOfStep('symmAnti','%s in (anti)symm parts of %s'%(fitvar,genvar))) ] for org in orgMuEl]
            templatess = [[ BV( SA(step[fitvar+'_anti'][org.indexOfSampleWithName(sample)])[1]) for sample in antisamples ] for step,org in zip(steps,orgMuEl)]
            basess = [ [ BV( SA(step[fitvar+'_symm'][org.indexOfSampleWithName(sample)])[0]) for sample in antisamples ] +
                       [ BV( step[fitvar][i]) for i,ss in enumerate(org.samples) if ss['name'] not in antisamples+omitSamples ]
                       for step,org in zip(steps,orgMuEl)]
            observeds = [ BV( step[fitvar][org.indexOfSampleWithName("top.Data 2012")] ) for step,org in zip(steps,orgMuEl) ]

            def stack(Lists) : return [sum(lists,[]) for lists in zip(*Lists)]
            templates = stack(templatess)
            bases = stack(basess)
            observed = sum(observeds,[])

            from supy.utils.fractions import componentSolver,drawComponentSolver
            if sumTemplates : templates = [np.sum(templates, axis=0)]
            cs = componentSolver(observed, templates, 1e4, base = np.sum(bases, axis=0) , normalize = False )
            stuff = drawComponentSolver(cs, canvas, distName = fitvar, showDifference = True,
                                        templateNames = ["anti %s-->t#bar{t}"%('q*' if sumTemplates else "qg" if 'QG' in t else "q#bar{q}" if 'QQ' in t else '') for t in antisamples][:1 if sumTemplates else None])
            supy.utils.tCanvasPrintPdf( canvas, maFileName, verbose = False)
            with open(maFileName+'.txt','a') as file : print >> file, "\n",fitvar+"\n", cs
            return steps,cs

        samples = ['top.ttj_%s.%s.%s'%(tt,w,rw) for w in ['wQQ','wQG','wAG']]
        measure('fitTopCosPhiBoost','genTopCosPhiBoost', samples[1:], sumTemplates=True)
        measure('fitTopDeltaBetazRel','genTopDeltaBetazRel', samples, sumTemplates=True)
        #measure('fitTopDeltaBetazRel','genTopDeltaBetazRel', samples)
        measure('fitTopCosThetaBoostAlt','genTopCosThetaBoostAlt', samples, sumTemplates=True)
        #measure('fitTopCosThetaBoostAlt','genTopCosThetaBoostAlt', samples)
        supy.utils.tCanvasPrintPdf( canvas, maFileName, option = ']')
