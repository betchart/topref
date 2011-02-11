#!/usr/bin/env python

import os,analysis,steps,calculables,samples,organizer,plotter,utils
import ROOT as r

class topAsymm(analysis.analysis) :
    def parameters(self) :
        objects = {}
        fields =                           [ "jet",            "met",            "muon",
                                             "electron",        "photon",
                                             "rechit", "muonsInJets"]
        
        objects["calo"] = dict(zip(fields, [("xcak5Jet","Pat"),"metP4AK5TypeII",("muon","Pat"),
                                            ("electron","Pat"),("photon","Pat"),
                                            "Calo",     False]))
        
        return { "objects": objects,
                 "lepton" : "muon",
                 "nJets" :  [{"min":3,"max":None}],
                 "triggerList" : ("HLT_Mu9","HLT_Mu15_v1"),#required to be a sorted tuple
                 }

    def listOfCalculables(self, params) :
        obj = params["objects"]
        lepton = obj[params["lepton"]]
        outList  = calculables.zeroArgs()
        outList += calculables.fromCollections(calculables.Muon, [obj["muon"]])
        outList += calculables.fromCollections(calculables.Electron, [obj["electron"]])
        outList += calculables.fromCollections(calculables.Photon, [obj["photon"]])
        outList += calculables.fromCollections(calculables.Jet, [obj["jet"]])
        outList += [
            calculables.Jet.IndicesBtagged(obj["jet"],"TrkCountingHighEffBJetTags",2.3), #Best guess at TCHEL
            calculables.Jet.Indices(      obj["jet"],      ptMin = 30, etaMax = 3.0, flagName = "JetIDloose"),
            calculables.Muon.Indices(     obj["muon"],     ptMin = 10, combinedRelIsoMax = 0.15),
            calculables.Electron.Indices( obj["electron"], ptMin = 10, simpleEleID = "95", useCombinedIso = True),
            calculables.Photon.photonIndicesPat(           ptMin = 25, flagName = "photonIDLooseFromTwikiPat"),

            calculables.Jet.SumP4(obj["jet"]),
            calculables.Jet.MhtOverMet(obj["jet"], obj["met"]),
            calculables.Jet.deadEcalDR(obj["jet"], minNXtals = 10),

            calculables.XClean.IndicesUnmatched(collection = obj["photon"], xcjets = obj["jet"], DR = 0.5),
            calculables.XClean.IndicesUnmatched(collection = obj["electron"], xcjets = obj["jet"], DR = 0.5),
            calculables.XClean.xcJet(obj["jet"], applyResidualCorrectionsToData = True,
                                     gamma    = obj["photon"],      gammaDR = 0.5,
                                     electron = obj["electron"], electronDR = 0.5,
                                     muon     = obj["muon"],         muonDR = 0.5, correctForMuons = not obj["muonsInJets"]),
            calculables.XClean.SumP4(obj["jet"], obj["photon"], obj["electron"], obj["muon"]),

            calculables.Vertex.ID(),
            calculables.Vertex.Indices(),
            calculables.Other.lowestUnPrescaledTrigger(params["triggerList"]),

            calculables.Other.SemileptonicTopIndex(lepton),
            calculables.Other.NeutrinoPz(lepton,"xcSumP4"),
            calculables.Other.NeutrinoP4P(lepton,"xcSumP4"),
            calculables.Other.NeutrinoP4M(lepton,"xcSumP4"),
            calculables.Other.SumP4NuP(lepton,"xcSumP4"),
            calculables.Other.SumP4NuM(lepton,"xcSumP4"),
            calculables.Other.SignedRapidity(lepton,"xcSumP4"),
            calculables.Other.RelativeRapidity(lepton,"xcSumP4"),
            calculables.Other.RelativeRapidity(lepton,"xcSumP4NuP"),
            calculables.Other.RelativeRapidity(lepton,"xcSumP4NuM"),

            calculables.Other.Mt(lepton,"%sNeutrinoP4P%s"%lepton),
            ]
        return outList
    
    def listOfSteps(self, pars) :
        obj = pars["objects"]
        _jet  = obj["jet"]
        _electron = obj["electron"]
        _muon = obj["muon"]
        _photon = obj["photon"]
        _met  = obj["met"]
        lepton = obj[pars["lepton"]]
        lPtMin = 20

        return [
            steps.Print.progressPrinter(),
            steps.Other.histogrammer("genpthat",200,0,1000,title=";#hat{p_{T}} (GeV);events / bin"),
            steps.Filter.pt("%sP4%s"%lepton, min = 20, indices = "%sIndicesAnyIso%s"%lepton, index = 0),
            steps.Filter.multiplicity("vertexIndices",min=1),
            steps.Filter.monster(),
            steps.Filter.hbheNoise(),
            steps.Trigger.techBitFilter([0],True),
            steps.Trigger.physicsDeclared(),            
            steps.Trigger.lowestUnPrescaledTrigger(),
            ]+[
            steps.Filter.multiplicity(s, max = 0) for s in ["%sIndicesOther%s"%_jet,
                                                            "%sIndicesOther%s"%_muon,
                                                            "%sIndicesUnmatched%s"%_electron,
                                                            "%sIndicesUnmatched%s"%_photon,
                                                            "%sIndices%s"%_electron,
                                                            "%sIndices%s"%_photon]
            ]+[
            steps.Jet.uniquelyMatchedNonisoMuons(_jet),

            steps.Histos.generic(("xcSumP4","metP4PF"),(50,50),(0,0),(500,500),
                                 title=";xcSumP4.pt; pfmet;events",funcString="lambda x: (x[0].pt(),x[1].pt())"),
            steps.Histos.pt("xcSumP4",50,0,500),
            steps.Histos.pt("metP4PF",50,0,500),
            steps.Filter.pt("xcSumP4",min=20),

            steps.Histos.multiplicity("%sIndices%s"%_jet),
            steps.Filter.multiplicity("%sIndices%s"%_jet, min=pars["nJets"]["min"], max=pars["nJets"]["max"]),

            steps.Histos.value("%sCombinedRelativeIso%s"%lepton, 50,0,1, indices = "%sIndicesAnyIso%s"%lepton),
            steps.Histos.multiplicity("%sIndices%s"%lepton),
            steps.Histos.multiplicity("%sIndicesNonIso%s"%lepton),

            #invert isolation requirement here
            #steps.Filter.multiplicity("%sIndices%s"%lepton, max = 0),
            steps.Filter.multiplicity("%sIndices%s"%lepton, min = 1, max = 1),
            steps.Filter.pt("%sP4%s"%lepton, min = lPtMin, indices = "%sIndices%s"%lepton, index = 0),

            steps.Histos.value("%sTrkCountingHighEffBJetTags%s"%calculables.Jet.xcStrip(_jet), 60,0,30, indices = "%sIndices%s"%_jet),
            steps.Histos.multiplicity("%sIndicesBtagged%s"%_jet),
            steps.Filter.multiplicity("%sIndicesBtagged%s"%_jet, min = 2),
            #steps.Filter.multiplicity("%sIndicesBtagged%s"%_jet, max = 0),
            
            steps.Histos.pt("xcSumP4",50,0,500),
            steps.Histos.pz("xcSumP4",50,-1500,1500),
            steps.Histos.generic(("xcSumP4NuM","xcSumP4NuP"),(75,75),(-1500,-1500),(1500,1500),
                                 funcString = "lambda x:(x[0].pz(),x[1].pz())", title = ";xcSumP4NuM.pz;xcSumP4NuP.pz;events / bin"),

            steps.Histos.value("%sSignedRapidity%s"%lepton, 31,-5,5),
            steps.Histos.value("%s%s"%lepton+"RelativeRapidityxcSumP4", 31,-5,5, xtitle = "#Delta y"),
            steps.Histos.value("%sMt%s"%lepton+"%sNeutrinoP4P%s"%lepton, 30,0,180, xtitle = "M_{T}"),

            steps.Histos.generic(("%s%s"%lepton+"RelativeRapidityxcSumP4NuM","%s%s"%lepton+"RelativeRapidityxcSumP4NuP"),
                                 (101,101), (-5,-5), (5,5), title = ";#Delta y #nu_{-};#Delta y #nu_{+};events / bin",
                                 funcString = "lambda x: (x[0],x[1])"),

            steps.Histos.generic(("%sNeutrinoPz%s"%lepton,"%sNeutrinoPz%s"%lepton),(100,100),(-1500,-500),(500,1500),
                                 title=";#nu_{-} p_{z} (GeV);#nu_{+} p_{z} (GeV);events / bin", funcString = "lambda x: (x[0][0],x[0][1])"),
            
            steps.Filter.multiplicity("%sIndices%s"%_jet, min=4, max=4),
            steps.Filter.multiplicity("%sIndicesBtagged%s"%_jet, max=2),

            steps.Histos.generic(("%sIndices%s"%_jet,"%sIndicesBtagged%s"%_jet,"%sCorrectedP4%s"%_jet),
                                 30,0,180, title=";M_{2-light};events / bin",
                                 funcString="lambda x: sum([x[2][i] for i in (set(x[0])-set(x[1]))],r.LorentzV()).M()")

            #steps.Other.productGreaterFilter(0,["%s%s"%lepton+"RelativeRapidityxcSumP4NuM","%s%s"%lepton+"RelativeRapidityxcSumP4NuP"]),
            #steps.Other.histogrammer("%sSignedRapidity%s"%lepton, 51, -5, 5, title = ";y_lep*q_lep*sign(boost);events / bin"),
            #steps.Other.histogrammer("%s%s"%lepton+"RelativeRapidityxcSumP4", 51, -5, 5, title = ";#Delta y;events / bin"),
            #steps.Other.histogrammer(("%s%s"%lepton+"RelativeRapidityxcSumP4NuM","%s%s"%lepton+"RelativeRapidityxcSumP4NuP"),
            #                         (101,101), (-5,-5), (5,5), title = ";#Delta y #nu_{-};#Delta y #nu_{+};events / bin",
            #                         funcString = "lambda x: (x[0],x[1])"),
            #steps.Other.histogrammer("%sMt%s"%lepton+"%sNeutrinoP4P%s"%lepton,50,0,200, title = ";M_{T};events / bin"),
            #
            #steps.Other.histogrammer(("xcSumP4NuM","xcSumP4NuP"),(75,75),(-1500,-1500),(1500,1500),
            #                         funcString = "lambda x:(x[0].pz(),x[1].pz())", title = ";xcSumP4NuM.pz;xcSumP4NuP.pz;events / bin")
            
            ]
    
    def listOfSampleDictionaries(self) :
        return [samples.mc, samples.muon]

    def listOfSamples(self,pars) :
        from samples import specify
        def data() :
            return specify( #nFilesMax = 1, nEventsMax = 40000,
                            names = [#"Mu.Run2010A-Nov4ReReco.RECO.Jad",
                                     #"Mu.Run2010B-Nov4ReReco.RECO.Jad",
                                     "Mu.Run2010A_skim",
                                     "Mu.Run2010B_skim"
                                     ])
        def qcd_py6(eL) :
            q6 = [0,5,15,30,50,80,120,170,300,470,600,800,1000,1400,1800]
            iCut = q6.index(50)
            return specify( effectiveLumi = eL, color = r.kBlue,
                            names = [("qcd_py6_pt_%dto%d"%t)[:None if t[1] else -3] for t in zip(q6,q6[1:]+[0])[iCut:]] )
        def qcd_py8(eL) :
            q8 = [0,15,30,50,80,120,170,300,470,600,800,1000,1400,1800]
            iCut = q8.index(50)
            return specify( effectiveLumi = eL, color = r.kBlue,
                            names = [("qcd_py8_pt%dto%d"%t)[:None if t[1] else -3] for t in zip(q8,q8[1:]+[0])[iCut:]] )
        def qcd_mg(eL) :
            qM = ["%d"%t for t in [50,100,250,500,1000]]
            return specify( effectiveLumi = eL, color = r.kBlue,
                            names = ["v12_qcd_mg_ht_%s_%s"%t for t in zip(qM,qM[1:]+["inf"])])
        def ttbar_mg(eL) :
            return specify( names = "tt_tauola_mg_v12", effectiveLumi = eL, color = r.kOrange)

        def ewk(eL) :
            return ( specify( names = "w_munu", effectiveLumi = eL, color = 28)  +
                     specify( names = "w_taunu", effectiveLumi = eL, color = r.kYellow-3) )

        eL = 400 # 1/pb
        return  ( data() +
                  qcd_py6(eL) +
                  ttbar_mg(eL) +
                  ewk(eL)[0:1]
                  )

    def conclude(self) :
        for tag in self.sideBySideAnalysisTags() :
            #organize
            org=organizer.organizer( self.sampleSpecs(tag) )
            org.mergeSamples(targetSpec = {"name":"Mu 2010", "color":r.kBlack, "markerStyle":20}, allWithPrefix="Mu.Run2010")
            org.mergeSamples(targetSpec = {"name":"qcd_py6", "color":r.kBlue}, allWithPrefix="qcd_py6")
            org.mergeSamples(targetSpec = {"name":"standard_model", "color":r.kGreen+3},
                             sources = ["qcd_py6","w_munu","w_taunu","tt_tauola_mg_v12"], keepSources = True)
            org.scale()
            
            #plot
            pl = plotter.plotter(org,
                                 psFileName = self.psFileName(tag),
                                 samplesForRatios = ("Mu 2010","standard_model"),
                                 sampleLabelsForRatios = ("data","s.m."),
                                 #whiteList = ["lowestUnPrescaledTrigger"],
                                 #doLog = False,
                                 #compactOutput = True,
                                 #noSci = True,
                                 #pegMinimum = 0.1,
                                 blackList = ["lumiHisto","xsHisto","nJobsHisto"],
                                 )
            pl.plotAll()
