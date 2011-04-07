#!/usr/bin/env python

import os,analysis,steps,calculables,samples,organizer,plotter,utils
import ROOT as r

class hadronicLook(analysis.analysis) :
    def parameters(self) :
        objects = {}
        fields =                                                [ "jet",            "met",            "muon",        "electron",        "photon",
                                                                  "compJet",    "compMet",
                                                                  "rechit", "muonsInJets", "jetPtMin", "jetId"]

        objects["caloAK5JetMet_recoLepPhot"] = dict(zip(fields, [("xcak5Jet","Pat"),"metP4AK5TypeII",("muon","Pat"),("electron","Pat"),("photon","Pat"),
                                                                 ("xcak5JetPF","Pat"),"metP4PF",
                                                                 "Calo",     False,         50.0,      "JetIDloose"]))
        
        objects["pfAK5JetMet_recoLepPhot"]   = dict(zip(fields, [("xcak5JetPF","Pat"), "metP4PF",    ("muon","Pat"),("electron","Pat"),("photon","Pat"),
                                                                 ("xcak5Jet","Pat"),"metP4AK5TypeII",
                                                                 "PF",        True,         50.0,      "JetIDtight"]))

        #objects["pfAK5JetMetLep_recoPhot"]   = dict(zip(fields, [("xcak5JetPF","Pat"), "metP4PF",    ("muon","PF"),("electron","PF"), ("photon","Pat"),
        #                                                         None, None,
        #                                                         "PF",        True,         50.0]))

        return { "objects": objects,
                 "nJetsMinMax" :      dict([ ("ge2",(2,None)),  ("2",(2,2)),  ("ge3",(3,None)) ]       [0:1] ),
                 "mcSoup" :           dict([ ("pythia6","py6"), ("pythia8","py8"), ("madgraph","mg") ] [0:1] ),
                 "etRatherThanPt" : [True,False]        [0],
                 "lowPtThreshold" : 30.0,
                 "lowPtName" : "lowPt",
                 #required to be a sorted tuple with length>1
                 #"triggerList" : ("HLT_HT100U","HLT_HT100U_v3","HLT_HT120U","HLT_HT140U","HLT_HT150U_v3"), #2010
                 "triggerList" : ("HLT_HT300_v3", "HLT_HT350_v2"),#early 2011
                 #"triggerList" : ("HLT_HT350_AlphaT0p51_v1", "HLT_HT350_AlphaT0p53_v1"), #mid 2011
                 }

    def calcListJet(self, obj, etRatherThanPt, lowPtThreshold, lowPtName) :
        def calcList(jet, met, photon, muon, electron, muonsInJets, jetPtMin, jetIdFlag) :
            outList = [
                calculables.XClean.xcJet(jet,
                                         applyResidualCorrectionsToData = False,
                                         gamma = photon,
                                         gammaDR = 0.5,
                                         muon = muon,
                                         muonDR = 0.5,
                                         correctForMuons = not muonsInJets,
                                         electron = electron,
                                         electronDR = 0.5),
                calculables.Jet.Indices( jet, jetPtMin, etaMax = 3.0, flagName = jetIdFlag),
                calculables.Jet.Indices( jet, lowPtThreshold, etaMax = 3.0, flagName = jetIdFlag, extraName = lowPtName),
                
                calculables.Jet.SumP4(jet),
                calculables.Jet.SumP4(jet, extraName = lowPtName),
                calculables.Jet.DeltaPhiStar(jet, extraName = lowPtName),
                calculables.Jet.DeltaPseudoJet(jet, etRatherThanPt),
                calculables.Jet.AlphaT(jet, etRatherThanPt),
                calculables.Jet.AlphaTMet(jet, etRatherThanPt, met),
                calculables.Jet.MhtOverMet(jet, met),
                calculables.Jet.deadEcalDR(jet, extraName = lowPtName, minNXtals = 10),
                ]
            return outList+calculables.fromCollections(calculables.Jet, [jet])

        outList = calcList(obj["jet"], obj["met"], obj["photon"], obj["muon"], obj["electron"], obj["muonsInJets"], obj["jetPtMin"], obj["jetId"])
        if obj["compJet"]!=None and obj["compMet"]!=None :
            outList += calcList(obj["compJet"], obj["compMet"], obj["photon"], obj["muon"], obj["electron"], obj["muonsInJets"], obj["jetPtMin"], obj["jetId"])
        return outList

    def calcListOther(self, obj, triggers) :
        return [
            calculables.XClean.IndicesUnmatched(collection = obj["photon"], xcjets = obj["jet"], DR = 0.5),
            calculables.XClean.IndicesUnmatched(collection = obj["electron"], xcjets = obj["jet"], DR = 0.5),

            calculables.Muon.Indices( obj["muon"], ptMin = 10, combinedRelIsoMax = 0.15),
            calculables.Electron.Indices( obj["electron"], ptMin = 10, simpleEleID = "95", useCombinedIso = True),
            calculables.Photon.photonIndicesPat(  ptMin = 25, flagName = "photonIDLooseFromTwikiPat"),
            #calculables.Photon.photonIndicesPat(  ptMin = 25, flagName = "photonIDTightFromTwikiPat"),
            
            calculables.Vertex.ID(),
            calculables.Vertex.Indices(),
            calculables.Other.lowestUnPrescaledTrigger(triggers),
            ]
    
    def listOfCalculables(self, params) :
        obj = params["objects"]
        outList  = calculables.zeroArgs()
        outList += calculables.fromCollections(calculables.Muon, [obj["muon"]])
        outList += calculables.fromCollections(calculables.Electron, [obj["electron"]])
        outList += calculables.fromCollections(calculables.Photon, [obj["photon"]])
        outList += self.calcListOther(obj, params["triggerList"])
        outList += self.calcListJet(obj, params["etRatherThanPt"], params["lowPtThreshold"], params["lowPtName"])
        return outList
    
    def listOfSteps(self, params) :
        _jet  = params["objects"]["jet"]
        _electron = params["objects"]["electron"]
        _muon = params["objects"]["muon"]
        _photon = params["objects"]["photon"]
        _met  = params["objects"]["met"]
        _etRatherThanPt = params["etRatherThanPt"]
        _et = "Et" if _etRatherThanPt else "Pt"

        return [
            steps.Print.progressPrinter(),
            steps.Other.histogrammer("genpthat",200,0,1000,title=";#hat{p_{T}} (GeV);events / bin"),
            steps.Jet.jetPtSelector(_jet, 100.0, 0),
            steps.Jet.jetPtSelector(_jet, 100.0, 1),
            steps.Jet.jetEtaSelector(_jet,2.5,0),
            steps.Trigger.lowestUnPrescaledTrigger(),
            steps.Other.vertexRequirementFilter(),
            steps.Trigger.l1Filter("L1Tech_BPTX_plus_AND_minus.v0"),
            steps.Trigger.physicsDeclared(),
            steps.Other.monsterEventFilter(),
            #steps.Other.cutSorter([
            steps.Other.hbheNoiseFilter(),
            
            steps.Trigger.hltPrescaleHistogrammer(params["triggerList"]),
            #steps.iterHistogrammer("ecalDeadTowerTrigPrimP4", 256, 0.0, 128.0, title=";E_{T} of ECAL TP in each dead region (GeV);TPs / bin", funcString="lambda x:x.Et()"),
            ]+(
            steps.Other.multiplicityPlotFilter("%sIndices%s"%_electron,          nMax = 0, xlabel = "N electrons") +
            steps.Other.multiplicityPlotFilter("%sIndices%s"%_muon,              nMax = 0, xlabel = "N muons") +
            steps.Other.multiplicityPlotFilter("%sIndices%s"%_photon,            nMax = 0, xlabel = "N photons") +
            steps.Other.multiplicityPlotFilter("%sIndicesOther%s"%_jet,          nMax = 0, xlabel = "number of %s%s above p_{T}#semicolon failing ID or #eta"%_jet) +
            steps.Other.multiplicityPlotFilter("%sIndicesOther%s"%_muon,         nMax = 0, xlabel = "number of %s%s above p_{T}#semicolon failing ID or #eta"%_muon) +
            steps.Other.multiplicityPlotFilter("%sIndicesUnmatched%s"%_electron, nMax = 0, xlabel = "N electrons unmatched") +
            steps.Other.multiplicityPlotFilter("%sIndicesUnmatched%s"%_photon,   nMax = 0, xlabel = "N photons unmatched") +
            steps.Other.multiplicityPlotFilter("%sIndices%s"%_jet, nMin=params["nJetsMinMax"][0], nMax=params["nJetsMinMax"][1], xlabel="number of %s%s passing ID#semicolon p_{T}#semicolon #eta cuts"%_jet)
            )+[
            steps.Jet.uniquelyMatchedNonisoMuons(_jet), 

            steps.Other.histogrammer("%sSum%s%s"%(_jet[0], _et, _jet[1]), 50, 0, 1500, title = ";H_{T} (GeV) from %s%s %ss;events / bin"%(_jet[0], _jet[1], _et)),
            steps.Other.variableGreaterFilter(350.0,"%sSum%s%s"%(_jet[0], _et, _jet[1]), suffix = "GeV"),
            
            #steps.Trigger.hltTurnOnHistogrammer( "%sSumEt%s"%_jet,    (60, 200.0, 500.0), "HLT_HT350_v2", ["HLT_HT300_v3"]),
            ##steps.Trigger.hltTurnOnHistogrammer( "%sSumEt%s"%_jet,    (60, 200.0, 500.0), "HLT_HT350_v2", ["HLT_HT250_v2"]),
            ##steps.Trigger.hltTurnOnHistogrammer( "%sSumEt%s"%_jet,    (60, 200.0, 500.0), "HLT_HT350_v2", ["HLT_HT200_v2"]),
            ##steps.Trigger.hltTurnOnHistogrammer( "%sSumEt%s"%_jet,    (60, 200.0, 500.0), "HLT_HT350_v2", ["HLT_HT150_v2"]),
            #
            #steps.Trigger.hltTurnOnHistogrammer( "%sSumEt%s"%_jet,    (60, 150.0, 450.0), "HLT_HT300_v3", ["HLT_HT250_v2"], permissive = True),
            ##steps.Trigger.hltTurnOnHistogrammer( "%sSumEt%s"%_jet,    (60, 200.0, 500.0), "HLT_HT300_v3", ["HLT_HT200_v2"]),
            ##steps.Trigger.hltTurnOnHistogrammer( "%sSumEt%s"%_jet,    (60, 200.0, 500.0), "HLT_HT300_v3", ["HLT_HT150_v2"]),
            #
            #steps.Trigger.hltTurnOnHistogrammer( "%sSumEt%s"%_jet,    (60, 100.0, 400.0), "HLT_HT250_v2", ["HLT_HT200_v2"], permissive = True),
            ##steps.Trigger.hltTurnOnHistogrammer( "%sSumEt%s"%_jet,    (60, 200.0, 500.0), "HLT_HT250_v2", ["HLT_HT150_v2"]),
            #
            #steps.Trigger.hltTurnOnHistogrammer( "%sSumEt%s"%_jet,    (60,  50.0, 350.0), "HLT_HT200_v2", ["HLT_HT150_v2"], permissive = True),
            #
            #steps.Other.variableGreaterFilter(150.0,"%sSum%s%s"%(_jet[0], _et, _jet[1]), suffix = "GeV"),
            #steps.Trigger.hltTurnOnHistogrammer( "%sAlphaTEt%s"%_jet, (70,   0.4,   0.75), "HLT_HT150_AlphaT0p60_v1", ["HLT_HT150_v2"]),
            #steps.Trigger.hltTurnOnHistogrammer( "%sAlphaTEt%s"%_jet, (70,   0.4,   0.75), "HLT_HT150_AlphaT0p70_v1", ["HLT_HT150_v2"]),
            #
            #steps.Other.variableGreaterFilter(200.0,"%sSum%s%s"%(_jet[0], _et, _jet[1]), suffix = "GeV"),
            #steps.Trigger.hltTurnOnHistogrammer( "%sAlphaTEt%s"%_jet, (60,   0.4,   0.7 ), "HLT_HT200_AlphaT0p60_v1", ["HLT_HT200_v2"]),
            #steps.Trigger.hltTurnOnHistogrammer( "%sAlphaTEt%s"%_jet, (60,   0.4,   0.7 ), "HLT_HT200_AlphaT0p65_v1", ["HLT_HT200_v2"]),
            #
            #steps.Other.variableGreaterFilter(250.0,"%sSum%s%s"%(_jet[0], _et, _jet[1]), suffix = "GeV"),
            #steps.Trigger.hltTurnOnHistogrammer( "%sAlphaTEt%s"%_jet, (60,   0.4,   0.7 ), "HLT_HT250_AlphaT0p55_v1", ["HLT_HT250_v2"]),
            #steps.Trigger.hltTurnOnHistogrammer( "%sAlphaTEt%s"%_jet, (60,   0.4,   0.7 ), "HLT_HT250_AlphaT0p62_v1", ["HLT_HT250_v2"]),
            #
            #steps.Other.variableGreaterFilter(300.0,"%sSum%s%s"%(_jet[0], _et, _jet[1]), suffix = "GeV"),
            #steps.Trigger.hltTurnOnHistogrammer( "%sAlphaTEt%s"%_jet, (80,   0.4,   0.6 ), "HLT_HT300_AlphaT0p52_v1", ["HLT_HT300_v3"]),
            #steps.Trigger.hltTurnOnHistogrammer( "%sAlphaTEt%s"%_jet, (80,   0.4,   0.6 ), "HLT_HT300_AlphaT0p54_v1", ["HLT_HT300_v3"]),
            #
            #steps.Other.variableGreaterFilter(350.0,"%sSum%s%s"%(_jet[0], _et, _jet[1]), suffix = "GeV"),
            #steps.Trigger.hltTurnOnHistogrammer( "%sAlphaTEt%s"%_jet, (80,   0.4,   0.6 ), "HLT_HT350_AlphaT0p51_v1", ["HLT_HT350_v2"]),
            #steps.Trigger.hltTurnOnHistogrammer( "%sAlphaTEt%s"%_jet, (80,   0.4,   0.6 ), "HLT_HT350_AlphaT0p53_v1", ["HLT_HT350_v2"]),
            
            
            #many plots
            steps.Trigger.lowestUnPrescaledTriggerHistogrammer(),
            steps.Other.passFilter("singleJetPlots1"),
            steps.Jet.singleJetHistogrammer(_jet),
            steps.Other.passFilter("jetSumPlots1"), 
            steps.Jet.cleanJetHtMhtHistogrammer(_jet,_etRatherThanPt),
            steps.Other.histogrammer(_met,100,0.0,500.0,title=";"+_met+" (GeV);events / bin", funcString = "lambda x: x.pt()"),
            steps.Other.passFilter("kinematicPlots1"), 
            steps.Jet.alphaHistogrammer(cs = _jet, deltaPhiStarExtraName = params["lowPtName"], etRatherThanPt = _etRatherThanPt),
            steps.Jet.alphaMetHistogrammer(cs = _jet, deltaPhiStarExtraName = params["lowPtName"], etRatherThanPt = _etRatherThanPt, metName = _met),
            
            #signal selection
            #steps.Other.variablePtGreaterFilter(140.0,"%sSumP4%s"%_jet,"GeV"),
            steps.Other.variableGreaterFilter(0.55,"%sAlphaT%s%s"%(_jet[0],"Et" if _etRatherThanPt else "Pt",_jet[1])),
            #]), #end cutSorter
            steps.Other.histogrammer("%sMht%sOver%s"%(_jet[0],_jet[1],_met), 100, 0.0, 3.0, title = ";MHT %s%s / %s;events / bin"%(_jet[0],_jet[1],_met)),
            steps.Other.variableLessFilter(1.25,"%sMht%sOver%s"%(_jet[0],_jet[1],_met)),
            steps.Other.deadEcalFilter(jets = _jet, extraName = params["lowPtName"], dR = 0.3, dPhiStarCut = 0.5),
            
            #steps.Other.skimmer(),
            #steps.Other.cutBitHistogrammer(self.togglePfJet(_jet), self.togglePfMet(_met)),
            #steps.Print.eventPrinter(),
            #steps.Print.jetPrinter(_jet),
            #steps.Print.particleP4Printer(_muon),
            #steps.Print.particleP4Printer(_photon),
            #steps.Print.recHitPrinter("clusterPF","Ecal"),
            #steps.Print.htMhtPrinter(_jet),
            #steps.Print.alphaTPrinter(_jet,_etRatherThanPt),
            #steps.Gen.genParticlePrinter(minPt=10.0,minStatus=3),
            #       
            #steps.Other.pickEventSpecMaker(),
            #steps.Displayer.displayer(jets = _jet,
            #                          muons = _muon,
            #                          met       = params["objects"]["met"],
            #                          electrons = params["objects"]["electron"],
            #                          photons   = params["objects"]["photon"],                            
            #                          recHits   = params["objects"]["rechit"],recHitPtThreshold=1.0,#GeV
            #                          scale = 400.0,#GeV
            #                          etRatherThanPt = _etRatherThanPt,
            #                          deltaPhiStarExtraName = params["lowPtName"],
            #                          deltaPhiStarCut = 0.5,
            #                          deltaPhiStarDR = 0.3,
            #                          printOtherJetAlgoQuantities = False,
            #                          jetsOtherAlgo = params["objects"]["compJet"],
            #                          metOtherAlgo  = params["objects"]["compMet"],
            #                          markusMode = False,
            #                          ),
            ]
    
    def listOfSampleDictionaries(self) :
        return [samples.mc, samples.jetmet, samples.signalSkim]

    def listOfSamples(self,params) :
        from samples import specify
        def data() : return specify( #nFilesMax = 4, nEventsMax = 2000,
                                     names = [#"Nov4_MJ_skim","Nov4_J_skim","Nov4_J_skim2","Nov4_JM_skim","Nov4_JMT_skim","Nov4_JMT_skim2",
                                              #"HT.Run2011A-PromptReco-v1.AOD.Bryn",
                                              "HT.Run2011A-PromptReco-v1.AOD.Henning",
                                              ])
        

        def qcd_py6(eL) :
            q6 = [0,5,15,30,50,80,120,170,300,470,600,800,1000,1400,1800]
            iCut = q6.index(80)
            return specify( effectiveLumi = eL, color = r.kBlue,
                            names = [("v14_qcd_py6_pt_%dto%d"%t)[:None if t[1] else -3] for t in zip(q6,q6[1:]+[0])[iCut:]] )

        def g_jets_py6(eL) :
            return specify( effectiveLumi = eL, color = r.kGreen,
                            names = ["v12_g_jets_py6_pt%d"%t for t in [30,80,170]] )

        def qcd_py8(eL) :
            q8 = [0,15,30,50,80,120,170,300,470,600,800,1000,1400,1800]
            iCut = q8.index(50)
            return specify( effectiveLumi = eL, color = r.kBlue,
                            names = [("v14_qcd_py8_pt%dto%d"%t)[:None if t[1] else -3] for t in zip(q8,q8[1:]+[0])[iCut:]] )

        def qcd_mg(eL) :
            qM = ["%d"%t for t in [50,100,250,500,1000]]
            return specify( effectiveLumi = eL, color = r.kBlue,
                            names = ["v12_qcd_mg_ht_%s_%s"%t for t in zip(qM,qM[1:]+["inf"])])

        def g_jets_mg(eL) :
            gM = [40,100,200]
            return specify( effectiveLumi = eL, color = r.kGreen,
                            names = [("v12_g_jets_mg_pt%d_%d")[:None if t[1] else -2] for t in zip(gM,gM[1:]+[0])] )

        def ttbar_mg(eL) :
            return specify( names = "tt_tauola_mg_v12", effectiveLumi = eL, color = r.kOrange)
        
        def ewk(eL) :
            return ( specify(names = "z_inv_mg_v12_skim",  effectiveLumi = eL, color = r.kMagenta ) +
                     specify(names = "z_jets_mg_v12_skim", effectiveLumi = eL, color = r.kYellow-3) +
                     specify(names = "w_jets_mg_v12_skim", effectiveLumi = eL, color = 28         ) )

        def susy(eL) :
            return ( specify(names = "lm0_v12", effectiveLumi = eL, color = r.kRed   ) +
                     specify(names = "lm1_v12", effectiveLumi = eL, color = r.kRed+1 ) )

        ##specify(name = "2010_data_calo_skim",       nFilesMax = -1, color = r.kBlack   , markerStyle = 20),            
        ##specify(name = "2010_data_pf_skim",         nFilesMax = -1, color = r.kBlack   , markerStyle = 20),
        ##specify(name = "test",                      nFilesMax = -1, color = r.kBlack   , markerStyle = 20),
        #
        #caloSkims = [
        #    specify(names = "2010_data_calo_skim",        nFilesMax = -1, color = r.kBlack   , markerStyle = 20),
        #    specify(names = "v12_qcd_py6_pt300_caloSkim", nFilesMax = -1, color = r.kBlue    ),
        #    specify(names = "tt_tauola_mg_v12_caloSkim",  nFilesMax =  3, color = r.kOrange  ),
        #    specify(names = "w_jets_mg_v12_skim_caloSkim",nFilesMax = -1, color = 28         ),
        #    specify(names = "z_inv_mg_v12_skim_caloSkim", nFilesMax = -1, color = r.kMagenta ),
        #    ]
        #pfSkims = [
        #    specify(names = "2010_data_pf_skim",          nFilesMax = -1, color = r.kBlack   , markerStyle = 20),
        #    specify(names = "v12_qcd_py6_pt300_pfSkim",   nFilesMax = -1, color = r.kBlue    ),
        #    specify(names = "tt_tauola_mg_v12_pfSkim",    nFilesMax =  3, color = r.kOrange  ),
        #    specify(names = "w_jets_mg_v12_skim_pfSkim",  nFilesMax = -1, color = 28         ),
        #    specify(names = "z_inv_mg_v12_skim_pfSkim",   nFilesMax = -1, color = r.kMagenta ),
        #    ]
        #
        #outList = []
        #outList = caloSkims
        #self.skimString = "_caloSkim"

        #outList = pfSkims
        #self.skimString = "_pfSkim"

        qcd_func,g_jets_func = {"py6": (qcd_py6,g_jets_py6),
                                "py8": (qcd_py8,g_jets_py6), # no g_jets_py8 available
                                "mg" : (qcd_mg, g_jets_mg ) }[params["mcSoup"]]
        eL = 50 # 1/pb
        #return data()
        return ( data() +
                 qcd_func(eL) + #g_jets_func(eL) +
                 ttbar_mg(eL) + ewk(eL) + susy(eL)
                 )

    def mergeSamples(self, org, tag) :
        def py6(org, smSources) :
            org.mergeSamples(targetSpec = {"name":"qcd_py6", "color":r.kBlue}, allWithPrefix="v14_qcd_py6")
            #org.mergeSamples(targetSpec = {"name":"g_jets_py6_v12", "color":r.kGreen}, allWithPrefix="v12_g_jets_py6")
            smSources.append("qcd_py6")
            #smSources.append("g_jets_py6_v12")

        def py8(org, smSources) :
            org.mergeSamples(targetSpec = {"name":"qcd_py8", "color":r.kBlue}, allWithPrefix="qcd_py8")
            org.mergeSamples(targetSpec = {"name":"g_jets_py6_v12", "color":r.kGreen}, allWithPrefix="v12_g_jets_py6")
            smSources.append("qcd_py8")
            smSources.append("g_jets_py6_v12")

        def mg(org, smSources) :
            org.mergeSamples(targetSpec = {"name":"qcd_mg_v12", "color":r.kBlue}, allWithPrefix="v12_qcd_mg")
            org.mergeSamples(targetSpec = {"name":"g_jets_mg_v12", "color":r.kGreen}, allWithPrefix="v12_g_jets_mg")
            smSources.append("qcd_mg_v12")
            smSources.append("g_jets_mg_v12")

        smSources = ["tt_tauola_mg_v12", "z_inv_mg_v12_skim", "z_jets_mg_v12_skim", "w_jets_mg_v12_skim"]
        for i in range(len(smSources)) :
            smSources[i] = smSources[i]+(self.skimString if hasattr(self,"skimString") else "")
            
        if "pythia6"  in tag : py6(org, smSources)
        if "pythia8"  in tag : py8(org, smSources)
        if "madgraph" in tag : mg (org, smSources)
        org.mergeSamples(targetSpec = {"name":"standard_model", "color":r.kGreen+3}, sources = smSources, keepSources = True)
        #org.mergeSamples(targetSpec = {"name":"2010 Data", "color":r.kBlack, "markerStyle":20}, allWithPrefix="Nov4")
        org.mergeSamples(targetSpec = {"name":"2011 Data", "color":r.kBlack, "markerStyle":20}, allWithPrefix="HT.Run2011A")
        
    def conclude(self) :
        for tag in self.sideBySideAnalysisTags() :
            ##for skimming only
            #org = organizer.organizer( self.sampleSpecs(tag) )
            #utils.printSkimResults(org)            
            
            #organize
            org=organizer.organizer( self.sampleSpecs(tag) )
            self.mergeSamples(org, tag)
            org.scale()
            
            #plot
            pl = plotter.plotter(org,
                                 psFileName = self.psFileName(tag),
                                 samplesForRatios = ("2011 Data","standard_model"),
                                 sampleLabelsForRatios = ("data","s.m."),
                                 #whiteList = ["lowestUnPrescaledTrigger"],
                                 #doLog = False,
                                 #compactOutput = True,
                                 #noSci = True,
                                 #pegMinimum = 0.1,
                                 blackList = ["lumiHisto","xsHisto","nJobsHisto"],
                                 )
            pl.plotAll()

            ##working on manipulation
            #dataIndex = org.indexOfSampleWithName("2010 Data")
            #csIndex = org.indicesOfSelectionsWithKey("cutSorterConfigurationCounts")[0]
            #csTriplet = [org.selections[csIndex][key][dataIndex] for key in ["cutSorterConfigurationCounts","cutSorterNames","cutSorterMoreNames"]]
            #cutSpecs = [(csTriplet[1].GetXaxis().GetBinLabel(i+1),csTriplet[2].GetXaxis().GetBinLabel(i+1)) for i in range(csTriplet[1].GetNbinsX())]

            
            
