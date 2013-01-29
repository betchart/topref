import supy,steps,samples,calculables

class jsonMaker(supy.analysis) :
    def parameters(self) :
        jf = [
            'json/Cert_190456-196531_8TeV_13Jul2012ReReco_Collisions12_JSON.txt',             #A1
            'json/Cert_190782-190949_8TeV_06Aug2012ReReco_Collisions12_JSON.txt',             #A2
            'json/Cert_190456-196531_8TeV_13Jul2012ReReco_Collisions12_JSON.txt',             #B1
            'json/Cert_198022-198523_8TeV_24Aug2012ReReco_Collisions12_JSON.txt',             #C1
            'json/Cert_190456-208686_8TeV_PromptReco-NoReprocessing_Collisions12CD_JSON.txt', #C2
            'json/Cert_201191-201191_8TeV_11Dec2012ReReco-recover_Collisions12_JSON.txt',     #C3
            'json/Cert_190456-208686_8TeV_PromptReco-NoReprocessing_Collisions12CD_JSON.txt', #D1
            ]

        group = self.vary()

        group['SingleMu'] = zip(["Mu.A.1",
                                 "Mu.A.2",
                                 "Mu.B.1",
                                 "Mu.C.1",
                                 "Mu.C.2",
                                 "Mu.C.3",
                                 "Mu.D.1"
                                 ], jf)

        group['SingleEl'] = zip(["El.A.1",
                                 "El.A.2",
                                 "El.B.1",
                                 "El.C.1",
                                 "El.C.2",
                                 "El.C.3",
                                 "El.D.1"
                                 ], jf)

        return {'group':group}

    def listOfSteps(self,pars) :
        return [ supy.steps.printer.progressPrinter(2,300),
                 steps.other.jsonMaker(),
                 ]

    def listOfCalculables(self,pars) :
        return supy.calculables.zeroArgs(supy.calculables)

    def listOfSamples(self,pars) :
        return sum([supy.samples.specify(names = samps, weights = (calculables.other.jw(jf) if jf else [])) for samps,jf in pars['group']],[])

    def listOfSampleDictionaries(self) :
        return [samples.lepton119]

    def mainTree(self) :
        return ("lumiTree","tree")

    def otherTreesToKeepWhenSkimming(self) :
        return []
