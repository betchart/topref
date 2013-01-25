import supy,steps,samples,calculables

class jsonMaker(supy.analysis) :
    def parameters(self) :
        jf = ['json/Cert_190456-196531_8TeV_13Jul2012ReReco_Collisions12_JSON.txt',
              'json/Cert_190782-190949_8TeV_06Aug2012ReReco_Collisions12_JSON.txt',
              'json/Cert_190456-196531_8TeV_13Jul2012ReReco_Collisions12_JSON.txt',
              'json/Cert_190456-208686_8TeV_PromptReco-NoReprocessing_Collisions12CD_JSON.txt',
              'json/Cert_198022-198523_8TeV_24Aug2012ReReco_Collisions12_JSON.txt',
              'json/Cert_190456-208686_8TeV_PromptReco-NoReprocessing_Collisions12CD_JSON.txt']

        group = self.vary()

        group['SingleMu'] = zip(["MuHad.2012A_1",
                                 "MuHad.2012A_2",
                                 "SingleMu.2012B",
                                 "SingleMu.2012C",
                                 "SingleMu.2012C_r",
                                 "SingleMu.2012D"
                                 ], jf)

        group['SingleEl'] = zip(["ElHad.2012A_1",
                                 "ElHad.2012A_2",
                                 "SingleEl.2012B",
                                 "SingleEl.2012C",
                                 "SingleEl.2012C_r",
                                 "SingleEl.2012D"
                                 ], jf)

        return {'group':group}

    def listOfSteps(self,pars) :
        return [ supy.steps.printer.progressPrinter(2,300),
                 steps.other.jsonMaker(),
                 ]

    def listOfCalculables(self,pars) :
        return supy.calculables.zeroArgs(supy.calculables)

    def listOfSamples(self,pars) :
        return sum([supy.samples.specify(names = samps, weights = calculables.other.jw(jf)) for samps,jf in pars['group']],[])

    def listOfSampleDictionaries(self) :
        return [samples.lepton118]

    def mainTree(self) :
        return ("lumiTree","tree")

    def otherTreesToKeepWhenSkimming(self) :
        return []
