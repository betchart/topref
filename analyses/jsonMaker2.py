import supy,steps,samples,calculables

class jsonMaker2(supy.analysis) :
    def parameters(self) :
        jf = [
            'lumi/json/A_JSON.txt',
            'lumi/json/B_JSON.txt',
            'lumi/json/C_JSON.txt',
            'lumi/json/D_JSON.txt',
            ]

        group = self.vary()

        group['SingleMu'] = zip(["Mu.A.1",
                                 "Mu.B.1",
                                 "Mu.C.1",
                                 "Mu.D.1"
                                 ], jf)

        group['SingleEl'] = zip(["El.A.1",
                                 "El.B.1",
                                 "El.C.1",
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
        return [samples.lepton406_fnal]

    def mainTree(self) :
        return ("lumiTree","tree")

    def otherTreesToKeepWhenSkimming(self) :
        return []
