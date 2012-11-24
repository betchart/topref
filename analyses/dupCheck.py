import supy,steps,samples,calculables

class dupCheck(supy.analysis) :
    def parameters(self) :
        group = self.vary()
        group['SingleMu'] = ["MuHad.2012A_1",
                             "MuHad.2012A_2",
                             "SingleMu.2012B",
                             "SingleMu.2012C"]
        group['SingleEl'] = ["ElHad.2012A_1",
                             "ElHad.2012A_2",
                             "SingleEl.2012B",
                             "SingleEl.2012C"]
        return {'group':group}
    
    def listOfSteps(self,pars) :
        return [supy.steps.printer.progressPrinter(2,300),
                steps.other.duplicateEventCheck()]

    def listOfCalculables(self,pars) :
        return supy.calculables.zeroArgs(supy.calculables)

    def listOfSamples(self,pars) :
        return sum([supy.samples.specify(names = samps) for samps in pars['group']],[])

    def listOfSampleDictionaries(self) :
        return [samples.lepton112]
