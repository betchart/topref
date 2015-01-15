import supy,steps,samples,calculables

class dupCheck(supy.analysis) :

    def parameters(self):
        return {'lepton':self.vary({'el':{'name':'el'},'mu':{'name':'mu'}})}

    @staticmethod
    def muons(suffix="") :     return [f+suffix for f in ['Mu.A.1','Mu.A.2','Mu.B.1','Mu.C.1','Mu.C.2','Mu.C.3','Mu.D.1']]
    @staticmethod
    def electrons(suffix="") : return [f+suffix for f in ['El.A.1','El.A.2','El.B.1','El.C.1','El.C.2','El.C.3','El.D.1']]
    @staticmethod
    def jsonFiles() :
        return ['lumi/json/'+f for f in ['Cert_190456-196531_8TeV_13Jul2012ReReco_Collisions12_JSON.txt',
                                         'Cert_190782-190949_8TeV_06Aug2012ReReco_Collisions12_JSON.txt',
                                         'Cert_190456-196531_8TeV_13Jul2012ReReco_Collisions12_JSON.txt',
                                         'Cert_198022-198523_8TeV_24Aug2012ReReco_Collisions12_JSON.txt',
                                         'Cert_190456-208686_8TeV_PromptReco-NoReprocessing_Collisions12CD_JSON.txt',
                                         'Cert_201191-201191_8TeV_11Dec2012ReReco-recover_Collisions12_JSON.txt',
                                         'Cert_190456-208686_8TeV_PromptReco-NoReprocessing_Collisions12CD_JSON.txt']]

    
    def listOfSteps(self,pars) :
        return [supy.steps.printer.progressPrinter(2,300),
                steps.other.duplicateEventCheck()]

    def listOfCalculables(self,pars) :
        return supy.calculables.zeroArgs(supy.calculables)

    def listOfSamples(self,pars) :
        return sum( [supy.samples.specify( names = ds, weights = calculables.other.jw(jfn))
                     for ds,jfn in zip({"mu":self.muons(),
                                        "el":self.electrons()}[pars['lepton']['name']],
                                       self.jsonFiles())],[])

    def listOfSampleDictionaries(self) :
        return [samples.lepton119_fnal]


    def conclude(self,pars):
        org = self.organizer(pars, verbose = True )
        supy.plotter(org, pdfFileName = self.pdfFileName(org.tag), doLog=False ).plotAll()
        
