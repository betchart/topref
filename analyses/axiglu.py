import supy,steps,samples,calculables
import ROOT as r

class axiglu(supy.analysis):
    def parameters(self):
        return {}
        
    def listOfSampleDictionaries(self): return [samples.top119_fnal]
        
    def listOfSamples(self,pars):
        suf = ['L200','A200','R200','A2K','R2K','ZP']
        col = [r.kGreen, r.kBlue, r.kRed, r.kCyan, r.kOrange, r.kBlack]
        return sum([supy.samples.specify(names = 'calib_'+s, color = c)
                    for s,c in zip(suf,col)], [])

    def listOfCalculables(self,pars):
        return sum( [supy.calculables.zeroArgs(module) for module in [calculables, supy.calculables]], [])

    def listOfSteps(self,pars):
        return [
            supy.steps.printer.progressPrinter(),
            supy.steps.histos.mass('genTopTTbarP4',100, 350, 1000)
            ]

    def conclude(self,pars) :
        org = self.organizer(pars)
        #org.mergeSamples(targetSpec = {"name":"SingleMu", "color":r.kBlack}, allWithPrefix="SingleMu")
        #org.scale()
        
        supy.plotter(org,
                     pdfFileName = self.pdfFileName(org.tag),
                     #samplesForRatios = ("2010 Data","standard_model"),
                     #sampleLabelsForRatios = ("data","s.m."),
                     #whiteList = ["lowestUnPrescaledTrigger"],
                     #doLog = False,
                     #compactOutput = True,
                     #noSci = True,
                     #pegMinimum = 0.1,
                     blackList = ["lumiHisto","xsHisto","nJobsHisto"],
                     ).plotAll()
