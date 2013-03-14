import supy,steps,samples,calculables

class topSync(supy.analysis) :
    def parameters(self) :
        lepton = {'name' : ['el', 'mu'],
                  'other' : ['mu', 'el'],
                  'relIso' : [0.10, 0.12]}
        lepton['fix'] = [(l,'') for l in lepton['name']]
        
        return {'CSVM':0.679,
                'smear': self.vary({'sm':'Smear',}),#'smU':'SmearUp','smD':'SmearDown'}),
                'lepton': self.vary([(lep,dict((key,val[i]) for key,val in lepton.items()))
                                     for i,lep in enumerate(lepton['name'])])}

    def listOfSampleDictionaries(self) : return [samples.top119]
        
    def listOfSamples(self,pars) : return supy.samples.specify(names = 'ttj_ph', effectiveLumi=10)

    def listOfCalculables(self,pars) :
        jet = ('jet','')
        mu = ('mu','')
        el = ('el','')
        return (supy.calculables.zeroArgs(supy.calculables) +
                [calculables.muon.Indices(mu),
                 calculables.muon.SignalID(mu),
                 calculables.electron.Indices(el),
                 calculables.electron.SignalID(el),
                 calculables.jet.AdjustedP4(jet, pars['smear']),
                 calculables.jet.Pt(jet),
                 calculables.jet.Indices(jet, ptMin = 20),
                 calculables.jet.IndicesBtagged(jet, 'CSV'),
                 supy.calculables.other.abbreviation('combinedSecondaryVertex','CSV',jet)])

    def listOfSteps(self,pars) :
        lep = pars['lepton']
        lepix = lep['fix']
        olepix = (lep['other'],'')
        return [
            supy.steps.printer.progressPrinter(),
            supy.steps.filters.value('mvaTrigV0Exists', min=True),
            supy.steps.filters.value('RelIso'.join(lepix), max=lep['relIso'], indices='Indices'.join(lepix), index=0),
            supy.steps.filters.multiplicity('muP4', max={'mu':1,'el':0}[lep['name']]),
            supy.steps.filters.multiplicity('elP4', max={'mu':0,'el':1}[lep['name']]),
            supy.steps.filters.value('PassConversionVeto'.join(lepix), min=True, indices='Indices'.join(lepix), index=0) if lep['name']=='el' else supy.steps.filters.label('empty'),
            supy.steps.filters.value('jetPt', min = 55, indices = 'jetIndices', index = 0),
            supy.steps.filters.value('jetPt', min = 45, indices = 'jetIndices', index = 1),
            supy.steps.filters.value('jetPt', min = 35, indices = 'jetIndices', index = 2),
            supy.steps.filters.value('jetPt', min = 20, indices = 'jetIndices', index = 3),
            supy.steps.filters.value('jetCSV', min = pars['CSVM'], indices = 'jetIndicesBtagged', index=0),
            steps.other.pickEventSpecMaker()
            ]
