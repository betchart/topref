import configuration
from supy import wrappedChain,utils
from calculables.other import ScaleFactors
##############################
barrelEtaMax = configuration.detectorSpecs()["cms"]["barrelEtaMax"]
endcapEtaMin = configuration.detectorSpecs()["cms"]["endcapEtaMin"]
##############################


class Indices(wrappedChain.calculable):
    def __init__(self, collection=None, ptMin=30, absEtaMax=2.5, ) :
        self.fixes = collection
        for item in ['ptMin','absEtaMax'] : setattr(self,item,eval(item))
        self.stashed = ['P4','SignalID']
        self.stash(self.stashed)
        self.moreName = 'pt>%.1f; |eta|<%.1f; SignalID' % (ptMin, absEtaMax)

    def accept(self, p4, id) : return all([id,
                                           p4.pt() > self.ptMin,
                                           abs(p4.eta()) < self.absEtaMax])

    def update(self,_) :
        self.value = [i for i,a in enumerate(utils.hackMap(self.accept,
                                                           *[self.source[getattr(self,item)]
                                                             for item in self.stashed]))
                      if a]


class SignalID(wrappedChain.calculable):
    def __init__(self, collection=None):
        self.fixes = collection
        self.stashed = ['Dxy',
                        'SuperClusterEta',
                        'GsfTrackInnerHits',
                        'mvaTrigV0']
        self.stash(self.stashed)

    @staticmethod
    def passID(dxy, scEta, mHits, mva):
        return all([dxy < 0.02,  mHits <= 0, mva > 0.5,
                    not (barrelEtaMax < abs(scEta) < endcapEtaMin ) ])

    def update(self,_) :
        self.value = utils.hackMap(self.passID,
                                   *[self.source[getattr(self,item)] for item in self.stashed])


class TriggerScaleFactors(ScaleFactors):
    '''Scale factors (eff_data / eff_mc) for HLT_Ele27_WP80.

    As found in AN-12-429, Table 12.'''
    def __init__(self, collection=None):
        self.fixes = collection
        self.stash(['P4','Indices'])

        self.rows = [(0,0.800), (0.800,1.478), (1.478,2.5)] # |eta|
        self.columns = [(30,40),(40,50),(50,200)] # pt

        self.central = [[0.987, 0.997, 0.998],
                        [0.964, 0.980, 0.988],
                        [1.004, 1.033, 0.976]]

        self.deltaUp = [[0.012, 0.001, 0.002],
                        [0.002, 0.001, 0.002],
                        [0.006, 0.007, 0.015]]

        self.deltaDn = [[-0.017, -0.001, -0.002],
                        [-0.001, -0.001, -0.002],
                        [-0.006, -0.007, -0.012]]


class SelectionScaleFactors(ScaleFactors):
    '''Scale factors (eff_data / eff_mc) for single electron selection (MVA>0.5).

    As found in AN-12-429, Table 7.'''
    def __init__(self, collection=None):
        self.fixes = collection
        self.stash(['P4','Indices'])

        self.rows = [(0,0.800), (0.800,1.478), (1.478,2.5)] # |eta|
        self.columns = [(20,30),(30,40),(40,50),(50,150)] # pt

        self.central = [[0.962, 0.948, 0.961, 0.962],
                        [0.940, 0.930, 0.965, 0.961],
                        [0.933, 0.924, 0.962, 0.960]]

        self.deltaUp = [[0.005, 0.002, 0.001, 0.001],
                        [0.009, 0.001, 0.002, 0.002],
                        [0.006, 0.003, 0.004, 0.006]]

        self.deltaDn = [[-0.003, -0.001, -0.001, -0.001],
                        [-0.010, -0.001, -0.001, -0.002],
                        [-0.017, -0.002, -0.003, -0.006]]


class SF(wrappedChain.calculable):
    def __init__(self, collection=None):
        self.fixes = collection
        self.stash(['TriggerScaleFactors','SelectionScaleFactors'])

    def update(self,_):
        self.value = self.source[self.TriggerScaleFactors][0] * self.source[self.SelectionScaleFactors][0]


class Reweights(wrappedChain.calculable):
    def __init__(self, collection=None):
        self.fixes = collection
        self.stash(['TriggerScaleFactors','SelectionScaleFactors'])

    def update(self,_):
        tSF = self.source[self.TriggerScaleFactors]
        sSF = self.source[self.SelectionScaleFactors]
        self.value = [sf/tSF[0] for sf in tSF[1:]] + [sf/sSF[0] for sf in sSF[1:]]
