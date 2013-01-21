import configuration
from supy import wrappedChain,utils
##############################
barrelEtaMax = configuration.detectorSpecs()["cms"]["barrelEtaMax"]
endcapEtaMin = configuration.detectorSpecs()["cms"]["endcapEtaMin"]
##############################

class Indices(wrappedChain.calculable) :
    def __init__(self, collection = None, ptMin = 30, absEtaMax = 2.5, ) :
        self.fixes = collection
        for item in ['ptMin','absEtaMax'] : setattr(self,item,eval(item))
        self.stashed = ['P4','SignalID']
        self.stash(self.stashed)
        self.moreName = 'pt>%.1f; |eta|<%.1f; SignalID'%(ptMin, absEtaMax)

    def accept(self, p4, id) : return all([id,
                                           p4.pt() > self.ptMin,
                                           abs(p4.eta()) < self.absEtaMax])

    def update(self,_) :
        self.value = [i for i,a in enumerate(utils.hackMap(self.accept,
                                                           *[self.source[getattr(self,item)] for item in self.stashed]))
                      if a]


class SignalID(wrappedChain.calculable) :
    def __init__(self, collection = None ) :
        self.fixes = collection
        self.stashed = ['Dxy',
                        'SuperClusterEta',
                        'GsfTrackInnerHits',
                        'mvaTrigV0']
        self.stash(self.stashed)

    @staticmethod
    def passID(dxy, scEta, mHits, mva) :
        return all([dxy < 0.02,  mHits <= 0, mva > 0.5,
                    not (barrelEtaMax < abs(scEta) < endcapEtaMin ) ])

    def update(self,ignored) :
        self.value = utils.hackMap(self.passID,
                                   *[self.source[getattr(self,item)] for item in self.stashed])
