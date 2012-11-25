from supy import wrappedChain,utils

class Indices(wrappedChain.calculable) :
    def __init__(self, collection = None, ptMin = 26, absEtaMax = 2.1, ) :
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
        self.stashed = ['IsPFMuon',       'IsGlobalMuon',  'NormChi2',
                        'TrackingLayers', 'ValidMuonHits', 'Db',#Dxy
                        'Dz',             'ValidPixelHits','MatchedStations']
        self.stash(self.stashed)

    @staticmethod
    def passID(isPF, isGlobal, chi2n, ntrack, nmu, db, dz, npx, nstat ):
        return all([isPF, isGlobal, chi2n < 10, ntrack > 5, nmu > 0,
                    abs(db) < 0.2, dz < 0.5, npx > 0, nstat > 1])

    def update(self,_) :
        self.value = utils.hackMap(self.passID,
                                   *[self.source[getattr(self,item)] for item in self.stashed])
