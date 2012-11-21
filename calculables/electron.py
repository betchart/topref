import configuration
from supy import wrappedChain
##############################
barrelEtaMax = configuration.detectorSpecs()["cms"]["barrelEtaMax"]
endcapEtaMin = configuration.detectorSpecs()["cms"]["endcapEtaMin"]
##############################
class SignalID(wrappedChain.calculable) :
    def update(self,ignored) :
        pass
