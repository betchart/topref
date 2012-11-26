from supy.defaults import *
from supy import whereami

def experiment() :
    return "cms"

def mainTree() :
    return ("topRef","tree")

def otherTreesToKeepWhenSkimming() :
    return [("lumiTree","tree")]

def trace() :
    return True

def useCachedFileLists() :
    return True

def cppFiles() :
    return ["cpp/linkdef.cxx"]

def hadd() :
    return ['hadd', whereami()+'/bin/phaddy'][1]

def cppROOTDictionariesToGenerate() :
    return [("vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > >", "vector;Math/LorentzVector.h")]

def detectorSpecs() :
    return {
        "cms": {"barrelEtaMax": 1.4442,
                "endcapEtaMin": 1.5660,
                },
        }

def LorentzVectorType() :
    return ('PtEtaPhiM4D','double')
