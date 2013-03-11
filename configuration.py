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
    return [("vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<Double32_t> > >", "vector;Math/LorentzVector.h")]

def detectorSpecs() :
    return {
        "cms": {"barrelEtaMax": 1.4442,
                "endcapEtaMin": 1.5660,
                "etaBE": 1.479, #from CMS PAS EGM-10-005
                },
        }

def LorentzVectorType() :
    return ('PtEtaPhiM4D','Double32_t')

def computeEntriesForReport() :
    return False
