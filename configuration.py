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
    return []

def hadd() :
    return ['hadd', whereami()+'/bin/phaddy'][1]

def cppROOTDictionariesToGenerate() :
    return [("vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > >", "vector;Math/LorentzVector.h")]

def detectorSpecs() :
    return {
        "cms": {"etaBE": 1.479, #from CMS PAS EGM-10-005
                "barrelEtaMax": 1.4442,
                "endcapEtaMin": 1.560,
                "CaloSubdetectors": ["Eb", "Ee", "Hbhe", "Hf"],
                "PFSubdetectors": ["Ecal", "Hcal", "Hfem", "Hfhad", "Ps"],
                "CaloRecHitCollections": [""],
                "PFRecHitCollections": ["","cluster"],
                },
        }
