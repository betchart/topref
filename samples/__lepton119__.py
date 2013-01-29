from supy.samples import SampleHolder
from supy.sites import srm
import collections
lepton119 = SampleHolder()

srm_burt = srm + "/bbetchar/TOP/automated"

# Mu lumi via 'recorded "HLT_IsoMu24_eta2p1_v*", is same as overview recorded (no prescale).
# El lumi via 'recorded "HLT_Ele27_WP80_v*", likewise unprescaled.

elMuLoc = [
    ({"Mu": 809.3793,"El": 809.3793}, "2013_01_25_04_54_06/Single%s.Run2012A-13Jul2012-v1.AOD"),
    ({"Mu":  82.1357,"El":  82.1357}, "2013_01_25_05_04_20/Single%s.Run2012A-recover-06Aug2012-v1.AOD"),
    ({"Mu":4403.7631,"El":4396.6948}, "2013_01_25_04_46_28/Single%s.Run2012B-13Jul2012-v1.AOD"),
    ({"Mu": 490.5101,"El": 481.9663}, "2013_01_25_04_59_14/Single%s.Run2012C-PromptReco-v1.AOD"),
    ({"Mu":6396.9964,"El":6400.9605}, "2013_01_25_04_39_43/Single%s.Run2012C-PromptReco-v2.AOD"),
    ({"Mu": 134.2416,"El": 134.2416}, "2013_01_25_05_08_46/Single%s.Run2012C-EcalRecover_11Dec2012-v1.AOD"),
    ({"Mu":7273.7043,"El":7273.1408}, "2013_01_25_04_27_54/Single%s.Run2012D-PromptReco-v1.AOD")
]

counts = collections.defaultdict(int)
for elmu,loc in elMuLoc :
    for lep,lumi in elmu.items() :
        run = loc[loc.index("2012"):][4:5]
        counts[lep+run]+=1
        lepton119.add("%s.%s.%d"%(lep,run,counts[lep+run]),
                      '%s/%s")'%(srm_burt, loc%({"El":"Electron","Mu":"Mu"}[lep])),
                      lumi = lumi)
