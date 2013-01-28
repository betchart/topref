from supy.samples import SampleHolder
from supy.sites import srm
import collections
lepton119 = SampleHolder()

srm_burt = srm + "/bbetchar/TOP/automated"

elMuLoc = [
    ({"Mu":9999,"El":9999}, "2013_01_25_04_54_06/Single%s.Run2012A-13Jul2012-v1.AOD"),
    ({"Mu":9999,"El":9999}, "2013_01_25_05_04_20/Single%s.Run2012A-recover-06Aug2012-v1.AOD"),
    ({"Mu":9999,"El":9999}, "2013_01_25_05_13_48/Single%s.Run2012A-20Nov2012-v2.AOD"),
    ({"Mu":9999,"El":9999}, "2013_01_25_04_46_28/Single%s.Run2012B-13Jul2012-v1.AOD"),
    ({"Mu":9999,"El":9999}, "2013_01_25_05_20_58/Single%s.Run2012B-20Nov2012-v2.AOD"),
    ({"Mu":9999,"El":9999}, "2013_01_25_04_59_14/Single%s.Run2012C-PromptReco-v1.AOD"),
    ({"Mu":9999,"El":9999}, "2013_01_25_04_39_43/Single%s.Run2012C-PromptReco-v2.AOD"),
    ({"Mu":9999,"El":9999}, "2013_01_25_05_08_46/Single%s.Run2012C-EcalRecover_11Dec2012-v1.AOD"),
    ({"Mu":9999,"El":9999}, "2013_01_25_04_27_54/Single%s.Run2012D-PromptReco-v1.AOD")
]

counts = collections.defaultdict(int)
for elmu,loc in elMuLoc :
    for lep,lumi in elmu.items() :
        run = loc[loc.index("2012"):][:5]
        counts[lep+run]+=1
        lepton119.add("Single%s.%s.%d"%(lep,run,counts[lep+run]),
                      '%s/%s")'%(srm_burt, loc%({"El":"Electron","Mu":"Mu"}[lep])),
                      lumi = lumi)
