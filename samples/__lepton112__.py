from supy.samples import SampleHolder
from supy.sites import srm
lepton112 = SampleHolder()

srm_burt = srm + "/bbetchar/TOP/automated"


loc = "2012_11_20_00_05_22/%sHad.Run2012A-13Jul2012-v1.AOD"
lepton112.add("ElHad.2012A_1", '%s/%s")'%(srm_burt, loc%"Electron"), lumi = 1)
lepton112.add("MuHad.2012A_1", '%s/%s")'%(srm_burt, loc%"Mu"),       lumi = 1)

loc = "2012_11_19_23_44_29//%sHad.Run2012A-recover-06Aug2012-v1.AOD"
lepton112.add("ElHad.2012A_2", '%s/%s")'%(srm_burt, loc%"Electron"), lumi = 1)
lepton112.add("MuHad.2012A_2", '%s/%s")'%(srm_burt, loc%"Mu"),       lumi = 1)


loc = "2012_11_20_00_32_39/Single%s.Run2012B-TOP%sPlusJets-13Jul2012-v1.AOD"
lepton112.add("SingleEl.2012B", '%s/%s")'%(srm_burt, loc%('Electron','Ele')), lumi = 1)
lepton112.add("SingleMu.2012B", '%s/%s")'%(srm_burt, loc%('Mu','Mu')),        lumi = 1)


loc = "2012_11_20_00_38_01/Single%s.Run2012C-TOP%sPlusJets-24Aug2012-v1.AOD"
lepton112.add("SingleEl.2012C", '%s/%s")'%(srm_burt, loc%('Electron','Ele')), lumi = 1)
lepton112.add("SingleMu.2012C", '%s/%s")'%(srm_burt, loc%('Mu','Mu')),        lumi = 1)
