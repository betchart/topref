from supy.samples import SampleHolder
from supy.sites import srm
lepton118 = SampleHolder()

srm_burt = srm + "/bbetchar/TOP/automated"


loc = "2013_01_15_22_16_08/%sHad.Run2012A-13Jul2012-v1.AOD"
lepton118.add("ElHad.2012A_1", '%s/%s")'%(srm_burt, loc%"Electron"), lumi = 809.3793 )
lepton118.add("MuHad.2012A_1", '%s/%s")'%(srm_burt, loc%"Mu"),       lumi = 809.3793 )

loc = "2013_01_15_22_22_12/%sHad.Run2012A-recover-06Aug2012-v1.AOD"
lepton118.add("ElHad.2012A_2", '%s/%s")'%(srm_burt, loc%"Electron"), lumi = 82.0896 )
lepton118.add("MuHad.2012A_2", '%s/%s")'%(srm_burt, loc%"Mu"),       lumi = 82.1357 )

loc = "2013_01_15_22_26_19/Single%s.Run2012B-13Jul2012-v1.AOD"
lepton118.add("SingleEl.2012B", '%s/%s")'%(srm_burt, loc%'Electron'),                         lumi = 4286.2901 )
lepton118.add("SingleMu.2012B", '%s/%s", alwaysUseLastAttempt = True)'%(srm_burt, loc%'Mu'),  lumi = 4403.7631 )

loc = "2013_01_15_22_31_59/Single%s.Run2012C-TOP%sPlusJets-19Dec2012-v1.AOD"
lepton118.add("SingleEl.2012C", '%s/%s")'%(srm_burt, loc%('Electron','Ele')), lumi = 5113.1059 )
lepton118.add("SingleMu.2012C", '%s/%s")'%(srm_burt, loc%('Mu','Mu')),        lumi = 5105.2774 )

loc = "2013_01_22_05_46_25//Single%s.Run2012C-TOP%sPlusJets-24Aug2012-v1.AOD"
lepton118.add("SingleEl.2012C_r", '%s/%s")'%(srm_burt, loc%('Electron','Ele')), lumi = 495.0029 )
lepton118.add("SingleMu.2012C_r", '%s/%s")'%(srm_burt, loc%('Mu','Mu')),        lumi = 495.0029 )

loc = "2013_01_16_17_48_02/Single%s.Run2012D-TOP%sPlusJets-19Dec2012-v1.AOD"
lepton118.add("SingleEl.2012D", '%s/%s")'%(srm_burt, loc%('Electron','Ele')), lumi = 7273.7043 )
lepton118.add("SingleMu.2012D", '%s/%s")'%(srm_burt, loc%('Mu','Mu')),        lumi = 7273.1408 )
