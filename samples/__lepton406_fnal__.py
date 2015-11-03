import supy
import collections

lepton406_fnal = supy.samples.SampleHolder()

eos = 'utils.fileListFromDisk(location="/eos/uscms/store/user/bbetchar/TOP/automated/'

# Mu lumi via 'recorded "HLT_IsoMu24_eta2p1_v*", is same as overview recorded (no prescale).
# El lumi via 'recorded "HLT_Ele27_WP80_v*", likewise unprescaled.

#!!!! lumi not checked for this run (406)

elMuLoc = [
    ({"Mu":( 809.3793, 19785316),"El":( 809.3793,  21229150)}, "2015_10_26_22_00_44/Single%s.Run2012A-22Jan2013-v1.AOD"),
    ({"Mu":(4403.7631, 59538958),"El":(4396.6948,  67940401)}, "2015_10_26_23_04_30/Single%s.Run2012B-22Jan2013-v1.AOD"),
    ({"Mu":(7021.7481, 87683348),"El":(7027.7248, 101883743)}, "2015_10_26_23_15_37/Single%s.Run2012C-22Jan2013-v1.AOD"),
    ({"Mu":(7273.7043, 95078490),"El":(7273.1408, 106555796)}, "2015_10_26_23_31_27/Single%s.Run2012D-22Jan2013-v1.AOD")
]

counts = collections.defaultdict(int)
for elmu,loc in elMuLoc :
    for lep,(lumi,N) in elmu.items() :
        run = loc[loc.index("2012"):][4:5]
        counts[lep+run]+=1
        lepton406_fnal.add("%s.%s.%d"%(lep,run,counts[lep+run]),
                           '%s/%s")'%(eos, loc%({"El":"Electron","Mu":"Mu"}[lep])),
                           lumi = lumi,
                           nCheck = N)
