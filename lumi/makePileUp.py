import sys,os,ROOT as r

# References
# https://twiki.cern.ch/twiki/bin/view/CMS/PileupSystematicErrors
# https://twiki.cern.ch/twiki/bin/view/CMS/PileupJSONFileforData#Calculating_Your_Pileup_Distribu
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupMCReweightingUtilities

if not len(sys.argv)>1 :
    print "Usage: makePileUp.py <jsonsDirectory>"

json = {}
for fname in os.listdir(sys.argv[1]) :
    with open(sys.argv[1]+'/'+fname) as infile:
        json.update(eval(''.join(infile.readlines())))

outJsonName = sys.argv[1].rstrip('/')+'_total.json'
with open(outJsonName,'w') as outfile:
    print >> outfile, str(json).replace("'",'"')

xsInelastic = 69300
systematics = 0.05
currentLumiJSON = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/PileUp/pileup_JSON_DCSONLY_190389-208686_corr.txt'


os.system('which pileupCalc.py')
for label,factor in [('',0),('_down',-systematics),('_up',systematics)] :
    rootName = sys.argv[1].rstrip('/')+label+'_pileup.root'
    print
    print xsInelastic, 1.0+factor
    cmd = ' '.join(["pileupCalc.py",
                    "-i %s"%outJsonName,
                    "--inputLumiJSON %s"%currentLumiJSON,
                    "--calcMode true",
                    "--minBiasXsec %d"%(xsInelastic * (1.0+factor)),
                    "--maxPileupBin 60",
                    "--numPileupBins 120",
                    rootName])
    print cmd
    os.system(cmd)
