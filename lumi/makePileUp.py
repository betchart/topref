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


hists = []
os.system('which pileupCalc.py')
for label,factor in [('',0),('_down',-systematics),('_up',systematics)] :
    rootName = sys.argv[1].rstrip('/')+label+'_pileup.root'
    os.system(' '.join(["pileupCalc.py",
                        "-i %s"%outJsonName,
                        "--inputLumiJSON %s"%currentLumiJSON,
                        "--calcMode true",
                        "--minBiasXsec %f"%(xsInelastic * 1.0+factor),
                        "--maxPileupBin 50",
                        "--numPileupBins 100",
                        rootName])
              )

    tfile = r.TFile.Open(rootName,'READ')
    hists.append( tfile.Get('pileup').Clone(rootName) )
    hists[-1].SetDirectory(0)
    tfile.Close()

canvas = r.TCanvas()
for i,(hist,color) in enumerate(zip(hists,[r.kBlack,r.kBlue,r.kRed])) :
    hist.SetLineColor(color)
    hist.SetLineWidth(3 if not i else 1)
    hist.Draw('hist' if not i else 'histsame')
    print hist.GetMean(), hist.GetRMS()
canvas.Print(sys.argv[1].rstrip('/')+'_pileups.pdf')
