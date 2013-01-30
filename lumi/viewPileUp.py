import sys,ROOT as r

if not len(sys.argv)>1 :
    print "Usage : viewPileUp.py <jsonsDirectory> after running makePileUp.py"

labels = ['_down','','_up']
suffix = '_pileup.root'

files = [r.TFile.Open(sys.argv[1]+label+suffix) for label in labels]
hists = [f.Get('pileup') for f in files]
colors = [r.kRed,r.kBlack,r.kBlue]
canvas = r.TCanvas()

for i,(h,c) in enumerate(zip(hists,colors)) :
    h.SetLineColor(c)
    h.Draw('histsame' if i else 'hist')

canvas.Print(sys.argv[1]+'.pdf')
