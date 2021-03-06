import math,collections, ROOT as r
from supy import utils,analysisStep,calculables
#####################################
class topPrinter(analysisStep) :

    def uponAcceptance (self,ev) :
        ids = list(ev['genPdgId'])
        if not (6 in ids and -6 in ids) : return
        iTop = ids.index(6)
        iTbar= ids.index(-6)
        iQs = range(2,min(iTop,iTbar)) # interested in stuff listed between protons and tops
        mom = ev['genMotherIndex']
        p4s = ev['genP4']
        status = ev['genStatus']
        labels = {1:' d', -1:r'/d', 2:' u', -2:r'/u', 3:' s', -3:r'/s', 4:' c', -4:r'/c', 5:' b', -5:r'/b', 21:' g', 6:' t', -6:'/t'}

        iGs = filter(lambda i: ids[i]==21, range(max(iTop,iTbar)+1,len(ids)))
        iDaughters = filter(lambda i : mom[i]<min(iTop,iTbar), range(max(iTop,iTbar)+1,len(ids)))

        print '-'*50
        print '\t'.join(['item','mom','pt','m','eta','phi'])
        print
        for i in iQs+[iTop,iTbar]+iDaughters :
            print '\t'.join(str(d)[:5] for d in ["%d(%s)"%(i,labels[ids[i]] if ids[i] in labels else str(ids[i])), mom[i], p4s[i].pt(), '-', p4s[i].eta(), p4s[i].phi(), status[i]])
        print

        pD = sum([p4s[i] for i in iDaughters],utils.LorentzV())
        print '\t'.join(str(d)[:5] for d in ["Daught", '-', pD.pt(), pD.M(), pD.eta(), pD.phi()])
        p4 = p4s[2]+p4s[3]
        print '\t'.join(str(d)[:5] for d in ["[%s,%s]"%(2,3), '-', p4.pt(), p4.M(), p4.eta(), p4.phi()])
        p4 = p4s[4]+p4s[5]
        print '\t'.join(str(d)[:5] for d in ["[%s,%s]"%(4,5), '-', p4.pt(), p4.M(), p4.eta(), p4.phi()])
        tt = p4s[iTop] + p4s[iTbar]
        print '\t'.join(str(d)[:5] for d in ["ttbar", '-', tt.pt(), tt.M(), tt.eta(), tt.phi()])
        print
        print
        if abs(tt.E() - (p4s[4]+p4s[5]).E())>0.5 : print (50*' '), "2 -> 3+"
#####################################
class qRecoilKinematics(analysisStep) :
    def hackStep(self,ev):
        if not ev['wQQ'] and ev['genTopTTbar'] : return
        qDir = ev['qDir']
        tdy = math.tanh(ev['genTopDeltaYttbar'])
        qtdy = qDir* tdy
        tday = ev['genTopTanhDeltaAbsY']
        self.book.fill(tday, 'genTopTanhDeltaAbsY', 100, -1, 1, title=';genTopTanhDeltaAbsY;events / bin')
        self.book.fill(tdy,  'genTopTanhDeltaY', 100, -1, 1, title=';genTopTanhDeltaY;events / bin')
        self.book.fill(qtdy,  'qDir_genTopTanhDeltaY', 100, -1, 1, title=';qDir*genTopTanhDeltaY;events / bin')

    def uponAcceptance(self,ev) :
        self.hackStep(ev)
        index = ev["genIndexTtbarExtraJet"]
        if index==None : return
        p4 = ev['genP4'][index]
        id = ev['genPdgId'][index]
        label = ('g' if id==21 else 'q' if id>0 else 'qbar')
        self.book.fill(p4.pt(), label + 'RecoilPt', 100, 0, 200, title = ';recoil %s.pt;events/bin'%label)
        self.book.fill(p4.eta(),label + 'RecoilEta', 50, -5, 5, title = ';recoil %s.eta;events/bin'%label)
#####################################
class pdfWeightsPlotter(calculables.secondary) :
    def __init__(self,vars = [], mins = [], maxs = []) :
        self.vars = vars
        self.mins = mins
        self.maxs = maxs
    def uponAcceptance(self,ev) :
        weights = ev['genPdfWeights']
        N = len(weights)
        for i,w in enumerate(weights):
            self.book.fill(i,'weights',N,0,N, title=";i^{th} weight;sum weights(i)", w=w)
            for var,lo,hi in zip(self.vars,self.mins,self.maxs) :
                self.book.fill((i,ev[var]), '%s_weights'%var, (N,100), (0,lo), (N,hi), w=w, title=';i^{th} weight;%s;sum weights(i)'%var)

    def update(self,_) : pass

    def onlySamples(self) : return ['ttj_ph.wGG.pu.sf','ttj_ph.wQQ.pu.sf','ttj_ph.wQG.pu.sf','ttj_ph.wAG.pu.sf']

    def reportCache(self) : 
        fileName = '/'.join(self.outputFileName.split('/')[:-1]+[self.name])
        histNames = ['weights']+['%s_weights'%v for v in self.vars]
        tops = self.onlySamples()
        cache = self.fromCache(tops,histNames)

        c = r.TCanvas()
        c.Print(fileName+'.pdf[')
        c.Divide(2,2)
        for n in histNames :
            hists = [cache[tt][n] for tt in tops]
            tt = hists[0].Clone(hists[0].GetName()+'_sum')
            for h in hists[1:] : tt.Add(h)
            [h.SetTitle(hn.split('.')[1]) for h,hn in zip(hists,tops)]
            hxs = [utils.divideX(h,tt if n=='weights' else  None) for h in [tt]+hists]
            for i,h in enumerate(hists if n=='weights' else hxs[1:]) :
                c.cd(i+1)
                if n=='weights' :
                    h.Divide(tt)
                    h.GetYaxis().SetTitle("fraction")
                else:
                    h.SetTitle(h.GetTitle() + " shape")
                h.SetMinimum(0)
                h.Draw('colz')
            c.Print(fileName+'.pdf')

            if n=='weights' : continue

            diffs=[]
            noms=[]
            for i,h in enumerate(hxs) :
                c.cd(i)
                hy = h.ProjectionY('_py',1,1)
                diff = utils.subtractX(h,hy,ratherY=True)
                diff.SetMinimum(-diff.GetMaximum())
                diff.SetMaximum(diff.GetMaximum())
                diff.SetTitle(diff.GetTitle() + " difference from nominal")
                if i: diff.Draw('colz')
                diffs.append(diff)
                noms.append(hy)
            c.Print(fileName+'.pdf')

            for i,(d,n) in enumerate(zip(diffs,noms)) :
                if not i: c.Clear()
                c.cd(i)
                d.Multiply(d)
                dy = d.ProjectionY('_py',1,d.GetNbinsX())
                for j in range(2+dy.GetNbinsX()):
                    n.SetBinError(j,math.sqrt(dy.GetBinContent(j)/2))
                n.SetMinimum(0)

                empty=n.GetFillColor()
                n.SetFillColor(r.kRed)
                n.SetLineColor(r.kBlack)
                n.DrawCopy("e4")
                n.SetFillColor(empty)
                n.Draw("hist C same")
                if not i:
                    c.Print(fileName+'.pdf')
                    c.Clear()
                    c.Divide(2,2)
            c.Print(fileName+'.pdf')
            
        c.Print(fileName+'.pdf]')
        print "Wrote %s.pdf"%fileName
