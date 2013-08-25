from supy import wrappedChain,utils,calculables
import math
import ROOT as r
try: import lhapdf
except: lhapdf=None
##############################
class qPtMin(wrappedChain.calculable) :
    def __init__(self,ptMin) : self.value = ptMin
    def update(self,_) : pass
##############################
class wGG(wrappedChain.calculable) :
    def update(self,_) :
        wOther = [self.source[g] for g in ['wQQ','wQG','wAG']]
        self.value = None if any(wOther) else 1
        if not self.value : assert wOther.count(None)==2
class wQQ(wrappedChain.calculable) :
    def update(self,_) :
        self.value = ( 1 if self.source['genQQbar'] else
                       None if any(self.source[g] for g in ['genGG','wAG','wQG']) else
                       self.resumToQQ() )
    def resumToQQ(self) :
        p4 = self.source['genP4']
        iEx = self.source['genIndexTtbarExtraJet']
        dPq,dPg = [ (p4[i]-p4[iEx]).P() for i in max(self.source[g] for g in ['genQG','genAG'])]
        return 1 if dPg<dPq else None

class w_G(wrappedChain.calculable) :
    def update(self,_) :
        self.value = ( 1 if self.source[self.qtype]
                       and self.source['genP4'][self.source['genIndexTtbarExtraJet']].pt() > self.source['qPtMin']
                       else None )
class wQG(w_G) : qtype = 'genQG'
class wAG(w_G) : qtype = 'genAG'
##############################
class genIndicesHardPartons(wrappedChain.calculable) :
    def __init__(self, ttType ) :
        self.value = {'POWHEG':(4,5),
                      'MC@NLO':(0,1)}[ttType]
    def update(self,_) : pass
##############################
class genIndexTtbarExtraJet(wrappedChain.calculable) :
    def __init__(self, ttType ) : self.ttType = ttType
    def powhegValue(self) : return 8 if self.source['genMotherIndex'][8]==4 else None
    def mcanloValue(self) : return 2 if abs(self.source['genPdgId'][2])!=6 else None
    def update(self,_) :
        self.value = {'POWHEG':self.powhegValue,
                      'MC@NLO':self.mcanloValue}[self.ttType]()
##############################
class genQQbar(wrappedChain.calculable) :
    '''Indices of quark and antiquark in hard collision'''
    def update(self,_) :
        if self.source['isRealData'] : self.value = (); return
        ids = list(self.source['genPdgId'])
        iQQ = tuple(sorted(self.source['genIndicesHardPartons'],key = ids.__getitem__,reverse = True))
        self.value = iQQ if not sum(ids[i] for i in iQQ) else tuple()
##############################
class gen_G(wrappedChain.calculable) :
    '''Indices of (anti)quark and gluon in hard collision.'''
    def update(self,_) :
        if self.source['isRealData'] : self.value = (); return
        ids = self.source['genPdgId']
        iHard = self.source['genIndicesHardPartons']
        iQg = (next((i for i in iHard if self.lo <= ids[i] <= self.up ), None),
               next((i for i in iHard if ids[i]==21), None))
        self.value = iQg if None not in iQg else ()
class genQG(gen_G) : lo,up = 1,6
class genAG(gen_G) : lo,up = -6,-1
class genGG(wrappedChain.calculable) :
    def update(self,_) :
        self.value = self.source['genIndicesHardPartons'] if not any(self.source[g] for g in ['genQQbar','genQG','genAG']) else tuple()
##############################
class qDir(wrappedChain.calculable) :
    def update(self,_) :
        iQ = next(iter(max(self.source[g] for g in ['genQQbar','genQG','genAG'])),None)
        self.value = (1 if self.source['genP4'][iQ].pz() > 0 else -1) if iQ!=None else None
##############################
class genSumP4(wrappedChain.calculable) :
    def update(self,_) :
        genP4 = self.source['genP4']
        iHard = self.source['genIndicesHardPartons']
        self.value = genP4.at(iHard[0]) + genP4.at(iHard[1])
##############################
class genIndicesB(wrappedChain.calculable) :
    def update(self,_) :
        ids = self.source['genPdgId']
        self.value = filter(lambda i: abs(ids[i]) is 5, range(len(ids)))
##############################
class genIndicesWqq(wrappedChain.calculable) :
    def update(self,_) :
        ids = self.source['genPdgId']
        mom = self.source['genMotherPdgId']
        self.value = filter(lambda i: abs(mom[i]) is 24 and abs(ids[i]) < 5, range(len(ids)))
##############################
class genPdfWeights(wrappedChain.calculable) :
    def __init__(self, pdfset, alphasset=None, alphasmems = (3,7), QisMt=True) :
        for item in ['pdfset','alphasset','QisMt','alphasmems'] : setattr(self,item,eval(item))
        self.moreName = pdfset + (alphasset if alphasset else "") + ("; Q=m_t" if QisMt else "")
        lhapdf.initPDFSet(0, pdfset)
        if alphasset:
            lhapdf.initPDFSet(1,alphasset)

    def update(self,_) :
        x1 = self.source['genx1']
        x2 = self.source['genx2']
        id1 = self.source['genid1']
        id2 = self.source['genid2']
        Q = 175.2 if self.QisMt else self.source['genQ']
        xpdfs = [self.xpdf((0,i),x1,x2,id1,id2,Q) for i in range(1+lhapdf.numberPDF(0))]
        xpdfs.extend([self.xpdf((0,0),x1,x2,id1,id2,q) for q in [Q/2,Q*2]])
        if self.alphasset:
            xpdfs.extend([self.xpdf((1,mem),x1,x2,id1,id2,Q) for mem in self.alphasmems])
        self.value = [xpdf/xpdfs[0] for xpdf in xpdfs]

    def xpdf(self, mem, x1, x2, id1, id2, Q) :
        lhapdf.usePDFMember(*mem)
        return lhapdf.xfx(mem[0], x1,Q,id1) * lhapdf.xfx(mem[0], x2,Q,id2)
##############################
class genPtWeights(wrappedChain.calculable):
    '''https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting'''
    def __init__(self, sqrts, channel):
        self.ab = {(7, 'all'):(0.199,-0.00166),
                   (7,'semi'):(0.174,-0.00137),
                   (7,  'di'):(0.222,-0.00197),
                   (8, 'all'):(0.156,-0.00137),
                   (8,'semi'):(0.159,-0.00141),
                   (8,  'di'):(0.148,-0.00129)
                   }[(sqrts,channel)]
        self.moreName = "Top Pt Reweight: %d TeV, %s"%(sqrts,channel)

    def weight(self,pt,ptbar):
        a,b = self.ab
        onent = 0.5 * (a+b*pt) * (a+b*ptbar)
        return math.exp(onent)

    def update(self,_):
        p4 = self.source['genP4']
        iT,iTbar = self.source['genTopTTbar']
        w = self.weight(p4[iT].pt(), p4[iTbar].pt())
        self.value = [1,w,w*w]
##############################
class qDirExpectation(calculables.secondary) :
    var = ""
    limit = 1
    tag = ""
    sample = ""
    p = None

    def onlySamples(self) : return [self.sample]

    def setup(self,*_) :
        import numpy as np
        orig = self.fromCache( [self.sample], [self.var], tag = self.tag)[self.sample][self.var]
        if not orig :
            print "No cache: %s; %s"%(self.sample,str(self))
            self.value = lambda *_ : 1
            return
        values = [orig.GetBinContent(i) for i in range(orig.GetNbinsX()+2)]
        neighbors = 10
        for i in range(neighbors,len(values)-neighbors) : orig.SetBinContent(i, sum(values[i-neighbors:i+neighbors+1]) / (1+2*neighbors))

        edges = utils.edgesRebinned(orig, targetUncRel = 0.1)
        hist = orig.Rebin(len(edges)-1, orig.GetName()+"_rebinned", edges)
        vals  = [hist.GetBinContent(i) for i in range(1,len(edges))]
        del hist
        iZero = edges.index(0)
        R = np.array(vals[iZero:])
        L = np.array(vals[:iZero])[::-1]
        p = (R-L) / ( R + L )

        self.p = r.TH1D(self.name, ";|%s|;|<qDir>|"%self.var, len(edges[iZero:])-1, edges[iZero:])
        for i in range(len(p)) : self.p.SetBinContent(i+1,p[i])
        self.p.SetBinContent(len(edges[iZero:])+2, edges[-1])

        widths = np.array([high-low for low,high in zip(edges[iZero:-1],edges[iZero+1:])])
        q = (R+L)/(widths * (sum(R)+sum(L)))
        self.q = self.p.Clone(self.name+"_pdist")
        self.q.Reset()
        self.q.SetTitle(";|%s|;probability of |%s|"%(self.var,self.var))
        for i in range(len(q)) : self.q.SetBinContent(i+1,q[i])

        self.mean = sum(a*b for a,b in zip(p,(R+L))) / (sum(R)+sum(L))
        self.value = self.calculate

    def reportCache(self) :
        fileName = '/'.join(self.outputFileName.split('/')[:-1]+[self.name])
        self.setup()
        if not self.p : return
        c = r.TCanvas()
        self.p.SetLineWidth(2)
        self.p.SetMaximum(1)
        self.p.SetMinimum(0)
        self.p.Draw('hist')
        self.q.SetLineWidth(2)
        self.q.SetLineColor(r.kRed)
        self.q.Draw('hist same')
        mean = self.q.Clone('mean'); mean.Reset();
        for i in range(mean.GetNbinsX()) : mean.SetBinContent(i+1,self.mean)
        mean.SetTitle("%.2f"%self.mean)
        mean.SetLineColor(r.kBlue)
        mean.Draw('hist same')
        utils.tCanvasPrintPdf(c,fileName)
        del mean
        del c

    def update(self,_) : pass

    def calculate(self, top, tbar) :
        var = self.varFunction(top,tbar)
        p = max(0,self.p.GetBinContent(self.p.FindFixBin(abs(var)))) if self.p else 0
        return p if var > 0 else -p

    def uponAcceptance(self,ev) :
        if ev['isRealData'] : return
        qdir = ev['qDir']
        iTT = ev['genTopTTbar']
        if not iTT or qdir==None : return
        p4 = ev['genP4']
        var = self.varFunction(p4[iTT[0]],p4[iTT[1]])
        self.book.fill(qdir*var, self.var, 1000, -self.limit, self.limit, title = ";qdir * %s;events / bin"%self.var )

    def varFunction(self,top,tbar) : return


class qDirExpectation_(qDirExpectation) :
    def __init__(self, var, limit, tag, sample) :
        for item in ['var','tag','sample', 'limit'] : setattr(self,item,eval(item))
        self.fixes = ('',var)

    def varFunction(self,top,tbar) : return self.source[self.var]

class qDirExpectation_SumRapidities(qDirExpectation) :
    def varFunction(self,top,tbar) : return top.Rapidity() + tbar.Rapidity()
    def __init__(self, tag, sample) :
        self.var = "SumRapidities"
        self.limit = 5
        self.tag = tag
        self.sample = sample

class qDirExpectation_EtaSum(qDirExpectation) :
    def varFunction(self,top,tbar) : return (top + tbar).Eta()
    def __init__(self, tag, sample) :
        self.var = "EtaSum"
        self.limit = 8
        self.tag = tag
        self.sample = sample

class qDirExpectation_RapiditySum(qDirExpectation) :
    def varFunction(self,top,tbar) : return (top + tbar).Rapidity()
    def __init__(self, tag, sample) :
        self.var = "RapiditySum"
        self.limit = 2.5
        self.tag = tag
        self.sample = sample

