from supy import wrappedChain,utils,calculables
import ROOT as r
##############################
class wGG(wrappedChain.calculable) :
    def update(self,_) :
        self.value = None if any(self.source[g] for g in ['genQG','genQQbar','genAG']) else 1
##############################
class wQQ(wrappedChain.calculable) :
    def update(self,_) :
        self.value = 1 if self.source['genQQbar'] else None
##############################
class wQG(wrappedChain.calculable) :
    def update(self,_) :
        self.value = 1 if self.source['genQG'] else None
##############################
class wAG(wrappedChain.calculable) :
    def update(self,_) :
        self.value = 1 if self.source['genAG'] else None
##############################
class genIndicesHardPartons(wrappedChain.calculable) :
    def __init__(self,indices = (4,5)) : self.value = indices
    def update(self,_) : pass
##############################
class genQQbar(wrappedChain.calculable) :
    '''Index of quark and of antiquark with hard collision'''
    def update(self,_) :
        if self.source['isRealData'] : self.value = (); return
        ids = list(self.source['genPdgId'])
        iHard = self.source['genIndicesHardPartons']
        self.value = tuple(sorted(iHard,key = ids.__getitem__,reverse = True)) \
                     if not sum([ids[i] for i in iHard]) else tuple()
##############################
class genQG(wrappedChain.calculable) :
    '''Index of (anti)quark and of gluon with hard collision.'''
    def update(self,_) :
        if self.source['isRealData'] : self.value = (); return
        ids = self.source['genPdgId']
        iHard = self.source['genIndicesHardPartons']
        iQg = (next((i for i in iHard if ids[i] in range(1,7)), None),
               next((i for i in iHard if ids[i]==21), None))
        self.value = iQg if None not in iQg else ()
##############################
class genAG(wrappedChain.calculable) :
    '''Index of (anti)quark and of gluon with hard collision.'''
    def update(self,_) :
        if self.source['isRealData'] : self.value = (); return
        ids = self.source['genPdgId']
        iHard = self.source['genIndicesHardPartons']
        iAg = (next((i for i in iHard if ids[i] in range(-6,0)), None),
               next((i for i in iHard if ids[i]==21), None))
        self.value = iAg if None not in iAg else ()
##############################
class genQuark(wrappedChain.calculable) :
    '''Indices of non-top quarks resulting from hard interaction.'''
    def update(self,_) :
        iHard = self.source["genIndicesHardPartons"]
        moms = self.source['genMotherIndex']
        iHardMax = max(iHard)
        self.value = sorted([i for i,(id,imom) in enumerate(zip(self.source['genPdgId'],
                                                                self.source['genMotherIndex']))
                             if iHardMax<i and imom in iHard and abs(id) in range(1,6)])
##############################
class genGlu(wrappedChain.calculable) :
    '''Indices of gluons resulting from hard interaction.'''
    def update(self,_) :
        iHardMax = max(self.source["genIndicesHardPartons"])
        self.value = sorted([i for i,(id,s) in enumerate(zip(self.source['genPdgId'],
                                                             self.source['genStatus'])) if id==21 and s==3 and iHardMax<i],
                            reverse=True,
                            key = lambda i: self.source['genP4'].at(i).Pt() )
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
class genIndices(wrappedChain.calculable) :
    @property
    def name(self) : return "genIndices" + self.label

    def __init__(self, pdgs = [], label = None, status = [], motherPdgs = []) :
        self.label = label
        self.PDGs = frozenset(pdgs)
        self.status = frozenset(status)
        self.motherPdgs = frozenset(motherPdgs)
        self.moreName = "; ".join(["pdgId in %s" %str(list(self.PDGs)),
                                   "status in %s"%str(list(self.status)),
                                   "motherPdg in %s"%str(list(self.motherPdgs))
                                   ])

    def update(self,_) :
        pdg = self.source["genPdgId"]
        status = self.source["genStatus"]
        motherPdg = self.source["genMotherPdgId"]
        self.value = filter( lambda i: ( (not self.PDGs) or (pdg.at(i) in self.PDGs) ) and \
                                 ( (not self.status) or (status.at(i) in self.status) ) and \
                                 ( (not self.motherPdgs) or (motherPdg.at(i) in self.motherPdgs) ),
                             range(pdg.size()) )

class genIndicesPtSorted(wrappedChain.calculable) :
    @property
    def name(self) :
        return "%sPtSorted"%self.label

    def __init__(self, label = "") :
        self.label = "genIndices"+label

    def update(self,_) :
        p4 = self.source["genP4"]
        self.value = sorted(self.source[self.label], key = lambda i:p4.at(i).pt(), reverse = True)

class genRootSHat(wrappedChain.calculable) :
    def update(self,_) :
        iHard = self.source["genIndicesHardPartons"]
        p4s = self.source["genP4"]
        self.value = None if not iHard else (p4s.at(iHard[0])+p4s.at(iHard[1])).mass()

class genSumPt(wrappedChain.calculable) :
    @property
    def name(self) :
        return "_".join(["genSumPt"]+self.indexLabels)

    def __init__(self, indexLabels = []) :
        self.indexLabels = map(lambda s:s.replace("genIndices",""), indexLabels)

    def update(self,_) :
        indices = []
        for label in self.indexLabels :
            indices += self.source["genIndices"+label]
        indices = set(indices)

        self.value = 0.0
        p4 = self.source["genP4"]
        for i in indices :
            self.value += p4.at(i).pt()

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
class genIndicesStatus3NoStatus3Daughter(wrappedChain.calculable) :
    def update(self,_) :
        status = self.source["genStatus"]
        mother = self.source["genMotherIndex"]

        status3List = filter( lambda i: status.at(i)==3, range(status.size()) )
        motherIndices = set([mother[i] for i in status3List])
        self.value = filter( lambda i: i not in motherIndices, status3List )
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

