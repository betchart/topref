import collections, ROOT as r
from supy import utils,analysisStep
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
    def uponAcceptance(self,ev) :
        index = ev["genIndexTtbarExtraJet"]
        if index==None : return
        p4 = ev['genP4'][index]
        id = ev['genPdgId'][index]
        label = ('g' if id==21 else 'q' if id>0 else 'qbar')
        self.book.fill(p4.pt(), label + 'RecoilPt', 100, 0, 200, title = ';recoil %s.pt;events/bin'%label)
        self.book.fill(p4.eta(),label + 'RecoilEta', 50, -5, 5, title = ';recoil %s.eta;events/bin'%label)
