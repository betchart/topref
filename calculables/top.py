from supy import wrappedChain,utils,calculables
import math,operator,itertools,collections,ROOT as r
try: import numpy as np
except: pass

######################################
class TopJets(wrappedChain.calculable) :
    def __init__(self, jets ) : self.value = jets
    def update(self,_): pass

class TopLeptons(wrappedChain.calculable) :
    def __init__(self, leptons ) : self.value = leptons
    def update(self,_): pass

class SemileptonicTopIndex(wrappedChain.calculable) :
    def update(self,_) :
        self.value = next( iter(self.source["Indices".join(self.source["TopLeptons"])]), None )

class fitTopRecoIndex(wrappedChain.calculable) :
    value = 0
    def update(self,_) : pass
class genTopRecoIndex(wrappedChain.calculable) :
    def update(self,_) :
        iPQHL = self.source['IndicesGenTopPQHL']
        self.value = next((i for i,TR in enumerate(self.source['TopReconstruction']) if TR['iPQHL']==iPQHL), -1)

class IndicesGenTopPQHL(wrappedChain.calculable) :
    rMax = 0.6
    def update(self,_) :
        if not self.source['genTopTTbar'] :
            self.value = (None,)*4
            return
        genP4 = self.source['genP4']
        iTT = self.source['genTTbarIndices']

        PQHL = [genP4[i] if i!=None else None for i in ( (iTT['q'][:2] if iTT['q'] else 2*[None]) + [ iTT['bhad'],iTT['blep'] ] )]

        jets = self.source['TopJets']
        indices = self.source['Indices'.join(jets)]
        p4 = self.source['AdjustedP4'.join(jets)]

        dRIs = [ min( (r.Math.VectorUtil.DeltaR( p4[i], gen ), i) for i in indices ) if gen else
                 (None,None)
                 for gen in PQHL ]

        PQHL = [i if dR<self.rMax else None for dR,i in dRIs ]
        self.value = tuple( sorted(PQHL[:2]) + PQHL[2:] )

class IndicesGenTopPQH(wrappedChain.calculable) :
    def update(self,_) : self.value = self.source['IndicesGenTopPQHL'][:3]
class IndexGenTopL(wrappedChain.calculable) :
    def update(self,_) : self.value = self.source['IndicesGenTopPQHL'][3]

class IndicesGenTopQQBB(wrappedChain.calculable) :
    def update(self,_) :
        iPQHL = self.source['IndicesGenTopPQHL']
        self.value =  iPQHL[:2]+tuple(sorted(iPQHL[2:]))

class IndicesGenTopExtra(wrappedChain.calculable) :
    rMax = 0.6
    def update(self,_) :
        imom = self.source['genMotherIndex']
        p4 = self.source['genP4']
        pdg = self.source['genPdgId']
        status = self.source['genStatus']

        extraP4 = [p4[i] for i in range(8,len(imom)) if 2<imom[i]<6 and abs(pdg[i]) in [1,2,3,4,5,11,13,15,21] and status[i]==3]
        jets = self.source['TopJets']
        indices = self.source['Indices'.join(jets)]
        jet = self.source['AdjustedP4'.join(jets)]
        self.value = [j for j in indices if any( self.rMax > r.Math.VectorUtil.DeltaR(jet[j],gen) for gen in extraP4 ) ]


######################################
class TopP4Calculable(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(['P4'])
######################################
class RawHadWmass(TopP4Calculable) :
    '''m_{W}^{had,raw}'''
    def update(self,_) : self.value = self.source[self.P4]['rawW'].M()
######################################
class Key(TopP4Calculable) :
    def update(self,_) : self.value = self.source[self.P4]['key']
######################################
class Chi2(TopP4Calculable) :
    '''#chi^{2}'''
    def update(self,_) : self.value = self.source[self.P4]['chi2']
######################################
class HadChi2(TopP4Calculable) :
    '''#chi^{2}_{had}'''
    def update(self,_) : self.value = self.source[self.P4]['hadChi2']
######################################
class SumP4(TopP4Calculable) :
    def update(self,_) : self.value = self.source[self.P4]['t'] + self.source[self.P4]['tbar']
######################################
class SumPt(TopP4Calculable) :
    '''t.pt+#bar{t}.pt'''
    def update(self,_) : self.value = self.source[self.P4]['t'].pt() + self.source[self.P4]['tbar'].pt()
######################################
class PtSum(TopP4Calculable) :
    '''t#bar{t}.pt'''
    def update(self,_) : self.value = self.source["SumP4".join(self.fixes)].pt()
######################################
class PtOverSumPt(TopP4Calculable) :
    '''t#bar{t}.pt / ( t.pt + #bar{t}.pt )'''
    def update(self,_) : self.value = self.source['PtSum'.join(self.fixes)]/self.source['SumPt'.join(self.fixes)]
######################################
class MassSum(TopP4Calculable) :
    '''(t+#bar{t}).m'''
    def update(self,_) : self.value = self.source["SumP4".join(self.fixes)].mass()
######################################
class RapiditySum(TopP4Calculable) :
    '''t#bar{t}.y'''
    def update(self,_) : self.value = abs(self.source["SumP4".join(self.fixes)].Rapidity())
######################################
class TanhRapiditySum(TopP4Calculable) :
    '''tanh(t#bar{t}.y)'''
    def update(self,_) : self.value = math.tanh( self.source['RapiditySum'.join(self.fixes)])
######################################
class SumRapidities(TopP4Calculable) :
    '''|t.y+#bar{t}.y|'''
    def update(self,_) : self.value = abs( self.source[self.P4]['t'].Rapidity() +
                                           self.source[self.P4]['tbar'].Rapidity() )
######################################
class TanhAvgRapidity(TopP4Calculable) :
    '''tanh(|t.y+#bar{t}.y|/2)'''
    def update(self,_) : self.value = math.tanh( 0.5*self.source['SumRapidities'.join(self.fixes)] )
######################################
class SumAbsRapidities(TopP4Calculable) :
    '''|t.y|+|#bar{t}.y|'''
    def update(self,_) : self.value = ( abs( self.source[self.P4]['t'].Rapidity() ) +
                                        abs( self.source[self.P4]['tbar'].Rapidity() ) )
######################################
class TanhAvgAbsRapidity(TopP4Calculable) :
    '''tanh((|t.y|+|#bar{t}.y|)/2)'''
    def update(self,_) : self.value = math.tanh( 0.5*self.source['SumAbsRapidities'.join(self.fixes)] )
######################################
class TtxMass(TopP4Calculable) :
    def update(self,_) : self.value = self.source[self.P4]['ttx'].mass()
######################################
class TtxPt(TopP4Calculable) :
    def update(self,_) : self.value = self.source[self.P4]['ttx'].pt()
######################################
class TtxPz(TopP4Calculable) :
    def update(self,_) : self.value = self.source[self.P4]['ttx'].z()
######################################
class TanhDeltaAbsY(TopP4Calculable) :
    '''tanh(|t.y|-|#bar{t}.y|)'''
    def update(self,_) :
        self.value = math.tanh(self.source['DeltaAbsYttbar'.join(self.fixes)])
######################################
class DeltaAbsYttbar(TopP4Calculable) :
    '''|t.y|-|#bar{t}.y|'''
    def update(self,_) :
        self.value = abs(self.source[self.P4]['t'].Rapidity()) - abs(self.source[self.P4]['tbar'].Rapidity())
######################################
class DeltaYttbar(TopP4Calculable) :
    '''t.y - #bar{t}.y'''
    def update(self,_) :
        self.value = self.source[self.P4]['t'].Rapidity() - self.source[self.P4]['tbar'].Rapidity()
######################################
class DPtDPhi(TopP4Calculable):
    def update(self,_):
        p4 = self.source[self.P4]
        ttbar = p4['t'] + p4['tbar']
        diff =  p4['tbar'] - p4['t']
        dphi = math.acos(math.cos(ttbar.phi() - diff.phi()))
        self.value = 1 - 2 * dphi / math.pi
#######################################
class TanhDirectedDeltaYttbar(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(['DeltaYttbar','SignQuarkZ'])
    def update(self,_) :
        self.value = math.tanh( self.source[self.SignQuarkZ] * self.source[self.DeltaYttbar] )
######################################
class SignQuarkZ(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
    def update(self,_) :
        self.value = 0 if not self.source['wQQ'] else self.source['qDir']
######################################
class FifthJetIndex(TopP4Calculable) :
    def update(self,_) :
        self.value = next(iter(self.source[self.P4]['iOtherJets']),None)
######################################

class genTopP4(wrappedChain.calculable) :
    def update(self,_) :
        indices = self.source['genTTbarIndices']
        p4 = self.source['genP4']
        qg = max(self.source['genQG'],self.source['genAG'])
        self.value = { 't':p4[indices['t']],
                       'tbar':p4[indices['tbar']],
                       'lepton': p4[indices['lplus']] if indices['lplus'] else p4[indices['lminus']] if indices['lminus'] else None,
                       'neutrino': None,
                       'p' : p4[indices['q'][0]] if indices['q'] else None,
                       'q' : p4[indices['q'][1]] if len(indices['q'])>1 else None,
                       'rawW': None,
                       'sumP4': None,
                       'key': None,
                       'chi2': None,
                       'hadChi2':None,
                       'ttx': None,
                       'fifth':None,
                       'ntracks':None
                       }

class fitTopP4(wrappedChain.calculable) :
    def update(self,_) :
        reco = self.source["TopReconstruction"][0]
        t = reco['top']
        tbar = reco['tbar']
        q_z = 0.5*(t+tbar).z()
        self.value = {'t':t,
                      'tbar':tbar,
                      'quark': utils.LorentzV().SetPxPyPzE(0,0,q_z,abs(q_z)),
                      'lepton': reco['lep'],
                      'neutrino' : reco['nu'],
                      'p' : reco['hadP'],
                      'q' : reco['hadQ'],
                      'rawW': reco['hadWraw'],
                      'metP4':reco['metP4'],
                      'key': reco['key'],
                      'chi2': reco['chi2'],
                      'hadChi2': reco['hadChi2'],
                      'ttx': reco['ttx'],
                      'fifth': reco['iX']!=None,
                      'iPQHL' : reco['iPQHL'],
                      'iOtherJets': [i for i in self.source['Indices'.join(self.source['TopJets'])] if i not in reco['iPQHL']]
                      }
######################################
######################################
######################################

class genTopTTbar(wrappedChain.calculable) :
    def update(self,_) :
        ids = [] if self.source['isRealData'] else list(self.source['genPdgId'])
        self.value = tuple(ids.index(i) for i in [6,-6]) if all([id in ids for id in [-6,6]]) else ()
######################################
class ttDecayMode(wrappedChain.calculable) :
    def update(self,_) :
        pdg = self.source['genPdgId']
        mom = self.source['genMotherPdgId']
        debris = [abs(pdg[i]) for i in range(len(pdg)) if abs(mom[i])==24]
        self.value = ('' if not self.source['genTopTTbar'] else
                      'ee' if debris.count(11) == 2 else
                      'mm' if debris.count(13) == 2 else
                      'tt' if debris.count(15) == 2 else
                      'em' if 11 in debris and 13 in debris else
                      'et' if 11 in debris and 15 in debris else
                      'mt' if 13 in debris and 15 in debris else
                      'ej' if debris.count(11) else
                      'mj' if debris.count(13) else
                      'tj' if debris.count(15) else
                      'jj')
######################################
class genTTbarIndices(wrappedChain.calculable) :
    def update(self,_) :
        ids = [i for i in self.source['genPdgId']]
        mom = self.source['genMotherPdgId']
        self.value = dict([(name, ids.index(i)) for name,i in [('t',6),
                                                               ('tbar',-6),
                                                               ('wplus',24),
                                                               ('wminus',-24)
                                                               ]])
        self.value.update(dict([ (w+"Child",filter(lambda i: mom[i]==ids[self.value[w]],range(len(ids)))) for w in ['wplus','wminus','t','tbar']]))
        self.value['b'] = next(i for i in self.value['tChild'] if abs(ids[i])!=24)
        self.value['bbar'] = next(i for i in self.value['tbarChild'] if abs(ids[i])!=24)
        self.value['lplus'] = next((i for i in self.value['wplusChild'] if ids[i] in [-11,-13]),None)
        self.value['lminus'] = next((i for i in self.value['wminusChild'] if ids[i] in [11,13]),None)
        self.value['q'] = ((self.value['wplusChild'] if not self.value['lplus'] else []) +
                            (self.value['wminusChild'] if not self.value['lminus'] else []))
        self.value['nu'] = next((i for i in (self.value['wplusChild']+self.value['wminusChild']) if abs(ids[i]) in [12,14]),None)
        self.value['semi'] = (self.value['lplus'] is None)^(self.value['lminus'] is None)
        self.value['blep'] = None if not self.value['semi'] else self.value['b'] if self.value['lminus']==None else self.value['bbar']
        self.value['bhad'] = None if not self.value['semi'] else self.value['bbar'] if self.value['lminus']==None else self.value['b']
######################################
class TopReconstruction(wrappedChain.calculable) :
    def __init__(self, eCoupling = 0.55, v2had = False) :
        self.epsilon = 1e-7
        for item in ['eCoupling', # percentage of jet resolution used to sharpen MET resolution
                     'v2had',     # v2 (1parameter,3residuals) is twice as fast as v1 (3parameters,5residuals) but 5% less accurate
                     ] : setattr(self,item,eval(item))
        self.moreName = "eCoupl:%.2f; v%dhad; v2lep"%(eCoupling,v2had+1)
        self.maxFits = 1

    def update(self,_) :
        
        jets = dict( (item, self.source[item.join(self.source["TopJets"])] )
                     for item in ["AdjustedP4","Indices","Resolution","CovariantResolution2","ScalingBQN"] )

        lepton = dict( (item, self.source[item.join(self.source['TopLeptons'])][self.source["SemileptonicTopIndex"]])
                       for item in ["Charge","P4"])

        topP = self.source["TopComboQQBBProbability"]
        TCL = self.source['TopCandidateLikelihood']
        
        recos = []
        for iPQH,i4s in itertools.groupby( sorted(sorted( TCL,
                                                          key = TCL.__getitem__,
                                                          reverse = True)[:self.maxFits] ),
                                           key = lambda x: x[:3]) :

            hadFit = utils.fitKinematic.leastsqHadronicTop2(*zip(*((jets["AdjustedP4"][i]*jets["ScalingBQN"]['B' if t==2 else 'Q'][i], jets["Resolution"][i]) for t,i in enumerate(iPQH))) ) if self.v2had else \
                     utils.fitKinematic.leastsqHadronicTop( *zip(*((jets["AdjustedP4"][i]*jets["ScalingBQN"]['B' if t==2 else 'Q'][i], jets["Resolution"][i]) for t,i in enumerate(iPQH))), widthW = 4./2 ) #tuned w width

            metP4 = self.source["metAdjustedP4"] + hadFit.rawT - hadFit.fitT
            nuXY = np.array([metP4.x(), metP4.y()])
            nuErr2 = sum([-self.eCoupling*jets["CovariantResolution2"][i] for i in iPQH], self.source["metCovariance"])

            for iPQHL in i4s :
                iL = iPQHL[3]
                iQQBB = iPQHL[:2]+tuple(sorted(iPQHL[2:]))
                b = jets["AdjustedP4"][iL]
                bscale = jets['ScalingBQN']['B'][iL]
                nuXY_b = nuXY - (bscale - 1)*np.array([b.y(),b.y()])
                nuErr2_b = nuErr2-self.eCoupling*jets["CovariantResolution2"][iL]
                lepFit = utils.fitKinematic.leastsqLeptonicTop2( b*bscale, jets["Resolution"][iL], lepton["P4"], nuXY_b, nuErr2_b)
                tt = hadFit.fitT + lepFit.fitT
                iX,ttx = min( [(None,tt)]+[(i,tt+jets["AdjustedP4"][i]) for i in jets["Indices"] if i not in iPQHL], key = lambda lv : lv[1].pt() )
                recos.append( {"nu"   : lepFit.fitNu,       "hadP" : hadFit.fitJ[0],
                               "lep"  : lepFit.mu,          "hadQ" : hadFit.fitJ[1],
                               "lepB" : lepFit.fitB,        "hadB" : hadFit.fitJ[2],
                               "lepW" : lepFit.fitW,        "hadW" : hadFit.fitW,
                               "lepTopP4" : lepFit.fitT,    "hadTopP4": hadFit.fitT,
                               "lepChi2" : lepFit.chi2,     "hadChi2" : hadFit.chi2,
                               "chi2" : hadFit.chi2 + lepFit.chi2,
                               "probability" : max(self.epsilon,topP[iQQBB]),

                               "top"  : lepFit.fitT if lepton["Charge"] > 0 else hadFit.fitT,
                               "tbar" : hadFit.fitT if lepton["Charge"] > 0 else lepFit.fitT,
                               "ttx" : ttx, "iX" : iX,

                               "iPQHL": iPQHL,
                               "lepCharge": lepton["Charge"], "hadTraw" : hadFit.rawT, "lepTraw" : lepFit.rawT,
                               "lepBound" : lepFit.bound,     "hadWraw" : hadFit.rawW, "lepWraw" : lepFit.rawW,
                               "metP4": metP4 + b - lepFit.fitB,
                               "nuErr2":nuErr2_b,
                               "residuals" : dict( zip(["lep"+i for i in "BSLWT"],  lepFit.residualsBSLWT ) +
                                                   zip(["had"+i for i in "PQBWT"], hadFit.residualsPQBWT ) ),
                               "nuEllipse"       : lepFit.Ellipse,
                               "nuSolutions"     : lepFit.solutions,
                               "nuChi2Matrix"    : lepFit.X
                               })
                recos[-1]["key"] = recos[-1]['chi2'] - 2*math.log(recos[-1]['probability'])

        self.value = sorted( recos,  key = lambda x: x["key"] )

######################################
class TTbarSignExpectation(wrappedChain.calculable) :

    def __init__(self, nSamples = 16, qDirFunc = None ) :
        self.nSamples = nSamples
        self.qDirFunc = qDirFunc
        self.moreName = "%d samples; %s"%(nSamples,qDirFunc)
        self.nu = utils.LorentzV()

    def signExpect(self,topReco) :
        samples = []
        lep = topReco['lepTopP4']
        had = topReco['hadTopP4']
        had_y = had.Rapidity()
        hadIsTop = topReco['lepCharge'] < 0
        bmu = topReco['lep'] + topReco['lepB']
        qDirFunc = self.source[self.qDirFunc] if self.qDirFunc else (lambda L,H : 1 if L.Rapidity()+H.Rapidity() > 0 else -1)
        
        c,s = topReco['nuSolutions'][0][:2]
        Ellipse = topReco['nuEllipse']
        M = topReco['nuChi2Matrix']
        if not (s or c) : return qDirFunc(had,lep) * (-1)**(hadIsTop^( lep.Rapidity() < had_y))
        tau_0 = math.atan2(s,c)
        for tau in np.arange(tau_0, tau_0 + 2*math.pi, 2*math.pi/self.nSamples)[::-1] :
            sol = np.array([math.cos(tau),math.sin(tau),1])
            chi2 = sol.T.dot(M.dot(sol))
            x,y,z = Ellipse.dot(sol)
            self.nu.SetPxPyPzE(x,y,z,0); self.nu.SetM(0)
            lep = bmu + self.nu
            samples.append( (math.exp(-0.5*chi2),
                             qDirFunc(had,lep) * (-1)**(hadIsTop^( lep.Rapidity() < had_y ) ) ) )

        xw = sum(p*sdy for p,sdy in samples)
        w = sum(p for p,sdy in samples)
        return xw / w if w else 0

    def update(self,_) :
        signs = [ (math.exp(-0.5*reco['key']),
                   self.signExpect(reco)) for reco in self.source["TopReconstruction"]]
        xw = sum(p*sign for p,sign in signs)
        w = sum(p for p,sign in signs)
        self.value =  xw / w if w else 0
######################################
class kinfitFailureModes(wrappedChain.calculable) :
    def update(self,_) :
        reco =  self.source["TopReconstruction"][0]
        genP4 = self.source["genP4"]
        pdg = self.source["genPdgId"]
        mom = self.source["genMotherPdgId"]
        status = self.source["genStatus"]
        jets = self.source["TopJets"]
        p4 = self.source["AdjustedP4".join(jets)]
        indices = self.source["Indices".join(jets)]
        igenTT = self.source["genTTbarIndices"]
        iTop,iTbar = self.source["genTopTTbar"]

        self.value = {}
        if self.source["ttDecayMode"] not in ["ej","mj"] : return
        self.value["met"] = (abs(r.Math.VectorUtil.DeltaPhi( self.source["metAdjustedP4"], genP4[igenTT['nu']])) < 0.7)
        self.value["nu"] = r.Math.VectorUtil.DeltaR( genP4[igenTT['nu']], reco['nu']) < 0.7

        igens = tuple(igenTT['q'])+(igenTT['bhad'],igenTT['blep'])
        self.value["jet"] = ( all( any(r.Math.VectorUtil.DeltaR( p4[iJet],genP4[igen]) < 0.5 for iJet in indices ) for igen in igens)
                              and not any(r.Math.VectorUtil.DeltaR( genP4[igen], genP4[jgen]) < 0.5 for igen,jgen in itertools.combinations(igens,2) ) )
        
        self.value["had"] = all( any(r.Math.VectorUtil.DeltaR( p4[iJet],genP4[igen]) < 0.5 for iJet in reco['iPQHL'][:3]) for igen in igens[:3] )
        self.value["blep"] = r.Math.VectorUtil.DeltaR( genP4[igenTT['blep']],p4[reco['iPQHL'][3]]) < 0.5
        self.value["bhad"] = r.Math.VectorUtil.DeltaR( genP4[igenTT['bhad']],p4[reco['iPQHL'][2]]) < 0.5
        glu = next( (genP4[i] for i in range(6,len(genP4)) if pdg[i]==21 and status[i]==3), None)
        if glu or reco['iX']!=None :
            self.value["glu"] = bool(glu) and (reco['iX'] != None) and (r.Math.VectorUtil.DeltaR( glu, p4[reco['iX']] ) < 0.5)
        self.value["t"] = r.Math.VectorUtil.DeltaR(genP4[iTop],reco['top']) < 0.5
        self.value["/t"] = r.Math.VectorUtil.DeltaR(genP4[iTbar],reco['tbar']) < 0.5
######################################
class lepDeltaRTopRecoGen(wrappedChain.calculable) :
    def update(self,_):
        indices = self.source['genTTbarIndices']
        genLep = self.source['genP4'][max(indices['lplus'],indices['lminus'])]
        self.value = [r.Math.VectorUtil.DeltaR(genLep,reco['lep']) for reco in self.source['TopReconstruction']]
class nuDeltaRTopRecoGen(wrappedChain.calculable) :
    def update(self,_):
        genNu = self.source['genP4'][self.source['genTTbarIndices']['nu']]
        self.value = [r.Math.VectorUtil.DeltaR(genNu,reco['nu']) for reco in self.source['TopReconstruction']]
class bLepDeltaRTopRecoGen(wrappedChain.calculable) :
    def update(self,_):
        indices = self.source['genTTbarIndices']
        genLepB = self.source['genP4'][indices['blep']]
        self.value = [r.Math.VectorUtil.DeltaR(genLepB,reco['lepB']) for reco in self.source['TopReconstruction']]
class bHadDeltaRTopRecoGen(wrappedChain.calculable) :
    def update(self,_):
        indices = self.source['genTTbarIndices']
        genHadB = self.source['genP4'][indices['bhad']]
        self.value = [r.Math.VectorUtil.DeltaR(genHadB,reco['hadB']) for reco in self.source['TopReconstruction']]
class pqDeltaRTopRecoGen(wrappedChain.calculable) :
    def update(self,_):
        PQ = tuple([self.source['genP4'][self.source['genTTbarIndices']['q'][i]] for i in range(2)])
        self.value= [min([ tuple(sorted([r.Math.VectorUtil.DeltaR(*t) for t in [(reco['hadP'],PQ[i]),(reco['hadQ'],PQ[j])]])) for i,j in [(0,1),(1,0)]])\
                     for reco in self.source['TopReconstruction']]
class qDeltaRTopRecoGen(wrappedChain.calculable) :
    def update(self,_): self.value = [pq[1] for pq in self.source['pqDeltaRTopRecoGen']]

######################################
class extraJetMoments2Sum(wrappedChain.calculable) :
    def update(self,_):
        jets = self.source['TopJets']
        i5 = self.source['fitTopFifthJetIndex']
        self.value = self.source['Moments2Sum'.join(jets)][i5] if i5!=None else 0
######################################

class ttSymmAnti(calculables.secondary) :
    def __init__(self, thisSample, samples, inspect=False, weights=[], altTag=None) :
        collection = ('genTop','')
        self.varX = 'DPtDPhi'.join(collection)
        self.varY = 'TanhDeltaAbsY'.join(collection)
        for item in ['thisSample','inspect','weights'] : setattr(self,item,eval(item))
        self.__symm, self.__anti = None,None
        self.samples = samples
        self.altTag = altTag

    def baseSamples(self): return self.samples

    def uponAcceptance(self,ev) :
        w = reduce(operator.mul, [ev[W] for W in self.weights], 1)

        x,y = ev[self.varX],ev[self.varY]
        self.book.fill(x, self.varX, 100, -1, 1, title = ";%s"%self.varX, w = w)
        self.book.fill(y, self.varY, 100, -1, 1, title = ";%s"%self.varY, w = w)
        self.book.fill((x,y), '2_x_y', (100,100), (-1,-1), (1,1), title=";%s;%s"%(self.varX,self.varY), w = w)
        if not self.inspect : return
        symmanti = ev[self.name]
        if not symmanti : return
        sumsymmanti = sum(symmanti)
        symm,anti = symmanti

        wsymm = w * symm / sumsymmanti
        wanti = w * anti / sumsymmanti
        wflat = w * 0.5 / max(1e-6, sumsymmanti)

        self.book.fill(x, 'x_0symm', 50, -1, 1, w = wsymm, title = ';(symm) %s;events / bin'%self.varX)
        self.book.fill(x, 'x_1anti', 50, -1, 1, w = wanti, title = ';(anti) %s;events / bin'%self.varX)
        self.book.fill(x, 'x_2flat', 50, -1, 1, w = wflat, title = ';(flat) %s;events / bin'%self.varX)

        self.book.fill(y, 'y_0symm', 100, -1, 1, w = wsymm, title = ';(symm) %s;events / bin'%self.varY)
        self.book.fill(y, 'y_1anti', 100, -1, 1, w = wanti, title = ';(anti) %s;events / bin'%self.varY)
        self.book.fill(y, 'y_2flat', 100, -1, 1, w = wflat, title = ';(flat) %s;events / bin'%self.varY)

        self.book.fill((x,y), '2_x_y_0symm', (50,100), (-1,-1), (1,1), w = wsymm, title = 'symm;%s;%s;'%(self.varX,self.varY))
        self.book.fill((x,y), '2_x_y_1anti', (50,100), (-1,-1), (1,1), w = wanti, title = 'anti;%s;%s;'%(self.varX,self.varY))
        self.book.fill((x,y), '2_x_y_2flat', (50,100), (-1,-1), (1,1), w = wflat, title = 'flat;%s;%s;'%(self.varX,self.varY))

    def update(self,_) :
        self.value = ()
        if not self.__symm : return
        x = np.array([self.source[self.varX]],'d')
        y = self.source[self.varY]

        anti = self.__anti[0].EvalPar(x, np.array([f.Eval(y) for f in self.__anti[1:]],'d'))
        symm = max(self.__symm[0].EvalPar(x, np.array([f.Eval(y) for f in self.__symm[1:]],'d')), 1.1*abs(anti))
        self.value = (symm,anti) if symm else (1,0)

    def setup(self,*_) :
        hist = self.fromCache([self.thisSample], ['2_x_y'], self.altTag)[self.thisSample]['2_x_y']
        if not hist : print "cannot find cache:", self.name ; return
        self.__stashsymmanti = self.prep(hist)
        self.__symm,self.__anti = [[next(iter(h.GetListOfFunctions())) for h in hs] for hs in self.__stashsymmanti]

    @staticmethod
    def prep(hist) :
        r.gStyle.SetLineStyleString(11,"0 20000"); # "invisible"
        hist.Scale(1./hist.Integral())
        symm,anti = utils.symmAnti2(hist)
        symmSliceX = r.TF1('func','[0] + [1]*x*(1-x**2/3) + [2]*x**2/2 * (1-x**2/2)',-1,1)
        symmSliceX.SetParNames("A^{+}","B^{+}","C^{+}")
        symmSliceX.SetLineWidth(0)
        symmSliceX.SetLineStyle(11)
        symm.FitSlicesX(symmSliceX)
        symm.GetListOfFunctions().Add(symmSliceX)
        symmHists = [symm] + [r.gDirectory.Get("%s_%d"%(symm.GetName(),i)) for i in range(3)]
        symmHists[1].Fit('++'.join(['(1-abs(x))']+['x**%d'%d for d in [0,2,4,6,8,10,12,14,16,18]]),'Q')
        symmHists[2].Fit('++'.join('x**%d'%d for d in [1,3,5,7,9,11,13,15,17]),'Q')
        symmHists[3].Fit('++'.join('x**%d'%d for d in [0,2,4,6,8,10,12,14,16]),'Q')
        antiSliceX = r.TF1('func','[0] + [1]*x*(1-x**2/3)',-1,1)
        antiSliceX.SetParNames("A^{-}","B^{-}")
        antiSliceX.SetLineWidth(0)
        antiSliceX.SetLineStyle(11)
        anti.FitSlicesX(antiSliceX)
        anti.GetListOfFunctions().Add(antiSliceX)
        antiHists = [anti] + [r.gDirectory.Get("%s_%d"%(anti.GetName(),i)) for i in range(2)]
        antiHists[1].Fit('++'.join('x**%d'%d for d in [1,3,5,7,9,11]),'Q')
        antiHists[2].Fit('++'.join('x**%d'%d for d in [0,2,4,6,8,10]),'Q')
        return symmHists,antiHists
        
    def reportCache(self) :
        fileName = '/'.join(self.outputFileName.split('/')[:-1]+[self.name])
        optstat = r.gStyle.GetOptStat()
        from supy import whereami
        r.gROOT.ProcessLine(".L %s/cpp/tdrstyle.C"%whereami())
        r.setTDRStyle()
        r.tdrStyle.SetOptStat(0)
        r.tdrStyle.SetOptFit(0)
        r.tdrStyle.SetPalette(1)
        r.tdrStyle.SetMarkerSize(0.5)
        r.tdrStyle.SetTitleH(0.05)
        r.tdrStyle.SetTitleW(0.4)
        r.tdrStyle.SetTitleX(0.4)
        r.tdrStyle.SetTitleY(1.0)
        r.tdrStyle.SetTitleBorderSize(0)

        samples = self.baseSamples()
        names = ['%s#rightarrow^{}t#bar{t} '%i for i in ['gg','qg','q#bar{q}','#bar{q}g']]
        colors = [r.kBlack,r.kBlue,r.kRed,r.kGreen]

        hists = [self.fromCache(samples, ['2_x_y'],self.altTag)[s]['2_x_y'] for s in samples]
        for h in hists:
            h.UseCurrentStyle()
            h.SetTitle(';X_{T};X_{L}')
        symmHists,antiHists = zip(*[self.prep(h) for h in hists])
        symm,symm0,symm1,symm2 = zip(*symmHists)
        anti,anti0,anti1 = zip(*antiHists)

        for h,n in zip(symm,names) : h.SetTitle(n + ' (symmetric)')
        for h,n in zip(anti,names) : h.SetTitle(n + ' (antisymmetric)')

        height = max(h.GetBinContent(h.GetMaximumBin()) for h in symm)
        c = r.TCanvas()
        c.Print(fileName+'.pdf[')
        c.Divide(2,2)
        r.tdrStyle.SetOptTitle(1)
        for i,h in enumerate(symm) :
            c.cd(i+1)
            h.SetMaximum(height)
            h.SetMinimum(0)
            h.Draw('colz')
        c.Print(fileName+'.pdf')
        #for i,h in enumerate(symm) :
        #    c.cd(i+1)
        #    h.Draw('surf2')
        #c.Print(fileName+'.pdf')
        r.tdrStyle.SetOptTitle(0)

        c.Clear()
        c.cd(0)

        for par in range(3) :
            hists = eval('symm%d'%par)
            height = 1.05*max(h.GetBinContent(h.GetMaximumBin()) for h in hists)
            for i,h in enumerate(hists) :
                h.GetYaxis().SetTitle('Fitted value of '+h.GetTitle().split('=')[-1])
                fit = next(iter(h.GetListOfFunctions()),None)
                if fit :
                    fit.SetLineColor(colors[i])
                    fit.SetLineWidth(1)
                h.SetLineColor(colors[i])
                h.SetMarkerColor(colors[i])
                h.SetMaximum(height)
                h.SetMinimum(-height if h.GetMinimum()<0 else 0)
                h.Draw('same' if i else '')
            c.Print(fileName+'.pdf')


        c.Clear()
        c.Divide(2,2)
        r.tdrStyle.SetOptTitle(1)
        height = max(h.GetBinContent(h.GetMaximumBin()) for h in anti)
        for i,h in enumerate(anti) :
            c.cd(i+1)
            h.SetMaximum(height)
            h.SetMinimum(-height)
            h.Draw('colz')
        c.Print(fileName+'.pdf')
        #for i,h in enumerate(anti) :
        #    c.cd(i+1)
        #    h.Draw('surf2')
        #c.Print(fileName+'.pdf')
        r.tdrStyle.SetOptTitle(0)

        c.Clear()
        c.cd(0)

        for par in range(2) :
            hists = eval('anti%d'%par)
            height = 1.05*max([h.GetBinContent(h.GetMaximumBin()) for h in hists]+[abs(h.GetMinimum()) for h in hists])
            for i,h in enumerate(hists) :
                h.GetYaxis().SetTitle('Fitted value of '+h.GetTitle().split('=')[-1])
                fit = next(iter(h.GetListOfFunctions()),None)
                if fit :
                    fit.SetLineColor(colors[i])
                    fit.SetLineWidth(1)
                h.SetLineColor(colors[i])
                h.SetMarkerColor(colors[i])
                h.SetMaximum(height)
                h.SetMinimum(-height)
                h.Draw('same' if i else '')
            c.Print(fileName+'.pdf')

        c.Clear()

        height = 1.1*max(h.GetMaximum() for h in symm)
        funcs = [next(iter(s.GetListOfFunctions()),None) for s in symm]
        [(func.SetLineStyle(1),func.SetLineWidth(1),func.SetLineColor(colors[i])) for i,func in enumerate(funcs)]
        for iY in range(1,1+symm[0].GetNbinsY()) :
            projs = []
            for i,(s,func) in enumerate(zip(symm,funcs)) :
                proj = s.ProjectionX('%s_%d'%(s.GetName(),iY),iY,iY)
                proj.SetLineColor(colors[i])
                proj.SetMarkerColor(colors[i])
                proj.SetMinimum(0)
                proj.SetMaximum(height)
                proj.Fit(func,'QR')
                projs.append(proj)
            for p in projs : p.Draw('same' if i else '')
            c.Print(fileName+'.pdf')

        c.Clear()


        height = 1.1*max([h.GetMaximum() for h in anti]+[abs(h.GetMinimum()) for y in anti])
        funcs = [next(iter(s.GetListOfFunctions()),None) for s in anti]
        [(func.SetLineStyle(1),func.SetLineWidth(1),func.SetLineColor(colors[i])) for i,func in enumerate(funcs)]
        for iY in range(1,1+anti[0].GetNbinsY()) :
            projs = []
            for i,(s,func) in enumerate(zip(anti,funcs)) :
                proj = s.ProjectionX('%s_%d'%(s.GetName(),iY),iY,iY)
                proj.SetLineColor(colors[i])
                proj.SetMinimum(-height)
                proj.SetMaximum(height)
                proj.Fit(func,'QR')
                projs.append(proj)
            for p in projs : p.Draw('same' if i else '')
            c.Print(fileName+'.pdf')

        c.Print(fileName+'.pdf]')
        r.gStyle.SetOptStat(optstat)
        print "Wrote file: %s.pdf"%fileName

class ttAltSymmAnti(ttSymmAnti):
    @property
    def name(self): return 'tt%sSymmAnti'%self.abbr

    @property
    def nameForCache(self): return 'ttSymmAnti'

    def __init__(self, sample, tag, abbr="Alt"):
        for item in ['sample','tag','abbr']: setattr(self,item,eval(item))
        super(ttAltSymmAnti,self).__init__(sample, [], altTag=tag)

    def uponAcceptance(self,_): pass
    def baseSamples(self): return ['bogusSampleDontCacheAnything']

class ttAltSymmAntiWeight(wrappedChain.calculable):
    @property
    def name(self): return 'tt%sSymmAntiWeight'%self.abbr

    def __init__(self, abbr):
        self.abbr = abbr
        self.var = 'tt%sSymmAnti'%self.abbr

    def update(self,_):
        self.value = [sum(self.source[self.var]) / sum(self.source['ttSymmAnti'])]
