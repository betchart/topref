from supy import wrappedChain,utils
from jet import xcStrip
import math,operator,itertools,ROOT as r
try: import numpy as np
except: pass

class TopJets(wrappedChain.calculable) :
    def __init__(self, jets ) : self.value = {"fixes":jets,"fixesStripped":xcStrip(jets)}
    def update(self,_): pass
######################################
class TopLeptons(wrappedChain.calculable) :
    def __init__(self, leptons ) : self.value = leptons
    def update(self,_): pass
######################################
class TopP4Calculable(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(['P4'])
######################################
class DeltaPhiLNu(TopP4Calculable) :
    def update(self,_) : self.value = abs(r.Math.VectorUtil.DeltaPhi(self.source[self.P4]['lepton'],self.source[self.P4]['sumP4']))
######################################
class LeptonPt(TopP4Calculable) :
    def update(self,_) : self.value = self.source[self.P4]['lepton'].pt()
######################################
class NuPt(TopP4Calculable) :
    def update(self,_) : self.value = self.source[self.P4]['neutrino'].pt()
######################################
class MET(TopP4Calculable) :
    def update(self,_) : self.value = self.source[self.P4]['sumP4'].pt()
######################################
class RawHadWmass(TopP4Calculable) :
    def update(self,_) : self.value = self.source[self.P4]['rawW'].M()
######################################
class Key(TopP4Calculable) :
    def update(self,_) : self.value = self.source[self.P4]['key']
######################################
class Chi2(TopP4Calculable) :
    def update(self,_) : self.value = self.source[self.P4]['chi2']
######################################
class HadChi2(TopP4Calculable) :
    def update(self,_) : self.value = self.source[self.P4]['hadChi2']
######################################
class SumP4(TopP4Calculable) :
    def update(self,_) : self.value = self.source[self.P4]['t'] + self.source[self.P4]['tbar']
######################################
class SumPt(TopP4Calculable) :
    def update(self,_) : self.value = self.source[self.P4]['t'].pt() + self.source[self.P4]['tbar'].pt()
######################################
class PtSum(TopP4Calculable) :
    def update(self,_) : self.value = self.source["SumP4".join(self.fixes)].pt()
######################################
class MassSum(TopP4Calculable) :
    def update(self,_) : self.value = self.source["SumP4".join(self.fixes)].mass()
######################################
class RapiditySum(TopP4Calculable) :
    def update(self,_) : self.value = abs(self.source["SumP4".join(self.fixes)].Rapidity())
######################################
class AbsSumRapidities(TopP4Calculable) :
    def update(self,_) : self.value = abs( self.source[self.P4]['t'].Rapidity() +
                                           self.source[self.P4]['tbar'].Rapidity() )
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
class FifthJet(TopP4Calculable) :
    def update(self,_) : self.value = self.source[self.P4]['fifth']
######################################
class Ntracks(TopP4Calculable) :
    def update(self,_) : self.value = self.source[self.P4]['ntracks']
######################################
class NtracksExtra(TopP4Calculable) :
    def update(self,_) : self.value = self.source['tracksCountwithPrimaryHighPurityTracks'] - self.source[self.P4]['ntracks']
######################################
class JetAbsEtaMax(wrappedChain.calculable) :
    def __init__(self,collection = None) :
        self.fixes = collection
    def update(self,_) :
        iPQHL = self.source["TopReconstruction"][0]['iPQHL']
        p4 = self.source["CorrectedP4".join(self.source["TopJets"]["fixes"])]
        self.value = max([abs(p4[i].eta()) for i in iPQHL])
######################################
class JetPtMin(wrappedChain.calculable) :
    def __init__(self,collection = None) :
        self.fixes = collection
    def update(self,_) :
        iPQHL = self.source["TopReconstruction"][0]['iPQHL']
        p4 = self.source["CorrectedP4".join(self.source["TopJets"]["fixes"])]
        self.value = min([abs(p4[i].pt()) for i in iPQHL])
######################################
class PartonFractions(TopP4Calculable) :
    def update(self,_) :
        sump4 = self.source['SumP4'.join(self.fixes)]
        m = sump4.mass() / 7000
        z = sump4.z() / 7000
        self.value = sorted([m+z,m-z])
######################################
class PartonLo(TopP4Calculable) :
    def update(self,_) : self.value = self.source['PartonFractions'.join(self.fixes)][0]
######################################
class PartonHi(TopP4Calculable) :
    def update(self,_) : self.value = self.source['PartonFractions'.join(self.fixes)][1]
######################################
class PartonXplusminus(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(['TtxMass','TtxPz'])
    def update(self,_) :
        self.value = ( (self.source[self.TtxPz] + self.source[self.TtxMass] ) / 7000 ,
                       (self.source[self.TtxPz] - self.source[self.TtxMass] ) / 7000 )
######################################
class PartonXhi(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(['PartonXplusminus'])
    def update(self,_) : self.value = max(self.source[self.PartonXplusminus], key = lambda i: abs(i))
######################################
class PartonXlo(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(['PartonXplusminus'])
    def update(self,_) : self.value = min(self.source[self.PartonXplusminus], key = lambda i: abs(i))
######################################
class Pt(wrappedChain.calculable) :
    def __init__(self,collection = None) :
        self.fixes = collection
        self.stash(['SumP4'])
    def update(self,_) : self.value = self.source[self.SumP4].pt()
######################################
class PtOverSumPt(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(['Pt','SumPt'])
    def update(self,_) : self.value = self.source[self.Pt] / self.source[self.SumPt]
######################################
class SqrtPtOverSumPt(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(['PtOverSumPt'])
    def update(self,_) : self.value = math.sqrt(self.source[self.PtOverSumPt])
######################################
class PtPlusSumPt(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(['Pt','SumPt'])
    def update(self,_) : self.value = self.source[self.Pt] + self.source[self.SumPt]
######################################
class SumP4Eta(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(['SumP4'])
    def update(self,_) : self.value = self.source[self.SumP4].eta()
######################################
class SumP4AbsEta(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(['SumP4'])
    def update(self,_) : self.value = abs(self.source[self.SumP4].eta())
######################################
class BoostZ(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(['SumP4'])
    def update(self,_) :
        self.value = r.Math.BoostZ( self.source[self.SumP4].BoostToCM().z())
######################################
class BoostZAlt(TopP4Calculable) :
    def update(self,_) :
        t,tbar = tuple(self.source[self.P4][i] for i in ['t','tbar'])
        pt = t.pt()
        ptbar = tbar.pt()
        self.value = r.Math.BoostZ( - (t.z()/pt + tbar.z()/ptbar) / (t.E()/pt + tbar.E()/ptbar) )
######################################
class DeltaPhi(TopP4Calculable) :
    def update(self,_) :
        self.value = r.Math.VectorUtil.DeltaPhi(self.source[self.P4]['t'], self.source[self.P4]['tbar'])
######################################
class WqqDeltaR(TopP4Calculable) :
    def update(self,_) :
        self.value = r.Math.VectorUtil.DeltaR(self.source[self.P4]['p'],self.source[self.P4]['q'])
######################################
class DeltaBetazRel(TopP4Calculable) :
    def update(self,_) :
        betazTop,betazTbar = [abs(self.source[self.P4][item].BoostToCM().z()) for item in ['t','tbar']]
        self.value = (betazTop - betazTbar) / (1-betazTop*betazTbar)
######################################
class DeltaAbsYttbar(TopP4Calculable) :
    def update(self,_) :
        self.value = abs(self.source[self.P4]['t'].Rapidity()) - abs(self.source[self.P4]['tbar'].Rapidity())
######################################
class DeltaYttbar(TopP4Calculable) :
    def update(self,_) :
        self.value = self.source[self.P4]['t'].Rapidity() - self.source[self.P4]['tbar'].Rapidity()
######################################
class CosEqualThetaZ(TopP4Calculable) :
    def update(self,_) :
        boost = self.source["BoostZAlt".join(self.fixes)]
        t = self.source[self.P4]['t']
        self.value = math.cos(boost(t).theta()) * (1 if boost.BetaVector().z()<0 else -1)
######################################
class CosEqualThetaZBar(TopP4Calculable) :
    def update(self,_) :
        boost = self.source["BoostZAlt".join(self.fixes)]
        t = self.source[self.P4]['tbar']
        self.value = math.cos(boost(t).theta()) * (1 if boost.BetaVector().z()<0 else -1)
######################################
class CosPhiBoost(TopP4Calculable) :
    def update(self,_) :
        p4 = self.source[self.P4]
        t = p4['t']
        ttbar = t + p4['tbar']
        beta = ttbar.BoostToCM()
        betaXYZM = beta.x(), beta.y(), beta.z(), 0
        boost = r.Math.Boost(*betaXYZM[:-1])

        self.value = math.cos(boost(p4['t']).phi() - beta.phi())
######################################
class CosThetaBoost(TopP4Calculable) :
    def update(self,_) :
        p4 = self.source[self.P4]
        t = p4['t']
        beta = (t+p4['tbar']).BoostToCM()
        betaXYZM = beta.x(),beta.y(),beta.z(),0
        boost = r.Math.Boost(*betaXYZM[:-1])
        bv = utils.LorentzV(); bv.SetPxPyPzE(*betaXYZM)

        self.value = -r.Math.VectorUtil.CosTheta(boost(t),bv)
######################################
class CosThetaBoostAlt(TopP4Calculable) :
    def update(self,_) :
        p4 = self.source[self.P4]
        t = p4['t']
        beta = (t+p4['tbar']).BoostToCM()
        betaXYZM = beta.x(),beta.y(),beta.z(),0
        boost = r.Math.Boost(*betaXYZM[:-1])
        bv = utils.LorentzV(); bv.SetPxPyPzE(-betaXYZM[0],
                                             -betaXYZM[1],
                                             +betaXYZM[2],
                                             +betaXYZM[3])

        self.value = -r.Math.VectorUtil.CosTheta(boost(t),bv)
######################################
class DirectedDeltaYttbar(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(['DeltaYttbar','SignQuarkZ'])
    def update(self,_) :
        self.value = self.source[self.SignQuarkZ] * self.source[self.DeltaYttbar]
######################################
class genttDDeltaY(wrappedChain.calculable) :
    def update(self,_) :
        iTTbar = self.source['genTopTTbar']
        if not iTTbar : self.value = None; return
        p4 = self.source['genP4']
        self.value = self.source['qDir']* ( p4[iTTbar[0]].Rapidity() - p4[iTTbar[1]].Rapidity())
######################################
class DeltaYLHadt(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(['P4','LeptonCharge'])
    def update(self,_) :
        p4 = self.source[self.P4]
        qLep = self.source[self.LeptonCharge]
        self.value =  qLep * (p4['lepton'].Rapidity() - p4['t' if qLep<0 else 'tbar'].Rapidity())
######################################
class DirectedDeltaYLHadt(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(['SignQuarkZ','DeltaYLHadt'])
    def update(self,_) : self.value = self.source[self.SignQuarkZ] * self.source[self.DeltaYLHadt]
######################################
class SignedLeptonRapidity(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(['P4',"LeptonCharge","SignQuarkZ"])
    def update(self,_) :
        self.value = (self.source[self.SignQuarkZ] *
                      self.source[self.LeptonCharge] *
                      self.source[self.P4]['lepton'].Rapidity() )
#####################################
class RelativeLeptonRapidity(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(["SumP4","P4","LeptonCharge"])
        self.moreName = "sign(y_sum) * q_lep * (y_lep-y_sum); %s%s; %s"%(collection+(SumP4,))

    def update(self,_) :
        q_lep = self.source[self.LeptonCharge]
        y_lep = self.source[self.P4]['lepton'].Rapidity()
        y_sum = self.source[self.SumP4].Rapidity()
        
        self.value = (1 if y_sum>0 else -1) * q_lep * (y_lep-y_sum)
#####################################
class DirectedLTopRapidity(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(['P4','SignQuarkZ','LeptonCharge'])
    def update(self,_) :
        lepQ = self.source[self.LeptonCharge]
        top = self.source[self.P4]['t' if lepQ>0 else 'tbar']
        self.value = self.source[self.SignQuarkZ] * top.Rapidity()
######################################
class DirectedHTopRapidity(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(['P4','SignQuarkZ','LeptonCharge'])
    def update(self,_) :
        lepQ = self.source[self.LeptonCharge]
        top = self.source[self.P4]['t' if lepQ<0 else 'tbar']
        self.value = self.source[self.SignQuarkZ] * top.Rapidity()
######################################
class SignQuarkZ(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(['P4'])
    def update(self,_) :
        self.value = 1 if self.source[self.P4]['quark'].z() > 0 else -1
#####################################
class genttOm(wrappedChain.calculable) :
    def update(self,_) :
        id = self.source['genPdgId']
        p4 = self.source['genP4']
        iHard = self.source['genIndicesHardPartons']
        iTTbar = self.source['genTopTTbar']
        beta = (p4[iTTbar[0]]+p4[iTTbar[1]]).BoostToCM()
        boost = r.Math.Boost(beta.x(),beta.y(),beta.z())
        iQs = [i for i in iHard if id[i]!=21]
        iQ = max(iQs,key=lambda i:id[i]) if iQs else iHard[0]
        iO = next(i for i in iHard if i!=iQ)
        iT = iTTbar[0]
        self.value = (r.Math.VectorUtil.CosTheta(boost(p4[iT]), boost(p4[iQ])) -
                      r.Math.VectorUtil.CosTheta(boost(p4[iT]), boost(p4[iO]))) * (1 if id[iQ] > 0 else -1)

class __genCosThetaStar__(wrappedChain.calculable) :
    def update(self,_) :
        '''Cosine angle between colliding and resulting particle in center of mass frame, three cases:

        case qqbar -> ttbarX : between light quark and top quark
        case gq(bar) -> ttbarX : between light (anti)quark and top (anti)quark
        case gg -> ttbarX : between arbitrary gluon and top quark
        '''
        id = self.source['genPdgId']
        p4 = self.source['genP4']
        iHard = self.source['genIndicesHardPartons']
        iTTbar = self.source['genTopTTbar']
        iCOM = self.iCOM()
        beta = (p4[iCOM[0]]+p4[iCOM[1]]).BoostToCM()
        boost = r.Math.Boost(beta.x(),beta.y(),beta.z())

        iQs = [i for i in iHard if id[i]!=21]
        iQ = max(iQs,key=lambda i:id[i]) if iQs else iHard[0]
        iT = self.topOrBar(*iTTbar)

        self.value = r.Math.VectorUtil.CosTheta( boost(p4[iT]), boost(p4[iQ]) ) * (1 if id[iQ]>0 else -1)
class genCosThetaStar(__genCosThetaStar__) :
    @staticmethod
    def topOrBar(iTop,iTbar) : return iTop
    def iCOM(self) : return self.source['genIndicesHardPartons']
class genCosThetaStarBar(__genCosThetaStar__) :
    @staticmethod
    def topOrBar(iTop,iTbar) : return iTbar
    def iCOM(self) : return self.source['genIndicesHardPartons']
class genttCosThetaStar(__genCosThetaStar__) :
    @staticmethod
    def topOrBar(iTop,iTbar) : return iTop
    def iCOM(self) : return self.source['genTopTTbar']
class genttCosThetaStarBar(__genCosThetaStar__) :
    @staticmethod
    def topOrBar(iTop,iTbar) : return iTbar
    def iCOM(self) : return self.source['genTopTTbar']
######################################
class __CosThetaStar__(wrappedChain.calculable) :
    def __init__(self, collection = None, topKey = 't', boostz = "BoostZ") :
        self.fixes = collection
        self.stash(['P4'])
        self.boostz = boostz.join(collection)
        self.TopKey = topKey
    def update(self,_) :
        p4 = self.source[self.P4] 
        sign = ( 1 if self.TopKey=='t' else -1)
        self.value = sign * r.Math.VectorUtil.CosTheta( self.source[self.boostz](p4[self.TopKey]),  p4['quark'] ) 
######################################
class CosThetaStarAlt(__CosThetaStar__) :
    def __init__(self, collection = None) : super(CosThetaStarAlt,self).__init__(collection,'t','BoostZAlt')
class CosThetaStar(__CosThetaStar__) :
    def __init__(self, collection = None) : super(CosThetaStar,self).__init__(collection, 't')
class CosThetaStarBar(__CosThetaStar__) :
    def __init__(self, collection = None) : super(CosThetaStarBar,self).__init__(collection, 'tbar')
class CosThetaStarAvg(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(['CosThetaStar','CosThetaStarBar'])
    def update(self,_) :
        star = self.source[self.CosThetaStar]
        bar =  self.source[self.CosThetaStarBar]
        sign = 1 if star>0 else -1
        self.value = sign - sign*math.sqrt((star-sign)*(bar-sign))
class CosThetaStarAngle(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(['CosThetaStar','CosThetaStarBar'])
    def update(self,_) :
        self.value = math.atan2(abs(self.source[self.CosThetaStar]),abs(self.source[self.CosThetaStarBar])) # on [0:pi/2]
######################################
######################################
class CosHelicityThetaL(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(['BoostZAlt','P4','LeptonCharge'])
    def update(self,_) :
        lepQ = self.source[self.LeptonCharge]
        p4 = self.source[self.P4]
        boost1 = self.source[self.BoostZAlt]
        top = boost1(p4['t' if lepQ>0 else 'tbar'])
        beta = top.BoostToCM()
        boost2 = r.Math.Boost(beta.x(), beta.y(), beta.z())
        self.value = r.Math.VectorUtil.CosTheta( top, boost2(boost1(p4['lepton'])))
######################################
class CosHelicityThetaQ(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(['BoostZAlt','P4','LeptonCharge'])
    def update(self,_) :
        lepQ = self.source[self.LeptonCharge]
        p4 = self.source[self.P4]
        boost1 = self.source[self.BoostZAlt]
        top = boost1(p4['t' if lepQ<0 else 'tbar'])
        beta = top.BoostToCM()
        boost2 = r.Math.Boost(beta.x(), beta.y(), beta.z())
        self.value = r.Math.VectorUtil.CosTheta( top, boost2(boost1(p4['q'])))
######################################
class CosThetaDaggerTT(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(['BoostZAlt','P4'])
    def update(self,_) :
        boost = self.source[self.BoostZAlt]
        p4 = self.source[self.P4]
        self.value = r.Math.VectorUtil.CosTheta(boost(p4['t']),boost(p4['tbar']))
######################################
class RadiativeCoherence(wrappedChain.calculable) :
    def __init__(self, collection = None ) :
        self.fixes = collection
        self.stash(['BoostZAlt','RecoIndex','SignQuarkZ'])

    def update(self,_) :
        jets = self.source["TopJets"]["fixes"]
        topReco = self.source["TopReconstruction"][self.source[self.RecoIndex]]
        boost = self.source[self.BoostZAlt]
        thetaTop = boost(topReco['top' if self.source[self.SignQuarkZ]>0 else 'tbar']).theta()
        p4 = self.source["CorrectedP4".join(jets)]
        extraPtTheta = [ ( p4[i].pt(), boost(p4[i]).theta() )  for i in (set(self.source["Indices".join(jets)]) - set(topReco['iPQHL']))]
        sumPtExtra = sum( pt for pt,theta in extraPtTheta )
        sumPtExtraInCones = (sum( pt for pt,theta in extraPtTheta if theta < thetaTop ) +
                             sum( pt for pt,theta in extraPtTheta if theta > math.pi - thetaTop ) )
        self.value = None if not sumPtExtra else sumPtExtraInCones/sumPtExtra - 1
######################################
######################################

class genTopP4(wrappedChain.calculable) :
    def update(self,_) :
        indices = self.source['genTTbarIndices']
        p4 = self.source['genP4']
        qqbar = self.source['genQQbar']
        qg = max(self.source['genQG'],self.source['genAG'])
        self.value = { 't':p4[indices['t']],
                       'tbar':p4[indices['tbar']],
                       'quark':p4[qqbar[0] if qqbar else qg[0] if qg else self.source['genIndicesHardPartons'][0]],
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
class genTopLeptonCharge(wrappedChain.calculable) :
    def update(self,_) : self.value = (1 if self.source['genTTbarIndices']['lplus'] else \
                                             -1 if self.source['genTTbarIndices']['lminus'] else 0)
        
class fitTopP4(wrappedChain.calculable) :
    def update(self,_) :
        reco = self.source["TopReconstruction"][0]
        tracks = self.source["CountwithPrimaryHighPurityTracks".join(self.source["TopJets"]['fixesStripped'])]
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
                      'sumP4':reco['sumP4'],
                      'key': reco['key'],
                      'chi2': reco['chi2'],
                      'hadChi2': reco['hadChi2'],
                      'ttx': reco['ttx'],
                      'fifth': reco['iX']!=None,
                      'ntracks': sum(tracks[i] for i in reco['iPQHL'])
                      }
class fitTopLeptonCharge(wrappedChain.calculable) :
    def update(self,_) :
        self.value = self.source["Charge".join(self.source["TopLeptons"])][self.source["SemileptonicTopIndex"]]

######################################
######################################
######################################

class genTopTTbar(wrappedChain.calculable) :
    def update(self,_) :
        ids = [] if self.source['isRealData'] else list(self.source['genPdgId'])
        self.value = tuple(ids.index(i) for i in [6,-6]) if all([id in ids for id in [-6,6]]) else ()
class genTopIndicesX(wrappedChain.calculable) :
    def update(self,_) :
        moms = self.source['genMotherIndex']
        ids = self.source['genPdgId']
        self.value = [i for i in range(len(moms)) if abs(ids[i])!=6 and moms[i]==4 ]
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
        mom = self.source['genMotherIndex']
        self.value = dict([(name, ids.index(i)) for name,i in [('t',6),
                                                               ('tbar',-6),
                                                               ('wplus',24),
                                                               ('wminus',-24)
                                                               ]])
        self.value.update(dict([ (w+"Child",filter(lambda i: mom[i]==self.value[w],range(len(ids)))) for w in ['wplus','wminus','t','tbar']]))
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
class genTopSemiLeptonicWithinAcceptance(wrappedChain.calculable) :
    def __init__(self, jetPtMin = None, jetAbsEtaMax = None, lepPtMin = None, lepAbsEtaMax = None) :
        for item in ['jetPtMin','jetAbsEtaMax','lepPtMin','lepAbsEtaMax'] : setattr(self,item,eval(item))
        self.moreName = 'jetPt>%0.1f; jet|eta|<%0.1f; lepPt>%0.1f; lep|eta|<%0.1f'%(jetPtMin,jetAbsEtaMax,lepPtMin,lepAbsEtaMax)
    def update(self,_) :
        self.value = False
        if not self.source["genTopTTbar"] : return
        indices = self.source['genTTbarIndices']
        if not indices['semi'] : return
        iLep = max(indices['lplus'],indices['lminus'])
        genP4 = self.source["genP4"]
        if genP4[iLep].pt() < self.lepPtMin : return
        if abs(genP4[iLep].eta()) > self.lepAbsEtaMax : return
        for iJet in ( indices['q'] + [indices['b'],indices['bbar']] ) :
            if genP4[iJet].pt() < self.jetPtMin : return
            if abs(genP4[iJet].eta()) > self.jetAbsEtaMax : return
        self.value = True
######################################
class genTopSemiMu(wrappedChain.calculable) :
    def update(self,_) :
        ids = self.source['genPdgId']
        iTT = self.source['genTTbarIndices']
        self.value = iTT['semi'] and abs(ids[max(iTT['lplus'],iTT['lminus'])])==13
######################################
class genTopSemiLeptonicAccepted(wrappedChain.calculable) :
    def update(self,_) :
        jets = self.source["TopJets"]["fixes"]
        self.value = len(self.source["IndicesGenB".join(jets)]) is 2 is len(self.source["IndicesGenWqq".join(jets)])
######################################
class mixedSumP4(wrappedChain.calculable) :
    def __init__(self, transverse = None, longitudinal = None) :
        self.trans = transverse
        self.longi = longitudinal
        self.moreName = "transverse: %s ; longitudinal: %s" % (transverse,longitudinal)
        self.value = utils.LorentzV()
    def update(self,_) :
        trans = self.source[self.trans]
        longi = self.source[self.longi]
        f = trans.pt() / longi.pt()
        self.value.SetPxPyPzE(-trans.px(),-trans.py(), f*longi.pz(), f*longi.E())
#####################################
class SemileptonicTopIndex(wrappedChain.calculable) :
    def update(self,_) :
        self.value = next( iter(self.source["IndicesAnyIsoIsoOrder".join(self.source["TopLeptons"])]), None )
#####################################
class TopReconstruction(wrappedChain.calculable) :
    def __init__(self, bscale = 1.1, eCoupling = 0.55, v2had = False, v2lep = True ) :
        theta = math.pi/6
        self.ellipseR = np.array([[math.cos(theta),-math.sin(theta)],[math.sin(theta), math.cos(theta)]])
        self.epsilon = 1e-7
        for item in ['bscale',    # factor by which to scale hypothesized b jets
                     'eCoupling', # percentage of jet resolution used to sharpen MET resolution
                     'v2had',     # v2 (1parameter,3residuals) is twice as fast as v1 (3parameters,5residuals) but 5% less accurate
                     'v2lep'      # no comment: read the codes
                     ] : setattr(self,item,eval(item))
        self.moreName = "bscale:%.1f; eCoupl:%.2f; v%dhad; v%dlep"%(bscale,eCoupling,v2had+1,v2lep+1)

    def update(self,_) :
        
        jets = dict( (item, self.source[item.join(self.source["TopJets"]["fixes"])] )
                     for item in ["CorrectedP4","IndicesBtagged","Indices","Resolution","CovariantResolution2","ComboPQBDeltaRawMassWTop"] )

        lepton = dict( (item, self.source[item.join(self.source['TopLeptons'])][self.source["SemileptonicTopIndex"]])
                       for item in ["Charge","P4"])

        topP = self.source["TopComboQQBBProbability"]
        bIndices = jets["IndicesBtagged"][:5] #consider only the first few b-tagged jets as possible b-candidates
        
        recos = []
        for iPQH in itertools.permutations(jets["Indices"],3) :
            if iPQH[0]>iPQH[1] : continue
            if iPQH[2] not in bIndices : continue
            if np.dot(*(2*[self.ellipseR.dot(jets["ComboPQBDeltaRawMassWTop"][iPQH]) / [35,70]])) > 1 : continue # elliptical window on raw masses

            hadFit = utils.fitKinematic.leastsqHadronicTop2(*zip(*((jets["CorrectedP4"][i]*(self.bscale if i==2 else 1), jets["Resolution"][i]) for i in iPQH)) ) if self.v2had else \
                     utils.fitKinematic.leastsqHadronicTop( *zip(*((jets["CorrectedP4"][i]*(self.bscale if i==2 else 1), jets["Resolution"][i]) for i in iPQH)), widthW = 4./2 ) #tuned w width

            sumP4 = self.source["mixedSumP4"] - hadFit.rawT + hadFit.fitT
            nuXY = -np.array([sumP4.x(), sumP4.y()])
            nuErr2 = sum([-self.eCoupling*jets["CovariantResolution2"][i] for i in iPQH], self.source["metCovariancePF"])

            for iL in set(bIndices)-set(iPQH) :
                iPQHL = iPQH+(iL,)
                iQQBB = iPQHL[:2]+tuple(sorted(iPQHL[2:]))
                b = jets["CorrectedP4"][iL]
                nuXY_b = nuXY - (self.bscale - 1)*np.array([b.y(),b.y()])
                nuErr2_b = nuErr2-self.eCoupling*jets["CovariantResolution2"][iL]
                lepFit = utils.fitKinematic.leastsqLeptonicTop2( b*self.bscale, jets["Resolution"][iL], lepton["P4"], nuXY_b, nuErr2_b) if self.v2lep else \
                         min( utils.fitKinematic.leastsqLeptonicTop( b*self.bscale, jets["Resolution"][iL], lepton["P4"], nuXY_b, nuErr2_b, zPlus = True ),
                              utils.fitKinematic.leastsqLeptonicTop( b*self.bscale, jets["Resolution"][iL], lepton["P4"], nuXY_b, nuErr2_b, zPlus = False ),
                              key = lambda x: x.chi2 )
                tt = hadFit.fitT + lepFit.fitT
                iX,ttx = min( [(None,tt)]+[(i,tt+jets["CorrectedP4"][i]) for i in jets["Indices"] if i not in iPQHL], key = lambda lv : lv[1].pt() )
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
                               "sumP4": sumP4 - b + lepFit.fitB,
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
        jets = self.source["TopJets"]['fixes']
        p4 = self.source["CorrectedP4".join(jets)]
        indices = self.source["Indices".join(jets)]
        igenTT = self.source["genTTbarIndices"]
        iTop,iTbar = self.source["genTopTTbar"]

        self.value = {}
        if self.source["ttDecayMode"] not in ["ej","mj"] : return
        self.value["met"] = (abs(r.Math.VectorUtil.DeltaPhi( -self.source["mixedSumP4"], genP4[igenTT['nu']])) < 0.7)
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
class fitTopRecoIndex(wrappedChain.calculable) :
    value = 0
    def update(self,_) : pass
class genTopRecoIndex(wrappedChain.calculable) :
    def __init__(self,rMax = 0.6, rMaxNu = 10.0) :
        for item in ['lep','bLep','bHad','q'] : setattr(self,"rMax"+item, rMax )
        self.rMaxnu = rMaxNu
        self.moreName = "deltaR[lep,b,b,q,q] <%0.1f; deltaRnum<%0.1f"%(rMax,rMaxNu)
    def update(self,_) :
        self.value = -1
        if not self.source['genTopTTbar'] : return
        iPQHL = self.source['IndicesGenTopPQHL'.join(self.source['TopJets']['fixes'])]
        iPass = [ i for i,R in enumerate(self.source['TopReconstruction'])
                  if R['iPQHL']==iPQHL and all( self.source["%sDeltaRTopRecoGen"%s][i]<getattr(self,"rMax%s"%s)
                                                for s in ['lep','nu'] ) ]

        if len(iPass) :
            self.value = sorted( iPass, key = lambda i: sum([self.source['%sDeltaRTopRecoGen'%s][i] for s in ['lep','nu']]))[0]

class IndicesGenTopPQHL(wrappedChain.calculable) :
    def __init__(self, jets=None, rMax = 0.6 ) :
        self.rMax = rMax
        self.fixes = jets
        self.stash( ['Indices','CorrectedP4'] )

    def update(self,_) :
        genP4 = self.source['genP4']
        iTT = self.source['genTTbarIndices']

        PQHL = [genP4[i] if i!=None else None for i in ( (iTT['q'][:2] if iTT['q'] else 2*[None]) + [ iTT['bhad'],iTT['blep'] ] )]

        indices = self.source[self.Indices]
        p4 = self.source[self.CorrectedP4]

        dRIs = [ min( (r.Math.VectorUtil.DeltaR( p4[i], gen ), i) for i in indices ) if gen else
                 (None,None)
                 for gen in PQHL ]

        PQHL = [i if dR<self.rMax else None for dR,i in dRIs ]
        self.value = tuple( sorted(PQHL[:2]) + PQHL[2:] )

class IndicesGenTopExtra(wrappedChain.calculable) :
    def __init__(self, jets=None, rMax = 0.6) :
        self.rMax = rMax
        self.fixes = jets
        self.stash(['Indices','CorrectedP4'])

    def update(self,_) :
        imom = self.source['genMotherIndex']
        p4 = self.source['genP4']
        pdg = self.source['genPdgId']
        status = self.source['genStatus']

        extraP4 = [p4[i] for i in range(8,len(imom)) if 2<imom[i]<6 and abs(pdg[i]) in [1,2,3,4,5,11,13,15,21] and status[i]==3]

        indices = self.source[self.Indices]
        jet = self.source[self.CorrectedP4]
        self.value = [j for j in indices if any( self.rMax > r.Math.VectorUtil.DeltaR(jet[j],gen) for gen in extraP4 ) ]

######################################
class __HardAsym__(wrappedChain.calculable) :
    def update(self,_) :
        self.value = None
        boxes,Is = zip(*[(self.source[item],i) for item,i in self.indices])
        if not all(boxes) : return
        indices = [box[i] for box,i in zip(boxes,Is)]
        p4s = self.source['genP4']
        self.hard.setMomenta([p4s.at(i) for i in indices], check=False)
        self.value = self.hard
class __wHardAsym__(wrappedChain.calculable) :
    def __init__(self,target,nominal) :
        self.fixes = ("","%+.2fon%+.2f"%(target,nominal))
        self.target = target
        self.nominal = nominal
    def update(self,_) :
        h = self.source[self.hardname()]
        self.value = 1 if not h else h.weight(self.target,self.nominal)
class QQbarHardAsym(__HardAsym__) :
    def __init__(self) :
        self.hard = utils.asymmWeighting.Asymm_qqbar_hard()
        self.indices = [('genQQbar',0),('genQQbar',1),('genTopTTbar',0),('genTopTTbar',1),('genGlu',0)]
class GqHardAsym(__HardAsym__) :
    def __init__(self) :
        self.hard = utils.asymmWeighting.Asymm_gq_hard()
        self.indices = [('genQG',1),('genQG',0),('genTopTTbar',0),('genTopTTbar',1),('genQuark',0)] # broken sincy intro of wAG
class GqbarHardAsym(__HardAsym__) :
    def __init__(self) :
        self.hard = utils.asymmWeighting.Asymm_gqbar_hard()
        self.indices = [('genQG',1),('genQG',0),('genTopTTbar',0),('genTopTTbar',1),('genQuark',0)] # broken sincy intro of wAG
class wQQbarHardAsym(__wHardAsym__) :
    def hardname(self) : return 'QQbarHardAsym'
class wQgHardAsym(__wHardAsym__) :
    def hardname(self) :
        return 'GqHardAsym' if self.source['genQG'][0]>0 else 'GqbarHardAsym' # broken sincy intro of wAG
class hardSymmAnti(wrappedChain.calculable) :
    def update(self,_) :
        qq = self.source['genQQbar']
        qg = self.source['genQG'] # broken sincy intro of wAG
        if not (qq or qg) : self.value = () ; return
        h = self.source['QQbarHardAsym' if qq else 'GqHardAsym' if qg[0]>0 else 'GqbarHardAsym']
        self.value = () if not h else (h.symm,h.anti)
######################################
class TopComboQQBBLikelihood(wrappedChain.calculable) :
    def __init__(self, tag = None) :
        self.tagProbabilityGivenBQN = tag+'ProbabilityGivenBQN'

    def update(self,_) :
        jets = self.source["TopJets"]["fixes"]
        indices = self.source["Indices".join(jets)]
        B,Q,N = zip(*self.source[self.tagProbabilityGivenBQN.join(jets)])
        self.value = {}
        for iPQHL in itertools.permutations(indices,4) :
            if iPQHL[0]>iPQHL[1] : continue
            if iPQHL[2]>iPQHL[3] : continue
            self.value[iPQHL] = reduce(operator.mul, ([Q[i] for i in iPQHL[:2]] +
                                                      [B[i] for i in iPQHL[2:]]  +
                                                      [N[k] for k in indices if k not in iPQHL]) )
######################################
class TopComboQQBBProbability(wrappedChain.calculable) :
    def update(self,_) :
        likelihoods = self.source['TopComboQQBBLikelihood']
        sumL = max(1e-20,sum(likelihoods.values()))
        self.value = dict([(key,val/sumL) for key,val in likelihoods.iteritems()])
######################################
class TopComboQQBBMaxProbability(wrappedChain.calculable) :
    def update(self,_) : self.value = max(self.source["TopComboQQBBProbability"].values())
######################################
class OtherJetsLikelihood(wrappedChain.calculable) :
    def __init__(self, tag = None) :
        self.tagProbabilityGivenBQN = tag+'ProbabilityGivenBQN'

    def update(self,_) :
        jets = self.source["TopJets"]["fixes"]
        indices = self.source["Indices".join(jets)]
        B,Q,N = zip(*self.source[self.tagProbabilityGivenBQN.join(jets)])
        self.value = reduce(operator.mul, [N[k] for k in indices])
######################################
class TopRatherThanWProbability(wrappedChain.calculable) :
    def __init__(self, priorTop = 0.05) :
        self.priorTop = priorTop
        self.invPriorTopMinusOne =  ( 1.0 / priorTop  - 1)
        self.moreName = "priorTop = %0.3f"%priorTop
        
    def update(self,_) :
        topLikes = self.source["TopComboQQBBLikelihood"]
        if not topLikes : self.value = self.priorTop; return
        topL = sum(topLikes.values()) / float(len(topLikes))
        wL = self.source["OtherJetsLikelihood"]
        denom = (topL + wL * self.invPriorTopMinusOne)
        self.value = (topL / denom) if denom else self.priorTop
######################################
class BMomentsSum2(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(["RecoIndex"])
    def update(self,_) :
        _,__,iH,iL = self.source['TopReconstruction'][self.source[self.RecoIndex]]['iPQHL']
        jets = self.source["TopJets"]["fixesStripped"]
        phi2 = self.source["Phi2Moment".join(jets)]
        eta2 = self.source["Eta2Moment".join(jets)]
        self.value = phi2[iH]+phi2[iL]+eta2[iH]+eta2[iL]
