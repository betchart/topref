import math,operator,itertools,ROOT as r
from supy import wrappedChain,calculables,utils
try: import numpy as np
except: np = None

#####################################
class TopCandidateIndices(wrappedChain.calculable) :
    '''4-tuples of jet indices.

    P,Q : indices of jets from W, from hadronically decaying top.
    H : index of b-jet from hadronically decaying top.
    L : index of b-jet from leptonically decaying top.
    '''

    def update(self,_) :
        indices = sorted(self.source['Indices'.join(self.source['TopJets'])])
        self.value = sorted([ PQ+HL
                              for PQ in itertools.combinations(indices,2)
                              for HL in itertools.permutations(indices,2)
                              if len(set(PQ+HL))==4 ])

class HTopCandidateIndices(wrappedChain.calculable) :
    '''3-tuples of jet indices from hadronically decaying top, last index the b-jet.'''

    def update(self,_) :
        self.value = sorted(set(iPQHL[:3] for iPQHL in self.source['TopCandidateIndices']))


#####################################
class TopFitLikelihoodCorrectIndex(wrappedChain.calculable) :
    def update(self,_):
        L = self.source['TopCandidateLikelihood']
        iPQHL = self.source['TopReconstruction'][0]['iPQHL']
        iPQHL_true = self.source['IndicesGenTopPQHL']
        self.value = sorted(L, key = L.__getitem__, reverse = True).index(iPQHL) if iPQHL in L and iPQHL==iPQHL_true else -1

class TopFitLikelihoodIndex(wrappedChain.calculable) :
    def update(self,_):
        L = self.source['TopCandidateLikelihood']
        iPQHL = self.source['TopReconstruction'][0]['iPQHL']
        self.value = sorted(L, key = L.__getitem__, reverse = True).index(iPQHL) if iPQHL in L else -1

class TopGenLikelihoodIndex(wrappedChain.calculable) :
    def update(self,_):
        L = self.source['TopCandidateLikelihood']
        iPQHL = self.source['IndicesGenTopPQHL']
        self.value = sorted(L, key = L.__getitem__, reverse = True).index(iPQHL) if iPQHL in L else -1


######################################
class TopCandidateLikelihood(wrappedChain.calculable) :
    def update(self,_) :
        self.value = {}
        qqbbL = self.source['TopComboQQBBLikelihood']
        sigmasLR = self.source['HTopSigmasPQBCombinationsLR']
        unfitX2LR = self.source['LTopUnfitSqrtChi2CombinationsLR']
        for iPQHL in self.source['TopCandidateIndices'] :
            self.value[iPQHL] = reduce(operator.mul, [qqbbL[ iPQHL[:2]+tuple(sorted(iPQHL[2:])) ],
                                                      sigmasLR[ iPQHL[:3] ],
                                                      unfitX2LR[ iPQHL[3] ]
                                                      ])


######################################
class LTopUnfitSqrtChi2(wrappedChain.calculable) :
    def update(self,_) :
        jets = self.source['TopJets']
        jP4s = self.source['AdjustedP4'.join(jets)]
        bscale = self.source['BScaling'.join(jets)]
        lP4 = self.source['P4'.join(self.source['TopLeptons'])][self.source['SemileptonicTopIndex']]
        met = self.source['metAdjustedP4']
        nuXY = np.array([met.px(),met.py()])
        nuErr2 = self.source['metCovariance']

        self.value = dict( (iL,
                            math.sqrt(utils.fitKinematic.leastsqLeptonicTop2( jP4s[iL]*bscale[iL], 0, lP4, nuXY, nuErr2 ).chi2) )
                           for iL in self.source['Indices'.join(jets)] )

#####################################
class RawMassWTopPQB(wrappedChain.calculable) :
    def update(self,_) :
        jets = self.source['TopJets']
        p4 = self.source['AdjustedP4'.join(jets)]
        bscale = self.source['BScaling'.join(jets)]
        self.value = {}
        for iPQB in self.source['HTopCandidateIndices'] :
            _,W,t = np.cumsum([p4[i]*(bscale[i] if i==2 else 1) for i in iPQB])
            self.value[iPQB] = (W.M(),t.M())

class RawMassWTopCorrectPQB(wrappedChain.calculable) :
    def update(self,_) :
        iPQB = self.source['IndicesGenTopPQH']
        self.value = self.source['RawMassWTopPQB'][iPQB] if len(set(iPQB))==3 and None not in iPQB else ()

class HTopSigmasPQB(wrappedChain.calculable) :
    matrix = None
    def update(self,_) :
        if self.matrix==None:
            self.matrix = dict.__getitem__(self.source.someDict if type(self.source).__name__=='keyTracer' else self.source,
                                           'RawMassWTopCorrectPQBTwoDChiSquared').matrix
        rawM = self.source['RawMassWTopPQB']
        self.value = dict( ( iPQH, math.sqrt(np.dot( v, np.dot(self.matrix, v) )) ) for iPQH,v in
                           [( i3, ms+(1,) ) for i3,ms in rawM.items()])

class ProbabilityHTopMasses(wrappedChain.calculable) :
    def update(self,_):
        LRs = self.source['HTopSigmasPQBCombinationsLR'].values()
        sumLRs = sum(LRs)
        self.value = sumLRs / (sumLRs + len(LRs) )


######################################
class TopComboQQBBLikelihood(wrappedChain.calculable) :
    def __init__(self, tag = None) :
        self.tagProbabilityGivenBQN = tag+'ProbabilityGivenBQN'

    def update(self,_) :
        self.value = {}
        jets = self.source["TopJets"]
        indices = self.source["Indices".join(jets)]
        B,Q,N = zip(*self.source[self.tagProbabilityGivenBQN.join(jets)])
        for iPQHL in self.source["TopCandidateIndices"] :
            if iPQHL[2] > iPQHL[3] : continue
            self.value[iPQHL] = reduce(operator.mul, ([Q[i] for i in iPQHL[:2]] +
                                                      [B[i] for i in iPQHL[2:]]  +
                                                      [N[k] for k in indices if k not in iPQHL]) )

class TopComboQQBBProbability(wrappedChain.calculable) :
    def update(self,_) :
        likelihoods = self.source['TopComboQQBBLikelihood']
        sumL = max(1e-20,sum(likelihoods.values()))
        self.value = dict([(key,val/sumL) for key,val in likelihoods.iteritems()])

class OtherJetsLikelihood(wrappedChain.calculable) :
    def __init__(self, tag = None) :
        self.tagProbabilityGivenBQN = tag+'ProbabilityGivenBQN'

    def update(self,_) :
        jets = self.source["TopJets"]
        indices = self.source["Indices".join(jets)]
        B,Q,N = zip(*self.source[self.tagProbabilityGivenBQN.join(jets)])
        self.value = reduce(operator.mul, [N[k] for k in indices])

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
