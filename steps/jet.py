import ROOT as r
from supy import analysisStep
#####################################
class singleJetHistogrammer(analysisStep) :

    def __init__(self,cs, maxIndex = 2) :
        self.cs = cs
        self.csbase = (cs[0].replace("xc",""),cs[1])
        self.maxIndex = maxIndex
        self.moreName="%s%s through index %d" % (self.cs+(maxIndex,))
        self.indicesName = "%sIndices%s" % self.cs
        self.p4sName = "%sCorrectedP4%s" % self.cs

    def uponAcceptance (self,eventVars) :
        p4s = eventVars[self.p4sName]
        cleanJetIndices = eventVars[self.indicesName]
        phi2mom = eventVars["%sPhi2Moment%s"%self.csbase]
        eta2mom = eventVars["%sEta2Moment%s"%self.csbase]

        self.book.fill( len(cleanJetIndices), "jetMultiplicity", 10, -0.5, 9.5,
                        title=";number of %s%s passing ID#semicolon p_{T}#semicolon #eta cuts;events / bin"%self.cs)
        
        for i,iJet in enumerate(cleanJetIndices) :
            jet = p4s.at(iJet)
            pt = jet.pt()
            eta = jet.eta()
            phi2 = phi2mom.at(iJet)
            eta2 = eta2mom.at(iJet)
            mom2Max = 0.1
            jetLabel = str(i+1) if i <= self.maxIndex else "_ge%d"%(self.maxIndex+2)

            self.book.fill(eta2,  "%s%s%sEta2mom" %(self.cs+(jetLabel,)), 50,  0.0, mom2Max, title=";jet%s #sigma_{#eta}^{2};events / bin"%jetLabel)
            self.book.fill(phi2,  "%s%s%sPhi2mom" %(self.cs+(jetLabel,)), 50,  0.0, mom2Max, title=";jet%s #sigma_{#phi}^{2};events / bin"%jetLabel)
            self.book.fill(pt,  "%s%s%sPt" %(self.cs+(jetLabel,)), 50,  0.0, 500.0, title=";jet%s p_{T} (GeV);events / bin"%jetLabel)
            self.book.fill(eta, "%s%s%seta"%(self.cs+(jetLabel,)), 50, -5.0,   5.0, title=";jet%s #eta;events / bin"%jetLabel)
            if i>self.maxIndex: continue
            for j,jJet in list(enumerate(cleanJetIndices))[i+1:self.maxIndex+1] :
                self.book.fill(abs(r.Math.VectorUtil.DeltaPhi(jet,p4s.at(jJet))), "%s%sdphi%d%d"%(self.cs+(i+1,j+1)), 50,0, r.TMath.Pi(),
                               title = ";#Delta#phi_{jet%d,jet%d};events / bin"%(i+1,j+1))
#####################################
