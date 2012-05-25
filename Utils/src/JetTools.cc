#include "MitPhysics/Utils/interface/JetTools.h"
#include <algorithm>
#include <vector>

ClassImp(mithep::JetTools)

using namespace mithep;

 
JetTools::JetTools()
{
   // Constructor
}

JetTools::~JetTools() 
{
  // Destructor.
}

//Remember to remove the signal from particles before inputting into the function
Double_t JetTools::NJettiness(const ParticleOArr *particles, const JetOArr *jets, double Q, double Y){
  if(particles->GetEntries() <= 0) return 0.0;

  Double_t fval = 0.0;
  Double_t fvalpart;

  for(int i=0;i<int(particles->GetEntries());i++){
    fvalpart = (particles->At(i)->Pt()) * TMath::Exp(-TMath::Abs(particles->At(i)->Eta()-Y)); 

    for(int j=0;j<int(jets->GetEntries());j++){
      fvalpart = TMath::Min(fvalpart,(jets->At(j)->Pt()) * 
                 (2 * TMath::CosH(TMath::Abs(jets->At(j)->Eta()-particles->At(i)->Eta()))
                  - 2 * TMath::Cos(fabs(MathUtils::DeltaPhi(jets->At(j)->Phi(),particles->At(i)->Phi())))));
    }
    fval = fval + fvalpart;
  }

  fval = fval / Q;

  return fval;
}

Double_t JetTools::NJettiness(const PFCandidateOArr *pfCandidates, const JetOArr *jets, double Q, double Y){
  if(pfCandidates->GetEntries() <= 0) return 0.0;

  Double_t fval = 0.0;
  Double_t fvalpart;

  for(int i=0;i<int(pfCandidates->GetEntries());i++){
    fvalpart = (pfCandidates->At(i)->Pt()) * TMath::Exp(-TMath::Abs(pfCandidates->At(i)->Eta()-Y)); 

    for(int j=0;j<int(jets->GetEntries());j++){
      fvalpart = TMath::Min(fvalpart,(jets->At(j)->Pt()) * 
                 (2 * TMath::CosH(TMath::Abs(jets->At(j)->Eta()-pfCandidates->At(i)->Eta()))
                  - 2 * TMath::Cos(fabs(MathUtils::DeltaPhi(jets->At(j)->Phi(),pfCandidates->At(i)->Phi())))));
    }
    fval = fval + fvalpart;
  }

  fval = fval / Q;

  return fval;
}

Double_t JetTools::NJettiness(const TrackOArr *tracks, const JetOArr *jets, double Q, double Y){
  if(tracks->GetEntries() <= 0) return 0.0;

  Double_t fval = 0.0;
  Double_t fvalpart;

  for(int i=0;i<int(tracks->GetEntries());i++){
    fvalpart = (tracks->At(i)->Pt()) * TMath::Exp(-TMath::Abs(tracks->At(i)->Eta()-Y));     

    for(int j=0;j<int(jets->GetEntries());j++){
      fvalpart = TMath::Min(fvalpart,(jets->At(j)->Pt()) * 
                 (2 * TMath::CosH(TMath::Abs(jets->At(j)->Eta()-tracks->At(i)->Eta()))
                  - 2 * TMath::Cos(fabs(MathUtils::DeltaPhi(jets->At(j)->Phi(),tracks->At(i)->Phi())))));
    }
    fval = fval + fvalpart;
  }

  fval = fval / Q;
  
  return fval;
}

Double_t JetTools::NJettiness(const JetOArr *jetsS, const JetOArr *jets, double Q, double Y){
  if(jetsS->GetEntries() <= 0) return 0.0;

  Double_t fval = 0.0;
  Double_t fvalpart;

  for(int i=0;i<int(jetsS->GetEntries());i++){
    fvalpart = (jetsS->At(i)->Pt()) * TMath::Exp(-TMath::Abs(jetsS->At(i)->Eta()-Y));     

    for(int j=0;j<int(jets->GetEntries());j++){
      fvalpart = TMath::Min(fvalpart,(jets->At(j)->Pt()) * 
                 (2 * TMath::CosH(TMath::Abs(jets->At(j)->Eta()-jetsS->At(i)->Eta()))
                  - 2 * TMath::Cos(fabs(MathUtils::DeltaPhi(jets->At(j)->Phi(),jetsS->At(i)->Phi())))));
    }
    fval = fval + fvalpart;
  }

  fval = fval / Q;
  
  return fval;
}

Double_t JetTools::NJettiness(const CaloTowerOArr *calos, const JetOArr *jets, double Q, double Y){
  if(calos->GetEntries() <= 0) return 0.0;

  Double_t fval = 0.0;
  Double_t fvalpart;

  for(int i=0;i<int(calos->GetEntries());i++){
    fvalpart = (calos->At(i)->Pt()) * TMath::Exp(-TMath::Abs(calos->At(i)->Eta()-Y));     

    for(int j=0;j<int(jets->GetEntries());j++){
      fvalpart = TMath::Min(fvalpart,(jets->At(j)->Pt()) * 
                 (2 * TMath::CosH(TMath::Abs(jets->At(j)->Eta()-calos->At(i)->Eta()))
                  - 2 * TMath::Cos(fabs(MathUtils::DeltaPhi(jets->At(j)->Phi(),calos->At(i)->Phi())))));
    }
    fval = fval + fvalpart;
  }

  fval = fval / Q;
  
  return fval;
}

//M_r
Double_t JetTools::M_r(const ParticleOArr *particles){

  if(particles->GetEntries() < 2) return -999.;

  Double_t E0  = particles->At(0)->E();
  Double_t E1  = particles->At(1)->E();
  Double_t Pz0 = particles->At(0)->Pz();
  Double_t Pz1 = particles->At(1)->Pz();

  Double_t den =  TMath::Power(Pz0-Pz1, 2) - TMath::Power(E0-E1,2);
  if(den <= 0) return -100.;

  return 2.0*TMath::Sqrt(TMath::Power(E0*Pz1 - E1*Pz0, 2)/den);
}

//Beta_r
Double_t JetTools::Beta_r(const ParticleOArr *particles){

  if(particles->GetEntries() < 2) return -999.;

  Double_t E0  = particles->At(0)->E();
  Double_t E1  = particles->At(1)->E();
  Double_t Pz0 = particles->At(0)->Pz();
  Double_t Pz1 = particles->At(1)->Pz();

  return (E0-E1)/(Pz0-Pz1);
}

//M_r_t
Double_t JetTools::M_r_t(const ParticleOArr *particles, const Met *met){

  if(particles->GetEntries() < 2) return -999.;

  Double_t Pt0    = particles->At(0)->Pt();
  Double_t Pt1    = particles->At(1)->Pt();
  Double_t etmiss = met->Pt();

  Double_t Px0    = particles->At(0)->Px();
  Double_t Px1    = particles->At(1)->Px();
  Double_t metx   = met->Px();
  Double_t Py0    = particles->At(0)->Py();
  Double_t Py1    = particles->At(1)->Py();
  Double_t mety   = met->Py();

  return TMath::Sqrt(0.5*etmiss*(Pt0 + Pt1) - 0.5*(metx*(Px0 + Px1) + mety*(Py0 + Py1)));
}

//Razor
Double_t JetTools::Razor(const ParticleOArr *particles, const Met *met){
  if(particles->GetEntries() < 2) return -999.;

  Double_t mr  = M_r(particles);
  Double_t mrt = M_r_t(particles,met);
  
  if(mr != 0) return mrt/mr;
  
  return -999.;
}

//Cosine Omega
Double_t JetTools::CosineOmega(const Particle *particles0, const Particle *particles1){

  TLorentzVector v_L1(particles0->Px(),particles0->Py(),particles0->Pz(),particles0->E());
  TLorentzVector v_L2(particles1->Px(),particles1->Py(),particles1->Pz(),particles1->E());

  Double_t beta = (v_L1.P()-v_L2.P())/(v_L1.Pz()-v_L2.Pz());

  TVector3 B;
  B.SetXYZ(0.0,0.0,-1.0*beta);

  v_L1.Boost(B);
  v_L2.Boost(B);

  Double_t cosomega = v_L1.Vect().Dot(v_L2.Vect())/(v_L1.P()*v_L2.P());

  return cosomega;
}

//Transverse Higgs mass
Double_t JetTools::MtHiggs(const ParticleOArr * leptons,
                           const Met *met, double metFraction[2], int nsel){
  if(leptons->Entries() < 2) return -999.0;

  double mtHiggs = -999.0;
  double enell = 0.0;
  double enenn = 0.0;
  double enex  = 0.0;
  double eney  = 0.0;
  double mll   = 0.0;
  double mnu   = 0.0;
  CompositeParticle *dilepton = new CompositeParticle();
  dilepton->AddDaughter(leptons->At(0));
  dilepton->AddDaughter(leptons->At(1));
  
  if     (nsel == 0){ // Use of Mt mass and mnu == mll
    enell = TMath::Sqrt(dilepton->Pt()*dilepton->Pt() + dilepton->Mt()*dilepton->Mt());
    enenn = TMath::Sqrt(met->Pt() *met->Pt()  + dilepton->Mt()*dilepton->Mt());
    enex  = dilepton->Px() + met->Px();
    eney  = dilepton->Py() + met->Py();
    mll   = dilepton->Mass();
    mnu   = mll;
  }
  else if(nsel == 1){ // Use of Mt mass and mnu == 0
    enell = TMath::Sqrt(dilepton->Pt()*dilepton->Pt() + dilepton->Mt()*dilepton->Mt());
    enenn = TMath::Sqrt(met->Pt() *met->Pt()  + 0.0*0.0);
    enex  = dilepton->Px() + met->Px();
    eney  = dilepton->Py() + met->Py();
    mll   = dilepton->Mass();
    mnu   = 0.0;
  }
  else if(nsel == 2){ // Use of M mass and mnu == mll
    enell = TMath::Sqrt(dilepton->Pt()*dilepton->Pt() + dilepton->Mass()*dilepton->Mass());
    enenn = TMath::Sqrt(met->Pt() *met->Pt()  + dilepton->Mass()*dilepton->Mass());
    enex  = dilepton->Px() + met->Px();
    eney  = dilepton->Py() + met->Py();
    mll   = dilepton->Mass();
    mnu   = mll;
  }
  else if(nsel == 3){ // Use of M mass and mnu == 0
    enell = TMath::Sqrt(dilepton->Pt()*dilepton->Pt() + dilepton->Mass()*dilepton->Mass());
    enenn = TMath::Sqrt(met->Pt() *met->Pt()  + 0.0*0.0);
    enex  = dilepton->Px() + met->Px();
    eney  = dilepton->Py() + met->Py();
    mll   = dilepton->Mass();
    mnu   = 0.0;
  }
  else if(nsel == 4){ // Use of Mt mass and replacing mnu using the met optimal
    enell = TMath::Sqrt(dilepton->Pt()*dilepton->Pt() + dilepton->Mt()*dilepton->Mt());
    enenn = TMath::Sqrt(met->Pt() *met->Pt()  + 0.0*0.0);
    enex  = dilepton->Px() + met->Px();
    eney  = dilepton->Py() + met->Py();
    mll   = dilepton->Mass();
    double metAuxPx[2] = {met->Px() * metFraction[0],
	    		  met->Px() * (1.0 - metFraction[0])};
    double metAuxPy[2] = {met->Py() * metFraction[1],
	 		  met->Py() * (1.0 - metFraction[1])};
    double ene = TMath::Sqrt(metAuxPx[0]*metAuxPx[0]+metAuxPy[0]*metAuxPy[0]) +
		 TMath::Sqrt(metAuxPx[1]*metAuxPx[1]+metAuxPy[1]*metAuxPy[1]);
    double px = metAuxPx[0] + metAuxPx[1];
    double py = metAuxPy[0] + metAuxPy[1];
    mnu = TMath::Sqrt(ene*ene - px*px - py*py);
  }
  else if(nsel == 5){ // Using the optimal met value
    double metAuxPx[2] = {met->Px() * metFraction[0],
	    		  met->Px() * (1.0 - metFraction[0])};
    double metAuxPy[2] = {met->Py() * metFraction[1],
	 		  met->Py() * (1.0 - metFraction[1])};
    double ene = leptons->At(0)->Pt() + leptons->At(1)->Pt() +
                 TMath::Sqrt(metAuxPx[0]*metAuxPx[0]+metAuxPy[0]*metAuxPy[0]) +
		 TMath::Sqrt(metAuxPx[1]*metAuxPx[1]+metAuxPy[1]*metAuxPy[1]);
    double px = leptons->At(0)->Px() + leptons->At(1)->Px() +
                metAuxPx[0] + metAuxPx[1];
    double py = leptons->At(0)->Py() + leptons->At(1)->Py() +
                metAuxPy[0] + metAuxPy[1];
    mtHiggs = ene*ene - px*px - py*py;
  }
  else if(nsel == 6){ // Use the formula from hep-ph:1006.4998
    mtHiggs = 2*leptons->At(0)->Pt()*leptons->At(0)->Pt() + 2*leptons->At(1)->Pt()*leptons->At(1)->Pt() + 3 * (
      leptons->At(0)->Pt()*leptons->At(1)->Pt() + met->Pt()*(leptons->At(0)->Pt()+leptons->At(1)->Pt())
      - met->Px()*dilepton->Px() - met->Py()*dilepton->Py()
      - leptons->At(0)->Px()*leptons->At(1)->Px() - leptons->At(0)->Py()*leptons->At(1)->Py());
  }
  else if(nsel == 7){ // Use of M mass and mnu == 0
    double deltaPhiDileptonMet = fabs(MathUtils::DeltaPhi(dilepton->Phi(), 
                                                          met->Phi()));
    mtHiggs = 2.0*dilepton->Pt()*met->Pt()*(1.0 - cos(deltaPhiDileptonMet));
  }

  if(nsel >= 0 && nsel <= 4){
    mtHiggs = mll*mll + mnu*mnu + 2.0*(enell*enenn - enex*enex - eney*eney);
  }

  if(mtHiggs <= 0) mtHiggs = 0.0;
  else             mtHiggs = TMath::Sqrt(mtHiggs);

  delete dilepton;

  return mtHiggs;
}

Double_t JetTools::Beta(const TrackCol *tracks, Jet *jet, const Vertex *vertex, Double_t  delta_z, Double_t delta_cone){  

  if(tracks->GetEntries() <= 0) return 1.0;

  double Pt_jets_X = 0. ;
  double Pt_jets_Y = 0. ;
  double Pt_jets_X_tot = 0. ;
  double Pt_jets_Y_tot = 0. ;

  for(int i=0;i<int(tracks->GetEntries());i++){
    if(MathUtils::DeltaR(tracks->At(i)->Mom(),jet->Mom()) < delta_cone){
      Pt_jets_X_tot += tracks->At(i)->Px();
      Pt_jets_Y_tot += tracks->At(i)->Py();  
      double pDz = TMath::Abs(tracks->At(i)->DzCorrected(*vertex));
      if(pDz < delta_z){
        Pt_jets_X += tracks->At(i)->Px();
        Pt_jets_Y += tracks->At(i)->Py();   
      }
    }
  }

  if(sqrt(Pt_jets_X_tot*Pt_jets_X_tot + Pt_jets_Y_tot*Pt_jets_Y_tot) > 0)
    return sqrt(Pt_jets_X*Pt_jets_X + Pt_jets_Y*Pt_jets_Y) / sqrt(Pt_jets_X_tot*Pt_jets_X_tot + Pt_jets_Y_tot*Pt_jets_Y_tot);

  return 1.0;
}


Double_t JetTools::Beta(const PFJet *jet, const Vertex *vertex, Double_t  delta_z){  
  double Pt_jets= 0. ;
  double Pt_jetsTot = 0. ;
  
  for(UInt_t i=0;i<jet->NPFCands();i++){
    if(jet->PFCand(i)->TrackerTrk()){ 
      Pt_jetsTot += jet->PFCand(i)->TrackerTrk()->Pt();
      double pDz = TMath::Abs(jet->PFCand(i)->TrackerTrk()->DzCorrected(*vertex));
      if(pDz < delta_z){
        Pt_jets += jet->PFCand(i)->TrackerTrk()->Pt();
      }
    }
  }

  Double_t beta = 0.;
  if (Pt_jetsTot > 0)
    beta = Pt_jets/Pt_jetsTot;

  return beta;
}

Double_t JetTools::Beta2(const PFJet *jet, const Vertex *vertex, Double_t  delta_z){  
  double Pt_jets= 0. ;
  double Pt_jetsTot = 0. ;

  for(UInt_t i=0;i<jet->NPFCands();i++){
    if(jet->PFCand(i)->BestTrk()){
      Pt_jetsTot += jet->PFCand(i)->BestTrk()->Pt()*jet->PFCand(i)->BestTrk()->Pt();
      double pDz = TMath::Abs(jet->PFCand(i)->BestTrk()->DzCorrected(*vertex));
      if(pDz < delta_z){
        Pt_jets += jet->PFCand(i)->BestTrk()->Pt()*jet->PFCand(i)->BestTrk()->Pt();
      }
    }
  }

  Double_t beta = 1.0;
  if (Pt_jetsTot > 0)
    beta = Pt_jets/Pt_jetsTot;
  return beta;
}


Bool_t  JetTools::PassBetaVertexAssociationCut(const PFJet *jet, const Vertex *referenceVertex, const VertexCol *vertices, Double_t delta_z) {

  Bool_t passBetaCut = kTRUE;
  if(vertices->GetEntries() > 0) {
    Double_t Beta = JetTools::Beta(jet, referenceVertex, 0.2);
    Double_t Beta_other = 0.0;
    for(UInt_t nv=0; nv<vertices->GetEntries(); nv++){
      if (referenceVertex == vertices->At(nv)) continue;
      Double_t BetaAux = JetTools::Beta(jet, vertices->At(nv), 0.2);
      if(BetaAux > Beta_other) Beta_other = BetaAux;
    }
    if(Beta_other > Beta) passBetaCut = kFALSE;
  }

  return passBetaCut; 

}

Bool_t  JetTools::PassBeta2VertexAssociationCut(const PFJet *jet, const Vertex *referenceVertex, const VertexCol *vertices, Double_t delta_z) {

  Bool_t passBetaCut = kTRUE;
  if(vertices->GetEntries() > 0) {
    Double_t Beta = JetTools::Beta2(jet, referenceVertex, 0.2);
    Double_t Beta_other = 0.0;
    for(UInt_t nv=0; nv<vertices->GetEntries(); nv++){
      if (referenceVertex == vertices->At(nv)) continue;
      Double_t BetaAux = JetTools::Beta2(jet, vertices->At(nv), 0.2);
      if(BetaAux > Beta_other) Beta_other = BetaAux;
    }
    if(Beta_other > Beta) passBetaCut = kFALSE;
  }

  return passBetaCut; 

}


Int_t JetTools::MaxBetaVertexIndex(const PFJet *jet, const VertexCol *vertices, Double_t  delta_z=0.2){  
  
  Int_t vertexIndex = -1;
  double beta = -0.1;
  for (UInt_t v=0; v < vertices->GetEntries(); v++){
    Double_t betaTmp = JetTools::Beta(jet, vertices->At(v), delta_z);
    if (betaTmp > beta) {
      beta = betaTmp;
      vertexIndex = v;
    }
  }
  return vertexIndex;

}

Int_t JetTools::MaxBeta2VertexIndex(const PFJet *jet, const VertexCol *vertices, Double_t  delta_z=0.2){  
  
  Int_t vertexIndex = -1;
  double beta = -0.1;
  for (UInt_t v=0; v < vertices->GetEntries(); v++){
    Double_t betaTmp = JetTools::Beta2(jet, vertices->At(v), delta_z);
    if (betaTmp > beta) {
      beta = betaTmp;
      vertexIndex = v;
    }
  }
  return vertexIndex;

}


Int_t JetTools::JetToPVAssociation(const PFJet *jet, const VertexCol *vertices, Double_t  delta_z=0.2){  

  std::vector<float> verticesPt2(vertices->GetEntries());
  for(UInt_t i=0;i<jet->NPFCands();i++){
    if(jet->PFCand(i)->BestTrk()){
      double minDZ = delta_z;
      int trackVertexIndex = -1;
      for (UInt_t v=0; v < vertices->GetEntries(); v++){
	if (minDZ > TMath::Abs(jet->PFCand(i)->BestTrk()->DzCorrected(*vertices->At(v)))) {
	  minDZ = TMath::Abs(jet->PFCand(i)->BestTrk()->DzCorrected(*vertices->At(v)));
	  trackVertexIndex = v;
	}
      }
      if (trackVertexIndex < 0) continue;
      verticesPt2[trackVertexIndex]+= jet->PFCand(i)->BestTrk()->Pt()*jet->PFCand(i)->BestTrk()->Pt();
    }
  }

  Int_t vertexIndex = 0;
  float pt2Max = 0; 
  for (uint i=0; i < verticesPt2.size(); ++i){
    if (pt2Max < verticesPt2[i]) {
      pt2Max = verticesPt2[i];
      vertexIndex = i;
    }
  }
  return vertexIndex;
}
const PFCandidate* JetTools::leadCand(const PFJet *iJet,int iPFType,bool i2nd) { 
  int lCount = 0;
  const PFCandidate *lCand = 0;
  for(UInt_t i0 = 0; i0 < iJet->NPFCands(); i0++) { 
    lCand = iJet->PFCand(i0);
    if(iPFType != -1 && lCand->PFType() != iPFType) continue;
    if(lCount == 0 && !i2nd) break;
    if(lCount >  0)          break;
    lCount++;
  }
  return lCand; 
}
Double_t JetTools::impactParameter(const PFJet *iJet,const Vertex *iVertex,bool iDZ) { 
  double lDZCorr = -1000;
  for(UInt_t i0 = 0; i0 < iJet->NPFCands(); i0++) { 
    const PFCandidate *pCand = iJet->PFCand(i0);
    if(pCand->TrackerTrk() == 0) continue;
    //if(pCand->Pt() < 1.) continue; => previous iterations
    if(iDZ)  lDZCorr = fabs(pCand->TrackerTrk()->DzCorrected(*iVertex));
    if(!iDZ) lDZCorr = fabs(pCand->TrackerTrk()->D0Corrected(*iVertex));
    break;
  }
  return lDZCorr;
}
Double_t JetTools::dRMean(const PFJet *iJet,int iPFType) { 
  double lDRMean = 0;
  for(UInt_t i0 = 0; i0 < iJet->NPFCands(); i0++) { 
    const PFCandidate *pCand = iJet->PFCand(i0);
    if(iPFType != -1 && pCand->PFType() != iPFType) continue;
    double pDR = MathUtils::DeltaR(iJet->Mom(),pCand->Mom());
    lDRMean    += pDR*(pCand->Pt())/iJet->RawMom().Pt();
  }
  return lDRMean;
}
Double_t JetTools::dR2Mean(const PFJet *iJet,int iPFType) { 
  double lDR2Mean = 0;
  double lSumPt2 = 0;
  for(UInt_t i0 = 0; i0 < iJet->NPFCands(); i0++) { 
    const PFCandidate *pCand = iJet->PFCand(i0);
    if(iPFType != -1 && pCand->PFType() != iPFType) continue;
    lSumPt2   += pCand->Pt()*pCand->Pt();
    double pDR = MathUtils::DeltaR(iJet->Mom(),pCand->Mom());
    lDR2Mean    += pDR*pDR*(pCand->Pt()*pCand->Pt());
  }
  lDR2Mean/=lSumPt2;
  return lDR2Mean;
}
Double_t JetTools::frac(const PFJet *iJet,Double_t iDRMax,Double_t iDRMin,Int_t iPFType) {
  double lFrac = 0;
  for(UInt_t i0 = 0; i0 < iJet->NPFCands(); i0++) {
    const PFCandidate *pCand = iJet->PFCand(i0);
    if(iPFType != -1 && pCand->PFType() != iPFType) continue;
    Double_t pDR = MathUtils::DeltaR(iJet->Mom(),pCand->Mom());
    if(pDR > iDRMax) continue;
    if(pDR < iDRMax-0.1) continue;
    lFrac += pCand->Pt()/iJet->RawMom().Pt();
  }
  return lFrac;
}
Double_t JetTools::betaStar(const PFJet *iJet,const Vertex *iVertex,const VertexCol* iVertices,Double_t iDZCut) {
  Double_t lTotal = 0;  
  Double_t lPileup = 0;
  for(UInt_t i0 = 0; i0 < iJet->NPFCands(); i0++) {
    const PFCandidate* pPF   = iJet->PFCand(i0);
    const Track* pTrack      = pPF->TrackerTrk();
    //if(pPF->GsfTrk()) pTrack = pPF->GsfTrk(); ==> not used in CMSSW
    if(pTrack == 0) continue;
    lTotal += pTrack->Pt();
    double pDZPV  = fabs(pTrack->DzCorrected(*iVertex));
    double pDZMin = pDZPV;
    for(unsigned int i1 = 0; i1 < iVertices->GetEntries(); i1++) {
      const Vertex *pV = iVertices->At(i1);
      if(pV->Ndof() < 4 ||
	 (pV->Position() - iVertex->Position()).R() < 0.02 ) continue;
      pDZMin = TMath::Min(pDZMin,fabs(pTrack->DzCorrected(*pV)));
    }
    if(pDZPV > 0.2 && pDZMin < 0.2) lPileup += pTrack->Pt(); 
  }
  if(lTotal == 0) lTotal = 1;
  return lPileup/(lTotal);
}
Double_t JetTools::betaStarClassic(const PFJet *iJet,const Vertex *iVertex,const VertexCol* iVertices) {
  Double_t lTotal = 0;  
  Double_t lPileup = 0;
  for(UInt_t i0 = 0; i0 < iJet->NPFCands(); i0++) {
    const PFCandidate* pPF   = iJet->PFCand(i0);
    const Track* pTrack      = pPF->TrackerTrk();
    //if(pPF->GsfTrk()) pTrack = pPF->GsfTrk(); ==> not used in CMSSW
    if(pTrack == 0) continue;
    lTotal += pTrack->Pt();
    bool isPV    = iVertex->HasTrack(pPF->TrackerTrk());
    bool isOtherV = false;
    for(unsigned int i1 = 0; i1 < iVertices->GetEntries(); i1++) {
      const Vertex *pV = iVertices->At(i1);
      if(isOtherV || isPV) continue;
      if(pV->Ndof() < 4 ||
	 (pV->Position() - iVertex->Position()).R() < 0.02 ) continue;
      isOtherV    = pV->HasTrack(pPF->TrackerTrk());
    }
    if(!isPV && isOtherV) lPileup += pTrack->Pt(); 
  }
  if(lTotal == 0) lTotal = 1;
  return lPileup/(lTotal);
}
Bool_t  JetTools::passPFLooseId(const PFJet *iJet) { 
  if(iJet->RawMom().E()                              == 0)       return false;
  if(iJet->NeutralHadronEnergy()/iJet->RawMom().E()  >  0.99)    return false;
  if(iJet->NeutralEmEnergy()/iJet->RawMom().E()      >  0.99)    return false;
  if(iJet->NConstituents()                           <  2)	 return false;
  if(iJet->ChargedHadronEnergy()/iJet->RawMom().E()  <= 0     && fabs(iJet->Eta()) < 2.4 ) return false;
  if(iJet->ChargedEmEnergy()/iJet->RawMom().E()      >  0.99  && fabs(iJet->Eta()) < 2.4 ) return false;
  if(iJet->ChargedMultiplicity()                     < 1      && fabs(iJet->Eta()) < 2.4 ) return false;
  if(fabs(iJet->Eta())                               > 4.99) return false;
  return true;
}

/*
double JetTools::genFrac(const PFJet *iJet) { 
  double lTrueFrac = 0;
  for(UInt_t i0 = 0; i0 < fParticles->GetEntries(); i0++) { 
    const MCParticle *p = fParticles->At(i0);
    if(p->Status() != 1) continue;
    double pDEta = iJet->Eta() - p->Eta(); 
    double pDPhi = fabs(iJet->Phi()-p->Phi()); if(pDPhi > 2.*TMath::Pi() - pDPhi) pDPhi =  2.*TMath::Pi() - pDPhi;
    double pDR   = sqrt(pDEta*pDEta + pDPhi*pDPhi);
    if(pDR > 0.5) continue;
    lTrueFrac += p->Pt();
  }
  lTrueFrac/=iJet->Pt();
  return lTrueFrac;
}
*/
