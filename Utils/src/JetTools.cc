#include "MitPhysics/Utils/interface/JetTools.h"

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
		- 2 * TMath::Cos(MathUtils::DeltaPhi(jets->At(j)->Phi(),particles->At(i)->Phi()))));
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
		- 2 * TMath::Cos(MathUtils::DeltaPhi(jets->At(j)->Phi(),pfCandidates->At(i)->Phi()))));
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
		- 2 * TMath::Cos(MathUtils::DeltaPhi(jets->At(j)->Phi(),tracks->At(i)->Phi()))));
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
		- 2 * TMath::Cos(MathUtils::DeltaPhi(jets->At(j)->Phi(),jetsS->At(i)->Phi()))));
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
		- 2 * TMath::Cos(MathUtils::DeltaPhi(jets->At(j)->Phi(),calos->At(i)->Phi()))));
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

  if(nsel >= 0 && nsel <= 4){
    mtHiggs = mll*mll + mnu*mnu + 2.0*(enell*enenn - enex*enex - eney*eney);
  }
  if(mtHiggs <= 0) mtHiggs = 0.0;
  else             mtHiggs = TMath::Sqrt(mtHiggs);

  delete dilepton;

  return mtHiggs;
}

void JetTools::Alpha(Double_t AlphaVar[2], const TrackCol *tracks, Jet *jet, const VertexCol *vertices, Double_t  delta_z, Double_t delta_cone){  
  AlphaVar[0] = -1.0;
  AlphaVar[1] = -1.0;
  if(tracks->GetEntries() <= 0) return;

  double Pt_jets_X = 0. ;
  double Pt_jets_Y = 0. ;
  double Pt_jets_X_tot = 0. ;
  double Pt_jets_Y_tot = 0. ;

  for(int i=0;i<int(tracks->GetEntries());i++){
    if(MathUtils::DeltaR(tracks->At(i)->Mom(),jet->Mom()) < delta_cone){
      Pt_jets_X_tot += tracks->At(i)->Px();
      Pt_jets_Y_tot += tracks->At(i)->Py();  
      double pDz = TMath::Abs(tracks->At(i)->DzCorrected(*vertices->At(0)));
      if(pDz < delta_z){
        Pt_jets_X += tracks->At(i)->Px();
        Pt_jets_Y += tracks->At(i)->Py();   
      }
    }
  }

  if(jet->Pt() > 0)
    AlphaVar[0] = sqrt(Pt_jets_X*Pt_jets_X + Pt_jets_Y*Pt_jets_Y) / jet->Pt();
  if(Pt_jets_X_tot > 0 || Pt_jets_Y_tot > 0)
    AlphaVar[1] = sqrt(Pt_jets_X*Pt_jets_X + Pt_jets_Y*Pt_jets_Y) / sqrt(Pt_jets_X_tot*Pt_jets_X_tot + Pt_jets_Y_tot*Pt_jets_Y_tot);
}

