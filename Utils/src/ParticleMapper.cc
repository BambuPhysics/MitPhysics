#include "MitPhysics/Utils/interface/ParticleMapper.h"

ClassImp(mithep::ParticleMapper)

using namespace mithep;

//--------------------------------------------------------------------------------------------------
ParticleMapper::ParticleMapper() :
  fNumParticles(0),
  fNumEtaBins(0),
  fNumPhiBins(0),
  fNumTotBins(0),
  fParticleLocation(0),
  fBinContents(0) { }

//--------------------------------------------------------------------------------------------------
ParticleMapper::~ParticleMapper(){
  delete fParticleLocation;
  delete fBinContents;
}

//--------------------------------------------------------------------------------------------------
void
ParticleMapper::Initialize( const PFCandidateCol &Particles, Double_t DeltaEta, Double_t DeltaPhi, Double_t EtaMax){

  fNumParticles = Particles.GetEntries();
  fNumEtaBins   = 2*ceil(EtaMax/DeltaEta);
  fNumPhiBins   = ceil(2*(TMath::Pi()/DeltaPhi));
  fNumTotBins   = fNumEtaBins * fNumPhiBins;

  Int_t AllocateParticleLocation[fNumParticles];
  fParticleLocation = AllocateParticleLocation;

  std::vector<Int_t> AllocateBinContents[fNumTotBins];
  fBinContents = AllocateBinContents;

  for(Int_t i0 = 0; i0 < fNumParticles; i0++){
    Int_t etaBin;
    Int_t phiBin;
    Double_t eta = Particles.At(i0)->Eta();
    if(abs(eta) > EtaMax) continue;
    Double_t phi = Particles.At(i0)->Phi();
    if(phi > TMath::Pi()) phi = phi - 2.0*(TMath::Pi());
    etaBin = floor(eta/DeltaEta) + fNumEtaBins/2;
    phiBin = floor(phi/DeltaPhi) + fNumPhiBins/2;

    Int_t finalBin = etaBin + phiBin*fNumEtaBins;
    fParticleLocation[i0] = finalBin;
    fBinContents[finalBin].push_back(i0);
  }

}

//--------------------------------------------------------------------------------------------------
std::vector<Int_t>
ParticleMapper::GetInBin(Int_t index){

  return fBinContents[fParticleLocation[index]];
  
}

//--------------------------------------------------------------------------------------------------
std::vector<Int_t>
ParticleMapper::GetSurrounding(Int_t index){

  Int_t bin = fParticleLocation[index];
  Int_t etaBin = bin % fNumEtaBins;
  Int_t phiBin = bin/fNumEtaBins;

  std::vector<Int_t> tempVec;

  Int_t tempEtaBin;
  Int_t tempPhiBin;
  Int_t tempBin;

  for(Int_t i0 = -1; i0 < 2; i0++){
    tempEtaBin = etaBin + i0;
    if(tempEtaBin < 0 || tempEtaBin >= fNumEtaBins) continue;
    for(Int_t i1 = -1; i1 < 2; i1++){
      tempPhiBin = phiBin + i0;
      if(tempPhiBin < 0 || tempPhiBin >= fNumPhiBins) continue;
      tempBin = tempEtaBin + tempPhiBin*fNumEtaBins;
      tempVec.insert(tempVec.end(),fBinContents[tempBin].begin(),fBinContents[tempBin].end());
    }
  }

  return tempVec;

}
