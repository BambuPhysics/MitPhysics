#include "MitPhysics/Utils/interface/ParticleMapper.h"

ClassImp(mithep::ParticleMapper)

using namespace mithep;

//--------------------------------------------------------------------------------------------------
ParticleMapper::ParticleMapper() :
  fDeltaEta(0),
  fDeltaPhi(0),
  fNumParticles(0),
  fNumEtaBins(0),
  fNumPhiBins(0),
  fNumTotBins(0),
  fParticleLocation(0),
  fBinContents(0)
{ }

//--------------------------------------------------------------------------------------------------
ParticleMapper::~ParticleMapper()
{
  delete fParticleLocation;
  delete fBinContents;
}

//--------------------------------------------------------------------------------------------------
void
ParticleMapper::Initialize(const PFCandidateCol &Particles, Double_t DeltaEta, Double_t DeltaPhi, Double_t EtaMax)
{
  fDeltaEta = DeltaEta;
  fDeltaPhi = DeltaPhi;

  fNumParticles = Particles.GetEntries();
  fNumEtaBins   = 2*ceil(EtaMax/fDeltaEta);
  fNumPhiBins   = floor(2*(TMath::Pi()/fDeltaPhi));  // The last bin will be larger than fDeltaPhi
  fNumTotBins   = fNumEtaBins * fNumPhiBins;

  fParticleLocation = new Int_t[fNumParticles];

  fBinContents = new std::vector<Int_t>[fNumTotBins];

  for (Int_t i0 = 0; i0 < fNumTotBins; i0++)
    fBinContents[i0].resize(0);

  for (Int_t i0 = 0; i0 < fNumParticles; i0++) {
    Int_t etaBin;
    Int_t phiBin;
    Double_t eta = Particles.At(i0)->Eta();
    if (fabs(eta) > EtaMax) {
      fParticleLocation[i0] = -1;
      continue;
    }
    Double_t phi = Particles.At(i0)->Phi();
    if (phi < 0)
      phi = phi + 2.0*(TMath::Pi());                // This way's easier so that there's only one bin with weird resolution
    etaBin = floor(eta/fDeltaEta) + fNumEtaBins/2;
    phiBin = floor(phi/fDeltaPhi);
    if (phiBin == fNumPhiBins)
      phiBin = phiBin - 1;                          // Sticks overflow into last bin
    Int_t finalBin = etaBin + phiBin*fNumEtaBins;
    fParticleLocation[i0] = finalBin;
    fBinContents[finalBin].push_back(i0);
  }
}

//--------------------------------------------------------------------------------------------------
std::vector<Int_t>
ParticleMapper::GetSurrounding(Int_t index)
{
  if (fParticleLocation[index] < 0) {
    std::vector<Int_t> blankVector;
    blankVector.resize(0);
    return blankVector;
  }

  Int_t bin = fParticleLocation[index];
  Int_t etaBin = bin % fNumEtaBins;
  Int_t phiBin = bin/fNumEtaBins;

  return ReturnNear(etaBin, phiBin);
}

//--------------------------------------------------------------------------------------------------
std::vector<Int_t>
ParticleMapper::GetNearEtaPhi(Double_t eta, Double_t phi)
{
  if (fabs(eta) > fDeltaEta*(fNumEtaBins/2)) {
    std::vector<Int_t> blankVector;
    blankVector.resize(0);
    return blankVector;
  }

  while (phi > 2*(TMath::Pi()))
    phi = phi - 2*(TMath::Pi());
  while (phi < 0) 
    phi = phi + 2*(TMath::Pi());

  Int_t etaBin = floor(eta/fDeltaEta) + fNumEtaBins/2;
  Int_t phiBin = floor(phi/fDeltaPhi);
  if (phiBin == fNumPhiBins) 
    phiBin = phiBin - 1;

  return ReturnNear(etaBin, phiBin);
}

//--------------------------------------------------------------------------------------------------
std::vector<Int_t>
ParticleMapper::ReturnNear(Int_t etaBin, Int_t phiBin)
{
  std::vector<Int_t> tempVec;
  tempVec.resize(0);

  Int_t tempEtaBin;
  Int_t tempPhiBin;
  Int_t tempBin;

  for (Int_t i0 = -1; i0 < 2; i0++) {
    tempEtaBin = etaBin + i0;
    if (tempEtaBin < 0 || tempEtaBin >= fNumEtaBins) 
      continue;
    for (Int_t i1 = -1; i1 < 2; i1++) {
      tempPhiBin = phiBin + i1;
      if (tempPhiBin == -1) 
        tempPhiBin = fNumPhiBins - 1;               // This allows wrapping
      if (tempPhiBin == fNumPhiBins)
        tempPhiBin = 0;
      tempBin = tempEtaBin + tempPhiBin*fNumEtaBins;
      tempVec.insert(tempVec.end(),fBinContents[tempBin].begin(),fBinContents[tempBin].end());
    }
  }
  return tempVec;
}
