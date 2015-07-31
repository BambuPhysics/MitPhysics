#include "MitPhysics/Utils/interface/ParticleMapper.h"

#include "TVector2.h"

ClassImp(mithep::ParticleMapper)

using namespace mithep;

//--------------------------------------------------------------------------------------------------
ParticleMapper::ParticleMapper(Double_t DeltaEta/* = 0.3*/, Double_t DeltaPhi/* = 0.3*/, Double_t EtaMax/* = 5.*/) :
  fDeltaEta(DeltaEta),
  fDeltaPhi(DeltaPhi),
  fEtaMax(EtaMax),
  fNumEtaBins(2*ceil(EtaMax/fDeltaEta)),
  fNumPhiBins(floor(2*(TMath::Pi()/fDeltaPhi))),  // The last bin will be larger than fDeltaPhi
  fParticleLocation(),
  fBinContents(fNumEtaBins * fNumPhiBins)
{
}

//--------------------------------------------------------------------------------------------------
ParticleMapper::~ParticleMapper()
{
}

//--------------------------------------------------------------------------------------------------
void
ParticleMapper::InitEvent(const PFCandidateCol &Particles)
{
  fParticleLocation.assign(Particles.GetEntries(), -1);

  for (auto& cont : fBinContents)
    cont.clear();

  for (UInt_t i0 = 0; i0 < fParticleLocation.size(); i0++) {
    Int_t etaBin;
    Int_t phiBin;
    Double_t eta = Particles.At(i0)->Eta();
    if (fabs(eta) > fEtaMax)
      continue;

    // Usng 0-2pi is easier so that there's only one bin with weird resolution
    Double_t phi = TVector2::Phi_0_2pi(Particles.At(i0)->Phi());

    etaBin = floor(eta/fDeltaEta) + fNumEtaBins/2;
    phiBin = floor(phi/fDeltaPhi);
    if (phiBin == Int_t(fNumPhiBins))
      phiBin = phiBin - 1;                          // Sticks overflow into last bin

    Int_t finalBin = etaBin + phiBin*fNumEtaBins;
    fParticleLocation[i0] = finalBin;
    fBinContents[finalBin].push_back(i0);
  }
}

//--------------------------------------------------------------------------------------------------
std::vector<UInt_t>
ParticleMapper::GetSurrounding(UInt_t index) const
{
  if (index >= fParticleLocation.size() || fParticleLocation[index] < 0)
    return {};

  Int_t bin = fParticleLocation[index];
  Int_t etaBin = bin % fNumEtaBins;
  Int_t phiBin = bin / fNumEtaBins;

  return ReturnNear(etaBin, phiBin);
}

//--------------------------------------------------------------------------------------------------
std::vector<UInt_t>
ParticleMapper::GetNearEtaPhi(Double_t eta, Double_t phi) const
{
  if (fabs(eta) > fDeltaEta*(fNumEtaBins/2))
    return {};

  phi = TVector2::Phi_0_2pi(phi);

  Int_t etaBin = floor(eta/fDeltaEta) + fNumEtaBins/2;
  Int_t phiBin = floor(phi/fDeltaPhi);
  if (phiBin == Int_t(fNumPhiBins))
    phiBin = phiBin - 1;

  return ReturnNear(etaBin, phiBin);
}

//--------------------------------------------------------------------------------------------------
std::vector<UInt_t>
ParticleMapper::ReturnNear(Int_t etaBin, Int_t phiBin) const
{
  std::vector<UInt_t> tempVec;

  Int_t tempEtaBin;
  Int_t tempPhiBin;
  Int_t tempBin;

  for (Int_t i0 = -1; i0 < 2; i0++) {
    tempEtaBin = etaBin + i0;
    if (tempEtaBin < 0 || tempEtaBin >= Int_t(fNumEtaBins))
      continue;

    for (Int_t i1 = -1; i1 < 2; i1++) {
      tempPhiBin = phiBin + i1;
      if (tempPhiBin == -1)
        tempPhiBin = fNumPhiBins - 1;               // This allows wrapping
      if (tempPhiBin == Int_t(fNumPhiBins))
        tempPhiBin = 0;
      tempBin = tempEtaBin + tempPhiBin*fNumEtaBins;
      tempVec.insert(tempVec.end(),fBinContents[tempBin].begin(),fBinContents[tempBin].end());
    }
  }
  return tempVec;
}
