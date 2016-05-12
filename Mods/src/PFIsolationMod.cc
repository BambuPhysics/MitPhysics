#include "MitPhysics/Mods/interface/PFIsolationMod.h"

#include "MitAna/DataTree/interface/ParticleCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"

#include "MitPhysics/Utils/interface/IsolationTools.h"

ClassImp(mithep::PFIsolationMod)

mithep::PFIsolationMod::PFIsolationMod(char const* name, char const* title) :
  BaseMod(name, title)
{
  for (unsigned iT = 0; iT != nIsolationTypes; ++iT) {
    fIsoSum[iT].Resize(32);

    fMinDR[iT] = 1.e-5;
    fMaxDR[iT] = 0.3;
    fMaxDRho[iT] = 1.;
    fVetoDRBarrel[iT] = 0.;
  }

  fVetoDREndcap[kChargedHadron] = 0.015;
  fVetoDREndcap[kNeutralHadron] = 0.;
  fVetoDREndcap[kPhoton] = 0.08;
}

void
mithep::PFIsolationMod::SlaveBegin()
{
  fIsoSum[kChargedHadron].SetName(fOutputName + "CH");
  fIsoSum[kNeutralHadron].SetName(fOutputName + "NH");
  fIsoSum[kPhoton].SetName(fOutputName + "Ph");

  for (unsigned iT = 0; iT != nIsolationTypes; ++iT) {
    if (!PublishObj(fIsoSum + iT))
      SendError(kAbortAnalysis, "SlaveBegin", "Cannot publish iso sum");
  }

  if (GetFillHist()) {
    for (unsigned iT = 0; iT != nIsolationTypes; ++iT)
      AddTH1(fHIso[iT], fIsoSum[iT].GetName(), "isolation", 1000, 0., 50.);
  }
}

void
mithep::PFIsolationMod::Process()
{
  auto& candidates = *GetObject<ParticleCol>(fCandsName);

  auto& pfs = *GetObject<PFCandidateCol>(fPFCandidatesName);
  auto& pfnopus = *GetObject<PFCandidateCol>(fPFNoPileupCandidatesName);

  for (unsigned iT = 0; iT != nIsolationTypes; ++iT)
    fIsoSum[iT].Resize(candidates.GetEntries());

  for (unsigned iC = 0; iC != candidates.GetEntries(); ++iC) {
    auto& candidate = *candidates.At(iC);

    std::map<Double_t, Double_t> deposits[nIsolationTypes] = {
      IsolationTools::PFIsoDeposit(candidate, pfnopus, fMinDR[kChargedHadron], fMaxDR[kChargedHadron], fMaxDRho[kChargedHadron]),
      IsolationTools::PFIsoDeposit(candidate, pfs, fMinDR[kNeutralHadron], fMaxDR[kNeutralHadron], fMaxDRho[kNeutralHadron], PFCandidate::eNeutralHadron),
      IsolationTools::PFIsoDeposit(candidate, pfs, fMinDR[kPhoton], fMaxDR[kPhoton], fMaxDRho[kPhoton], PFCandidate::eGamma)
    };

    for (unsigned iT = 0; iT != nIsolationTypes; ++iT) {
      fIsoSum[iT].At(iC) = IsolationTools::IsoFromDeposit(candidate, deposits[iT], fMaxDR[iT], fVetoDRBarrel[iT], fVetoDREndcap[iT]);
      if (GetFillHist())
        fHIso[iT]->Fill(fIsoSum[iT].At(iC));
    }
  }
}

void
mithep::PFIsolationMod::SlaveTerminate()
{
  for (unsigned iT = 0; iT != nIsolationTypes; ++iT)
    RetractObj(fIsoSum[iT].GetName());
}
