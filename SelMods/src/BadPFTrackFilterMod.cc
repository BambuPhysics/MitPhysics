#include "MitPhysics/SelMods/interface/BadPFTrackFilterMod.h"

#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/Track.h"
#include "MitCommon/MathTools/interface/MathUtils.h"

ClassImp(mithep::BadPFTrackFilterMod)

char const*
mithep::BadPFTrackFilterMod::GetOutputName(Int_t i/* = -1*/) const
{
if (i == 0)
    return "BadPFMuonFilter";
  else if (i == 1)
    return "BadChargedCandidateFilter";
  else
    return fOutput.GetName();
}

void
mithep::BadPFTrackFilterMod::SlaveBegin()
{
  if (fTaggingMode) {
    fOutput.Resize(2);
    fOutput.Resize(2); // bug in MitAna/DataCont up to bambu 044
    PublishObj(&fOutput);
  }
}

void
mithep::BadPFTrackFilterMod::SlaveTerminate()
{
  if (fTaggingMode)
    RetractObj(fOutput.GetName());
}

void
mithep::BadPFTrackFilterMod::Process()
{
  if (fTaggingMode) {
    // will set flag to true if bad event
    fOutput.At(0) = false;
    fOutput.At(1) = false;
  }

  auto* muons = GetObject<MuonCol>(fMuonsName);
  auto* pfs = GetObject<PFCandidateCol>(fPFCandidatesName);

  PFCandidateOArr pfMuons;
  for (unsigned iPF = 0; iPF != pfs->GetEntries(); ++iPF) {
    auto& cand(*pfs->At(iPF));
    if (cand.PFType() == PFCandidate::eMuon && cand.Pt() > fMinPt)
      pfMuons.Add(&cand);
  }

  PFCandidateOArr pfCHs;
  for (unsigned iPF = 0; iPF != pfs->GetEntries(); ++iPF) {
    auto& cand(*pfs->At(iPF));
    if (cand.PFType() == PFCandidate::eHadron)
      pfCHs.Add(&cand);
  }

  // BadPFMuonFilter

  unsigned iMu = 0;
  for (; iMu != muons->GetEntries(); ++iMu) {
    auto& muon(*muons->At(iMu));

    auto* innerTrack = muon.TrackerTrk();

    if (!innerTrack)
      continue;

    if (innerTrack->Pt() < fMinPt)
      continue;

    if (innerTrack->Quality().Quality(TrackQuality::highPurity))
      continue;

    if (innerTrack->PtErr() / innerTrack->Pt() < fMinTrackRelErr)
      continue;

    // CMSSW implementation also checks "originalAlgo" which we don't have access to
    // 14 = Track::muonSeededStepOutIn
    //    if (innerTrack->Algo() != Track::muonSeededStepOutIn)
    if (innerTrack->Algo() != 14)
      continue;

    unsigned iPF = 0;
    for (; iPF != pfMuons.GetEntries(); ++iPF) {
      auto& pfMuon(*pfMuons.At(iPF));
      if (MathUtils::DeltaR(pfMuon, muon) < 0.001)
        break;
    }

    if (iPF != pfMuons.GetEntries())
      break;
  }

  if (iMu != muons->GetEntries()) {
    if (fTaggingMode)
      fOutput.At(0) = true;
    else {
      SkipEvent();
      return;
    }
  }

  // BadChargedCandidateFilter

  iMu = 0;
  for (; iMu != muons->GetEntries(); ++iMu) {
    auto& muon(*muons->At(iMu));

    if (muon.Pt() < fMinPt)
      continue;

    auto* innerTrack = muon.TrackerTrk();

    if (!innerTrack)
      continue;

    if (innerTrack->Quality().Quality(TrackQuality::highPurity))
      continue;

    if (innerTrack->PtErr() / innerTrack->Pt() < fMinTrackRelErr)
      continue;


    unsigned iPF = 0;
    for (; iPF != pfCHs.GetEntries(); ++iPF) {
      auto& ch(*pfCHs.At(iPF));

      if (MathUtils::DeltaR(ch, muon) > 0.001)
        continue;

      double pfPt = ch.Pt();
      double inPt = innerTrack->Pt();

      if ((pfPt - inPt) / ((pfPt + inPt) * 0.5) > fMinRelDPt)
        break;
    }

    if (iPF != pfCHs.GetEntries())
      break;
  }

  if (iMu != muons->GetEntries()) {
    if (fTaggingMode)
      fOutput.At(1) = true;
    else {
      SkipEvent();
      return;
    }
  }
}
