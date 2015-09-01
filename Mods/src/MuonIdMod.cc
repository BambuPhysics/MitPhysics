#include "MitPhysics/Mods/interface/MuonIdMod.h"
#include "MitPhysics/Utils/interface/MuonTools.h"
#include "MitAna/DataTree/interface/Names.h"

ClassImp(mithep::MuonIdMod)

// -------------------------------------------------------------------------
mithep::MuonIdMod::MuonIdMod(char const* name/* = "MuonIdMod"*/, char const* title/* = "Muon Identification"*/) :
  IdMod<mithep::Muon>(name, title)
{
  fInputName = mithep::Names::gkMuonBrn;
}

// -------------------------------------------------------------------------
mithep::MuonIdMod::~MuonIdMod()
{
}

// -------------------------------------------------------------------------
Bool_t
mithep::MuonIdMod::IsGood(mithep::Muon const& muon)
{
  fCutFlow->Fill(cAll);

  // Using bool PassClass
  if (!MuonTools::PassClass(&muon, MuonTools::EMuClassType(fMuonClassType), GetVertices()))
    return false;

  fCutFlow->Fill(cMuonClass);

  double pt = 0.;
  double absEta = 0.;
  MuonTools::MuonPtEta(&muon, MuonTools::EMuClassType(fMuonClassType), pt, absEta);

  if (pt < fPtMin)
    return false;

  fCutFlow->Fill(cPt);

  if (absEta > fEtaMax)
    return false;

  fCutFlow->Fill(cEta);

  // Apply d0 cut
  if (fApplyD0Cut) {
    if (fWhichVertex >= -1) {
      if (!MuonTools::PassD0Cut(&muon, GetVertices(), MuonTools::EMuIdType(fIdType), fWhichVertex))
        return false;
    }
    else {
      if (!MuonTools::PassD0Cut(&muon, GetBeamSpot(), MuonTools::EMuIdType(fIdType)))
        return false;
    }
  }

  fCutFlow->Fill(cD0);

  // Apply dz cut
  if (fApplyDZCut) {
    if (!MuonTools::PassDZCut(&muon, GetVertices(), MuonTools::EMuIdType(fIdType), fWhichVertex))
      return false;
  }

  fCutFlow->Fill(cDZ);

  // calling MuonTools::PassId directly (no wrapping / forking necessary at the moment)
  if (!MuonTools::PassId(&muon, MuonTools::EMuIdType(fIdType)))
    return false;

  fCutFlow->Fill(cId);
    
  if (!PassIsolation(muon))
    return false;

  fCutFlow->Fill(cIsolation);

  return true;
}

// -------------------------------------------------------------------------
void
mithep::MuonIdMod::IdBegin()
{
  if (fIsoType == MuonTools::kCustomIso)
    SendError(kAbortAnalysis, "IdBegin", "Custom muon isolation is not yet implemented.");

  fCutFlow->SetBins(nCuts, 0., double(nCuts));
  TAxis* xaxis = fCutFlow->GetXaxis();
  xaxis->SetBinLabel(cAll + 1, "All");
  xaxis->SetBinLabel(cMuonClass + 1, "MuonClass");
  xaxis->SetBinLabel(cPt + 1, "Pt");
  xaxis->SetBinLabel(cEta + 1, "Eta");
  xaxis->SetBinLabel(cD0 + 1, "D0");
  xaxis->SetBinLabel(cDZ + 1, "DZ");
  xaxis->SetBinLabel(cId + 1, "Id");
  xaxis->SetBinLabel(cIsolation + 1, "Isolation");
}

Bool_t
mithep::MuonIdMod::PassIsolation(mithep::Muon const& muon)
{
  switch (fIsoType) {
  case MuonTools::kTrackCaloSliding:
    return MuonTools::PassIsoRhoCorr(&muon, MuonTools::kTrackCaloSliding, GetPileupEnergyDensity()->At(0)->Rho(fRhoAlgo));

  case MuonTools::kPFIso:
  case MuonTools::kMVAIso_BDTG_IDIso:
    return MuonTools::PassPFIso(&muon, MuonTools::EMuIsoType(fIsoType), GetPFCandidates(), GetVertices()->At(0));

  case MuonTools::kPFRadialIso:
    return MuonTools::PassPFIso(&muon, MuonTools::kPFRadialIso, GetPFNoPileupCandidates());

  case MuonTools::kPFIsoBetaPUCorrected:
  case MuonTools::kPFIsoBetaPUCorrectedTight:
    return MuonTools::PassPFIso(&muon, MuonTools::EMuIsoType(fIsoType), GetPFNoPileupCandidates(), 0, GetPFPileupCandidates());

  case MuonTools::kPFIsoNoL:
    return MuonTools::PassPFIso(&muon, MuonTools::kPFIsoNoL, GetPFCandidates(), GetVertices()->At(0), 0, GetNonIsolatedElectrons(), GetNonIsolatedMuons());

  case MuonTools::kNoIso:
    return kTRUE;

  default:
    return MuonTools::PassIso(&muon, MuonTools::EMuIsoType(fIsoType));
  }
}
