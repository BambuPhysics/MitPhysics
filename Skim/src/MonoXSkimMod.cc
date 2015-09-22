#include "MitPhysics/Skim/interface/MonoXSkimMod.h"

#include "MitAna/DataTree/interface/MetCol.h"
#include "MitAna/DataTree/interface/ParticleCol.h"
#include "MitAna/DataTree/interface/TriggerTable.h"
#include "MitAna/DataTree/interface/TriggerMask.h"
#include "MitAna/DataTree/interface/Names.h"
#include "MitAna/TreeMod/interface/HLTFwkMod.h"
#include "MitCommon/MathTools/interface/MathUtils.h"

#include <limits>

using namespace mithep;
ClassImp(mithep::MonoXSkimMod)

TString MonoXSkimMod::fgEventCategories[nEventCategories] = {
  "Met",
  "Dielectron",
  "Dimuon",
  "SingleElectron",
  "SingleMuon",
  "Photon"
};

//--------------------------------------------------------------------------------------------------
MonoXSkimMod::MonoXSkimMod(const char *name, const char *title) :
  BaseMod(name, title),
  fCategoryFlags(nEventCategories, "EventCategories")
{
  fCategoryFlags.Resize(nEventCategories);

  std::fill_n(fCategoryActive, nEventCategories, true);
  std::fill_n(fMinNumJets, nEventCategories, 0);
  std::fill_n(fMaxNumJets, nEventCategories, 0xffffffff);
}

//--------------------------------------------------------------------------------------------------
void
MonoXSkimMod::SlaveBegin()
{
  // Selection Histograms
  for(unsigned iCat = 0; iCat != nEventCategories; ++iCat){
    if (!fCategoryActive[iCat])
      continue;

    // Histograms after preselection
    AddTH1(fCutflow[iCat], "hCutflow" + fgEventCategories[iCat], ";;Events", nCuts, 0., nCuts);
    auto* xaxis = fCutflow[iCat]->GetXaxis();
    xaxis->SetBinLabel(cMonoX + 1, "MonoX");
    xaxis->SetBinLabel(cTrigger + 1, "Trigger");
    xaxis->SetBinLabel(cNElectrons + 1, "NElectrons");
    xaxis->SetBinLabel(cNMuons + 1, "NMuons");
    xaxis->SetBinLabel(cNPhotons + 1, "NPhotons");
    xaxis->SetBinLabel(cNTaus + 1, "NTaus");
    xaxis->SetBinLabel(cNJets + 1, "NJets");
    xaxis->SetBinLabel(cMet + 1, "Met");
  }

  PublishObj(&fCategoryFlags);

  fNEventsSelected = 0;
}

void
MonoXSkimMod::BeginRun()
{
  if (!HasHLTInfo() || !GetHltFwkMod()->HasData())
    return;

  auto* hltTable = GetHLTTable();
  for (unsigned iCat = 0; iCat != nEventCategories; ++iCat) {
    for (auto trigger : fTriggerNames[iCat]) {
      if (trigger.EndsWith("*"))
        trigger.ReplaceAll("*", "");
      auto* triggerName = hltTable->GetWildcard(trigger);
      if (triggerName)
        fTriggerIds[iCat].push_back(triggerName->Id());
    }
  }
}

//--------------------------------------------------------------------------------------------------
void
MonoXSkimMod::Process()
{
  auto getObject([this](ParticleCol*& objects, TString const& colName) {
      if (colName.Length() != 0) {
        objects = this->GetObject<ParticleCol>(colName);
        if (!objects)
          this->SendError(kAbortAnalysis, "Process", "Could not find " + colName);
      }
    });

  // Retrieve objects

  auto* mets = GetObject<MetCol>(fMetName);
  if (!mets)
    SendError(kAbortAnalysis, "Process", "Could not find " + fMetName);
  auto& met(*mets->At(0));

  ParticleCol* jets = 0;
  getObject(jets, fJetsName);

  ParticleCol* vetoElectrons = 0;
  getObject(vetoElectrons, fVetoElectronsName);

  if (fVetoElectrons && !vetoElectrons)
    SendError(kAbortAnalysis, "Process", "Electron veto used but no electron collection name set");

  ParticleCol* goodElectrons = 0;
  getObject(goodElectrons, fGoodElectronsName);

  if ((fCategoryActive[kDielectron] || fCategoryActive[kSingleElectron]) && !goodElectrons)
    SendError(kAbortAnalysis, "Process", "Electron categories used but no electron collection name set");

  ParticleCol* vetoMuons = 0;
  getObject(vetoMuons, fVetoMuonsName);

  if (fVetoMuons && !vetoMuons)
    SendError(kAbortAnalysis, "Process", "Muon veto used but no muon collection name set");

  ParticleCol* goodMuons = 0;
  getObject(goodMuons, fGoodMuonsName);

  if ((fCategoryActive[kDimuon] || fCategoryActive[kSingleMuon]) && !goodMuons)
    SendError(kAbortAnalysis, "Process", "Muon categories used but no muon collection name set");

  ParticleCol* taus = 0;
  getObject(taus, fVetoTausName);

  if (fVetoTaus && !taus)
    SendError(kAbortAnalysis, "Process", "Tau veto used but no tau collection name set");

  ParticleCol* vetoPhotons = 0;
  getObject(vetoPhotons, fVetoPhotonsName);

  if (fVetoPhotons && !vetoPhotons)
    SendError(kAbortAnalysis, "Process", "Photon veto used but no photon collection name set");

  ParticleCol* goodPhotons = 0;
  getObject(goodPhotons, fGoodPhotonsName);

  if (fCategoryActive[kPhoton] && !goodPhotons)
    SendError(kAbortAnalysis, "Process", "Photon category used but no photon collection name set");

  TriggerMask* triggerMask = 0;
  if (HasHLTInfo() && GetHltFwkMod()->HasData()) {
    triggerMask = GetObject<TriggerMask>(Names::gkHltBitBrn);
    if (!triggerMask)
      SendError(kAbortAnalysis, "Process", "Could not find %s", Names::gkHltBitBrn);
  }

  // start event selection

  for (unsigned iCat = 0; iCat != nEventCategories; ++iCat)
    fCategoryFlags.At(iCat) = false;

  // cut on main object (overlap removal is not considered in the skim)
  ParticleCol* monoX = 0;  
  switch (fAnalysisType) {
  case kMonoJet:
    monoX = jets;
    break;

  case kMonoPhoton:
    monoX = goodPhotons;
    break;

  default:
    SendError(kAbortAnalysis, "Process", "No analysis type set!");
  }

  if (!monoX)
    SendError(kAbortAnalysis, "Process", "Main object collection is NULL");

  unsigned iMono = 0;
  for (; iMono != monoX->GetEntries(); ++iMono) {
    if (monoX->At(iMono)->Pt() > fMinMonoXPt)
      break;
  }

  if (iMono == monoX->GetEntries()) {
    SkipEvent();
    return;
  }

  // Met calculators for control regions

  auto calcMet1([&met](ParticleCol const* goodObjs)->double {
      double maxMet2 = 0.;
      for (unsigned iO = 0; iO != goodObjs->GetEntries(); ++iO) {
        double mex = met.Px() + goodObjs->At(iO)->Px();
        double mey = met.Py() + goodObjs->At(iO)->Py();
        double met2 = mex * mex + mey * mey;
        if (met2 > maxMet2)
          maxMet2 = met2;
      }
      return maxMet2;
    });

  auto calcMet2([&met](ParticleCol const* goodObjs)->double {
      double maxMet2 = 0.;
      for (unsigned iO0 = 0; iO0 != goodObjs->GetEntries(); ++iO0) {
        for (unsigned iO1 = iO0 + 1; iO1 != goodObjs->GetEntries(); ++iO1) {
          double mex = met.Px() + goodObjs->At(iO0)->Px() + goodObjs->At(iO1)->Px();
          double mey = met.Py() + goodObjs->At(iO0)->Py() + goodObjs->At(iO1)->Py();
          double met2 = mex * mex + mey * mey;
          if (met2 > maxMet2)
            maxMet2 = met2;
        }
      }
      return maxMet2;
    });

  bool eventSelected = false;

  for (unsigned iCat = 0; iCat != nEventCategories; ++iCat) {
    if (!fCategoryActive[iCat])
      continue;

    fCutflow[iCat]->Fill(cMonoX);

    if (triggerMask && fTriggerIds[iCat].size() != 0) {
      unsigned iT = 0;
      for (; iT != fTriggerIds[iCat].size(); ++iT) {
        if (triggerMask->At(fTriggerIds[iCat].at(iT)))
          break;
      }
      // if at least one trigger fired, fill the cutflow
      if (iT != fTriggerIds[iCat].size()) {
        fCutflow[iCat]->Fill(cTrigger);
      }
      else if (!fIgnoreTrigger) {
        // abort if not ignoring the trigger
        continue;
      }
    }

    // Cuts based on object multiplicities

    switch (iCat) {
    case kDielectron:
      if (goodElectrons->GetEntries() < 2)
        continue;
      break;
    case kSingleElectron:
      if (goodElectrons->GetEntries() == 0)
        continue;
      break;
    default:
      if (fVetoElectrons && vetoElectrons->GetEntries() != 0)
        continue;
      break;
    }

    fCutflow[iCat]->Fill(cNElectrons);

    switch (iCat) {
    case kDimuon:
      if (goodMuons->GetEntries() < 2)
        continue;
      break;
    case kSingleMuon:
      if (goodMuons->GetEntries() == 0)
        continue;
      break;
    default:
      if (fVetoMuons && vetoMuons->GetEntries() != 0)
        continue;
      break;
    }

    fCutflow[iCat]->Fill(cNMuons);

    if (fVetoTaus && taus->GetEntries() != 0)
      continue;

    fCutflow[iCat]->Fill(cNTaus);

    switch (iCat) {
    case kPhoton:
      if (goodPhotons->GetEntries() == 0)
        continue;
      break;
    default:
      if (fVetoPhotons && vetoPhotons->GetEntries() != 0)
        continue;
      break;
    };

    fCutflow[iCat]->Fill(cNPhotons);

    if (jets && (jets->GetEntries() < fMinNumJets[iCat] || jets->GetEntries() > fMaxNumJets[iCat]))
      continue;

    fCutflow[iCat]->Fill(cNJets);

    // Cut on met
    // Object definitions in skim can be looser than in the final event selection.
    // Need to try all combinations and take the maximum MET to cut on.

    double maxMet2 = 0.;
    switch (iCat) {
    case kDielectron:
      maxMet2 = calcMet2(goodElectrons);
      break;

    case kSingleElectron:
      maxMet2 = calcMet1(goodElectrons);
      break;

    case kDimuon:
      maxMet2 = calcMet2(goodMuons);
      break;

    case kSingleMuon:
      maxMet2 = calcMet1(goodMuons);
      break;

    case kPhoton:
      maxMet2 = calcMet1(goodPhotons);
      break;

    default:
      maxMet2 = met.Mom().Perp2();
      break;
    }

    if (maxMet2 < fMinMetPt * fMinMetPt)
      continue;

    fCutflow[iCat]->Fill(cMet);

    // Passed category cuts
    fCategoryFlags.At(iCat) = true;
    eventSelected = true;
  }

  if (eventSelected)
    ++fNEventsSelected;
  else
    SkipEvent();
}

//--------------------------------------------------------------------------------------------------
void
MonoXSkimMod::SlaveTerminate()
{
  Info("SlaveTerminate", "Selected events on MonoXSkimMod: %d", fNEventsSelected);

  RetractObj(fCategoryFlags.GetName());
}
