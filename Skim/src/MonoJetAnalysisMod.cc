#include "MitPhysics/Skim/interface/MonoJetAnalysisMod.h"

#include "MitAna/DataTree/interface/MetCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/PFTauCol.h"
#include "MitAna/DataTree/interface/PhotonCol.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/PFJet.h"
#include "MitAna/DataTree/interface/TriggerTable.h"
#include "MitAna/DataTree/interface/TriggerMask.h"
#include "MitAna/TreeMod/interface/HLTFwkMod.h"
#include "MitCommon/MathTools/interface/MathUtils.h"

#include <limits>

using namespace mithep;
ClassImp(mithep::MonoJetAnalysisMod)

TString MonoJetAnalysisMod::fgMonoJetCategories[nMonoJetCategories] = {
  "Signal",
  "Dielectron",
  "Dimuon",
  "SingleElectron",
  "SingleMuon",
  "Photon"
};

//--------------------------------------------------------------------------------------------------
MonoJetAnalysisMod::MonoJetAnalysisMod(const char *name, const char *title) :
  BaseMod(name, title)
{
  // cuts
  double dbig = std::numeric_limits<double>::max();
  std::fill_n(fCategoryActive, nMonoJetCategories, false);
  std::fill_n(fMinNumJets, nMonoJetCategories, 0xffffffff);
  std::fill_n(fMaxNumJets, nMonoJetCategories, 0);
  std::fill_n(fMinLeadJetPt, nMonoJetCategories, dbig);
  std::fill_n(fMaxJetEta, nMonoJetCategories, 0.);
  std::fill_n(fMinMetPt, nMonoJetCategories, dbig);
  std::fill_n(fMinChargedHadronFrac, nMonoJetCategories, dbig);
  std::fill_n(fMaxNeutralHadronFrac, nMonoJetCategories, 0.);
  std::fill_n(fMaxNeutralEmFrac, nMonoJetCategories, 0.);

  fCategoryFlags.Resize(nMonoJetCategories);
}

//--------------------------------------------------------------------------------------------------
void
MonoJetAnalysisMod::SlaveBegin()
{
  // Selection Histograms
  for(unsigned iCat = 0; iCat != nMonoJetCategories; ++iCat){
    if (!fCategoryActive[iCat])
      continue;

    // Histograms after preselection
    AddTH1(fCutflow[iCat], "hCutflow" + fgMonoJetCategories[iCat], ";;Events", nCuts, 0., nCuts);
    auto* xaxis = fCutflow[iCat]->GetXaxis();
    xaxis->SetBinLabel(cAll + 1, "All");
    xaxis->SetBinLabel(cTrigger + 1, "Trigger");
    xaxis->SetBinLabel(cNElectrons + 1, "NElectrons");
    xaxis->SetBinLabel(cNMuons + 1, "NMuons");
    xaxis->SetBinLabel(cNPhotons + 1, "NPhotons");
    xaxis->SetBinLabel(cNTaus + 1, "NTaus");
    xaxis->SetBinLabel(cLeadJet + 1, "LeadJet");
    xaxis->SetBinLabel(cNJets + 1, "NJets");
    xaxis->SetBinLabel(cMet + 1, "Met");

    AddTH1(fRealMetPt[iCat],  "hRealMetPt" + fgMonoJetCategories[iCat], ";Real E_{T}^{miss} (GeV);Events / 4 GeV", 100, 0., 400.);
    AddTH1(fMetPt[iCat],  "hMetPt" + fgMonoJetCategories[iCat], ";Emul. E_{T}^{miss} (GeV);Events / 4 GeV", 100, 0., 400.);
    AddTH1(fJetPt[iCat],  "hJetPt" + fgMonoJetCategories[iCat], ";Jet p_{T} (GeV);Events / 4 GeV", 100, 0., 400.);
    AddTH1(fJetEta[iCat], "hJetEta" + fgMonoJetCategories[iCat],";Jet #eta;Events / 0.2", 50, -5., 5.);
  }

  PublishObj(&fCategoryFlags);

  fNEventsSelected = 0;
}

void
MonoJetAnalysisMod::BeginRun()
{
  if (!HasHLTInfo() || !GetHltFwkMod()->HasData())
    return;

  auto* hltTable = GetHLTTable();
  for (unsigned iCat = 0; iCat != nMonoJetCategories; ++iCat) {
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
MonoJetAnalysisMod::Process()
{
  auto* mets = GetObject<MetCol>(fMetName);
  if (!mets)
    SendError(kAbortAnalysis, "Process", "Could not find " + fMetName);
  auto* met = mets->At(0);

  auto* jets = GetObject<JetCol>(fJetsName);
  if (!jets)
    SendError(kAbortAnalysis, "Process", "Could not find " + fJetsName);

  auto* electrons = GetObject<ElectronCol>(fVetoElectronsName);
  if (!electrons)
    SendError(kAbortAnalysis, "Process", "Could not find " + fVetoElectronsName);
  NFArrBool* electronMask = 0;
  if (fCategoryActive[kDielectron] || fCategoryActive[kSingleElectron]) {
    electronMask = GetObject<NFArrBool>(fElectronMaskName);
    if (!electronMask)
      SendError(kAbortAnalysis, "Process", "Could not find " + fElectronMaskName);
  }

  auto* muons = GetObject<MuonCol>(fVetoMuonsName);
  if (!muons)
    SendError(kAbortAnalysis, "Process", "Could not find " + fVetoMuonsName);
  NFArrBool* muonMask = 0;
  if (fCategoryActive[kDimuon] || fCategoryActive[kSingleMuon]) {
    muonMask = GetObject<NFArrBool>(fMuonMaskName);
    if (!muonMask)
      SendError(kAbortAnalysis, "Process", "Could not find " + fMuonMaskName);
  }

  auto* taus = GetObject<PFTauCol>(fVetoTausName);
  if (!taus)
    SendError(kAbortAnalysis, "Process", "Could not find " + fVetoTausName);

  auto* photons = GetObject<PhotonCol>(fVetoPhotonsName);
  if (!photons)
    SendError(kAbortAnalysis, "Process", "Could not find " + fVetoPhotonsName);
  NFArrBool* photonMask = 0;
  if (fCategoryActive[kPhoton]) {
    photonMask = GetObject<NFArrBool>(fPhotonMaskName);
    if (!photonMask)
      SendError(kAbortAnalysis, "Process", "Could not find " + fPhotonMaskName);
  }

  TriggerMask* triggerMask = 0;
  if (HasHLTInfo() && GetHltFwkMod()->HasData()) {
    triggerMask = GetObject<TriggerMask>(fTriggerBitsName);
    if (!triggerMask)
      SendError(kAbortAnalysis, "Process", "Could not find " + fTriggerBitsName);
  }

  bool skipEvent = true;

  for (unsigned iCat = 0; iCat != nMonoJetCategories; ++iCat) {
    fCategoryFlags.At(iCat) = false;
    
    if (!fCategoryActive[iCat])
      continue;

    fCutflow[iCat]->Fill(cAll);

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
      // Mask applied on top of veto electrons
      // -> there should be only two entries and both must be true
      if (electronMask->GetEntries() != 2 || !(electronMask->At(0) && electronMask->At(1)))
        continue;
      break;
    case kSingleElectron:
      if (electronMask->GetEntries() != 1 || !electronMask->At(0))
        continue;
      break;
    default:
      if (electrons->GetEntries() != 0)
        continue;
      break;
    }

    fCutflow[iCat]->Fill(cNElectrons);

    switch (iCat) {
    case kDimuon:
      if (muonMask->GetEntries() != 2 || !(muonMask->At(0) && muonMask->At(1)))
        continue;
      break;
    case kSingleMuon:
      if (muonMask->GetEntries() != 1 || !muonMask->At(0))
        continue;
      break;
    default:
      if (muons->GetEntries() != 0)
        continue;
      break;
    }

    fCutflow[iCat]->Fill(cNMuons);

    if (taus->GetEntries() != 0)
      continue;

    fCutflow[iCat]->Fill(cNTaus);

    switch (iCat) {
    case kPhoton:
      if (photonMask->GetEntries() != 1 || !photonMask->At(0))
        continue;
      break;
    default:
      if (photons->GetEntries() != 0)
        continue;
      break;
    };

    fCutflow[iCat]->Fill(cNPhotons);
    
    // Cuts on the jets
    std::vector<Jet*> goodJets;

    for (unsigned iJ = 0; iJ < jets->GetEntries(); ++iJ) {
      auto& jet(*jets->At(iJ));
      if (jet.AbsEta() >= fMaxJetEta[iCat])
        continue;

      double rawE = jet.RawMom().E();
      auto& pfJet(static_cast<PFJet&>(jet));

      if (pfJet.ChargedHadronEnergy() / rawE < fMinChargedHadronFrac[iCat])
        continue;

      if (pfJet.NeutralHadronEnergy() / rawE > fMaxNeutralHadronFrac[iCat])
        continue;

      if (pfJet.NeutralEmEnergy() / rawE > fMaxNeutralEmFrac[iCat])
        continue;

      switch (iCat) {
      case kDielectron:
        if (MathUtils::DeltaR(jet, *electrons->At(1)) < 0.4)
          continue;
        //fallthrough
      case kSingleElectron:
        if (MathUtils::DeltaR(jet, *electrons->At(0)) < 0.4)
          continue;
        break;
      case kDimuon:
        if (MathUtils::DeltaR(jet, *muons->At(1)) < 0.4)
          continue;
        //fallthrough
      case kSingleMuon:
        if (MathUtils::DeltaR(jet, *muons->At(0)) < 0.4)
          continue;
        break;
      case kPhoton:
        if (MathUtils::DeltaR(jet, *photons->At(0)) < 0.4)
          continue;
        break;
      default:
        break;
      }

      goodJets.push_back(&jet);
    }

    if (fMinNumJets[iCat] != 0) {
      if (goodJets.size() < fMinNumJets[iCat] || goodJets[0]->Pt() < fMinLeadJetPt[iCat])
        continue;

      fCutflow[iCat]->Fill(cLeadJet);
    }

    if (goodJets.size() > fMaxNumJets[iCat])
      continue;

    fCutflow[iCat]->Fill(cNJets);

    // Cut on met
    double metPt = 0.;
    switch (iCat) {
    case kSignal:
      metPt = met->Pt();
      break;

    case kDielectron:
      {
        auto emulMet(met->Mom());
        emulMet += electrons->At(0)->Mom();
        emulMet += electrons->At(1)->Mom();
        metPt = emulMet.Pt();
      }
      break;

    case kDimuon:
      {
        auto emulMet(met->Mom());
        emulMet += muons->At(0)->Mom();
        emulMet += muons->At(1)->Mom();
        metPt = emulMet.Pt();
      }
      break;

    case kSingleElectron:
      {
        auto emulMet(met->Mom());
        emulMet += electrons->At(0)->Mom();
        metPt = emulMet.Pt();
      }
      break;

    case kSingleMuon:
      {
        auto emulMet(met->Mom());
        emulMet += muons->At(0)->Mom();
        metPt = emulMet.Pt();
      }
      break;

    case kPhoton:
      {
        auto& photon(*photons->At(0));
        auto emulMet(met->Mom());
        emulMet += photon.Mom();
        for (unsigned iF = 0; iF != photon.NFootprintCandidates(); ++iF)
          emulMet += photon.FootprintCandidate(iF)->Mom();
        metPt = emulMet.Pt();
      }
      break;

    default:
      break;
    }

    if (metPt < fMinMetPt[iCat])
      continue;

    fCutflow[iCat]->Fill(cMet);

    // Passed category cuts
    fCategoryFlags.At(iCat) = true;
    skipEvent = false;

    // Fill Preselection Histograms
    fRealMetPt[iCat]->Fill(met->Pt());
    fMetPt[iCat]->Fill(metPt);
    for (auto* jet : goodJets) {
      fJetPt[iCat]->Fill(jet->Pt());
      fJetEta[iCat]->Fill(jet->Eta());
    }
  }

  // Skip the event if does not passes the cuts
  if (skipEvent)
    SkipEvent();
  else
    ++fNEventsSelected;
}

//--------------------------------------------------------------------------------------------------
void
MonoJetAnalysisMod::SlaveTerminate()
{
  Info("SlaveTerminate", "Selected events on MonoJetAnalysisMod: %d", fNEventsSelected);

  RetractObj(fCategoryFlags.GetName());
}
