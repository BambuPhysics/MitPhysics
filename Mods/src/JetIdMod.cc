#include "MitPhysics/Mods/interface/JetIdMod.h"
#include "MitPhysics/Utils/interface/JetIDMVA.h"
#include "MitAna/DataTree/interface/CaloJet.h"
#include "MitAna/DataTree/interface/PFJet.h"
#include "MitAna/DataTree/interface/Names.h"
#include "MitCommon/MathTools/interface/MathUtils.h"

ClassImp(mithep::JetIdMod)

//--------------------------------------------------------------------------------------------------
mithep::JetIdMod::JetIdMod(const char *name, const char *title) :
  IdMod<mithep::Jet>(name, title)
{
  fInputName = mithep::Names::gkPFJetBrn;
}

//--------------------------------------------------------------------------------------------------
void
mithep::JetIdMod::IdBegin()
{
  // If we use MVA Id, need to load MVA weights
  if (fMVATrainingSet != JetIDMVA::nMVATypes && !fJetIDMVA) {
    fJetIDMVA = new JetIDMVA();
    fJetIDMVA->Initialize(JetIDMVA::CutType(fMVACutWP), JetIDMVA::MVAType(fMVATrainingSet), fMVAWeightsFile, fMVACutsFile);

    fOwnJetIDMVA = true;
  }

  if (fJetIDMVA && !fJetIDMVA->IsInitialized())
    SendError(kAbortModule, "SlaveBegin", "Jet ID MVA is not initialized.");

  if (fUseClassicBetaForMVA)
    fJetIDMVA->fDZCut = -1.; // negative dz cut triggers classic beta calculation

  fCutFlow->SetBins(nCuts, 0., double(nCuts));
  TAxis* xaxis = fCutFlow->GetXaxis();
  xaxis->SetBinLabel(cAll + 1, "All");
  xaxis->SetBinLabel(cPt + 1, "Pt");
  xaxis->SetBinLabel(cEta + 1, "Eta");
  xaxis->SetBinLabel(cEEMFraction + 1, "EEMFraction");
  xaxis->SetBinLabel(cChargedHFrac + 1, "chargedHadronFraction");
  xaxis->SetBinLabel(cNeutralHFrac + 1, "neutralHadronFraction");
  xaxis->SetBinLabel(cChargedEMFrac + 1, "chargedEMFraction");
  xaxis->SetBinLabel(cNeutralEMFrac + 1, "neutralEMFraction");
  xaxis->SetBinLabel(cMuonFrac + 1, "muonFraction");
  xaxis->SetBinLabel(cNPFCandidates + 1, "nPFCandidates");
  xaxis->SetBinLabel(cNChargedPFCandidates + 1, "nChargedPFCandidates");
  xaxis->SetBinLabel(cPFId + 1, "PFId");
  xaxis->SetBinLabel(cBeta + 1, "Beta");
  xaxis->SetBinLabel(cMVA + 1, "MVA");
}

void
mithep::JetIdMod::IdTerminate()
{
  if (fOwnJetIDMVA) {
    delete fJetIDMVA;
    fJetIDMVA = 0;
  }
}

//--------------------------------------------------------------------------------------------------
Bool_t
mithep::JetIdMod::IsGood(mithep::Jet const& jet)
{
  fCutFlow->Fill(cAll);

  if(jet.AbsEta() > fEtaMax)
    return false;
  fCutFlow->Fill(cEta);

  // Must use corrected jets as input
  if (jet.Pt() < fPtMin)
    return false;
  fCutFlow->Fill(cPt);

  if (fJetEEMFractionMinCut > 0) {
    mithep::CaloJet const* caloJet = dynamic_cast<mithep::CaloJet const*>(&jet); 
    // The 2.6 value is hardcoded, no reason to change that value in CMS
    if (caloJet && caloJet->AbsEta() < 2.6 &&
        caloJet->EnergyFractionEm() < fJetEEMFractionMinCut)
      return false;
  }
  fCutFlow->Fill(cEEMFraction);

  PFJet const* pfJet = dynamic_cast<mithep::PFJet const*>(&jet);

  if (pfJet) {
    double chargedHadronFraction = pfJet->ChargedHadronEnergy() / pfJet->E();
    if (pfJet->AbsEta()<2.4 && (chargedHadronFraction < fMinChargedHadronFraction || chargedHadronFraction > fMaxChargedHadronFraction))
      return false;
    fCutFlow->Fill(cChargedHFrac);

    double neutralHadronFraction = pfJet->NeutralHadronEnergy() / pfJet->E();
    if (neutralHadronFraction < fMinNeutralHadronFraction || neutralHadronFraction > fMaxNeutralHadronFraction)
      return false;
    fCutFlow->Fill(cNeutralHFrac);

    double chargedEMFraction = pfJet->ChargedEmEnergy() / pfJet->E();
    if (pfJet->AbsEta()<2.4 && (chargedEMFraction < fMinChargedEMFraction || chargedEMFraction > fMaxChargedEMFraction))
      return false;
    fCutFlow->Fill(cChargedEMFrac);

    double neutralEMFraction = pfJet->NeutralEmEnergy() / pfJet->E();
    if (neutralEMFraction < fMinNeutralEMFraction || neutralEMFraction > fMaxNeutralEMFraction)
      return false;
    fCutFlow->Fill(cNeutralEMFrac);

    double muonFraction = pfJet->MuonEnergy() / pfJet->E();
    if (muonFraction < fMinMuonFraction || muonFraction > fMaxMuonFraction)
      return false;
    fCutFlow->Fill(cMuonFrac);

    if (pfJet->NPFCands() < fMinNPFCandidates)
      return false;
    fCutFlow->Fill(cNPFCandidates);

    if (pfJet->AbsEta()<2.4 && pfJet->ChargedMultiplicity() < fMinNChargedPFCandidates)
      return false;
    fCutFlow->Fill(cNChargedPFCandidates);

    if (fPFId != JetTools::nPFIdWorkingPoints && !JetTools::passPFId(pfJet, JetTools::PFIdWorkingPoint(fPFId)))
      return false;
    fCutFlow->Fill(cPFId);

    if (fApplyBetaCut && !JetTools::PassBetaVertexAssociationCut(pfJet, GetVertices()->At(0), GetVertices(), 0.2))
      return false;
    fCutFlow->Fill(cBeta);

    if (fJetIDMVA && !fJetIDMVA->pass(pfJet, GetVertices()->At(0), GetVertices()))
      return false;
    fCutFlow->Fill(cMVA);
  }

  return true;
}
