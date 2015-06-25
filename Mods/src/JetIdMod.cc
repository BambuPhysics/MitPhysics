#include "MitPhysics/Mods/interface/JetIdMod.h"
#include "MitPhysics/Utils/interface/JetTools.h"
#include "MitPhysics/Utils/interface/JetIDMVA.h"
#include "MitAna/DataTree/interface/CaloJet.h"
#include "MitAna/DataTree/interface/PFJet.h"
#include "MitAna/DataTree/interface/Names.h"
#include "MitCommon/MathTools/interface/MathUtils.h"

#include "TSystem.h"

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
  if (fApplyMVACut && !fJetIDMVA) {
    fJetIDMVA = new JetIDMVA();
    TString dataDir(gSystem->Getenv("MIT_DATA"));
    if (dataDir.Length() == 0) {
      SendError(kAbortModule, "SlaveBegin", "MIT_DATA environment is not set.");
      return;
    }

    if (fApplyMVACHS)
      fJetIDMVA->Initialize(JetIDMVA::kLoose,
                            dataDir + "/TMVAClassificationCategory_JetID_53X_chs_Dec2012.weights.xml",
                            dataDir + "/TMVAClassificationCategory_JetID_53X_chs_Dec2012.weights.xml",
                            JetIDMVA::k53,
                            dataDir + "/JetIDMVA_JetIdParams.py");
    else
      fJetIDMVA->Initialize(JetIDMVA::kLoose,
                            dataDir + "/mva_JetID_lowpt.weights.xml",
                            dataDir + "/mva_JetID_highpt.weights.xml",
                            JetIDMVA::kBaseline,
                            dataDir + "/JetIDMVA_JetIdParams.py");
  }

  if (fJetIDMVA && !fJetIDMVA->IsInitialized())
    SendError(kAbortModule, "SlaveBegin", "Jet ID MVA is not initialized.");

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
  xaxis->SetBinLabel(cPFLooseId + 1, "PFLooseId");
  xaxis->SetBinLabel(cBeta + 1, "Beta");
  xaxis->SetBinLabel(cMVA + 1, "MVA");
}

//--------------------------------------------------------------------------------------------------
Bool_t
mithep::JetIdMod::IsGood(mithep::Jet const& jet)
{
  fCutFlow->Fill(cAll);

  if(jet.AbsEta() > fEtaMax)
    return false;
  fCutFlow->Fill(cEta);

  Double_t jetpt;
  if (fUseJetCorrection)
    jetpt = jet.Pt();
  else
    jetpt = jet.RawMom().Pt();

  if (jetpt < fPtMin)
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
    if (chargedHadronFraction < fMinChargedHadronFraction || chargedHadronFraction > fMaxChargedHadronFraction)
      return false;
    fCutFlow->Fill(cChargedHFrac);

    double neutralHadronFraction = pfJet->NeutralHadronEnergy() / pfJet->E();
    if (neutralHadronFraction < fMinNeutralHadronFraction || neutralHadronFraction > fMaxNeutralHadronFraction)
      return false;
    fCutFlow->Fill(cNeutralHFrac);

    double chargedEMFraction = pfJet->ChargedEmEnergy() / pfJet->E();
    if (chargedEMFraction < fMinChargedEMFraction || chargedEMFraction > fMaxChargedEMFraction)
      return false;
    fCutFlow->Fill(cChargedEMFrac);

    double neutralEMFraction = pfJet->NeutralEmEnergy() / pfJet->E();
    if (neutralEMFraction < fMinNeutralEMFraction || neutralEMFraction > fMaxNeutralEMFraction)
      return false;
    fCutFlow->Fill(cNeutralEMFrac);

    if (fApplyPFLooseId && !JetTools::passPFLooseId(pfJet))
      return false;
    fCutFlow->Fill(cPFLooseId);

    if (fApplyBetaCut && !JetTools::PassBetaVertexAssociationCut(pfJet, GetVertices()->At(0), GetVertices(), 0.2))
      return false;
    fCutFlow->Fill(cBeta);

    if (fApplyMVACut && !fJetIDMVA->pass(pfJet, GetVertices()->At(0), GetVertices()))
      return false;
    fCutFlow->Fill(cMVA);
  }

  return true;
}
