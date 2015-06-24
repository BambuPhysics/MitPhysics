#include "MitPhysics/Mods/interface/JetIdMod.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/CaloJetCol.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitPhysics/Utils/interface/JetTools.h"
#include "MitPhysics/Utils/interface/JetIDMVA.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataTree/interface/JetFwd.h"

#include "TSystem.h"
#include <algorithm>
#include <limits>

using namespace mithep;

ClassImp(mithep::JetIdMod)

template<class T>
void
JetIdMod::GetAuxInput(JetIdMod::AuxInput inputCol, TObject const** aux)
{
  aux[inputCol] = GetObject<T>(fAuxInputNames[inputCol], true);
  if (!aux[inputCol])
    SendError(kAbortAnalysis, "GetAuxInput", "Could not retrieve auxiliary input " + fAuxInputNames[inputCol]);
}

//--------------------------------------------------------------------------------------------------
JetIdMod::JetIdMod(const char *name, const char *title) :
  IdMod(name, title)
{
  fOutput = new JetOArr(32, TString(name) + "Output");
  SetInputName(ModNames::gkPubJetsName);
  fAuxInputNames[kVertices] = ModNames::gkGoodVertexesName;
}

//--------------------------------------------------------------------------------------------------
void
JetIdMod::IdBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we just request the jet collection branch.

  // If we use MVA Id, need to load MVA weights
  if (fApplyMVACut) {
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

  fCutFlow->SetBins(nCuts, 0., double(nCuts));
  TAxis* xaxis = fCutFlow->GetXaxis();
  xaxis->SetBinLabel(cAll + 1, "All");
  xaxis->SetBinLabel(cPt + 1, "Pt");
  xaxis->SetBinLabel(cEta + 1, "Eta");
  xaxis->SetBinLabel(cEEMFraction+1,"EEMFraction");
  xaxis->SetBinLabel(cBeta+1,"Beta");
  xaxis->SetBinLabel(cMVA+1,"MVA");
  xaxis->SetBinLabel(cPFLooseId+1,"PFLooseId");
  xaxis->SetBinLabel(cChargedHFrac,"chargedHadronFraction");
  xaxis->SetBinLabel(cNeutralHFrac,"neutralHadronFraction");
  xaxis->SetBinLabel(cChargedEMFrac,"chargedEMFraction");
  xaxis->SetBinLabel(cNeutralEMFrac,"neutralEMFraction");
}

//--------------------------------------------------------------------------------------------------
void
JetIdMod::Process()
{
  // Process entries of the tree.
  auto* jets = GetObject<mithep::JetCol>(fInputName, true);
  if (!jets)
    SendError(kAbortAnalysis, "Process", "Jets not found");

  TObject const* aux[nAuxInputs] = {};

  GetAuxInput<mithep::VertexCol>(kVertices,aux);
  if (!aux[kVertices])
    SendError(kAbortAnalysis,"Process","Vertices not found");

  mithep::JetOArr* goodJets = 0;
  if (fIsFilterMode) {
    goodJets = static_cast<mithep::JetOArr*>(fOutput);
    goodJets->Reset();
  }
  else {
    fFlags.Resize(jets->GetEntries());
    for (unsigned int i = 0; i < jets->GetEntries(); i++)
      fFlags.At(i) = false;
  }
  UInt_t nGoodJets = 0;
  for (UInt_t i = 0; i < jets->GetEntries(); i++) {
    Jet const& jet = *jets->At(i);

    fCutFlow->Fill(cAll);

    if(jet.AbsEta() > fEtaMax)
      continue;
    fCutFlow->Fill(cEta);

    Double_t jetpt;
    if (fUseJetCorrection)
      jetpt = jet.Pt();
    else
      jetpt = jet.RawMom().Pt();

    if (jetpt < fJetPtCut)
      continue;
    fCutFlow->Fill(cPt);

    if (fJetEEMFractionMinCut > 0) {
      mithep::CaloJet const* caloJet = dynamic_cast<mithep::CaloJet const*>(&jet); 
      // The 2.6 value is hardcoded, no reason to change that value in CMS
      if (caloJet && caloJet->AbsEta() < 2.6 &&
          caloJet->EnergyFractionEm() < fJetEEMFractionMinCut)
        continue;
    }
    fCutFlow->Fill(cEEMFraction);

    PFJet const* pfJet = dynamic_cast<mithep::PFJet const*>(&jet);

    if (!pfJet) continue;
  
    if (fApplyBetaCut && !JetTools::PassBetaVertexAssociationCut(pfJet, static_cast<mithep::VertexCol const*>(aux[kVertices])->At(0), static_cast<mithep::VertexCol const*>(aux[kVertices]), 0.2))
      continue;
    fCutFlow->Fill(cBeta);
    if (fApplyPFLooseId && !JetTools::passPFLooseId(pfJet))
      continue;
    fCutFlow->Fill(cPFLooseId);
    
    if (fApplyMVACut && !fJetIDMVA->pass(pfJet, static_cast<mithep::VertexCol const*>(aux[kVertices])->At(0), static_cast<mithep::VertexCol const*>(aux[kVertices])))
      continue;
    fCutFlow->Fill(cMVA);

    double chargedHadronFraction = pfJet->ChargedHadronEnergy() / pfJet->E();
    if (chargedHadronFraction < fMinChargedHadronFraction || chargedHadronFraction > fMaxChargedHadronFraction)
      continue;
    fCutFlow->Fill(cChargedHFrac);

    double neutralHadronFraction = pfJet->NeutralHadronEnergy() / pfJet->E();
    if (neutralHadronFraction < fMinNeutralHadronFraction || neutralHadronFraction > fMaxNeutralHadronFraction)
      continue;
    fCutFlow->Fill(cNeutralHFrac);

    double chargedEMFraction = pfJet->ChargedEmEnergy() / pfJet->E();
    if (chargedEMFraction < fMinChargedEMFraction || chargedEMFraction > fMaxChargedEMFraction)
      continue;
    fCutFlow->Fill(cChargedEMFrac);

    double neutralEMFraction = pfJet->NeutralEmEnergy() / pfJet->E();
    if (neutralEMFraction < fMinNeutralEMFraction || neutralEMFraction > fMaxNeutralEMFraction)
      continue;
    fCutFlow->Fill(cNeutralEMFrac);

    if (fIsFilterMode)
      goodJets->Add(&jet);
    else
      fFlags.At(i) = true;
    nGoodJets++;
  }

  if (nGoodJets < fMinNJets) {
    SkipEvent();
    return;
  }

  if (fIsFilterMode) {
    // sort according to pt
    goodJets->Sort();
  }
}
