#include "MitPhysics/Mods/interface/JetIDModRun1.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/CaloJetCol.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitPhysics/Utils/interface/JetTools.h"
#include "MitPhysics/Utils/interface/JetIDMVA.h"
#include "MitCommon/MathTools/interface/MathUtils.h"

#include "TSystem.h"

ClassImp(mithep::JetIDModRun1)

mithep::JetIDModRun1::JetIDModRun1(const char* name, const char* title) :
  BaseMod(name, title),
  fJetsName(ModNames::gkPubJetsName),
  fGoodJetsName(ModNames::gkGoodJetsName),  
  fVertexName(ModNames::gkGoodVertexesName),
  fUseJetCorrection(kTRUE),
  fJetPtCut(35.0),
  fJetEtaMaxCut(5.0),
  fJetEEMFractionMinCut(0.01),
  fMinChargedHadronFraction(0.),
  fMaxChargedHadronFraction(1.),
  fMinNeutralHadronFraction(0.),
  fMaxNeutralHadronFraction(1.),
  fMinChargedEMFraction(0.),
  fMaxChargedEMFraction(1.),
  fMinNeutralEMFraction(0.),
  fMaxNeutralEMFraction(1.),
  fApplyBetaCut(kFALSE),
  fApplyPFLooseId(kFALSE),
  fApplyMVACut(kFALSE),
  fApplyMVACHS(kFALSE)
{
}

void
mithep::JetIDModRun1::Process()
{
  // Process entries of the tree. 

  mithep::JetCol const* inJets = GetObjThisEvt<JetCol>(fJetsName);
  if (!inJets) {
    SendError(kAbortModule, "Process","Pointer to input jet collection %s null.", fJetsName.Data());
    return;
  }
  mithep::VertexCol const* inVertices = GetObjThisEvt<VertexOArr>(fVertexName);
  if (!inVertices) {
    SendError(kAbortModule, "Process","Pointer to input jet collection %s null.", fVertexName.Data());
    return;
  }

  mithep::JetOArr* GoodJets = new mithep::JetOArr; 
  GoodJets->SetName(fGoodJetsName);

  // loop over jets
  for (UInt_t i=0; i < inJets->GetEntries(); ++i) {
    mithep::Jet const* jet = inJets->At(i);

    if (jet->AbsEta() > fJetEtaMaxCut) 
      continue;

    Double_t jetpt;
    if (fUseJetCorrection)
      jetpt = jet->Pt();
    else
      jetpt = jet->RawMom().Pt();

    if (jetpt < fJetPtCut)
      continue;

    if (fJetEEMFractionMinCut > 0) {
      mithep::CaloJet const* caloJet = dynamic_cast<mithep::CaloJet const*>(jet); 
      // The 2.6 value is hardcoded, no reason to change that value in CMS
      if (caloJet && caloJet->AbsEta() < 2.6 &&
          caloJet->EnergyFractionEm() < fJetEEMFractionMinCut)
        continue;
    }

    PFJet const* pfJet = dynamic_cast<mithep::PFJet const*>(jet);

    // Jet to vertex association cut
    if (fApplyBetaCut || fApplyPFLooseId || fApplyMVACut) {
      if (!pfJet)
        continue;

      if (fApplyBetaCut && !JetTools::PassBetaVertexAssociationCut(pfJet, inVertices->At(0), inVertices, 0.2))
        continue;

      if (fApplyPFLooseId && !JetTools::passPFLooseId(pfJet))
        continue;
      
      if (fApplyMVACut && !fJetIDMVA->pass(pfJet, inVertices->At(0), inVertices))
        continue;
    }

    if (pfJet) {
      double chargedHadronFraction = pfJet->ChargedHadronEnergy() / pfJet->E();
      if (chargedHadronFraction < fMinChargedHadronFraction || chargedHadronFraction > fMaxChargedHadronFraction)
        continue;

      double neutralHadronFraction = pfJet->NeutralHadronEnergy() / pfJet->E();
      if (neutralHadronFraction < fMinNeutralHadronFraction || neutralHadronFraction > fMaxNeutralHadronFraction)
        continue;

      double chargedEMFraction = pfJet->ChargedEmEnergy() / pfJet->E();
      if (chargedEMFraction < fMinChargedEMFraction || chargedEMFraction > fMaxChargedEMFraction)
        continue;

      double neutralEMFraction = pfJet->NeutralEmEnergy() / pfJet->E();
      if (neutralEMFraction < fMinNeutralEMFraction || neutralEMFraction > fMaxNeutralEMFraction)
        continue;
    }

    // add good jet to collection
    GoodJets->Add(jet);
  }

  // sort according to pt
  GoodJets->Sort();

  if (GoodJets->GetEntries() < fMinNJets) {
    SkipEvent();
    return;
  }
  
  // add to event for other modules to use
  AddObjThisEvt(GoodJets);  
}

void
mithep::JetIDModRun1::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here, we typically
  // initialize histograms and other analysis objects and request branches. For this module, we
  // request a branch of the MitTree.
  
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
}

void
mithep::JetIDModRun1::SlaveTerminate()
{
  delete fJetIDMVA;
}
