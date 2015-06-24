#include "MitPhysics/Mods/interface/PFTauIDMod.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitAna/DataTree/interface/PFTauCol.h"

using namespace mithep;

ClassImp(mithep::PFTauIDMod)

//--------------------------------------------------------------------------------------------------
PFTauIDMod::PFTauIDMod(const char *name, const char *title) : 
BaseMod(name,title),
  fPFTausName(Names::gkPFTauBrn),
  fGoodPFTausName(ModNames::gkGoodPFTausName),
  fPtMin(15.0),
  fEtaMax(2.4),
  fPtLeadChargedHadronPFCandMin(5.0),
  fIsoChargedHadronPtSumMax(2.0),
  fIsoGammaEtSumMax(3.0),
  fSignalMassMin(0.13),
  fSignalMassMax(2.00),
  fIsLooseId(kTRUE),
  fIsHPSSel(kTRUE),
  fHPSIso("loose"),
  fPFTaus(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void PFTauIDMod::Process()
{
  // Process entries of the tree. 

  LoadBranch(fPFTausName);

  PFTauOArr *GoodTaus = new PFTauOArr;
  GoodTaus->SetName(fGoodPFTausName);

  for (UInt_t i=0; i<fPFTaus->GetEntries(); ++i) {
    const PFTau *tau = fPFTaus->At(i);

    // non-HPS
    if (fIsHPSSel == kFALSE) {
      if (tau->NSignalPFCands() == 0) continue;
      CompositeParticle tauSystem;
      CompositeParticle tauChargedSystem;
      UInt_t nTrk = 0;
      for (UInt_t j=0; j<tau->NSignalPFCands(); ++j) {
        tauSystem.AddDaughter(tau->SignalPFCand(j));
        if (tau->SignalPFCand(j) != 0 &&
            tau->SignalPFCand(j)->Charge() != 0){
          nTrk++;
          tauChargedSystem.AddDaughter(tau->SignalPFCand(j));
        }
      }
      if (nTrk != 1 && nTrk != 3) continue;
      if (TMath::Abs(tau->Charge()) - 1.0 > 0.0001) continue;
      if (tauSystem.Pt() <= fPtMin) continue;
      if (!tau->LeadChargedHadronPFCand()) continue;
      if (tau->LeadChargedHadronPFCand()->Pt() <= fPtLeadChargedHadronPFCandMin) continue;
      if (tau->IsoChargedHadronPtSum() >= fIsoChargedHadronPtSumMax) continue;
      if (tau->IsoGammaEtSum() >= fIsoGammaEtSumMax) continue;
      if (tauChargedSystem.Mass() <= fSignalMassMin || tauChargedSystem.Mass() >= fSignalMassMax) continue;
    }

    // HPS
    else { 
      // basic selections
      if (tau->Pt() < fPtMin) continue;
      if (fabs(tau->Eta()) > fEtaMax) continue;
      if (tau->DiscriminationByDecayModeFinding() < 0.5) continue;

      // default HPS
      if (!fIsLooseId) { 
        // reject leptons 
        if (tau->DiscriminationByLooseElectronRejection() < 0.5) continue;
        if (tau->DiscriminationByLooseMuonRejection() < 0.5) continue;
        // default isolation working points
        if (fHPSIso.Contains("loose",TString::kIgnoreCase)) {
          if (!tau->LooseCombinedIsolationDBSumPtCorr3Hits()) continue;
        }
        else if (fHPSIso.Contains("med",TString::kIgnoreCase)) {
          if (!tau->MediumCombinedIsolationDBSumPtCorr3Hits()) continue;
        }
        else if (fHPSIso.Contains("tight",TString::kIgnoreCase)) {
          if (!tau->TightCombinedIsolationDBSumPtCorr3Hits()) continue;
        }
        else {
          SendError(kWarning,"Process","ERROR: HPS Isolation not properly defined!");
        }
      }
      // loose Id
      else {
        if (tau->RawCombinedIsolationDBSumPtCorr3Hits() >5) continue;
      }
    }

    // add good tau to output collection
    tau->Mark();
    GoodTaus->Add(tau);
  }

  // sort according to pt
  GoodTaus->Sort();

  // add to event for other modules to use
  AddObjThisEvt(GoodTaus);  
}

//--------------------------------------------------------------------------------------------------
void PFTauIDMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here, we just request the
  // Tau collection branch.

  // note: for safety, this overrides SetPFTausName
  if (fIsHPSSel)
    fPFTausName = TString(Names::gkHPSTauBrn);
  else
    fPFTausName = TString(Names::gkPFTauBrn);

  ReqBranch(fPFTausName, fPFTaus);
}
