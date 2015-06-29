// $Id: JetCleaningMod.cc,v 1.17 2010/06/28 21:07:06 ceballos Exp $

#include "MitPhysics/Mods/interface/JetCleaningMod.h"
#include "MitAna/DataTree/interface/PhotonCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/TauCol.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitPhysics/Init/interface/ModNames.h"

using namespace mithep;

ClassImp(mithep::JetCleaningMod)

//--------------------------------------------------------------------------------------------------
JetCleaningMod::JetCleaningMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fCleanElectronsName(ModNames::gkCleanElectronsName),        
  fCleanMuonsName(ModNames::gkCleanMuonsName),        
  fCleanPhotonsName(ModNames::gkCleanPhotonsName),        
  fCleanTausName(ModNames::gkCleanTausName),        
  fGoodJetsName(ModNames::gkGoodJetsName),        
  fCleanJets(0),
  fMinDeltaRToElectron(0.3),
  fMinDeltaRToMuon(0.3),
  fMinDeltaRToPhoton(0.3),
  fMinDeltaRToTau(0.3),
  fApplyPhotonRemoval(kFALSE),
  fApplyTauRemoval(kFALSE)
{
  // Constructor.
  fCleanJets = new JetOArr;
  fCleanJets->SetName(ModNames::gkCleanJetsName);
}

//--------------------------------------------------------------------------------------------------
void JetCleaningMod::Process()
{
  // Process entries of the tree.

  // get input collections
  const JetCol  *GoodJets       = GetObject<JetCol>(fGoodJetsName);
  const ElectronCol *CleanElectrons = 0;
  if (!fCleanElectronsName.IsNull())
    CleanElectrons = GetObject<ElectronCol>(fCleanElectronsName);
  const MuonCol *CleanMuons = 0;
  if (!fCleanMuonsName.IsNull())
    CleanMuons = GetObject<MuonCol>(fCleanMuonsName);
  const PhotonCol   *CleanPhotons   = 0;
  if (!fCleanPhotonsName.IsNull())
    CleanPhotons    = GetObject<PhotonCol>(fCleanPhotonsName);
  const TauCol   *CleanTaus   = 0;
  if (fApplyTauRemoval && !fCleanTausName.IsNull())
    CleanTaus    = GetObject<TauCol>(fCleanTausName);

  // create output collection

  // remove any jet that overlaps in eta, phi with an isolated electron.    
  for (UInt_t i=0; i<GoodJets->GetEntries(); ++i) {
    const Jet *jet = GoodJets->At(i);        

    // check for overlap with an Electron
    Bool_t isElectronOverlap = kFALSE;
    if (CleanElectrons) {
      UInt_t n = CleanElectrons->GetEntries();
      for (UInt_t j=0; j<n; ++j) {
        Double_t deltaR = MathUtils::DeltaR(CleanElectrons->At(j)->Phi(),
                                            CleanElectrons->At(j)->Eta(), 
                                            jet->Phi(), jet->Eta());  
        if (deltaR < fMinDeltaRToElectron) {
          isElectronOverlap = kTRUE;
          break;	 	 
        }      
      }
    }

    if (isElectronOverlap) continue;

    // check for overlap with an Muon
    Bool_t isMuonOverlap = kFALSE;
    if (CleanMuons) {
      UInt_t n = CleanMuons->GetEntries();
      for (UInt_t j=0; j<n; ++j) {
        Double_t deltaR = MathUtils::DeltaR(CleanMuons->At(j)->Mom(),jet->Mom());  
        if (deltaR < fMinDeltaRToMuon) {
          isMuonOverlap = kTRUE;
          break;	 	 
        }      
      }
    }

    if (isMuonOverlap) continue;

    // check for overlap with a photon
    Bool_t isPhotonOverlap = kFALSE;
    if (fApplyPhotonRemoval && CleanPhotons) {
      UInt_t n = CleanPhotons->GetEntries();
      for (UInt_t j=0; j<n; ++j) {
        Double_t deltaR = MathUtils::DeltaR(CleanPhotons->At(j)->SCluster()->Phi(), 
                                            CleanPhotons->At(j)->SCluster()->Eta(),
                                            jet->Phi(), jet->Eta());  
        if (deltaR < fMinDeltaRToPhoton) {
          isPhotonOverlap = kTRUE;
          break;	 	 
        }      
      }
    }

    if (isPhotonOverlap)
      continue;

    // check for overlap with a tau
    Bool_t isTauOverlap = kFALSE;
    if (fApplyTauRemoval && CleanTaus) {
      UInt_t n = CleanTaus->GetEntries();
      for (UInt_t j=0; j<n; ++j) {
        Double_t deltaR = MathUtils::DeltaR(CleanTaus->At(j)->Mom(),jet->Mom());  
        if (deltaR < fMinDeltaRToTau) {
          isTauOverlap = kTRUE;
          break;	 	 
        }      
      }
    }

    if (isTauOverlap)
      continue;

    fCleanJets->Add(jet);
  }

  // sort according to pt
  fCleanJets->Sort();

}

void
JetCleaningMod::SlaveBegin () {
  PublishObj(fCleanJets);
}

void 
JetCleaningMod::SlaveEnd () {
  RetractObj(fCleanJets->GetName());
}

