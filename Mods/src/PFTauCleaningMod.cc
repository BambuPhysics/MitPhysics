#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Mods/interface/PFTauCleaningMod.h"

using namespace mithep;

ClassImp(mithep::PFTauCleaningMod)

//--------------------------------------------------------------------------------------------------
PFTauCleaningMod::PFTauCleaningMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fCleanElectronsName(ModNames::gkCleanElectronsName),        
  fCleanMuonsName(ModNames::gkCleanMuonsName),        
  fGoodPFTausName(ModNames::gkGoodPFTausName),        
  fMinDeltaRToElectron(0.3),
  fMinDeltaRToMuon(0.3)
{
  // Constructor.
  fCleanPFTaus = new PFTauOArr;
  SetOutputName(ModNames::gkCleanPFTausName);
}

//--------------------------------------------------------------------------------------------------
void PFTauCleaningMod::Process()
{
  // Process entries of the tree.

  // get input collections
  const PFTauCol *GoodPFTaus = GetObject<PFTauCol>(fGoodPFTausName);
  const ElectronCol *CleanElectrons = GetObject<ElectronCol>(fCleanElectronsName);
  const MuonCol *CleanMuons = GetObject<MuonCol>(fCleanMuonsName);


  // remove any Tau that overlaps in eta, phi with an isolated electron.
  UInt_t n = GoodPFTaus->GetEntries();
  for (UInt_t i=0; i<n; ++i) {
    const PFTau *tau = GoodPFTaus->At(i);        

    Bool_t isElectronOverlap =  false;
     
    // check for overlap with an electron
    if (CleanElectrons) {
      UInt_t n1 = CleanElectrons->GetEntries();
      for (UInt_t j=0; j<n1; j++) {
        Double_t deltaR = MathUtils::DeltaR(CleanElectrons->At(j)->Mom(),
	                                    tau->Mom());  
        if (deltaR < fMinDeltaRToElectron) {
          isElectronOverlap = kTRUE;
          break;	 	 
        }      
      }
    }

    if (isElectronOverlap)
      continue;

    Bool_t isMuonOverlap =  false;
     
    // check for overlap with an Muon
    if (CleanMuons) {
      UInt_t n2 = CleanMuons->GetEntries();
      for (UInt_t j=0; j<n2; j++) {
        Double_t deltaR = MathUtils::DeltaR(CleanMuons->At(j)->Mom(),
	                                    tau->Mom());  
        if (deltaR < fMinDeltaRToMuon) {
          isMuonOverlap = kTRUE;
          break;	 	 
        }      
      }
    }

    if (isMuonOverlap)
      continue;

    fCleanPFTaus->Add(tau);     
  }

  // sort according to pt
 fCleanPFTaus->Sort();

}

void
PFTauCleaningMod::SlaveBegin () {
  PublishObj(fCleanPFTaus);
}

void 
PFTauCleaningMod::SlaveEnd () {
	RetractObj(fCleanPFTaus->GetName());
}
