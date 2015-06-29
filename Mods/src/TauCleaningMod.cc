#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Mods/interface/TauCleaningMod.h"

using namespace mithep;

ClassImp(mithep::TauCleaningMod)

//--------------------------------------------------------------------------------------------------
TauCleaningMod::TauCleaningMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fCleanElectronsName(ModNames::gkCleanElectronsName),        
  fCleanMuonsName(ModNames::gkCleanMuonsName),        
  fGoodTausName(ModNames::gkGoodTausName),        
  fMinDeltaRToElectron(0.3),
  fMinDeltaRToMuon(0.3),
  fCleanCaloTaus(new CaloTauOArr())
{
  fCleanCaloTaus->SetName(ModNames::gkCleanTausName);
}

TauCleaningMod::~TauCleaningMod()
{
  delete fCleanCaloTaus;
}

//--------------------------------------------------------------------------------------------------
void TauCleaningMod::Process()
{
  // Process entries of the tree.

  fCleanCaloTaus->Reset();

  // get input collections
  auto* GoodTaus = GetObject<CaloTauCol>(fGoodTausName);
  auto* CleanElectrons = GetObject<ElectronCol>(fCleanElectronsName);
  auto* CleanMuons = GetObject<MuonCol>(fCleanMuonsName);

  // remove any Tau that overlaps in eta, phi with an isolated electron.
  UInt_t n = GoodTaus->GetEntries();
  for (UInt_t i=0; i<n; ++i) {
    const CaloTau *tau = GoodTaus->At(i);        

    Bool_t isElectronOverlap =  false;
     
    // check for overlap with an electron
    if (CleanElectrons) {
      UInt_t n1 = CleanElectrons->GetEntries();
      for (UInt_t j=0; j<n1; j++) {
        Double_t deltaR = MathUtils::DeltaR(CleanElectrons->At(j)->Mom(),
	                                    tau->SourceCaloJet()->Mom());  
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
	                                    tau->SourceCaloJet()->Mom());  
        if (deltaR < fMinDeltaRToMuon) {
          isMuonOverlap = kTRUE;
          break;	 	 
        }      
      }
    }

    if (isMuonOverlap)
      continue;

    fCleanCaloTaus->Add(tau);     
  }

  // sort according to pt
  fCleanCaloTaus->Sort();
}

void
TauCleaningMod::SlaveBegin()
{
  PublishObj(fCleanCaloTaus);
}

void
TauCleaningMod::SlaveTerminate()
{
  RetractObj(fCleanCaloTaus->GetName());
}
