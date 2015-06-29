#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Mods/interface/PhotonCleaningMod.h"

using namespace mithep;

ClassImp(mithep::PhotonCleaningMod)

//--------------------------------------------------------------------------------------------------
PhotonCleaningMod::PhotonCleaningMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fCleanElectronsName(ModNames::gkCleanElectronsName),        
  fGoodPhotonsName(ModNames::gkGoodPhotonsName),        
  fCleanPhotons(new PhotonOArr),
  fMinDeltaRToElectron(0.3)
{
  // Constructor.
  fCleanPhotons->SetName(ModNames::gkCleanPhotonsName);
}

PhotonCleaningMod::~PhotonCleaningMod()
{
  delete fCleanPhotons;
}

//--------------------------------------------------------------------------------------------------
void PhotonCleaningMod::Process()
{
  // Process entries of the tree.

  fCleanPhotons->Reset();

  // get input collections
  const PhotonCol   *GoodPhotons    = GetObject<PhotonCol>(fGoodPhotonsName);
  const ElectronCol *CleanElectrons = GetObject<ElectronCol>(fCleanElectronsName);

  // remove any photon that overlaps in eta, phi with an isolated electron.
  UInt_t n = GoodPhotons->GetEntries();
  for (UInt_t i=0; i<n; ++i) {
    const Photon *ph = GoodPhotons->At(i);        

    Bool_t isElectronOverlap =  false;
     
    // check for overlap with an electron
    if (CleanElectrons) {
      UInt_t n2 = CleanElectrons->GetEntries();
      for (UInt_t j=0; j<n2; j++) {
        Double_t deltaR = MathUtils::DeltaR(CleanElectrons->At(j)->SCluster()->Phi(), 
                                            CleanElectrons->At(j)->SCluster()->Eta(),
                                             ph->SCluster()->Phi(),ph->SCluster()->Eta());  
        if (deltaR < fMinDeltaRToElectron) {
          isElectronOverlap = kTRUE;
          break;	 	 
        }      
      }
    }

    if (isElectronOverlap)
      continue;

    fCleanPhotons->Add(ph);     
  }

  // sort according to pt
  fCleanPhotons->Sort();

}

void
PhotonCleaningMod::SlaveBegin ()
{
  PublishObj(fCleanPhotons);
}

void 
PhotonCleaningMod::SlaveEnd ()
{
  RetractObj(fCleanPhotons->GetName());
}


