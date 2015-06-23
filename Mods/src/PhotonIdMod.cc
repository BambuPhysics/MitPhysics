#include "MitPhysics/Mods/interface/PhotonIdMod.h"
#include "MitPhysics/Utils/interface/PhotonTools.h" 
#include "MitPhysics/Utils/interface/IsolationTools.h" 
#include "MitAna/DataTree/interface/StableData.h" 
#include "MitAna/DataTree/interface/TriggerObjectCol.h"  
#include "MitPhysics/Init/interface/ModNames.h"   
#include "MitPhysics/Init/interface/Constants.h"

using namespace mithep;
ClassImp(mithep::PhotonIdMod)

template<class T> 
void
PhotonIdMod::GetAuxInput(PhotonIdMod::AuxInput inputCol, TObject const** aux)
{

  aux[inputCol] = GetObject<T>(fAuxInputNames[inputCol], true);  
  if (!aux[inputCol])       
    SendError(kAbortAnalysis, "GetAuxInput", "Could not retrieve auxiliary input " + fAuxInputNames[inputCol]);    
}


mithep::PhotonIdMod::PhotonIdMod(char const* name/* = "PhotonIdMod"*/, char const* title/* = "Photon Identification"*/) :
  IdMod(name, title)
{
  fOutput = new PhotonOArr(32, TString(name) + "Output");
  fAuxInputNames[kPileupEnergyDensity] = Names::gkPileupEnergyDensityBrn; 

}

mithep::PhotonIdMod::~PhotonIdMod()
{
}

void
mithep::PhotonIdMod::Process()
{

  // Process entries of the tree.  
  auto* photons = GetObject<mithep::PhotonCol>(fInputName, true); 
  if (!photons) 
    SendError(kAbortAnalysis, "Process", "Photons not found");  

  TObject const* aux[nAuxInputs] = {};

  //get trigger object collection if trigger matching is enabled  
  if (fApplyTriggerMatching) {       
    aux[kTrigObjects] = GetHLTObjects(fAuxInputNames[kTrigObjects]);
    if (!aux[kTrigObjects])			
      SendError(kAbortAnalysis, "Process", "TrigObjects not found");   
  }

  mithep::PhotonOArr* goodPhotons = 0; 
  if (fIsFilterMode) { 
    goodPhotons = static_cast<mithep::PhotonOArr*>(fOutput);  
    goodPhotons->Reset();  
  }

  else{
    fFlags.Resize(photons->GetEntries());
    for (unsigned iP = 0; iP != photons->GetEntries(); ++iP) 
      fFlags.At(iP) = false; 
  }
 
  for (UInt_t iP = 0; iP != photons->GetEntries(); ++iP) { 
    Photon const& photon = *photons->At(iP);
    fCutFlow->Fill(cAll); 
  
    if (photon.Pt() < fPtMin)
      continue; 
    fCutFlow->Fill(cPt); 

    if (photon.AbsEta() > fEtaMax)
      continue;
  
    fCutFlow->Fill(cEta); 
    
    //apply trigger matching      
    if (fApplyTriggerMatching &&
	!PhotonTools::PassTriggerMatching(&photon, static_cast<mithep::TriggerObjectCol const*>(aux[kTrigObjects])))
      continue; 
    fCutFlow->Fill(cTriggerMatching);   

    //apply id cut 
    if (!PassIdCut(photon, aux))
      continue;         
    fCutFlow->Fill(cId); 

    //apply Isolation Cut 
    if (!PassIsolationCut(photon, aux))
      continue;
    fCutFlow->Fill(cIsolation); 

    if (fIsFilterMode)  
      goodPhotons->Add(&photon);
    else
      fFlags.At(iP) = true;
  }

  if(fIsFilterMode){
    // sort according to pt
    goodPhotons->Sort();
  }


}

void
mithep::PhotonIdMod::IdBegin()
{

  fCutFlow->SetBins(nCuts, 0., double(nCuts)); 
  TAxis* xaxis = fCutFlow->GetXaxis(); 
  xaxis->SetBinLabel(cAll + 1, "All");
  xaxis->SetBinLabel(cPt + 1,  "Pt");
  xaxis->SetBinLabel(cEta + 1, "Eta");
  xaxis->SetBinLabel(cTriggerMatching + 1, "TriggerMatching");
  xaxis->SetBinLabel(cId + 1,  "Id"); 
  xaxis->SetBinLabel(cIsolation + 1, "Isolation");

}

void
mithep::PhotonIdMod::IdTerminate()
{
}

//--------------------------------------------------------------------------------------------------
Bool_t 
mithep::PhotonIdMod::PassIdCut(Photon const& pho, TObject const**)
{
  switch(fIdType){

  case PhotonTools::kRun2Tight:
  case PhotonTools::kRun2Medium:
  case PhotonTools::kRun2Loose:
    return PhotonTools::PassID(&pho, PhotonTools::EPhIdType(fIdType));

  default:
    return false;
  }
}

//--------------------------------------------------------------------------------------------------
Bool_t
mithep::PhotonIdMod::PassIsolationCut(Photon const& pho, TObject const** aux)
{
  switch (fIsoType) { 
    
  case PhotonTools::kRun2LooseIso:
  case PhotonTools::kRun2MediumIso:
  case PhotonTools::kRun2TightIso:
    if (!aux[kPileupEnergyDensity])  
      GetAuxInput<mithep::PileupEnergyDensityCol>(kPileupEnergyDensity, aux);  
    return PhotonTools::PassIsoRhoCorr(&pho, PhotonTools::EPhIsoType(fIsoType),
				       static_cast<mithep::PileupEnergyDensityCol const*>(aux[kPileupEnergyDensity])->At(0)->Rho(fRhoAlgo));

  case PhotonTools::kNoIso:
    return true;

  default:
    return false;
  }
}
