#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataTree/interface/CaloJetCol.h"
#include "MitAna/DataTree/interface/CaloMetCol.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Mods/interface/CaloMetCorrectionMod.h"

using namespace mithep;

ClassImp(mithep::CaloMetCorrectionMod)

//--------------------------------------------------------------------------------------------------
CaloMetCorrectionMod::CaloMetCorrectionMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fMetName("CorMuonMet"),
  fCorrectedMetName("CorrectedMet"),
  fCorrectedJetsName(ModNames::gkCorrectedJetsName),  
  fMinJetPt(20.0),
  fMaxJetEMF(0.9)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void CaloMetCorrectionMod::Process()
{
  // Process entries of the tree. 

  const CaloMetCol *inMets = GetObject<CaloMetCol>(fMetName);
  if (!inMets) {
    SendError(kAbortModule, "Process", 
              "Pointer to input met collection %s is null.",
              fMetName.Data());
    return;
  }

  const CaloJetCol *inJets = GetObject<CaloJetCol>(fCorrectedJetsName);
  if (!inJets) {
    SendError(kAbortModule, "Process", 
              "Pointer to input jet collection %s is null.",
              fCorrectedJetsName.Data());
    return;
  }

  CaloMetOArr *CorrectedMet = new CaloMetOArr;
  CorrectedMet->SetOwner(kTRUE);
  CorrectedMet->SetName(fCorrectedMetName);

  Double_t deltaPx = 0.0;
  Double_t deltaPy = 0.0;
  Double_t deltaSumEt = 0.0;

  // loop over jets
  for (UInt_t i=0; i<inJets->GetEntries(); ++i) {
    const CaloJet *inJet = inJets->At(i);

    //cache uncorrected momentum
    const FourVectorM rawMom = inJet->RawMom();

    //check if jet passes cuts
    if (rawMom.Pt()>fMinJetPt && inJet->EnergyFractionEm()<fMaxJetEMF) {
      //get jet correction
      Double_t correction = inJet->CombinedCorrectionScale() - 1.0;
    
      //compute met correction factors
      deltaPx -= rawMom.Px()*correction;
      deltaPy -= rawMom.Py()*correction;
      deltaSumEt += rawMom.Et()*correction;
    }
          
  }

  //loop over met (should normally be only one)
  for (UInt_t i=0; i<inMets->GetEntries(); ++i) {
    const CaloMet *inMet = inMets->At(i);

    //copy input met
    CaloMet *met = new CaloMet(*inMet);

    //apply corrections computed above
    met->PushCorrectionX(deltaPx);
    met->PushCorrectionY(deltaPy);
    met->PushCorrectionSumEt(deltaSumEt);

    met->SetMex(inMet->Mex()+deltaPx);
    met->SetMey(inMet->Mey()+deltaPy);
    met->SetSumEt(inMet->SumEt()+deltaSumEt);

    if (0) {
      UInt_t corSize = inMet->Dmex().GetEntries();
      if (corSize>0) {
        printf("inMet Dmex = %5f, outMet Dmex= %5f\n",inMet->Dmex().At(corSize-1),deltaPx);
        printf("inMet Dmey = %5f, outMet Dmey= %5f\n",inMet->Dmey().At(corSize-1),deltaPy);
        printf("inMet DSumet = %5f, outMet DSumet= %5f\n",inMet->DSumEt().At(corSize-1),deltaSumEt);
      }
      printf("inMet Pt = %5f, outMet Pt = %5f\n",inMet->Pt(),met->Pt());
      printf("inMet SumEt = %5f, outMet SumEt = %5f\n",inMet->SumEt(),met->SumEt());
    }

    // add corrected met to collection
    CorrectedMet->AddOwned(met);             
  }
  
  // add to event for other modules to use
  AddObjThisEvt(CorrectedMet);
}


