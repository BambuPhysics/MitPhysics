#include <TH1D.h>
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataCont/interface/BaseCollection.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Mods/interface/EffMod.h"

using namespace mithep;

ClassImp(mithep::EffMod)

//--------------------------------------------------------------------------------------------------
EffMod::EffMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fCol1Name(ModNames::gkMCLeptonsName), 
  fCol2Name(ModNames::gkMergedLeptonsName),
  fMinPt(10),
  fMaxEta(2.4),
  fRadius(0.1),
  fPartType(MCParticle::kUnknown)
{
  // Constructor.
}


//--------------------------------------------------------------------------------------------------
void EffMod::Process()
{
  // Process entries of the tree.

  const BaseCollection *col1 = GetObject<BaseCollection>(fCol1Name);
  const BaseCollection *col2 = GetObject<BaseCollection>(fCol2Name);

  UInt_t ents1 = 0;
  if (col1)
    ents1 = col1->GetEntries();

  UInt_t ents2 = 0;
  if (col2)
    ents2 = col2->GetEntries();

  Bool_t *found = new Bool_t[ents2];
  for (UInt_t i=0; i<ents2; ++i) 
    found[i] = kFALSE;
  
  // find matches
  for (UInt_t j=0; j<ents1; ++j) {
    const Particle *p1 = dynamic_cast<const Particle*>(col1->ObjAt(j));
    if (!p1)
      continue;
    
    if (fPartType != MCParticle::kUnknown) {
      const MCParticle *mc = dynamic_cast<const MCParticle*>(p1);
      if (!mc || !mc->Is(fPartType))
        continue;
    }

    Double_t pt = p1->Pt();
    if (pt < fMinPt) 
      continue;
    if (p1->AbsEta() > fMaxEta) 
      continue;

    Double_t phi = p1->Phi();
    Double_t eta = p1->Eta();

    fCol1Pt->Fill(TMath::Min(pt, 299.999));
    fCol1Eta->Fill(eta);

    UInt_t foundInd  = ents2;
    Double_t foundD  = 1e12;

    for (UInt_t i=0; i<ents2; ++i) {
      if (found[i]) continue;
      const Particle *p2 = dynamic_cast<const Particle*>(col2->ObjAt(i));
      if (!p2)
        continue;

      if(MathUtils::DeltaR(phi, eta, p2->Phi(), p2->Eta()) < fRadius) {
        Double_t newDiff = TMath::Abs(pt - p2->Pt());
        if (newDiff < foundD) {
          foundInd = i;
          foundD = newDiff;
        } 
      }
    }
    if (foundInd < ents2) {
      fCol2Pt->Fill(TMath::Min(p1->Pt(),299.999));
      fCol2Eta->Fill(p1->Eta());
      found[foundInd] = 1;
    }
  }

  // find fakes
  for (UInt_t i=0; i<ents2; ++i) {
      const Particle *p2 = dynamic_cast<const Particle*>(col2->ObjAt(i));
      if (!p2)
        continue;
      if (found[i]) {
        fGoodPt->Fill(TMath::Min(p2->Pt(),299.999));
        fGoodEta->Fill(p2->Eta());
      } else {
        fFakePt->Fill(TMath::Min(p2->Pt(),299.999));
        fFakeEta->Fill(p2->Eta());
      }
  }        

}

//--------------------------------------------------------------------------------------------------
void EffMod::SlaveBegin()
{
  // Create and add histograms to the output list.

  AddTH1(fCol1Pt,"hCol1Pt",";p_{t} [GeV];#",100,0,300);
  AddTH1(fCol1Eta,"hCol1Eta",";#eta;#",50,-5,5);
  AddTH1(fCol2Pt,"hCol2Pt",";p_{t} [GeV];#",100,0,300);
  AddTH1(fCol2Eta,"hCol2Eta",";#eta;#",50,-5,5);
  AddTH1(fGoodPt,"hGoodPt",";p_{t} [GeV];#",100,0,300);
  AddTH1(fGoodEta,"hGoodEta",";#eta;#",50,-5,5);
  AddTH1(fFakePt,"hFakePt",";p_{t} [GeV];#",100,0,300);
  AddTH1(fFakeEta,"hFakeEta",";#eta;#",50,-5,5);
}

//--------------------------------------------------------------------------------------------------
void EffMod::Terminate()
{
  // Create and add ratio histograms to the output list.

  if (gROOT->IsBatch()) {
    TH1D *pteff = static_cast<TH1D*>(fCol2Pt->Clone("hPtEff"));
    pteff->Divide(fCol1Pt);
    AddOutput(pteff);

    TH1D *etaeff = static_cast<TH1D*>(fCol2Eta->Clone("hEtaEff"));
    etaeff->Divide(fCol1Eta);
    AddOutput(etaeff);

    TH1D *ptfrate = static_cast<TH1D*>(fFakePt->Clone("hPtFakeRate"));
    ptfrate->Divide(fGoodPt);
    AddOutput(ptfrate);

    TH1D *etafrate = static_cast<TH1D*>(fFakeEta->Clone("hEtaFakeRate"));
    etafrate->Divide(fGoodEta);
    AddOutput(etafrate);
  }
}
