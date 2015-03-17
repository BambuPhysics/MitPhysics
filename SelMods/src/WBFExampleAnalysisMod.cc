#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TParameter.h>
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/DataTree/interface/Names.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitAna/DataCont/interface/ObjArray.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataTree/interface/ParticleCol.h"
#include "MitPhysics/SelMods/interface/WBFExampleAnalysisMod.h"

using namespace mithep;
ClassImp(mithep::WBFExampleAnalysisMod)

//--------------------------------------------------------------------------------------------------
WBFExampleAnalysisMod::WBFExampleAnalysisMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fCleanJetsName("NoDefaultNameSet"),
  fJetPtMax(30),
  fJetPtMin(30),
  fDeltaEtaMin(4),
  fDiJetMassMin(500)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void WBFExampleAnalysisMod::Begin()
{
  // Run startup code on the client machine. For this module, we dont do
  // anything here.
}

//--------------------------------------------------------------------------------------------------
void WBFExampleAnalysisMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we typically initialize histograms and other analysis objects and request
  // branches. For this module, we request a branch of the MitTree.

  //Create your histograms here

  //*************************************************************************************************
  // Selection Histograms
  //*************************************************************************************************
  AddTH1(fWBFSelection,"fWBFSelection", ";Cut Number;Number of Events", 7, -1.5, 5.5);
  AddTH1(fWBFDeltaPhi, "fWBFDeltaPhi",  ";Delta Phi;Number of Events",   90,0.,180  );
 
  //***********************************************************************************************
  // N-1 Histograms
  //***********************************************************************************************
  //All events
  AddTH1(fWBFPtJetMax_NMinusOne, "fWBFPtJetMax_NMinusOne", ";Pt Jet Max;Number of Events",300,0.,300. );
  AddTH1(fWBFPtJetMin_NMinusOne, "fWBFPtJetMin_NMinusOne", ";Pt Jet Min;Number of Events",300,0.,300. );
  AddTH1(fWBFEta12_NMinusOne,    "fWBFEta12_NMinusOne",    ";eta1*eta2;Number of Events", 2,-1.5,1.5  );
  AddTH1(fWBFdeltaEta_NMinusOne, "fWBFdeltaEta_NMinusOne", ";Delta Eta;Number of Events", 100,0.,10.  );
  AddTH1(fWBFdijetMass_NMinusOne,"fWBFdijetMass_NMinusOne",";DiJet Mass;Number of Events",300,0.,3000.);
  AddTH1(fWBFZVar_NMinusOne,     "fWBFZVar_NMinusOne",     ";Z variable;Number of Events",50,0.,5.    );

  //***********************************************************************************************
  // After all cuts Histograms
  //***********************************************************************************************
  AddTH1(fWBFSSMass_afterCuts,"fWBFSSMass_afterCuts",";SS dilepton Mass;Number of Events",200,0.,400);
  AddTH1(fWBFSSDeltaPhi_afterCuts,"fWBFSSDeltaPhi_afterCuts",";SS DeltaPhi 2l;Number of Events",90,0.,180);
  AddTH1(fWBFOSMass_afterCuts,"fWBFOSMass_afterCuts",";OS dilepton Mass;Number of Events",200,0.,400);
  AddTH1(fWBFOSDeltaPhi_afterCuts,"fWBFOSDeltaPhi_afterCuts",";OS DeltaPhi 2l;Number of Events",90,0.,180);

}

//--------------------------------------------------------------------------------------------------
void WBFExampleAnalysisMod::Process()
{
  //Obtain all the good objects from the event cleaning module
  ParticleOArr *CleanLeptons = dynamic_cast<mithep::ParticleOArr*>
     (FindObjThisEvt(ModNames::gkMergedLeptonsName));
  ObjArray<Jet> *CleanJets = dynamic_cast<ObjArray<Jet>* >
    (FindObjThisEvt(fCleanJetsName.Data()));
  TParameter<Double_t> *NNLOWeight = GetObjThisEvt<TParameter<Double_t> >("NNLOWeight");

  if (CleanJets->GetEntries() < 2)
    return;

  //*********************************************************************************************
  //Define Cuts
  //*********************************************************************************************
  const int nCuts = 6;
  bool passCut[nCuts] = {false, false, false, false, false, false};
  
  // ptjet max cut
  if (CleanJets->At(0)->Pt() > fJetPtMax)
    passCut[0] = true;

  // ptjet min cut
  if (CleanJets->At(1)->Pt() > fJetPtMin)
    passCut[1] = true;
  
  // this cut is always applied
  if (CleanJets->At(0)->Eta()*CleanJets->At(1)->Eta() < 0)
    passCut[2] = true;
  
  // deltaEta cut
  double deltaEta = TMath::Abs(CleanJets->At(0)->Eta()-CleanJets->At(1)->Eta());
  if (deltaEta > fDeltaEtaMin)
    passCut[3] = true;
  
  // dijet mass cut
  CompositeParticle dijet;
  dijet.AddDaughter(CleanJets->At(0));
  dijet.AddDaughter(CleanJets->At(1));

  if (dijet.Mass() > fDiJetMassMin)
    passCut[4] = true;
  
  // jet veto cut, use of Zeppenfeld variable
  passCut[5] = true;
  double zVarMin = 30.;
  if (CleanJets->GetEntries() >=2) {
    for (UInt_t i=2; i<CleanJets->GetEntries(); i++) {
      if (CleanJets->At(i)->Pt() <= 30.0)
	return;
      double zVar = TMath::Abs(CleanJets->At(i)->Eta()-(CleanJets->At(0)->Eta()+CleanJets->At(1)->Eta())/2)/
                    TMath::Abs(CleanJets->At(0)->Eta()-CleanJets->At(1)->Eta());
      if (zVar < zVarMin)
	zVarMin = zVar;
    }
    if (zVarMin < 1)
      passCut[5] = false;
  }

  //*********************************************************************************************
  //Make Selection Histograms. Number of events passing each level of cut
  //*********************************************************************************************  
  bool passAllCuts = true;
  for(int c=0; c<nCuts; c++) passAllCuts = passAllCuts & passCut[c];
  
  if(passAllCuts){
    double deltaPhiJJ = MathUtils::DeltaPhi(CleanJets->At(0)->Phi(), 
                                            CleanJets->At(1)->Phi())* 180.0 / TMath::Pi();
    fWBFDeltaPhi->Fill(deltaPhiJJ,NNLOWeight->GetVal()); 
  }
  
  //Cut Selection Histograms
  fWBFSelection->Fill(-1,NNLOWeight->GetVal());
  for (int k=0;k<nCuts;k++) {
    bool pass = true;
    bool passPreviousCut = true;
    for (int p=0;p<=k;p++) {
      pass = (pass && passCut[p]);
      if (p<k)
        passPreviousCut = (passPreviousCut&& passCut[p]);
    }
    
    if (pass) {
      fWBFSelection->Fill(k,NNLOWeight->GetVal());
    }
  }
  
  //*********************************************************************************************
  // N-1 Histograms
  //*********************************************************************************************
  // N-1 ptjet max
  bool pass = true;
  for (int k=0;k<nCuts;k++) {
    if (k != 0) {
      pass = (pass && passCut[k]);      
    }
  }
  if (pass) {
    fWBFPtJetMax_NMinusOne->Fill(TMath::Min(CleanJets->At(0)->Pt(),299.999),NNLOWeight->GetVal());
  }     
  
  // N-1 ptjet min
  pass = true;
  for (int k=0;k<nCuts;k++) {
    if (k != 1) {
      pass = (pass && passCut[k]);      
    }
  }
  if (pass) {
    fWBFPtJetMin_NMinusOne->Fill(TMath::Min(CleanJets->At(1)->Pt(),299.999),NNLOWeight->GetVal());
  }     
  
  // N-1 eta1*eta2
  pass = true;
  for (int k=0;k<nCuts;k++) {
    if (k != 2) {
      pass = (pass && passCut[k]);      
    }
  }
  if (pass) {
    fWBFEta12_NMinusOne->Fill(CleanJets->At(0)->Eta()*CleanJets->At(1)->Eta()/
                   TMath::Abs(CleanJets->At(0)->Eta()*CleanJets->At(1)->Eta()),NNLOWeight->GetVal());
  }     
  
  // N-1 deltaEta
  pass = true;
  for (int k=0;k<nCuts;k++) {
    if (k != 3) {
      pass = (pass && passCut[k]);      
    }
  }
  if (pass) {
    fWBFdeltaEta_NMinusOne->Fill(TMath::Min(deltaEta,9.999),NNLOWeight->GetVal());
  }
  
  // N-1 dijet mass
  pass = true;
  for (int k=0;k<nCuts;k++) {
    if (k != 4) {
      pass = (pass && passCut[k]);      
    }
  }
  if (pass) {
    fWBFdijetMass_NMinusOne->Fill(TMath::Min(dijet.Mass(),2999.999),NNLOWeight->GetVal());
  }
  
  // N-1 dijet mass
  pass = true;
  for (int k=0;k<nCuts;k++) {
    if (k != 5) {
      pass = (pass && passCut[k]);      
    }
  }
  if (pass) {
    fWBFZVar_NMinusOne->Fill(TMath::Min(zVarMin,4.999),NNLOWeight->GetVal());
  }

  //*********************************************************************************************
  //Plots after all Cuts
  //*********************************************************************************************
  if (passAllCuts) {
    // Distributions for dileptons events
    if (CleanLeptons->GetEntries() >= 2 && 
        CleanLeptons->At(0)->Pt()  > 20 && CleanLeptons->At(1)->Pt() > 10){

      CompositeParticle dilepton;
      dilepton.AddDaughter(CleanLeptons->At(0));
      dilepton.AddDaughter(CleanLeptons->At(1));
      double deltaPhiLeptons = MathUtils::DeltaPhi(CleanLeptons->At(0)->Phi(), 
                                                   CleanLeptons->At(1)->Phi())* 180.0 / TMath::Pi();
      
      if (CleanLeptons->At(0)->Charge() * CleanLeptons->At(1)->Charge() > 0){  // same-sign
        fWBFSSMass_afterCuts->Fill(TMath::Min(dilepton.Mass(),399.999),NNLOWeight->GetVal());
        fWBFSSDeltaPhi_afterCuts->Fill(deltaPhiLeptons,NNLOWeight->GetVal());
      }
      else { // opposite-sign
        fWBFOSMass_afterCuts->Fill(TMath::Min(dilepton.Mass(),399.999),NNLOWeight->GetVal());
        fWBFOSDeltaPhi_afterCuts->Fill(deltaPhiLeptons,NNLOWeight->GetVal());
      }
    }
  }
  
  return;
}

//--------------------------------------------------------------------------------------------------
void WBFExampleAnalysisMod::SlaveTerminate()
{
  
  // Run finishing code on the computer (slave) that did the analysis. For this
  // module, we dont do anything here.

} 
//--------------------------------------------------------------------------------------------------
void WBFExampleAnalysisMod::Terminate()
{
  // Run finishing code on the client computer. For this module, we dont do
  // anything here.

}
