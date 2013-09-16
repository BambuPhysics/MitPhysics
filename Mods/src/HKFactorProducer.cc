// $Id: HKFactorProducer.cc,v 1.18 2012/09/24 12:23:15 ceballos Exp $

#include "MitPhysics/Mods/interface/HKFactorProducer.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataTree/interface/MCParticleCol.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TParameter.h>
#include <TTree.h>
#include <TFile.h>

using namespace mithep;

ClassImp(mithep::HKFactorProducer)

//--------------------------------------------------------------------------------------------------
HKFactorProducer::HKFactorProducer(const char *name, const char *title) : 
  BaseMod(name,title),
  fProcessID(102),
  fInputFileName("setme"),
  fMCBosonsName(ModNames::gkMCBosonsName),
  fMCEvInfoName(Names::gkMCEvtInfoBrn),
  fEmbedWeightName(Names::gkEmbedWeightBrn),
  fIsData(kFALSE),
  fMakePDFNtuple(kFALSE),
  fDoHiggsMhReweighting(kFALSE),
  fMh(0),
  fWidth(0),
  fBWflag(0),
  fMCEventInfo(0),
  fEmbedWeight(0),
  fOutputFile(0),
  fOutputName("ntuple.root")
{
  // Constructor
}

//--------------------------------------------------------------------------------------------------
HKFactorProducer::~HKFactorProducer()
{
}

//--------------------------------------------------------------------------------------------------
void HKFactorProducer::Process()
{
  // Process entries of the tree.

  Double_t theWeight = 1.0;

  if(fIsData == kFALSE){
    // get the bosons
    MCParticleCol *GenBosons = GetObjThisEvt<MCParticleCol>(fMCBosonsName);

    LoadBranch(fMCEvInfoName);

    // only accept the exact process id
    if (fMh > 0 && fWidth > 0 && fBWflag >=0) { 

      Double_t mH  = -1.0;
      for (UInt_t i=0; i<GenBosons->GetEntries(); ++i) {
	if(GenBosons->At(i)->PdgId() == MCParticle::kH) {
	  mH  = GenBosons->At(i)->Mass();
    	  break;
	}
      }

      if(mH >= 0) {
	Double_t MTop = 172.5;
	theWeight = fWeightAlgo.getweight(fMh,fWidth,MTop,mH,fBWflag);

	if (theWeight > 3.0) {
          cout << "MHweights: " << fMh << " " << fWidth  << " " << MTop      << " " 
	                        << mH  << " " << fBWflag << " " << theWeight << endl;
	}

	if (GetFillHist()) {
          hDHKFactor[0]->Fill(0.5);
          hDHKFactor[1]->Fill(TMath::Min(mH,1999.999));
          hDHKFactor[2]->Fill(TMath::Min(mH,1999.999),theWeight);
          hDHKFactor[3]->Fill(TMath::Max(TMath::Min(theWeight,5.999),-3.999));
	}
      }
    }
    else if (fProcessID == 999){ // for MCatNLO we care about positive or negative weights
      theWeight = fMCEventInfo->Weight();
      if (GetFillHist()) hDHKFactor[3]->Fill(TMath::Max(TMath::Min(theWeight,5.999),-3.999));
      theWeight = theWeight/fabs(theWeight);
      if (GetFillHist()) hDHKFactor[0]->Fill(0.5,theWeight);
    }
    else if (fProcessID == 998){ // for other samples we care about the actual weights
      theWeight = fMCEventInfo->Weight();
      if (GetFillHist()) hDHKFactor[3]->Fill(TMath::Max(TMath::Min(theWeight,5.999),-3.999));
      if (GetFillHist()) hDHKFactor[0]->Fill(0.5,theWeight);
    }
    else {
      if (GetFillHist()) hDHKFactor[0]->Fill(0.5,1.0);
    }
    // process id distribution
    if (GetFillHist()) hDHKFactor[4]->Fill(TMath::Min((Double_t)fMCEventInfo->ProcessId(),999.499));

    if (fMakePDFNtuple == kTRUE){
      fTreeVariables[0] = fMCEventInfo->Weight();
      fTreeVariables[1] = fMCEventInfo->Scale();
      fTreeVariables[2] = fMCEventInfo->Id1();
      fTreeVariables[3] = fMCEventInfo->X1();
      fTreeVariables[4] = fMCEventInfo->Pdf1();
      fTreeVariables[5] = fMCEventInfo->Id2();
      fTreeVariables[6] = fMCEventInfo->X2();
      fTreeVariables[7] = fMCEventInfo->Pdf2();
      fTree->Fill();
    }
  }
  else if (fProcessID == 997){ // for tau embedding samples
    LoadBranch(fEmbedWeightName);
    theWeight = fEmbedWeight->At(fEmbedWeight->GetEntries()-1)->GenWeight();
    if (GetFillHist()) hDHKFactor[3]->Fill(TMath::Max(TMath::Min(theWeight,5.999),-3.999));
    if (GetFillHist()) hDHKFactor[0]->Fill(0.5,theWeight);
  }

  TParameter<Double_t> *NNLOWeight = new TParameter<Double_t>("NNLOWeight", theWeight);
  AddObjThisEvt(NNLOWeight);
}

//--------------------------------------------------------------------------------------------------
void HKFactorProducer::SlaveBegin()
{
  // Book branch and histograms if wanted.

  if(fIsData == kFALSE){
    ReqBranch(fMCEvInfoName, fMCEventInfo);
  }
  if(fProcessID == 997){
    ReqBranch(fEmbedWeightName, fEmbedWeight);
  }

  if (GetFillHist()) {
    char sb[1024];
    sprintf(sb,"hDHKFactor_%d", 0);  hDHKFactor[0]  = new TH1D(sb,sb,1,0,1); 
    sprintf(sb,"hDHKFactor_%d", 1);  hDHKFactor[1]  = new TH1D(sb,sb,1000,0.,2000.); 
    sprintf(sb,"hDHKFactor_%d", 2);  hDHKFactor[2]  = new TH1D(sb,sb,1000,0.,2000.); 
    sprintf(sb,"hDHKFactor_%d", 3);  hDHKFactor[3]  = new TH1D(sb,sb,500,-4.0,6.0); 
    sprintf(sb,"hDHKFactor_%d", 4);  hDHKFactor[4]  = new TH1D(sb,sb,1000,-0.5,999.5); 
    for(Int_t i=0; i<5; i++) AddOutput(hDHKFactor[i]);
  }

  //***********************************************************************************************
  //Create Ntuple Tree  
  //***********************************************************************************************
  if (fMakePDFNtuple == kTRUE){
    printf("... init PDF ntuple ...\n");
    fOutputFile = new TFile(fOutputName, "RECREATE");
    fTree = new TTree("PDFTree", "PDFTree");
    const char* TreeFormat;
    TreeFormat = "weight/F:lq:lid1:lx1:lpdf1:lid2:lx2:lpdf2";
    fTree->Branch("H", &fTreeVariables,TreeFormat);    
  }
}

//--------------------------------------------------------------------------------------------------
void HKFactorProducer::SlaveTerminate(){
  if (fMakePDFNtuple == kTRUE){
    fOutputFile->cd();
    fTree->Write();
    fOutputFile->Close();
 }
}
