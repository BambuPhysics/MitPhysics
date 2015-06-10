#include "MitPhysics/Utils/interface/MVAMet.h"
#include "MitPhysics/Utils/interface/MetLeptonTools.h"
#include "MitPhysics/Utils/interface/JetTools.h"
#include "MitPhysics/Utils/interface/RecoilTools.h"
#include "MitAna/DataTree/interface/StableData.h"
#include "MitPhysics/Utils/interface/RecoilTools.h"
#include "MitPhysics/Utils/interface/MetLeptonTools.h"
#include <TFile.h>
#include <TRandom3.h>
#include <TSystem.h>
#include "CondFormats/EgammaObjects/interface/GBRForest.h"
// #include "Cintex/Cintex.h"
#include <utility>

ClassImp(mithep::MVAMet)

using namespace mithep;

//--------------------------------------------------------------------------------------------------
MVAMet::MVAMet() :
  fRecoilTools(0),
  fPhiMethodName     ("PhiCorrection"),
  fU1MethodName      ("U1Correction"),
  fCovU1MethodName ("CovU1"),
  fCovU2MethodName ("CovU2"),
  fIsInitialized(kFALSE),
  fU      (0),
  fUPhi   (0),
  fTKSumEt(0),
  fTKU    (0),
  fTKUPhi (0),
  fNPSumEt(0),
  fNPU    (0),
  fNPUPhi (0),
  fPUSumEt(0),
  fPUMet  (0),
  fPUMetPhi(0),
  fPCSumEt(0),
  fPCU    (0),
  fPCUPhi (0),
  fJSPt1  (0),
  fJSEta1 (0),
  fJSPhi1 (0),
  fJSPt2  (0),
  fJSEta2 (0),
  fJSPhi2 (0),
  fNJet   (0),
  fNAllJet(0),
  fNPV    (0),
  fPhiReader(0),
  fU1Reader(0),
  fCovU1Reader(0),
  fCovU2Reader(0),
  fMetLeptonTools(0) { }
//--------------------------------------------------------------------------------------------------
MVAMet::~MVAMet()
{
  delete fPhiReader;
  delete fU1Reader;
  delete fCovU1Reader;
  delete fCovU2Reader;

}
//--------------------------------------------------------------------------------------------------
void MVAMet::Initialize(TString iJetLowPtFile, 
                        TString iJetHighPtFile,
                        TString iJetCutFile,
                        TString iU1Weights, 
                        TString iPhiWeights, 
                        TString iCovU1Weights,
                        TString iCovU2Weights,
                        JetIDMVA::MVAType     iType,MVAMetType iMETType) { 

  if (iJetLowPtFile.Length() == 0 ||
      iJetHighPtFile.Length() == 0 ||
      iJetCutFile.Length() == 0 ||
      iU1Weights.Length() == 0 || 
      iPhiWeights.Length() == 0 || 
      iCovU1Weights.Length() == 0 ||
      iCovU2Weights.Length() == 0) {
    TString dataDir(gSystem->Getenv("MIT_DATA"));
    if (dataDir.Length() == 0)
      throw std::runtime_error("MIT_DATA environment is not set.");

    if (iJetLowPtFile.Length() == 0)
      iJetLowPtFile = dataDir + "/mva_RecoilPhiRegress_baseline.weights.xml";
    if (iJetHighPtFile.Length() == 0)
      iJetHighPtFile = dataDir + "/mva_RecoilPhiRegress_baseline.weights.xml";
    if (iJetCutFile.Length() == 0)
      iJetCutFile = dataDir + "/mva_RecoilPhiRegress_baseline.weights.xml";
    if (iU1Weights.Length() == 0)
      iU1Weights = dataDir + "/gbrmet.root";
    if (iPhiWeights.Length() == 0)
      iPhiWeights = dataDir + "/gbrmetphi.root";
    if (iCovU1Weights.Length() == 0)
      iCovU1Weights = dataDir + "/gbrcovu1_52.root";
    if (iCovU2Weights.Length() == 0)
      iCovU2Weights = dataDir + "/gbrcovu2_52.root";
  }

  fIsInitialized = kTRUE;
  fType          = iType;
  f42            = iU1Weights.Contains("42");
  fType1         = (iMETType == kUseType1 || iMETType == kUseType1Rho);
  fOld42         = (iMETType == kOld42);
  bool lUseRho   = true; if(iMETType == kUseType1Rho || iMETType == kUseRho) lUseRho = false; 
  fRecoilTools = new RecoilTools(iJetLowPtFile,iJetHighPtFile,iJetCutFile,(iMETType == kOld42),fType,lUseRho);
  if(f42) fRecoilTools->fJetIDMVA->fJetPtMin = 1.;

  // Is this necessary? (Y.I. Apr 2 2015)
  //  ROOT::Cintex::Cintex::Enable();   
  
  TFile *lPhiForest = new TFile(iPhiWeights,"READ");
  fPhiReader = (GBRForest*)lPhiForest->Get(fPhiMethodName);
  lPhiForest->Close();

  TFile *lU1Forest = new TFile(iU1Weights,"READ");
  fU1Reader  = (GBRForest*)lU1Forest->Get(fU1MethodName);
  lU1Forest->Close();

  TFile *lCovU1Forest = new TFile(iCovU1Weights,"READ");
  fCovU1Reader       = (GBRForest*)lCovU1Forest->Get(fCovU1MethodName);
  lCovU1Forest->Close();

  TFile *lCovU2Forest = new TFile(iCovU2Weights,"READ");
  fCovU2Reader       = (GBRForest*)lCovU2Forest->Get(fCovU2MethodName);
  lCovU2Forest->Close();
  
  fCov = new TMatrixD(2,2);
  fPhiVals = new Float_t[23];
  fU1Vals  = new Float_t[25];
  fCovVals = new Float_t[26];

  fMetLeptonTools = new MetLeptonTools();
}
//--------------------------------------------------------------------------------------------------
Float_t* MVAMet::getVals() { 
  fU1Vals[0]  =  fSumEt   ;
  fU1Vals[1]  =  fNPV     ;
  fU1Vals[2]  =  fU       ;
  fU1Vals[3]  =  fUPhi    ;
  fU1Vals[4]  =  fTKSumEt ;
  fU1Vals[5]  =  fTKU     ;
  fU1Vals[6]  =  fTKUPhi  ;
  fU1Vals[7]  =  fNPSumEt ;
  fU1Vals[8]  =  fNPU     ;
  fU1Vals[9]  =  fNPUPhi  ;
  fU1Vals[10] =  fPUSumEt ;
  fU1Vals[11] =  fPUMet   ;
  fU1Vals[12] =  fPUMetPhi;
  fU1Vals[13] =  fPCSumEt ;
  fU1Vals[14] =  fPCU     ;
  fU1Vals[15] =  fPCUPhi  ;
  fU1Vals[16] =  fJSPt1   ;
  fU1Vals[17] =  fJSEta1  ;
  fU1Vals[18] =  fJSPhi1  ;
  fU1Vals[19] =  fJSPt2   ;
  fU1Vals[20] =  fJSEta2  ;
  fU1Vals[21] =  fJSPhi2  ;
  fU1Vals[22] =  fNAllJet ;
  fU1Vals[23] =  fNJet    ;
  fU1Vals[24] =  fUPhiMVA ;
  return fU1Vals;
}
//--------------------------------------------------------------------------------------------------
Double_t MVAMet::evaluatePhi() { 
  fPhiVals[0]  =  fNPV     ;
  fPhiVals[1]  =  fU       ;
  fPhiVals[2]  =  fUPhi    ;
  fPhiVals[3]  =  fTKSumEt ;
  fPhiVals[4]  =  fTKU     ;
  fPhiVals[5]  =  fTKUPhi  ;
  fPhiVals[6]  =  fNPSumEt ;
  fPhiVals[7]  =  fNPU     ;
  fPhiVals[8]  =  fNPUPhi  ;
  fPhiVals[9]  =  fPUSumEt ;
  fPhiVals[10] =  fPUMet   ;
  fPhiVals[11] =  fPUMetPhi;
  fPhiVals[12] =  fPCSumEt ;
  fPhiVals[13] =  fPCU     ;
  fPhiVals[14] =  fPCUPhi  ;
  fPhiVals[15] =  fJSPt1   ;
  fPhiVals[16] =  fJSEta1  ;
  fPhiVals[17] =  fJSPhi1  ;
  fPhiVals[18] =  fJSPt2   ;
  fPhiVals[19] =  fJSEta2  ;
  fPhiVals[20] =  fJSPhi2  ;
  fPhiVals[21] =  fNAllJet ;
  fPhiVals[22] =  fNJet    ;
  return fPhiReader->GetResponse(fPhiVals);
}
//--------------------------------------------------------------------------------------------------
Double_t MVAMet::evaluateU1() { 
  fU1Vals[0]  =  fSumEt   ;
  fU1Vals[1]  =  fNPV     ;
  fU1Vals[2]  =  fU       ;
  fU1Vals[3]  =  fUPhi    ;
  fU1Vals[4]  =  fTKSumEt ;
  fU1Vals[5]  =  fTKU     ;
  fU1Vals[6]  =  fTKUPhi  ;
  fU1Vals[7]  =  fNPSumEt ;
  fU1Vals[8]  =  fNPU     ;
  fU1Vals[9]  =  fNPUPhi  ;
  fU1Vals[10] =  fPUSumEt ;
  fU1Vals[11] =  fPUMet   ;
  fU1Vals[12] =  fPUMetPhi;
  fU1Vals[13] =  fPCSumEt ;
  fU1Vals[14] =  fPCU     ;
  fU1Vals[15] =  fPCUPhi  ;
  fU1Vals[16] =  fJSPt1   ;
  fU1Vals[17] =  fJSEta1  ;
  fU1Vals[18] =  fJSPhi1  ;
  fU1Vals[19] =  fJSPt2   ;
  fU1Vals[20] =  fJSEta2  ;
  fU1Vals[21] =  fJSPhi2  ;
  fU1Vals[22] =  fNAllJet ;
  fU1Vals[23] =  fNJet    ;
  fU1Vals[24] =  fUPhiMVA ;
  return fU1Reader->GetResponse(fU1Vals);
}
//--------------------------------------------------------------------------------------------------
Double_t MVAMet::evaluateCovU1() {
  fCovVals[0]  =  fSumEt   ;
  fCovVals[1]  =  fNPV     ;
  fCovVals[2]  =  fU       ;
  fCovVals[3]  =  fUPhi    ;
  fCovVals[4]  =  fTKSumEt ;
  fCovVals[5]  =  fTKU     ;
  fCovVals[6]  =  fTKUPhi  ;
  fCovVals[7]  =  fNPSumEt ;
  fCovVals[8]  =  fNPU     ;
  fCovVals[9]  =  fNPUPhi  ;
  fCovVals[10]  =  fPUSumEt ;
  fCovVals[11] =  fPUMet   ;
  fCovVals[12] =  fPUMetPhi;
  fCovVals[13] =  fPCSumEt ;
  fCovVals[14] =  fPCU     ;
  fCovVals[15] =  fPCUPhi  ;
  fCovVals[16] =  fJSPt1   ;
  fCovVals[17] =  fJSEta1  ;
  fCovVals[18] =  fJSPhi1  ;
  fCovVals[19] =  fJSPt2   ;
  fCovVals[20] =  fJSEta2  ;
  fCovVals[21] =  fJSPhi2  ;
  fCovVals[22] =  fNAllJet ;
  fCovVals[23] =  fNJet    ;
  fCovVals[24] =  fUPhiMVA ;
  fCovVals[25] =  fUMVA    ;
  double lCovU1 = fCovU1Reader->GetResponse(fCovVals);
  if(!fOld42) lCovU1 = lCovU1*lCovU1*fUMVA*fUMVA; 
  return lCovU1;
}
//--------------------------------------------------------------------------------------------------
Double_t MVAMet::evaluateCovU2() {
  fCovVals[0]  =  fSumEt   ;
  fCovVals[1]  =  fNPV     ;
  fCovVals[2]  =  fU       ;
  fCovVals[3]  =  fUPhi    ;
  fCovVals[4]  =  fTKSumEt ;
  fCovVals[5]  =  fTKU     ;
  fCovVals[6]  =  fTKUPhi  ;
  fCovVals[7]  =  fNPSumEt ;
  fCovVals[8]  =  fNPU     ;
  fCovVals[9]  =  fNPUPhi  ;
  fCovVals[10] =  fPUSumEt ;
  fCovVals[11] =  fPUMet   ;
  fCovVals[12] =  fPUMetPhi;
  fCovVals[13] =  fPCSumEt ;
  fCovVals[14] =  fPCU     ;
  fCovVals[15] =  fPCUPhi  ;
  fCovVals[16] =  fJSPt1   ;
  fCovVals[17] =  fJSEta1  ;
  fCovVals[18] =  fJSPhi1  ;
  fCovVals[19] =  fJSPt2   ;
  fCovVals[20] =  fJSEta2  ;
  fCovVals[21] =  fJSPhi2  ;
  fCovVals[22] =  fNAllJet ;
  fCovVals[23] =  fNJet    ;
  fCovVals[24] =  fUPhiMVA ;
  fCovVals[25] =  fUMVA    ;
  double lCovU2 = fCovU2Reader->GetResponse(fCovVals);
  if(!fOld42) lCovU2 = lCovU2*lCovU2*fUMVA*fUMVA; 
  return lCovU2;
}
//-------------------------------------------------------------------------------------------------- 
Double_t MVAMet::MVAValue(  Bool_t iPhi, 
			    Float_t iPFSumEt, 
			    Float_t iU      ,
			    Float_t iUPhi   ,
			    Float_t iTKSumEt,
			    Float_t iTKU    ,
			    Float_t iTKUPhi ,
			    Float_t iNPSumEt,
			    Float_t iNPU    ,
			    Float_t iNPUPhi ,
			    Float_t iPUSumEt,
			    Float_t iPUMet  ,
			    Float_t iPUMetPhi,
			    Float_t iPCSumEt,
			    Float_t iPCU    ,
			    Float_t iPCUPhi ,
			    Float_t iJSPt1  ,
			    Float_t iJSEta1 ,
			    Float_t iJSPhi1 ,
			    Float_t iJSPt2  ,
			    Float_t iJSEta2 ,
			    Float_t iJSPhi2 ,
			    Float_t iNAllJet,
			    Float_t iNJet   ,
			    Float_t iNPV    
			    ){
  
  if (!fIsInitialized) { 
    std::cout << "Error: MVA Met not properly initialized.\n"; 
    return -9999;
  }
  
  fSumEt    = iPFSumEt;
  fU        = iU      ;
  fUPhi     = iUPhi   ;
  fTKSumEt  = iTKSumEt/iPFSumEt;
  fTKU      = iTKU    ;
  fTKUPhi   = iTKUPhi ;
  fNPSumEt  = iNPSumEt/iPFSumEt;
  fNPU      = iNPU    ;
  fNPUPhi   = iNPUPhi ;
  fPUSumEt  = iPUSumEt/iPFSumEt;
  fPUMet    = iPUMet  ;
  fPUMetPhi = iPUMetPhi;
  fPCSumEt  = iPCSumEt/iPFSumEt;
  fPCU      = iPCU    ;
  fPCUPhi   = iPCUPhi ;
  fJSPt1    = iJSPt1  ;
  fJSEta1   = iJSEta1 ;
  fJSPhi1   = iJSPhi1 ;
  fJSPt2    = iJSPt2  ;
  fJSEta2   = iJSEta2 ;
  fJSPhi2   = iJSPhi2 ;
  fNAllJet  = iNAllJet;
  fNJet     = iNJet   ;
  fNPV      = iNPV    ;
 
  Double_t lMVA = -9999;  
  lMVA = evaluatePhi();
  if(!iPhi) fUPhiMVA = iUPhi+lMVA;
  //Not no nice feature of the training
  //fTKSumEt  /= iPFSumEt;
  //fNPSumEt  /= iPFSumEt;
  //fPUSumEt  /= iPFSumEt;
  //fPCSumEt  /= iPFSumEt; 
  if(!iPhi) lMVA  = evaluateU1();
  return lMVA;
}
//--------------------------------------------------------------------------------------------------
//====> Please not that the jet collection must be cleaned => all jets near leptons must be removed
Met MVAMet::GetMet(	Bool_t iPhi,
			Float_t iPtVis,Float_t iPhiVis,Float_t iSumEtVis,
			Float_t iPtQ  ,Float_t iPhiQ  ,Float_t iSumEtQ,
			const PFMet            *iMet  ,
			const PFCandidateCol   *iCands,
			const Vertex *iVertex,const VertexCol *iVertices,
			const PFJetCol         *iJets ,
			FactorizedJetCorrector *iJetCorrector,
			const PileupEnergyDensityCol *iPUEnergyDensity,
			int iNPV,
      Bool_t hasPhotons,
			Bool_t printDebug) {

  int lNPV = 0;
  for(unsigned int i0 = 0; i0 < iVertices->GetEntries(); i0++) { 
    const Vertex *lPV = iVertices->At(i0);
    if(lPV->Ndof()           < 4.0)       continue;
    if(fabs(lPV->Z())        > 24.)       continue;
    if(lPV->Position().Rho() > 2.)        continue;
    lNPV++;
  }
  
  Met lPFRec = fRecoilTools->pfRecoil   (iPtVis,iPhiVis,iSumEtVis,iCands);
  Met lTKRec = fRecoilTools->trackRecoil(iPtQ  ,iPhiQ  ,iSumEtQ  ,iCands,iVertex); 
  Met lNPRec = fRecoilTools->NoPURecoil (iPtQ  ,iPhiQ  ,iSumEtQ  ,iJets,iJetCorrector,iPUEnergyDensity,iCands,iVertex,iVertices);
  Met lPCRec = fRecoilTools->PUCRecoil  (iPtVis,iPhiVis,iSumEtVis,iJets,iJetCorrector,iPUEnergyDensity,iCands,iVertex,iVertices,fType1);  
  Met lPUMet = fRecoilTools->PUMet      (iJets,iJetCorrector,iPUEnergyDensity,iCands,iVertex,iVertices);
  if(fType1) lPFRec = fRecoilTools->pfRecoilType1(iPtVis,iPhiVis,iSumEtVis,iCands,iJets,iJetCorrector,iPUEnergyDensity);

  Double_t lPt0 = 0; const PFJet *lLead = 0; 
  Double_t lPt1 = 0; const PFJet *l2nd  = 0; 
  int lNAllJet  = 0;
  int lNJet     = 0;
  for(unsigned int i0 = 0; i0 < iJets->GetEntries(); i0++) {
    const PFJet *pJet = iJets->At(i0);
    Double_t pPt = fRecoilTools->fJetIDMVA->correctedPt(pJet,iJetCorrector,iPUEnergyDensity);
    //if(pPt < 0) continue;
    if(!JetTools::passPFLooseId(pJet)) continue;
    if(f42 && pPt > 1.) lNAllJet++;
    if(!f42)            lNAllJet++;
    if(pPt  > 30)  lNJet++;
    if(lPt0 < pPt) {lPt1 = lPt0; lPt0 = pPt;   l2nd  = lLead; lLead = pJet; continue;}    
    if(lPt1 < pPt) {lPt1 = pPt;  l2nd  = pJet; continue;}
  }

  fSumEt    = lPFRec.SumEt();
  fU        = lPFRec.Pt();
  fUPhi     = lPFRec.Phi();
  fTKSumEt  = lTKRec.SumEt()/lPFRec.SumEt();
  fTKU      = lTKRec.Pt();
  fTKUPhi   = lTKRec.Phi();
  fNPSumEt  = lNPRec.SumEt()/lPFRec.SumEt();
  fNPU      = lNPRec.Pt();
  fNPUPhi   = lNPRec.Phi();
  fPUSumEt  = lPUMet.SumEt()/lPFRec.SumEt();
  fPUMet    = lPUMet.Pt();
  fPUMetPhi = lPUMet.Phi();
  fPCSumEt  = lPCRec.SumEt()/lPFRec.SumEt();
  fPCU      = lPCRec.Pt()    ;
  fPCUPhi   = lPCRec.Phi()   ;
  fJSPt1    = lPt0; 
  fJSEta1   = 0; if(lLead != 0) fJSEta1 = lLead->Eta();
  fJSPhi1   = 0; if(lLead != 0) fJSPhi1 = lLead->Phi();
  fJSPt2    = lPt1; 
  fJSEta2   = 0; if(l2nd  != 0) fJSEta2 = l2nd ->Eta();
  fJSPhi2   = 0; if(l2nd  != 0) fJSPhi2 = l2nd ->Phi();
  fNJet     = lNJet   ;
  fNAllJet  = lNAllJet;
  fNPV      = lNPV    ;
  Float_t lMVA = evaluatePhi();
  //Perform photon check and fix events with unphysical properties
  if (hasPhotons)
    if (fTKU/iPtVis < 0.1 || fNPU/iPtVis < 0.1) lMVA = 0.;
  if(!iPhi) fUPhiMVA = fUPhi + lMVA; 
  //Not no nice feature of teh training
  //fTKSumEt  /= lPFRec.SumEt();
  //fNPSumEt  /= lPFRec.SumEt();
  //fPUSumEt  /= lPFRec.SumEt();
  //fPCSumEt  /= lPFRec.SumEt();
  if(!iPhi) lMVA     = evaluateU1();//fU1Reader    ->EvaluateMVA( fU1MethodName );  
  //Perform photon check and fix events with unphysical properties
  if (hasPhotons)
    if (fTKU/iPtVis < 0.1 || fNPU/iPtVis < 0.1) lMVA = 1.;
  fUMVA              = fU*lMVA;

  TLorentzVector lUVec(0,0,0,0);   lUVec .SetPtEtaPhiM(fUMVA,0,fUPhiMVA,0);
  TLorentzVector lVVec(0,0,0,0);   lVVec .SetPtEtaPhiM(iPtVis ,0,iPhiVis ,0);
  if(lMVA < 0) lUVec .RotateZ(TMath::Pi());                                                   
  lUVec      -= lVVec;
  Met lMet(lUVec.Px(),lUVec.Py());
  
  //TMatrixD     lCov(2,2);
  //Covariance matrix perpendicular and parallel to the recoil (sort of correct)                                                                                                                            
  double lCovU1 = evaluateCovU1();
  double lCovU2 = evaluateCovU2();

  //Now Compute the covariance matrix in X and Y                                                                                                                                                           
  double lCos2 = cos(fUPhiMVA)*cos(fUPhiMVA);
  double lSin2 = sin(fUPhiMVA)*sin(fUPhiMVA);

  //Now Compute the covariance matrix in X and Y                                                                                                                                                           
  (*fCov)(0,0)   =  lCovU1*lCos2+lCovU2*lSin2;
  (*fCov)(1,0)   = -lCovU1*sin(fUPhiMVA)*cos(fUPhiMVA)+lCovU2*sin(fUPhiMVA)*cos(fUPhiMVA);
  (*fCov)(0,1)   =  (*fCov)(1,0);
  (*fCov)(1,1)   =  lCovU1*lSin2+lCovU2*lCos2;
  TMatrixD lInv(2,2); 
  lInv(0,0) = (*fCov)(0,0); lInv(1,1) = (*fCov)(1,1); lInv(1,0) = (*fCov)(1,0); lInv(0,1) = (*fCov)(0,1);
  if(lInv.Determinant() != 0) lInv.Invert();
  fSignificance   = TMath::Sqrt(lUVec.Px()*lUVec.Px()*(lInv)(0,0) + 2.*lUVec.Px()*lUVec.Py()*(lInv)(1,0)  + lUVec.Py()*lUVec.Py()*(lInv)(1,1)); 
  fUncertainty     = sqrt(lCovU1+lCovU2);
  
  if (printDebug == kTRUE) {
    std::cout << "Debug Met MVA: "
	      <<  fU        << " : "
	      <<  fUPhi     << " : "
	      <<  fTKSumEt  << " : "
	      <<  fTKU      << " : "
	      <<  fTKUPhi   << " : "
	      <<  fNPSumEt  << " : "
	      <<  fNPU      << " : "
	      <<  fNPUPhi   << " : "
	      <<  fPUSumEt  << " : "
	      <<  fPUMet    << " : "
	      <<  fPUMetPhi << " : "
	      <<  fPCUPhi   << " : "
	      <<  fPCSumEt  << " : "
	      <<  fPCU     << " : "
	      <<  fPCUPhi  << " : "
	      <<  fJSPt1   << " : "
	      <<  fJSEta1  << " : "
	      <<  fJSPhi1  << " : "
	      <<  fJSPt2   << " : "
	      <<  fJSEta2  << " : "
	      <<  fJSPhi2  << " : "
	      <<  fNJet    << " : "
	      <<  fNAllJet << " : "
	      <<  fNPV     << " : "
              << " === : === "
              << lMet.Pt()  << " : "
              << lMet.Phi() << " : "
              << std::endl;
  }

  return lMet;
}
//--------------------------------------------------------------------------------------------------
//====> Please not that the jet collection must be cleaned => all jets near leptons must be removed
//====> Corrected Jets
Met MVAMet::GetMet(	Bool_t iPhi,
			Float_t iPtVis,Float_t iPhiVis,Float_t iSumEtVis,
			Float_t iPtQ  ,Float_t iPhiQ  ,Float_t iSumEtQ  , //Charged components
			const PFMet            *iMet  ,
			const PFCandidateCol   *iCands,
			const Vertex *iVertex         ,const VertexCol *iVertices,Double_t iRho,
			const PFJetCol         *iJets ,
			int iNPV,
			Bool_t printDebug) {

  int lNPV = 0;
  for(unsigned int i0 = 0; i0 < iVertices->GetEntries(); i0++) { 
    const Vertex *lPV = iVertices->At(i0);
    if(lPV->Ndof()           < 4.0)       continue;
    if(fabs(lPV->Z())        > 24.)       continue;
    if(lPV->Position().Rho() > 2.)        continue;
    lNPV++;
  }
  Met lPFRec = fRecoilTools->pfRecoil   (iPtVis,iPhiVis,iSumEtVis,iCands);
  Met lTKRec = fRecoilTools->trackRecoil(iPtQ  ,iPhiQ  ,iSumEtQ  ,      iCands,iVertex); 
  Met lNPRec = fRecoilTools->NoPURecoil (iPtQ  ,iPhiQ  ,iSumEtQ  ,iJets,iCands,iVertex,iVertices,iRho);  
  Met lPCRec = fRecoilTools->PUCRecoil  (iPtVis,iPhiVis,iSumEtVis,iJets,iCands,iVertex,iVertices,iRho,fType1);
  Met lPUMet = fRecoilTools->PUMet      (                         iJets,iCands,iVertex,iVertices,iRho);
  //FIXME Not yet supported
  //if(fType1) lPFRec = fRecoilTools->pfRecoilType1(iPtVis,iPhiVis,iSumEtVis,iCands,iJets);

  Double_t lPt0 = 0; const PFJet *lLead = 0; 
  Double_t lPt1 = 0; const PFJet *l2nd  = 0; 
  int lNAllJet  = 0;
  int lNJet     = 0;
  for(unsigned int i0 = 0; i0 < iJets->GetEntries(); i0++) {
    const PFJet *pJet = iJets->At(i0);
    Double_t pPt = pJet->Pt();
    if(!JetTools::passPFLooseId(pJet)) continue;
    if(f42 && pPt > 1.) lNAllJet++;
    if(!f42)            lNAllJet++;
    if(pPt  > 30.)  lNJet++;
    if(lPt0 < pPt) {lPt1 = lPt0; lPt0  = pPt; l2nd  = lLead; lLead = pJet; continue;}    
    if(lPt1 < pPt) {lPt1 = pPt;  l2nd  = pJet; continue;}
  }
  
  fSumEt    = lPFRec.SumEt();
  fU        = lPFRec.Pt();
  fUPhi     = lPFRec.Phi();
  fTKSumEt  = lTKRec.SumEt()/lPFRec.SumEt();
  fTKU      = lTKRec.Pt();
  fTKUPhi   = lTKRec.Phi();
  fNPSumEt  = lNPRec.SumEt()/lPFRec.SumEt();
  fNPU      = lNPRec.Pt();
  fNPUPhi   = lNPRec.Phi();
  fPUSumEt  = lPUMet.SumEt()/lPFRec.SumEt();
  fPUMet    = lPUMet.Pt();
  fPUMetPhi = lPUMet.Phi();
  fPCSumEt  = lPCRec.SumEt()/lPFRec.SumEt();
  fPCU      = lPCRec.Pt()    ;
  fPCUPhi   = lPCRec.Phi()   ;
  fJSPt1    = lPt0; 
  fJSEta1   = 0; if(lLead != 0) fJSEta1 = lLead->Eta();
  fJSPhi1   = 0; if(lLead != 0) fJSPhi1 = lLead->Phi();
  fJSPt2    = lPt1; 
  fJSEta2   = 0; if(l2nd  != 0) fJSEta2 = l2nd ->Eta();
  fJSPhi2   = 0; if(l2nd  != 0) fJSPhi2 = l2nd ->Phi();
  fNJet     = lNJet   ;
  fNAllJet  = lNAllJet;
  fNPV      = lNPV    ;
  
  Float_t lMVA = evaluatePhi();
  
  if(!iPhi) fUPhiMVA = fUPhi + lMVA; 
  //Not no nice feature of teh training
  //fTKSumEt  /= lPFRec.SumEt();
  //fNPSumEt  /= lPFRec.SumEt();
  //fPUSumEt  /= lPFRec.SumEt();
  //fPCSumEt  /= lPFRec.SumEt();
  if(!iPhi) lMVA    = evaluateU1();//fU1Reader    ->EvaluateMVA( fU1MethodName );  
  fUMVA        = lMVA*fU;
  
  TLorentzVector lUVec(0,0,0,0);   lUVec .SetPtEtaPhiM(fUMVA  ,0,fUPhiMVA,0);
  TLorentzVector lVVec(0,0,0,0);   lVVec .SetPtEtaPhiM(iPtVis ,0,iPhiVis ,0);
  if(lMVA < 0) lUVec .RotateZ(TMath::Pi());                                                   
  lUVec      -= lVVec;
  Met lMet(lUVec.Px(),lUVec.Py());
  //Cov matrix
  double lCovU1 = evaluateCovU1();
  double lCovU2 = evaluateCovU2();

  double lCos2 = cos(fUPhiMVA)*cos(fUPhiMVA);
  double lSin2 = sin(fUPhiMVA)*sin(fUPhiMVA);
  
  //Now Compute teh covariance matrix in X and Y                                 
  (*fCov)(0,0)   =  lCovU1*lCos2+lCovU2*lSin2;
  (*fCov)(1,0)   = -lCovU1*sin(fUPhiMVA)*cos(fUPhiMVA)+lCovU2*sin(fUPhiMVA)*cos(fUPhiMVA);
  (*fCov)(0,1)   =  (*fCov)(1,0);
  (*fCov)(1,1)   =  lCovU1*lSin2+lCovU2*lCos2;
  TMatrixD lInv(2,2); 
  lInv(0,0) = (*fCov)(0,0); lInv(1,1) = (*fCov)(1,1); lInv(1,0) = (*fCov)(1,0); lInv(0,1) = (*fCov)(0,1);
  if(lInv.Determinant() != 0) lInv.Invert();
  fSignificance   = TMath::Sqrt(lUVec.Px()*lUVec.Px()*(lInv)(0,0) + 2.*lUVec.Px()*lUVec.Py()*(lInv)(1,0)  + lUVec.Py()*lUVec.Py()*(lInv)(1,1)); 
  fUncertainty     = sqrt(lCovU1+lCovU2);

  if (printDebug == kTRUE) {
    std::cout << "Debug Met MVA: "
	      <<  fU        << " : "
	      <<  fUPhi     << " : "
	      <<  fTKSumEt  << " : "
	      <<  fTKU      << " : "
	      <<  fTKUPhi   << " : "
	      <<  fNPSumEt  << " : "
	      <<  fNPU      << " : "
	      <<  fNPUPhi   << " : "
	      <<  fPUSumEt  << " : "
	      <<  fPUMet    << " : "
	      <<  fPUMetPhi << " : "
	      <<  fPCUPhi   << " : "
	      <<  fPCSumEt  << " : "
	      <<  fPCU     << " : "
	      <<  fPCUPhi  << " : "
	      <<  fJSPt1   << " : "
	      <<  fJSEta1  << " : "
	      <<  fJSPhi1  << " : "
	      <<  fJSPt2   << " : "
	      <<  fJSEta2  << " : "
	      <<  fJSPhi2  << " : "
	      <<  fNJet    << " : "
	      <<  fNAllJet << " : "
	      <<  fNPV     << " : "
              << " === : === "
              << lMet.Pt()  << " : "
              << lMet.Phi() << " : "
              << std::endl;
  }

  return lMet;
}
//--------------------------------------------------------------------------------------------------
//Interms of the two candidates
Met MVAMet::GetMet(	bool iPhi,
			Float_t iPt1,Float_t iPhi1,Float_t iEta1,Float_t iChargedFrac1,
			Float_t iPt2,Float_t iPhi2,Float_t iEta2,Float_t iChargedFrac2,
			const PFMet            *iMet  ,
			const PFCandidateCol   *iCands,
			const Vertex           *iVertex,const VertexCol *iVertices,
			const PFJetCol         *iJets ,
			FactorizedJetCorrector *iJetCorrector,
			const PileupEnergyDensityCol *iPUEnergyDensity,
			int iNPV,
			Bool_t printDebug) {

  int lNPV = 0;
  for(unsigned int i0 = 0; i0 < iVertices->GetEntries(); i0++) { 
    const Vertex *lPV = iVertices->At(i0);
    if(lPV->Ndof()           < 4.0)       continue;
    if(fabs(lPV->Z())        > 24.)       continue;
    if(lPV->Position().Rho() > 2.)        continue;
    lNPV++;
  }
  //Met lMVec1 = fRecoilTools->pfCone(iPhi1,iEta1,iCands);
  //Met lMVec2 = fRecoilTools->pfCone(iPhi2,iEta2,iCands);
  //double lChargedFrac1 = 1.;
  //double lChargedFrac2 = 1.;
  //if(iId == MVAMet::TauTau){
  //  lChargedFrac1 = chargedFracTau(
  // }
  
  TLorentzVector lVVec1   (0,0,0,0);   lVVec1   .SetPtEtaPhiM(iPt1,0              ,iPhi1,0);
  TLorentzVector lVVec2   (0,0,0,0);   lVVec2   .SetPtEtaPhiM(iPt2,0              ,iPhi2,0);
  TLorentzVector lVVec1Q  (0,0,0,0);   lVVec1Q  .SetPtEtaPhiM(iPt1*iChargedFrac1,0,iPhi1,0);
  TLorentzVector lVVec2Q  (0,0,0,0);   lVVec2Q  .SetPtEtaPhiM(iPt2*iChargedFrac2,0,iPhi2,0);

  if(iChargedFrac1 < 0 ) {
    Met lQ1Cone = fRecoilTools->pfCone(iPhi1,iEta1,iCands,iVertex,true);
    lVVec1Q.SetPtEtaPhiM(lQ1Cone.Pt(),0,lQ1Cone.Phi(),0);
  }
  if(iChargedFrac2 < 0 ) {
    Met lQ2Cone = fRecoilTools->pfCone(iPhi2,iEta2,iCands,iVertex,true);
    lVVec2Q.SetPtEtaPhiM(lQ2Cone.Pt(),0,lQ2Cone.Phi(),0);
  }
  //Met lPFRec = fRecoilTools->pfRecoil   (iPhi1,iEta1,iPhi2,iEta2,      iCands);
  //Met lTKRec = fRecoilTools->trackRecoil(iPhi1,iEta1,iPhi2,iEta2,      iCands,iVertex); 
  //Met lNPRec = fRecoilTools->NoPURecoil (iPhi1,iEta1,iPhi2,iEta2,iJets,iCands,iVertex,iVertices,iJetCorrector,iPUEnergyDensity); 
  //Met lPCRec = fRecoilTools->PUCRecoil  (iPhi1,iEta1,iPhi2,iEta2,iJets,iCands,iVertex,iVertices,iJetCorrector,iPUEnergyDensity);
  //Met lPUMet = fRecoilTools->PUMet      (iPhi1,iEta1,iPhi2,iEta2,iJets,iCands,iVertex,iVertices,iJetCorrector,iPUEnergyDensity);

  //Met lNPMet  = fRecoilTools->NoPUMet   (iJets,iJetCorrector,iPUEnergyDensity,iCands,iVertex,iVertices,iPhi1,iEta1,iPhi2,iEta2); 
  //std::cout << "===> " << lNPMet.Pt() << " -- " << lNPMet.Phi() << std::endl;
  //std::cout << "===> " << lVVec1.Pt() << " -- " << lVVec1.Phi() << " -- " << iPt1 << " -- " << iPhi1 << std::endl;
  //std::cout << "===> " << lVVec2.Pt() << " -- " << lVVec2.Phi() << " -- " << iPt2 << " -- " << iPhi2 << std::endl;
  //std::cout << "===> " << lNPRec.Pt() << " -- " << lNPRec.Phi() << std::endl;

  Float_t lSumEtVis  = lVVec1.Pt() + lVVec2.Pt();
  lVVec1+=lVVec2;
  Float_t lPtVis     = lVVec1.Pt();
  Float_t lPhiVis    = lVVec1.Phi();

  Float_t lSumEtQ    = lVVec1Q.Pt() + lVVec2Q.Pt();
  lVVec1Q+=lVVec2Q;
  Float_t lPtVisQ    = lVVec1Q.Pt();
  Float_t lPhiVisQ   = lVVec1Q.Phi();
 
  Met lPFRec = fRecoilTools->pfRecoil   (lPtVis       ,lPhiVis,       lSumEtVis,iCands);   
  Met lTKRec = fRecoilTools->trackRecoil(lPtVisQ      ,lPhiVisQ      ,lSumEtQ  ,iCands,iVertex);  
  Met lNPRec = fRecoilTools->NoPURecoil (lPtVisQ      ,lPhiVisQ      ,lSumEtQ  ,iJets,iJetCorrector,iPUEnergyDensity,iCands,iVertex,iVertices,       iPhi1,iEta1,iPhi2,iEta2);    
  Met lPCRec = fRecoilTools->PUCRecoil  (lPtVis       ,lPhiVis       ,lSumEtVis,iJets,iJetCorrector,iPUEnergyDensity,iCands,iVertex,iVertices,fType1,iPhi1,iEta1,iPhi2,iEta2);    
  Met lPUMet = fRecoilTools->PUMet      (                                       iJets,iJetCorrector,iPUEnergyDensity,iCands,iVertex,iVertices,       iPhi1,iEta1,iPhi2,iEta2);    
  if(fType1) lPFRec = fRecoilTools->pfRecoilType1(lPtVis,lPhiVis,lSumEtVis,iCands,iJets,iJetCorrector,iPUEnergyDensity,iPhi1,iEta1,iPhi2,iEta2);

  //Met lNPRec1 = fRecoilTools->NoPUMet (iJets,iJetCorrector,iPUEnergyDensity,iCands,iVertex,iVertices,0,0,0,0);
  //std::cout << " ====> PF => " << lPFRec.Pt() << " -- " << lPFRec1.Pt() << std::endl;
  //std::cout << " ====> TK => " << lTKRec.Pt() << " -- " << lTKRec1.Pt() << std::endl;
  //std::cout << " ====> NP =>  " << lNPRec.Pt() << " -- " << lNPRec1.Pt() << std::endl;
  //std::cout << " ====> PC => " << lPCRec.Pt() << " -- " << lPCRec1.Pt() << std::endl;
  
  Double_t lPt0 = 0; const PFJet *lLead = 0; 
  Double_t lPt1 = 0; const PFJet *l2nd  = 0; 
  int lNAllJet  = 0;
  int lNJet     = 0;
  for(unsigned int i0 = 0; i0 < iJets->GetEntries(); i0++) {
    const PFJet *pJet = iJets->At(i0);
    Double_t pPt = fRecoilTools->fJetIDMVA->correctedPt(pJet,iJetCorrector,iPUEnergyDensity);
    double pDEta1 = pJet->Eta() - iEta1;
    double pDPhi1 = fabs(pJet->Phi() - iPhi1); if(pDPhi1 > 2.*TMath::Pi()-pDPhi1) pDPhi1 = 2.*TMath::Pi()-pDPhi1;
    double pDR1   = sqrt(pDEta1*pDEta1 + pDPhi1*pDPhi1);
    if(pDR1 < 0.5) continue;
    double pDEta2 = pJet->Eta() - iEta2;
    double pDPhi2 = fabs(pJet->Phi() - iPhi2); if(pDPhi2 > 2.*TMath::Pi()-pDPhi2) pDPhi2 = 2.*TMath::Pi()-pDPhi2;
    double pDR2   = sqrt(pDEta2*pDEta2 + pDPhi2*pDPhi2);
    if(pDR2 < 0.5) continue;  
    if(!JetTools::passPFLooseId(pJet)) continue;
    if(f42 && pPt > 1.) lNAllJet++;
    if(!f42)            lNAllJet++;
    if(pPt  > 30.)  lNJet++;
    if(lPt0 < pPt) {lPt1 = lPt0; lPt0 = pPt; l2nd = lLead; lLead = pJet; continue;}    
    if(lPt1 < pPt) {lPt1 = pPt;  l2nd  = pJet; continue;}
  }
  
  fSumEt   = lPFRec.SumEt();
  fU       = lPFRec.Pt();
  fUPhi    = lPFRec.Phi();
  fTKSumEt = lTKRec.SumEt()/lPFRec.SumEt();
  fTKU     = lTKRec.Pt();
  fTKUPhi  = lTKRec.Phi();
  fNPSumEt = lNPRec.SumEt()/lPFRec.SumEt();
  fNPU     = lNPRec.Pt();
  fNPUPhi  = lNPRec.Phi();
  fPUSumEt = lPUMet.SumEt()/lPFRec.SumEt();
  fPUMet   = lPUMet.Pt();
  fPUMetPhi= lPUMet.Phi();
  fPCSumEt = lPCRec.SumEt()/lPFRec.SumEt();
  fPCU     = lPCRec.Pt()    ;
  fPCUPhi  = lPCRec.Phi()   ;
  fJSPt1   = lPt0; 
  fJSEta1  = 0; if(lLead != 0) fJSEta1 = lLead->Eta();
  fJSPhi1  = 0; if(lLead != 0) fJSPhi1 = lLead->Phi();
  fJSPt2   = lPt1; 
  fJSEta2  = 0; if(l2nd  != 0) fJSEta2 = l2nd ->Eta();
  fJSPhi2  = 0; if(l2nd  != 0) fJSPhi2 = l2nd ->Phi();
  fNJet    = lNJet   ;
  fNAllJet = lNAllJet;
  fNPV     = lNPV    ;
  
  Float_t lMVA = evaluatePhi();//fPhiReader                 ->EvaluateMVA( fPhiMethodName );
  if(!iPhi) fUPhiMVA = fUPhi + lMVA; 
  //fTKSumEt  /= lPFRec.SumEt();
  //fNPSumEt  /= lPFRec.SumEt();
  //fPUSumEt  /= lPFRec.SumEt();
  //fPCSumEt  /= lPFRec.SumEt();
  if(!iPhi) lMVA    = evaluateU1();//fU1Reader    ->EvaluateMVA( fU1MethodName );  
  fUMVA             = fU*lMVA;

  TLorentzVector lUVec (0,0,0,0);   lUVec .SetPtEtaPhiM(fUMVA  ,0,fUPhiMVA,0);
  TLorentzVector lVVec (0,0,0,0);   lVVec .SetPtEtaPhiM(lPtVis ,0,lPhiVis ,0);
  if(lMVA < 0) lUVec .RotateZ(TMath::Pi());                                                   
  lUVec      -= lVVec;
  Met lMet(lUVec.Px(),lUVec.Py());
  //Cov matrix                                                                                                                                                                                                
  double lCovU1 = evaluateCovU1();
  double lCovU2 = evaluateCovU2();

  double lCos2 = cos(fUPhiMVA)*cos(fUPhiMVA);
  double lSin2 = sin(fUPhiMVA)*sin(fUPhiMVA);

  //Now Compute teh covariance matrix in X and Y                                                                                                                                                              
  (*fCov)(0,0)   =  lCovU1*lCos2+lCovU2*lSin2;
  (*fCov)(1,0)   = -lCovU1*sin(fUPhiMVA)*cos(fUPhiMVA)+lCovU2*sin(fUPhiMVA)*cos(fUPhiMVA);
  (*fCov)(0,1)   =  (*fCov)(1,0);
  (*fCov)(1,1)   =  lCovU1*lSin2+lCovU2*lCos2;
  TMatrixD lInv(2,2); 
  lInv(0,0) = (*fCov)(0,0); lInv(1,1) = (*fCov)(1,1); lInv(1,0) = (*fCov)(1,0); lInv(0,1) = (*fCov)(0,1);
  if(lInv.Determinant() != 0) lInv.Invert();
  fSignificance   = TMath::Sqrt(lUVec.Px()*lUVec.Px()*(lInv)(0,0) + 2.*lUVec.Px()*lUVec.Py()*(lInv)(1,0)  + lUVec.Py()*lUVec.Py()*(lInv)(1,1)); 
  fUncertainty     = sqrt(lCovU1+lCovU2);


  if (printDebug == kTRUE) {
    std::cout << "Debug Met MVA: "
	      <<  fSumEt    << " : "
	      <<  fU        << " : "
	      <<  fUPhi     << " : "
	      <<  fTKSumEt  << " : "
	      <<  fTKU      << " : "
	      <<  fTKUPhi   << " : "
	      <<  fNPSumEt  << " : "
	      <<  fNPU      << " : "
	      <<  fNPUPhi   << " : "
	      <<  fPUSumEt  << " : "
	      <<  fPUMet    << " : "
	      <<  fPUMetPhi << " : "
	      <<  fPCSumEt  << " : "
	      <<  fPCU      << " : "
	      <<  fPCUPhi   << " : "
	      <<  fJSPt1    << " : "
	      <<  fJSEta1   << " : "
	      <<  fJSPhi1   << " : "
	      <<  fJSPt2    << " : "
	      <<  fJSEta2   << " : "
	      <<  fJSPhi2   << " : "
	      <<  fNJet     << " : "
	      <<  fNAllJet  << " : "
	      <<  fNPV      << " : "
              << " === : === "
              << lMet.Pt()  << " : "
              << lMet.Phi() << " : "
              << std::endl;
  }

  return lMet;
}
//--------------------------------------------------------------------------------------------------
//Interms of the two candidates => corrected Jets
Met MVAMet::GetMet(	Bool_t iPhi,
			Float_t iPt1,Float_t iPhi1,Float_t iEta1,Float_t iChargedFrac1,
			Float_t iPt2,Float_t iPhi2,Float_t iEta2,Float_t iChargedFrac2,
			const PFMet            *iMet  ,
			const PFCandidateCol   *iCands,
			const Vertex           *iVertex,const VertexCol *iVertices,Double_t iRho,
			const PFJetCol         *iJets ,
			int iNPV,
			Bool_t printDebug) {
  int lNPV = 0;
  for(unsigned int i0 = 0; i0 < iVertices->GetEntries(); i0++) { 
    const Vertex *lPV = iVertices->At(i0);
    if(lPV->Ndof()           < 4.0)       continue;
    if(fabs(lPV->Z())        > 24.)       continue;
    if(lPV->Position().Rho() > 2.)        continue;
    lNPV++;
  }
  
  TLorentzVector lVVec1(0,0,0,0);   lVVec1.SetPtEtaPhiM(iPt1,0,iPhi1 ,0);
  TLorentzVector lVVec2(0,0,0,0);   lVVec2.SetPtEtaPhiM(iPt2,0,iPhi2 ,0);
  TLorentzVector lVVec1Q  (0,0,0,0);   lVVec1Q  .SetPtEtaPhiM(iPt1*iChargedFrac1,0,iPhi1,0);
  TLorentzVector lVVec2Q  (0,0,0,0);   lVVec2Q  .SetPtEtaPhiM(iPt2*iChargedFrac2,0,iPhi2,0);

  if(iChargedFrac1 < 0 ) {
    Met lQ1Cone = fRecoilTools->pfCone(iPhi1,iEta1,iCands,iVertex,true);
    lVVec1Q.SetPtEtaPhiM(lQ1Cone.Pt(),0,lQ1Cone.Phi(),0);
   }
  if(iChargedFrac2 < 0 ) {
    Met lQ2Cone = fRecoilTools->pfCone(iPhi2,iEta2,iCands,iVertex,true);
    lVVec2Q.SetPtEtaPhiM(lQ2Cone.Pt(),0,lQ2Cone.Phi(),0);
  }

  Float_t lSumEtVis  = lVVec1.Pt() + lVVec2.Pt();
  lVVec1+=lVVec2;
  Float_t lPtVis     = lVVec1.Pt();
  Float_t lPhiVis    = lVVec1.Phi();

  Float_t lSumEtQ    = lVVec1Q.Pt() + lVVec2Q.Pt();
  lVVec1Q+=lVVec2Q;
  Float_t lPtVisQ    = lVVec1Q.Pt();
  Float_t lPhiVisQ   = lVVec1Q.Phi();

  Met lPFRec = fRecoilTools->pfRecoil   (lPtVis       ,lPhiVis,       lSumEtVis,iCands);   
  Met lTKRec = fRecoilTools->trackRecoil(lPtVisQ      ,lPhiVisQ      ,lSumEtQ  ,iCands,iVertex);  
  Met lNPRec = fRecoilTools->NoPURecoil (lPtVisQ      ,lPhiVisQ      ,lSumEtQ  ,iJets,iCands,iVertex,iVertices,1.,       iPhi1,iEta1,iPhi2,iEta2);    
  Met lPCRec = fRecoilTools->PUCRecoil  (lPtVis       ,lPhiVis       ,lSumEtVis,iJets,iCands,iVertex,iVertices,1.,fType1,iPhi1,iEta1,iPhi2,iEta2);    
  Met lPUMet = fRecoilTools->PUMet      (                                       iJets,iCands,iVertex,iVertices,1.,       iPhi1,iEta1,iPhi2,iEta2);    
  //FIXME Not yet supported
  //if(fType1) lPFRec = fRecoilTools->pfRecoilType1(iPtVis,iPhiVis,iSumEtVis,iCands);
  
  //Met lPFRec = fRecoilTools->pfRecoil   (iPhi1,iEta1,iPhi2,iEta2,      iCands);
  //Met lTKRec = fRecoilTools->trackRecoil(iPhi1,iEta1,iPhi2,iEta2,      iCands,iVertex); 
  //Met lNPRec = fRecoilTools->NoPURecoil (iPhi1,iEta1,iPhi2,iEta2,iJets,iCands,iVertex,iVertices);  //Not using rho
  //Met lPCRec = fRecoilTools->PUCRecoil  (iPhi1,iEta1,iPhi2,iEta2,iJets,iCands,iVertex,iVertices);
  //Met lPUMet = fRecoilTools->PUMet      (iPhi1,iEta1,iPhi2,iEta2,iJets,iCands,iVertex,iVertices);

  Double_t lPt0 = 0; const PFJet *lLead = 0; 
  Double_t lPt1 = 0; const PFJet *l2nd  = 0; 
  int lNAllJet  = 0;
  int lNJet     = 0;
  for(unsigned int i0 = 0; i0 < iJets->GetEntries(); i0++) {
    const PFJet *pJet = iJets->At(i0);
    Double_t pPt = pJet->Pt();
    double pDEta1 = pJet->Eta() - iEta1;
    double pDPhi1 = fabs(pJet->Phi() - iPhi1); if(pDPhi1 > 2.*TMath::Pi()-pDPhi1) pDPhi1 = 2.*TMath::Pi()-pDPhi1;
    double pDR1   = sqrt(pDEta1*pDEta1 + pDPhi1*pDPhi1);
    if(pDR1 < 0.5) continue;
    double pDEta2 = pJet->Eta() - iEta2;
    double pDPhi2 = fabs(pJet->Phi() - iPhi2); if(pDPhi2 > 2.*TMath::Pi()-pDPhi2) pDPhi2 = 2.*TMath::Pi()-pDPhi2;
    double pDR2   = sqrt(pDEta2*pDEta2 + pDPhi2*pDPhi2);
    if(pDR2 < 0.5) continue;  
    if(!JetTools::passPFLooseId(pJet)) continue;
    if(f42 && pPt > 1.) lNAllJet++;
    if(!f42)            lNAllJet++;
    if(pPt  > 30.)  lNJet++;
    if(lPt0 < pPt) {lPt1 = lPt0; lPt0  = pPt; l2nd  = lLead; lPt1 = l2nd->Pt(); lLead = pJet; continue;}        
    if(lPt1 < pPt) {lPt1 = pPt;  l2nd  = pJet; continue;}
  }
  fSumEt   = lPFRec.SumEt();
  fU       = lPFRec.Pt();
  fUPhi    = lPFRec.Phi();
  fTKSumEt = lTKRec.SumEt()/lPFRec.SumEt();
  fTKU     = lTKRec.Pt();
  fTKUPhi  = lTKRec.Phi();
  fNPSumEt = lNPRec.SumEt()/lPFRec.SumEt();
  fNPU     = lNPRec.Pt();
  fNPUPhi  = lNPRec.Phi();
  fPUSumEt = lPUMet.SumEt()/lPFRec.SumEt();
  fPUMet   = lPUMet.Pt();
  fPUMetPhi= lPUMet.Phi();
  fPCSumEt = lPCRec.SumEt()/lPFRec.SumEt();
  fPCU     = lPCRec.Pt()    ;
  fPCUPhi  = lPCRec.Phi()   ;
  fJSPt1   = lPt0; 
  fJSEta1  = 0; if(lLead != 0) fJSEta1 = lLead->Eta();
  fJSPhi1  = 0; if(lLead != 0) fJSPhi1 = lLead->Phi();
  fJSPt2   = lPt1; 
  fJSEta2  = 0; if(l2nd  != 0) fJSEta2 = l2nd ->Eta();
  fJSPhi2  = 0; if(l2nd  != 0) fJSPhi2 = l2nd ->Phi();
  fNJet    = lNJet   ;
  fNAllJet = lNAllJet;
  fNPV     = lNPV    ;

  Float_t lMVA = evaluatePhi();//fPhiReader                 ->EvaluateMVA( fPhiMethodName );
  
  if(!iPhi) fUPhiMVA = fUPhi + lMVA; 
  //fTKSumEt  /= lPFRec.SumEt();
  //fNPSumEt  /= lPFRec.SumEt();
  //fPUSumEt  /= lPFRec.SumEt();
  //fPCSumEt  /= lPFRec.SumEt();
  if(!iPhi) lMVA    = evaluateU1();//fU1Reader    ->EvaluateMVA( fU1MethodName );  
  fUMVA        = lMVA*fU;
  
  TLorentzVector lUVec (0,0,0,0);   lUVec .SetPtEtaPhiM(fUMVA  ,0,fUPhiMVA,0);
  TLorentzVector lVVec (0,0,0,0);   lVVec .SetPtEtaPhiM(lPtVis ,0,lPhiVis ,0);
  if(lMVA   < 0) lUVec .RotateZ(TMath::Pi());                                                   
  lUVec      -= lVVec;
  Met lMet(lUVec.Px(),lUVec.Py());
  //Cov matrix                                                                                                                                                                           
  double lCovU1 = evaluateCovU1();
  double lCovU2 = evaluateCovU2();
  
  double lCos2 = cos(fUPhiMVA)*cos(fUPhiMVA);
  double lSin2 = sin(fUPhiMVA)*sin(fUPhiMVA);
  //Now Compute teh covariance matrix in X and Y                                                                                                                               
  (*fCov)(0,0)   =  lCovU1*lCos2+lCovU2*lSin2;
  (*fCov)(1,0)   = -lCovU1*sin(fUPhiMVA)*cos(fUPhiMVA)+lCovU2*sin(fUPhiMVA)*cos(fUPhiMVA);
  (*fCov)(0,1)   =  (*fCov)(1,0);
  (*fCov)(1,1)   =  lCovU1*lSin2+lCovU2*lCos2;
  TMatrixD lInv(2,2); 
  lInv(0,0) = (*fCov)(0,0); lInv(1,1) = (*fCov)(1,1); lInv(1,0) = (*fCov)(1,0); lInv(0,1) = (*fCov)(0,1);
  if(lInv.Determinant() != 0) lInv.Invert();
  fSignificance   = TMath::Sqrt(lUVec.Px()*lUVec.Px()*(lInv)(0,0) + 2.*lUVec.Px()*lUVec.Py()*(lInv)(1,0)  + lUVec.Py()*lUVec.Py()*(lInv)(1,1)); 
  fUncertainty     = sqrt(lCovU1+lCovU2);  
  
  if (printDebug == kTRUE) {
    std::cout << "Debug Met MVA: "
	      <<  fU        << " : "
	      <<  fUPhi     << " : "
	      <<  fTKSumEt  << " : "
	      <<  fTKU      << " : "
	      <<  fTKUPhi   << " : "
	      <<  fNPSumEt  << " : "
	      <<  fNPU      << " : "
	      <<  fNPUPhi   << " : "
	      <<  fPUSumEt  << " : "
	      <<  fPUMet    << " : "
	      <<  fPUMetPhi << " : "
	      <<  fPCSumEt  << " : "
	      <<  fPCU      << " : "
	      <<  fPCUPhi   << " : "
	      <<  fJSPt1    << " : "
	      <<  fJSEta1   << " : "
	      <<  fJSPhi1   << " : "
	      <<  fJSPt2    << " : "
	      <<  fJSEta2   << " : "
	      <<  fJSPhi2   << " : "
	      <<  fNJet     << " : "
	      <<  fNAllJet  << " : "
	      <<  fNPV      << " : "
              << " === : === "
              << lMet.Pt()  << " : "
              << lMet.Phi() << " : "
              << std::endl;
  }

  return lMet;
}
//--------------------------------------------------------------------------------------------------
Met MVAMet::GetMet(const MuonCol *iMuons,const ElectronCol *iElectrons,const PFTauCol *iTaus,
                   const PFCandidateCol *iCands,const PFJetCol *iJets,const Vertex *iPV,const VertexCol *iVertices,const PFMetCol *iPFMet,
                   FactorizedJetCorrector *iJetCorrector,const PileupEnergyDensityCol* iPUEnergyDensity) {
  const Vertex *lPV = iPV; if(iPV == 0 && iVertices->GetEntries() > 0) lPV = iVertices->At(0); 
  FourVectorM lTotVec(0,0,0,0); double lTotSumEt = 0;
  FourVectorM lVisVec(0,0,0,0); double lVisSumEt = 0;
  std::vector<std::pair<FourVectorM,FourVectorM> > lDecay; std::vector<int> lId;
  fNMuons = 0;
  for(UInt_t i0 = 0; i0 < iMuons->GetEntries(); i0++) {
    const Muon *pMu = iMuons->At(i0);
    if(!MetLeptonTools::looseMuId(pMu,iCands,lPV,iVertices)) continue;
    std::pair<FourVectorM,FourVectorM> pVec(pMu->Mom(),pMu->Mom());
    lDecay  .push_back(pVec);
    lId     .push_back(0);
    fNMuons++;
  }
  fNElectrons=0;
  for(UInt_t i0 = 0; i0 < iElectrons->GetEntries(); i0++) {
    const Electron *pElectron = iElectrons->At(i0);
    if(!MetLeptonTools::looseEleId(pElectron,iPUEnergyDensity,iCands,lPV,iVertices)) continue;
    std::pair<FourVectorM,FourVectorM> pVec(pElectron->Mom(),pElectron->Mom());
    lDecay  .push_back(pVec);
    lId     .push_back(1);
    fNElectrons++;
  }
  fNTaus = 0;
  for(UInt_t i0 = 0; i0 < iTaus->GetEntries(); i0++) {
    const PFTau *pTau = iTaus->At(i0);
    if(!fMetLeptonTools->looseTauId(pTau,iPUEnergyDensity)) continue;
    FourVectorM pVis(0,0,0,0); pVis.SetCoordinates(pTau->Pt()*MetLeptonTools::vis(pTau),pTau->Eta(),pTau->Phi(),pTau->Mass()*MetLeptonTools::vis(pTau));
    std::pair<FourVectorM,FourVectorM> pVec(pTau->Mom(),pVis);
    lDecay  .push_back(pVec);
    lId     .push_back(2);
    fNTaus++;
  }
  std::vector<std::pair<FourVectorM,FourVectorM> > lFinalDecay;
  std::vector<int>                                 lFinalId;
  for(unsigned int i0 = 0; i0 < lDecay.size(); i0++) { 
    bool pAdd = true;
    for(unsigned int i1 = 0; i1 < lDecay.size(); i1++) { 
      if(i0 == i1) continue;
      if(MathUtils::DeltaR(lDecay[i0].first,lDecay[i1].first) < 0.5)                                                     pAdd = false;
      if(!pAdd  &&   lId[i0] != 2 && lId[i1] == 2)                                                                       pAdd = true;
      if(!pAdd  &&  ((lId[i0] != 2 && lId[i1] != 2) || (lId[i0] == 2 && lId[i1] == 2))
	        &&                  lDecay[i0].first.pt() >  lDecay[i1].first.pt())                                     pAdd = true;
      if(MathUtils::DeltaR(lDecay[i0].first,lDecay[i1].first) < 0.5 && lDecay[i0].first.pt() == lDecay[i1].first.pt()) { pAdd = true;
	for(unsigned int i2 = 0; i2 < lFinalDecay.size(); i2++) if(fabs(lFinalDecay[i2].first.pt() - lDecay[i0].first.pt()) < 0.1) pAdd = false; 
      }
      if(!pAdd) break;
      //if(!pAdd && lId[i2] == 2 && lId[i1] != 2) pAdd = true;
    }
    if(pAdd) lFinalDecay.push_back(lDecay[i0]);
    if(pAdd) lFinalId   .push_back(lId   [i0]);
  }
  for(unsigned int i0 = 0; i0 < lFinalDecay.size(); i0++) { 
    lTotVec   += lFinalDecay[i0].first;
    lVisVec   += lFinalDecay[i0].second;
    lTotSumEt += lFinalDecay[i0].first.pt();
    lVisSumEt += lFinalDecay[i0].second.pt();
  }
  PFJetOArr *lCleanJets = new PFJetOArr();
  for(UInt_t i0 = 0; i0 < iJets->GetEntries(); i0++) {
    const PFJet *pJet = iJets->At(i0);
    bool pClean = false;
    for(unsigned int i1 = 0; i1 < lFinalDecay.size(); i1++) {
      if(MathUtils::DeltaR(pJet->Mom(),lFinalDecay[i1].first) < 0.5) pClean = true;
    }
    if(!pClean) lCleanJets->Add(pJet);
  }
  //for(unsigned int i0 = 0; i0 < lFinalDecay.size(); i0++) std::cout << "----> " << lFinalDecay[i0].first.pt() << " -- " << lFinalDecay[i0].first.phi() << " -- " << lFinalDecay[i0].first.eta() << " -- " << lFinalId[i0] << " -- " << MathUtils::DeltaR(lFinalDecay[i0].first,lFinalDecay[0].first)<< std::endl;
  Met lMVAMet = GetMet(  false,
			 lTotVec.Pt(),lTotVec.Phi(),lTotSumEt,
			 lVisVec.Pt(),lVisVec.Phi(),lVisSumEt,
			 iPFMet->At(0),
			 iCands,lPV,iVertices,
			 lCleanJets,
			 iJetCorrector,
			 iPUEnergyDensity,
			 int(iVertices->GetEntries()));//,true);  
  return lMVAMet;
}
//--------------------------------------------------------------------------------------------------
Met MVAMet::GetMet(const PhotonCol *iPhotons,
                   const PFCandidateCol *iCands,const PFJetCol  *iJets,const Vertex *iPV,const VertexCol *iVertices,const PFMetCol *iPFMet,
                   FactorizedJetCorrector *iJetCorrector,const PileupEnergyDensityCol* iPUEnergyDensity) {
  const Vertex *lPV = iPV; if(iPV == 0 && iVertices->GetEntries() > 0) lPV = iVertices->At(0); 
  FourVectorM lTotVec(0,0,0,0); double lTotSumEt = 0;
  FourVectorM lVisVec(0,0,0,0); double lVisSumEt = 0;
  std::vector<std::pair<FourVectorM,FourVectorM> > lDecay; std::vector<int> lId;
  fNPhotons = 0;
  for(UInt_t i0 = 0; i0 < iPhotons->GetEntries(); i0++) {
    const Photon *pPhoton = iPhotons->At(i0);
    if(!MetLeptonTools::loosePhotonId(pPhoton)) continue;
    FourVectorM pVis(0,0,0,0); pVis.SetCoordinates(0.,pPhoton->Eta(),pPhoton->Phi(),0);
    std::pair<FourVectorM,FourVectorM> pVec(pPhoton->Mom(),pVis);
    lDecay  .push_back(pVec);
    lId     .push_back(0);
    fNPhotons++;
  }
  std::vector<std::pair<FourVectorM,FourVectorM> > lFinalDecay;
  std::vector<int>                                 lFinalId;
  for(unsigned int i0 = 0; i0 < lDecay.size(); i0++) { 
    bool pAdd = true;
    for(unsigned int i1 = 0; i1 < lDecay.size(); i1++) { 
      if(i0 == i1) continue;
      if(MathUtils::DeltaR(lDecay[i0].first,lDecay[i1].first) < 0.5)                                                     pAdd = false;
      if(!pAdd  &&   lId[i0] != 2 && lId[i1] == 2)                                                                       pAdd = true;
      if(!pAdd  &&  ((lId[i0] != 2 && lId[i1] != 2) || (lId[i0] == 2 && lId[i1] == 2))
	        &&                  lDecay[i0].first.pt() >  lDecay[i1].first.pt())                                     pAdd = true;
      if(MathUtils::DeltaR(lDecay[i0].first,lDecay[i1].first) < 0.5 && lDecay[i0].first.pt() == lDecay[i1].first.pt()) { pAdd = true;
	for(unsigned int i2 = 0; i2 < lFinalDecay.size(); i2++) if(fabs(lFinalDecay[i2].first.pt() - lDecay[i0].first.pt()) < 0.1) pAdd = false; 
      }
      if(!pAdd) break;
      //if(!pAdd && lId[i2] == 2 && lId[i1] != 2) pAdd = true;
    }
    if(pAdd) lFinalDecay.push_back(lDecay[i0]);
    if(pAdd) lFinalId   .push_back(lId   [i0]);
  }
  for(unsigned int i0 = 0; i0 < lFinalDecay.size(); i0++) { 
    lTotVec   += lFinalDecay[i0].first;
    lVisVec   += lFinalDecay[i0].second;
    lTotSumEt += lFinalDecay[i0].first.pt();
    lVisSumEt += lFinalDecay[i0].second.pt();
  }
  PFJetOArr *lCleanJets = new PFJetOArr();
  for(UInt_t i0 = 0; i0 < iJets->GetEntries(); i0++) {
    const PFJet *pJet = iJets->At(i0);
    bool pClean = false;
    for(unsigned int i1 = 0; i1 < lFinalDecay.size(); i1++) {
      if(MathUtils::DeltaR(pJet->Mom(),lFinalDecay[i1].first) < 0.5) pClean = true;
    }
    if(!pClean) lCleanJets->Add(pJet);
  }
  //for(unsigned int i0 = 0; i0 < lFinalDecay.size(); i0++) std::cout << "----> " << lFinalDecay[i0].first.pt() << " -- " << lFinalDecay[i0].first.phi() << " -- " << lFinalDecay[i0].first.eta() << " -- " << lFinalId[i0] << " -- " << MathUtils::DeltaR(lFinalDecay[i0].first,lFinalDecay[0].first)<< std::endl;
  Met lMVAMet = GetMet(  false,
			 lTotVec.Pt(),lTotVec.Phi(),lTotSumEt,
			 lVisVec.Pt(),lVisVec.Phi(),lVisSumEt,
			 iPFMet->At(0),
			 iCands,lPV,iVertices,
			 lCleanJets,
			 iJetCorrector,
			 iPUEnergyDensity,
			 int(iVertices->GetEntries()));//,true);  
  return lMVAMet;
}
//--------------------------------------------------------------------------------------------------
Met MVAMet::GetMet(const MuonCol *iMuons,const ElectronCol *iElectrons,const PFTauCol *iTaus,const PhotonCol *iPhotons,
                   const PFCandidateCol *iCands,const PFJetCol  *iJets,const Vertex *iPV,const VertexCol *iVertices,const PFMetCol *iPFMet,
                   FactorizedJetCorrector *iJetCorrector,const PileupEnergyDensityCol* iPUEnergyDensity, Bool_t printDebug) {
  const Vertex *lPV = iPV; if(iPV == 0 && iVertices->GetEntries() > 0) lPV = iVertices->At(0); 
  FourVectorM lTotVec(0,0,0,0); double lTotSumEt = 0;
  FourVectorM lVisVec(0,0,0,0); double lVisSumEt = 0;
  std::vector<std::pair<FourVectorM,FourVectorM> > lDecay; std::vector<int> lId;
  fNMuons = 0;
  for(UInt_t i0 = 0; i0 < iMuons->GetEntries(); i0++) {
    const Muon *pMu = iMuons->At(i0);
    if(!MetLeptonTools::looseMuId(pMu,iCands,lPV,iVertices)) continue;
    std::pair<FourVectorM,FourVectorM> pVec(pMu->Mom(),pMu->Mom());
    lDecay  .push_back(pVec);
    lId     .push_back(0);
    fNMuons++;
  }
  fNElectrons=0;
  for(UInt_t i0 = 0; i0 < iElectrons->GetEntries(); i0++) {
    const Electron *pElectron = iElectrons->At(i0);
    if(!MetLeptonTools::looseEleId(pElectron,iPUEnergyDensity,iCands,lPV,iVertices)) continue;
    std::pair<FourVectorM,FourVectorM> pVec(pElectron->Mom(),pElectron->Mom());
    lDecay  .push_back(pVec);
    lId     .push_back(1);
    fNElectrons++;
  }
  fNTaus = 0;
  for(UInt_t i0 = 0; i0 < iTaus->GetEntries(); i0++) {
    const PFTau *pTau = iTaus->At(i0);
    if(!fMetLeptonTools->looseTauId(pTau,iPUEnergyDensity)) continue;
    FourVectorM pVis(0,0,0,0); pVis.SetCoordinates(pTau->Pt()*MetLeptonTools::vis(pTau),pTau->Eta(),pTau->Phi(),pTau->Mass()*MetLeptonTools::vis(pTau));
    std::pair<FourVectorM,FourVectorM> pVec(pTau->Mom(),pVis);
    lDecay  .push_back(pVec);
    lId     .push_back(2);
    fNTaus++;
  }
  fNPhotons = 0;
  for(UInt_t i0 = 0; i0 < iPhotons->GetEntries(); i0++) {
    const Photon *pPhoton = iPhotons->At(i0);
    if(!MetLeptonTools::loosePhotonId(pPhoton)) continue;
    Float_t photonChargedFrac = MetLeptonTools::chargedFracInCone(pPhoton,iCands,lPV);
    FourVectorM pVis(0,0,0,0); pVis.SetCoordinates(pPhoton->Pt()*photonChargedFrac,pPhoton->Eta(),pPhoton->Phi(),0);
    std::pair<FourVectorM,FourVectorM> pVec(pPhoton->Mom(),pVis);
    lDecay  .push_back(pVec);
    lId     .push_back(3);
    fNPhotons++;
  }
  std::vector<std::pair<FourVectorM,FourVectorM> > lFinalDecay;
  std::vector<int>                                 lFinalId;
  bool hasPhotons = false;
  for(unsigned int i0 = 0; i0 < lDecay.size(); i0++) { 
    bool pAdd = true;
    for(unsigned int i1 = 0; i1 < lDecay.size(); i1++) { 
      if(i0 == i1) continue;
      if(MathUtils::DeltaR(lDecay[i0].first,lDecay[i1].first) < 0.5)                                                     pAdd = false;
      if(!pAdd  &&   lId[i0] != 2 && lId[i1] == 2)                                                                       pAdd = true;
      if(!pAdd  &&  ((lId[i0] != 2 && lId[i1] != 2) || (lId[i0] == 2 && lId[i1] == 2))
	        &&                  lDecay[i0].first.pt() >  lDecay[i1].first.pt())                                     pAdd = true;
      if(MathUtils::DeltaR(lDecay[i0].first,lDecay[i1].first) < 0.5 && lDecay[i0].first.pt() == lDecay[i1].first.pt()) { pAdd = true;
	for(unsigned int i2 = 0; i2 < lFinalDecay.size(); i2++) if(fabs(lFinalDecay[i2].first.pt() - lDecay[i0].first.pt()) < 0.1) pAdd = false; 
      }
      if(!pAdd) break;
      //if(!pAdd && lId[i2] == 2 && lId[i1] != 2) pAdd = true;
    }
    if(pAdd) lFinalDecay.push_back(lDecay[i0]);
    if(pAdd) lFinalId   .push_back(lId   [i0]);
    // check of "lepton" is a photon
    if (lId   [i0] == 3)
      hasPhotons = true;  
  }
  for(unsigned int i0 = 0; i0 < lFinalDecay.size(); i0++) { 
    lTotVec   += lFinalDecay[i0].first;
    lVisVec   += lFinalDecay[i0].second;
    lTotSumEt += lFinalDecay[i0].first.pt();
    lVisSumEt += lFinalDecay[i0].second.pt();
  }
  PFJetOArr *lCleanJets = new PFJetOArr();
  for(UInt_t i0 = 0; i0 < iJets->GetEntries(); i0++) {
    const PFJet *pJet = iJets->At(i0);
    bool pClean = false;
    for(unsigned int i1 = 0; i1 < lFinalDecay.size(); i1++) {
      if(MathUtils::DeltaR(pJet->Mom(),lFinalDecay[i1].first) < 0.5) pClean = true;
    }
    if(!pClean) lCleanJets->Add(pJet);
  }
  //for(unsigned int i0 = 0; i0 < lFinalDecay.size(); i0++) std::cout << "----> " << lFinalDecay[i0].first.pt() << " -- " << lFinalDecay[i0].first.phi() << " -- " << lFinalDecay[i0].first.eta() << " -- " << lFinalId[i0] << " -- " << MathUtils::DeltaR(lFinalDecay[i0].first,lFinalDecay[0].first)<< std::endl;
  Met lMVAMet = GetMet(  false,
			 lTotVec.Pt(),lTotVec.Phi(),lTotSumEt,
			 lVisVec.Pt(),lVisVec.Phi(),lVisSumEt,
			 iPFMet->At(0),
			 iCands,lPV,iVertices,
			 lCleanJets,
			 iJetCorrector,
			 iPUEnergyDensity,
			 int(iVertices->GetEntries()),
       hasPhotons,
       printDebug);  
  return lMVAMet;
}


