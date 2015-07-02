#include "MitPhysics/Utils/interface/MVAMet.h"
#include "MitPhysics/Utils/interface/MetLeptonTools.h"
#include "MitPhysics/Utils/interface/JetTools.h"
#include "MitPhysics/Utils/interface/RecoilTools.h"
#include "MitAna/DataTree/interface/StableData.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/PFTauCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/PhotonCol.h"
#include "MitAna/DataTree/interface/PFMetCol.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityCol.h"
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
