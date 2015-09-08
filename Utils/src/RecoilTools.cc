#include "MitPhysics/Utils/interface/RecoilTools.h"
#include "MitPhysics/Utils/interface/JetTools.h"
#include <TFile.h>
#include "TSystem.h"
#include <stdexcept>

ClassImp(mithep::RecoilTools)

using namespace mithep;

//--------------------------------------------------------------------------------------------------
RecoilTools::RecoilTools(TString iJetLowPtMVAFile, TString iJetHighPtMVAFile, TString iCutFile,
			 bool i42, JetIDMVA::MVAType iType, bool iUseRho)
{ 
  if (iJetLowPtMVAFile.Length() == 0 || iJetHighPtMVAFile.Length() == 0 || iCutFile.Length() == 0) {
    //default
    TString dataDir(gSystem->Getenv("MIT_DATA"));
    if (dataDir.Length() == 0)
      throw std::runtime_error("MIT_DATA environment is not set.");

    if (iJetLowPtMVAFile.Length() == 0)
      iJetLowPtMVAFile = dataDir + "/mva_JetID_lowpt.weights.xml";
    if (iJetHighPtMVAFile.Length() == 0)
      iJetHighPtMVAFile = dataDir + "/mva_JetID_highpt.weights.xml";
    if (iCutFile.Length() == 0)
      iCutFile = dataDir + "/JetIDMVA_JetIdParams.py";
  }

  fJetIDMVA = new JetIDMVA();
  fJetIDMVA->Initialize( JetIDMVA::kMET,iJetLowPtMVAFile,iJetHighPtMVAFile,iType,iCutFile);
  f42       = i42;
  fUseRho   = iUseRho;
}

//--------------------------------------------------------------------------------------------------
RecoilTools::~RecoilTools()
{
  delete fJetIDMVA;
}

//--------------------------------------------------------------------------------------------------
bool RecoilTools::filter(const PFJet *iJet, Double_t iPhi1, Double_t iEta1, Double_t iPhi2,
			 Double_t iEta2)
{
  double pDEta1 = iJet->Eta() - iEta1;
  double pDPhi1 = fabs(iJet->Phi() - iPhi1);
  if (pDPhi1 > 2.*TMath::Pi()-pDPhi1)
    pDPhi1 = 2.*TMath::Pi()-pDPhi1;
  double pDR1   = sqrt(pDEta1*pDEta1 + pDPhi1*pDPhi1);
  if (pDR1 < 0.5)
    return false;
  double pDEta2 = iJet->Eta() - iEta2;
  double pDPhi2 = fabs(iJet->Phi() - iPhi2);
  if (pDPhi2 > 2.*TMath::Pi()-pDPhi2)
    pDPhi2 = 2.*TMath::Pi()-pDPhi2;
  double pDR2   = sqrt(pDEta2*pDEta2 + pDPhi2*pDPhi2);
  if (pDR2 < 0.5)
    return false;
  return true;
}
//--------------------------------------------------------------------------------------------------
bool RecoilTools::filter(const PFCandidate *iCand,Double_t iPhi1,Double_t iEta1,Double_t iPhi2,Double_t iEta2) { 
  double pDEta1 = iCand->Eta() - iEta1;
  double pDPhi1 = fabs(iCand->Phi() - iPhi1);
  if (pDPhi1 > 2.*TMath::Pi()-pDPhi1)
    pDPhi1 = 2.*TMath::Pi()-pDPhi1;
  double pDR1   = sqrt(pDEta1*pDEta1 + pDPhi1*pDPhi1);
  if (pDR1 < 0.3)
    return false;
  double pDEta2 = iCand->Eta() - iEta2;
  double pDPhi2 = fabs(iCand->Phi() - iPhi2);
  if (pDPhi2 > 2.*TMath::Pi()-pDPhi2)
    pDPhi2 = 2.*TMath::Pi()-pDPhi2;
  double pDR2   = sqrt(pDEta2*pDEta2 + pDPhi2*pDPhi2);
  if (pDR2 < 0.3)
    return false;
  return true;
}
//--------------------------------------------------------------------------------------------------
Met RecoilTools::pfRecoil(Double_t iVisPt,Double_t iVisPhi,Double_t iVisSumEt,
			  const PFCandidateCol *iCands) { 
  double lSumEt = 0;
  FourVectorM lVec(0,0,0,0);
  for (UInt_t i0 = 0; i0 < iCands->GetEntries(); i0++) { 
    lVec -= iCands->At(i0)->Mom(); 
    lSumEt += iCands->At(i0)->Pt();
  }
  Met lPFMet(lVec.Px(),lVec.Py()); 
  lPFMet.SetSumEt(lSumEt);
  lPFMet.SetMex  (lPFMet.Mex()+iVisPt*cos(iVisPhi));  
  lPFMet.SetMey  (lPFMet.Mey()+iVisPt*sin(iVisPhi));
  lPFMet.SetSumEt(lPFMet.SumEt()-iVisSumEt);
  return lPFMet;
}
//--------------------------------------------------------------------------------------------------
Met RecoilTools::pfRecoil(double iPhi1,double iEta1,double iPhi2,double iEta2,
			  const PFCandidateCol *iCands) { 
  double lSumEt = 0;
  FourVectorM lVec(0,0,0,0);
  for (UInt_t i0 = 0; i0 < iCands->GetEntries(); i0++) { 
    const PFCandidate *pfcand = iCands->At(i0);
    if (!filter(pfcand,iPhi1,iEta1,iPhi2,iEta2))
      continue;
    lVec   -= pfcand->Mom(); 
    lSumEt += pfcand->Pt();
  }
  Met lPFMet(lVec.Px(),lVec.Py()); 
  lPFMet.SetSumEt(lSumEt);
  return lPFMet;

}

//--------------------------------------------------------------------------------------------------
Met RecoilTools::pfCone(double iPhi1, double iEta1, const PFCandidateCol *iCands,
			const Vertex *iVertex, bool iCharge, Double_t iDZCut) { 
  double lSumEt = 0;
  FourVectorM lVec(0,0,0,0);
  for (UInt_t i0 = 0; i0 < iCands->GetEntries(); i0++) { 
    const PFCandidate *pfcand = iCands->At(i0);
    if (filter(pfcand,iPhi1,iEta1,100,100)) continue;
    if (iCharge) { 
      double lDZ = -999;
      if      (pfcand->HasTrackerTrk())
	lDZ =  fabs(pfcand->TrackerTrk()->DzCorrected(*iVertex));//
      else if (pfcand->HasGsfTrk())
	lDZ =  fabs(pfcand->GsfTrk()    ->DzCorrected(*iVertex));
      if (fabs(lDZ) > iDZCut)
	continue;
    }
    lVec   += pfcand->Mom(); 
    lSumEt += pfcand->Pt();
  }
  Met lPFMet(lVec.Px(),lVec.Py()); 
  lPFMet.SetSumEt(lSumEt);

  return lPFMet;
}
//--------------------------------------------------------------------------------------------------
Met RecoilTools::trackMet(const PFCandidateCol *iCands,const Vertex *iVertex,Double_t iDZCut) { 
  double trkMetx  = 0;
  double trkMety  = 0;
  double trkSumEt = 0; 
  for (UInt_t i0=0; i0<iCands->GetEntries(); i0++) {
    const PFCandidate *pfcand = iCands->At(i0);
    double lDZ = -999;
    if      (pfcand->HasTrackerTrk())
      lDZ =  fabs(pfcand->TrackerTrk()->DzCorrected(*iVertex));//
    else if (pfcand->HasGsfTrk())
      lDZ =  fabs(pfcand->GsfTrk()    ->DzCorrected(*iVertex));
    if (fabs(lDZ) > iDZCut)
      continue;
    trkMetx  -= pfcand->Px();
    trkMety  -= pfcand->Py();
    trkSumEt += pfcand->Pt();
  }
  Met lMet(trkMetx,trkMety);
  lMet.SetSumEt(trkSumEt);

  return lMet;
}
//--------------------------------------------------------------------------------------------------
Met RecoilTools::trackRecoil(double iPhi1, double iEta1, double iPhi2, double iEta2,
			     const PFCandidateCol *iCands, const Vertex *iVertex, Double_t iDZCut)
{
  double trkMetx  = 0;
  double trkMety  = 0;
  double trkSumEt = 0; 
  for (UInt_t i0=0; i0<iCands->GetEntries(); i0++) {
    const PFCandidate *pfcand = iCands->At(i0);
    if (!filter(pfcand,iPhi1,iEta1,iPhi2,iEta2))
      continue;
    if ((pfcand->HasTrackerTrk() && (fabs(pfcand->TrackerTrk()->DzCorrected(*iVertex))<iDZCut)) ||
	(pfcand->HasGsfTrk()     && (fabs(pfcand->GsfTrk()->DzCorrected(*iVertex))    <iDZCut))  ) {
      trkMetx  -= pfcand->Px();
      trkMety  -= pfcand->Py();
      trkSumEt += pfcand->Pt();
    }
  }
  Met lMet(trkMetx,trkMety);
  lMet.SetSumEt(trkSumEt);

  return lMet;
}
//--------------------------------------------------------------------------------------------------
//Compute the recoil => here this requires the vector sum of the visible components
//VisPt    => Vector sum pT  of the visible non-recoiling components
//VisPhi   => Vector sum Phi of the visible non-recoiling components
//visSumEt => Vector sum Et  =f the visible non-recoiling components
Met RecoilTools::trackRecoil(Double_t iVisPt,Double_t iVisPhi,Double_t iVisSumEt,
			     const PFCandidateCol *iCands,const Vertex *iVertex,double iDZCut)
{
  Met lTrkMet = trackMet(iCands,iVertex,iDZCut);
  lTrkMet.SetMex  (lTrkMet.Mex()+iVisPt*cos(iVisPhi));  
  lTrkMet.SetMey  (lTrkMet.Mey()+iVisPt*sin(iVisPhi));
  lTrkMet.SetSumEt(lTrkMet.SumEt()-iVisSumEt);

  return lTrkMet;
}

//--------------------------------------------------------------------------------------------------
//Corrected Jets
void RecoilTools::addNeut(const PFJet *iJet,FourVectorM &iVec,Double_t &iSumEt,double iRho,int iSign) { 
  FourVectorM lVec(0,0,0,0);
  double lPt = iJet->Pt();
  //if (iJet->RawMom().Pt() < 10)
  //  lPt = TMath::Max(iJet->RawMom().Pt()-iJet->JetArea()*iRho,0.); //to be fixed
  //if (iJet->RawMom().Pt() < 10)
  //  lPt = iJet->RawMom().Pt()*iJet->L1OffsetCorrectionScale();
  lPt *= (iJet->NeutralEmEnergy()/iJet->RawMom().E() + iJet->NeutralHadronEnergy()/iJet->RawMom().E());
  lVec.SetPt(lPt); lVec.SetEta(iJet->Eta()); lVec.SetPhi(iJet->Phi()); lVec.SetM(iJet->Mass());
  if (iSign > 0)
    iVec   -= lVec;
  if (iSign < 0)
    iVec   += lVec;
  if (iSign > 0)
    iSumEt += lPt;
  if (iSign < 0)
    iSumEt -= lPt;
}

//-------------------------------------------------------------------------------------------------
//Corrected Jets
Met RecoilTools::NoPUMet(const PFJetCol *iJets, const PFCandidateCol *iCands,
			 const Vertex *iVertex, const VertexCol *iVertices, Double_t iRho,
			 Double_t iPhi1, Double_t iEta1, Double_t iPhi2, Double_t iEta2,
			 Double_t iDZCut)
{
  FourVectorM lVec(0,0,0,0);
  double      lSumEt = 0; 
  for (UInt_t i0 = 0; i0 < iCands->GetEntries(); i0++) { 
    const PFCandidate *pPF = iCands->At(i0);
    double lDZ = 999;
    if      (pPF->HasTrackerTrk())
      lDZ =  fabs(pPF->TrackerTrk()->DzCorrected(*iVertex));//
    else if (pPF->HasGsfTrk())
      lDZ =  fabs(pPF->GsfTrk()    ->DzCorrected(*iVertex));
    if (fabs(lDZ) > iDZCut)
      continue;
    lVec     -= pPF->Mom();
    lSumEt   += pPF->Pt();  
  }
  for (UInt_t i0 = 0; i0 < iJets->GetEntries(); i0++) {
    const PFJet *pJet = iJets->At(i0);
    if (!filter(pJet,iPhi1,iEta1,iPhi2,iEta2))
      continue;
    if (!fJetIDMVA->pass(pJet,iVertex,iVertices))
      continue;
    addNeut(pJet,lVec,lSumEt,iRho);
  }
  Met lMet(lVec.Px(),lVec.Py());
  lMet.SetSumEt(lSumEt);

  return lMet;
}

//--------------------------------------------------------------------------------------------------
//Corrected Jets
Met RecoilTools::NoPURecoil(Double_t iVisPt,Double_t iVisPhi,Double_t iVisSumEt,
			    const PFJetCol *iJets,const PFCandidateCol *iCands,
			    const Vertex *iVertex,const VertexCol *iVertices,Double_t iRho,
			    Double_t iPhi1,Double_t iEta1,Double_t iPhi2,Double_t iEta2,
			    Double_t iDZCut)
{
  Met lNoPUMet = NoPUMet(iJets,iCands,iVertex,iVertices,iRho,iPhi1,iEta1,iPhi2,iEta2,iDZCut);
  lNoPUMet.SetMex  (lNoPUMet.Mex()+iVisPt*cos(iVisPhi));  
  lNoPUMet.SetMey  (lNoPUMet.Mey()+iVisPt*sin(iVisPhi));
  lNoPUMet.SetSumEt(lNoPUMet.SumEt()-iVisSumEt);

  return lNoPUMet;
}

//--------------------------------------------------------------------------------------------------
//Corrected jets
Met RecoilTools::PUCMet( const PFJetCol       *iJets,const PFCandidateCol *iCands,
			 const Vertex *iVertex,const VertexCol *iVertices,Double_t iRho,bool iAddType1,
			 Double_t iPhi1,Double_t iEta1,Double_t iPhi2,Double_t iEta2,
			 Double_t iDZCut)
{
  FourVectorM lVec(0,0,0,0);
  double      lSumEt = 0; 
  for (UInt_t i0 = 0; i0 < iCands->GetEntries(); i0++) { 
    const PFCandidate *pPF = iCands->At(i0);
    const Track* pTrack = pPF->TrackerTrk();
    if (pPF->GsfTrk())
      pTrack = pPF->GsfTrk();
    if (pTrack == 0                                   && 
	(pPF->PFType() == PFCandidate::eGamma         || 
	 pPF->PFType() == PFCandidate::eEGammaHF      || 
	 pPF->PFType() == PFCandidate::eNeutralHadron || 
	 pPF->PFType() == PFCandidate::eHadronHF        )) {
      lVec -= pPF->Mom();
      lSumEt += pPF->Pt();
    }
    if (pTrack == 0)
      continue;
    double lDZ = 999;
    if      (pPF->HasTrackerTrk())
      lDZ =  fabs(pPF->TrackerTrk()->DzCorrected(*iVertex));//
    else if (pPF->HasGsfTrk())
      lDZ =  fabs(pPF->GsfTrk()    ->DzCorrected(*iVertex));
    if ( fabs(lDZ) > iDZCut)
      continue;
    lVec     -= pPF->Mom();
    lSumEt   += pPF->Pt();
  }

  //FIXME Not Yet Supported
  //if (iAddType1)
  //  addType1(lVec,lSumEt,iJets,iJetCorrector,iPUEnergyDensity,iPhi1,iEta1,iPhi2,iEta2);

  for (UInt_t i0 = 0; i0 < iJets->GetEntries(); i0++) {
    const PFJet *pJet = iJets->At(i0);
    if (!JetTools::passPFId(pJet, JetTools::kPFLoose))
      continue;
    if (fJetIDMVA->pass(pJet,iVertex,iVertices))
      continue;
    if (!filter(pJet,iPhi1,iEta1,iPhi2,iEta2))
      continue; // Quick cleaning ==> if not done already
    addNeut(pJet,lVec,lSumEt,iRho,-1);
  }
  Met lMet(lVec.Px(),lVec.Py());
  lMet.SetSumEt(lSumEt);

  return lMet;
}

//--------------------------------------------------------------------------------------------------
//Corrected Jets
Met RecoilTools::PUCRecoil(Double_t iVisPt,Double_t iVisPhi,Double_t iVisSumEt,
			   const PFJetCol       *iJets,
			   const PFCandidateCol *iCands,const Vertex *iVertex,
			   const VertexCol *iVertices,Double_t iRho,bool iAddType1,
			   Double_t iPhi1,Double_t iEta1,Double_t iPhi2,Double_t iEta2,
			   Double_t iDZCut)
{
  Met lPUCMet = PUCMet(iJets,iCands,iVertex,iVertices,iRho,iAddType1,iPhi1,iEta1,iPhi2,iEta2,iDZCut);
  lPUCMet.SetMex  (lPUCMet.Mex()+iVisPt*cos(iVisPhi));  
  lPUCMet.SetMey  (lPUCMet.Mey()+iVisPt*sin(iVisPhi));
  lPUCMet.SetSumEt(lPUCMet.SumEt()-iVisSumEt);

  return lPUCMet;
}

//--------------------------------------------------------------------------------------------------
//Corrected Jets
Met RecoilTools::PUMet(const PFJetCol *iJets, const PFCandidateCol *iCands,const Vertex *iVertex,
		       const VertexCol *iVertices,Double_t iRho, Double_t iPhi1,Double_t iEta1,
		       Double_t iPhi2,Double_t iEta2, Double_t iDZCut)
{ 
  FourVectorM lVec(0,0,0,0);
  double      lSumEt = 0; 
  for (UInt_t i0 = 0; i0 < iCands->GetEntries(); i0++) { 
    const PFCandidate *pPF = iCands->At(i0);
    const Track* pTrack = pPF->TrackerTrk();
    if (pPF->GsfTrk())
      pTrack = pPF->GsfTrk();
    if (pTrack == 0)
      continue;
    double lDZ = 999;
    if (pPF->HasTrackerTrk())
      lDZ =  fabs(pPF->TrackerTrk()->DzCorrected(*iVertex));//
    else if (pPF->HasGsfTrk())
      lDZ =  fabs(pPF->GsfTrk()    ->DzCorrected(*iVertex));
    if ( fabs(lDZ) < iDZCut)
      continue;
    lVec     -= pPF->Mom();
    lSumEt   += pPF->Pt();
  }
  for (UInt_t i0 = 0; i0 < iJets->GetEntries(); i0++) {
    const PFJet *pJet = iJets->At(i0);
    if (!JetTools::passPFId(pJet, JetTools::kPFLoose))
      continue;
    if (!filter(pJet,iPhi1,iEta1,iPhi2,iEta2))
      continue; //Quick cleaning
    if (fJetIDMVA->pass(pJet,iVertex,iVertices))
      continue;
    addNeut(pJet,lVec,lSumEt,iRho);
  }
  Met lMet(lVec.Px(),lVec.Py());
  lMet.SetSumEt(lSumEt);

  return lMet;
}

//--------------------------------------------------------------------------------------------------
//Corrected Jets
Met RecoilTools::PURecoil(Double_t iVisPt, Double_t iVisPhi, Double_t iVisSumEt,
			  const PFJetCol *iJets, const PFCandidateCol *iCands,
			  const Vertex *iVertex, const VertexCol *iVertices, Double_t iRho,
			  Double_t iPhi1,Double_t iEta1, Double_t iPhi2, Double_t iEta2,
			  Double_t iDZCut)
{
  Met lPUMet = PUMet(iJets,iCands,iVertex,iVertices,iRho,iPhi1,iEta1,iPhi2,iEta2,iDZCut);
  lPUMet.SetMex  (lPUMet.Mex()+iVisPt*cos(iVisPhi));  
  lPUMet.SetMey  (lPUMet.Mey()+iVisPt*sin(iVisPhi));
  lPUMet.SetSumEt(lPUMet.SumEt()-iVisSumEt);

  return lPUMet;
}
