#include "MitPhysics/Utils/interface/MetLeptonTools.h"

#include "MitAna/DataTree/interface/Muon.h"
#include "MitAna/DataTree/interface/PFTau.h"
#include "MitAna/DataTree/interface/Electron.h"
#include "MitAna/DataTree/interface/Photon.h"
#include "MitAna/DataTree/interface/PFMet.h"
#include "MitAna/DataTree/interface/ChargedParticle.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityCol.h"
#include "MitPhysics/Utils/interface/TauIsoMVA.h"
#include "MitCommon/MathTools/interface/MathUtils.h"

#include <algorithm>
#include <vector>
#include "TString.h"
#include "TSystem.h"

using namespace mithep;

ClassImp(mithep::MetLeptonTools)

MetLeptonTools::MetLeptonTools() { 
  TString dataDir(gSystem->Getenv("MIT_DATA"));
  if (dataDir.Length() == 0)
    throw std::runtime_error("MIT_DATA environment is not set.");

  //fTauIsoMVA = new TauIsoMVA();
  //fTauIsoMVA->Initialize(dataDir + "/SXIsoMVA_BDTG.weights.xml");
  fTauIsoMVA = new TauIsoMVA();
  fTauIsoMVA->InitializeGBR(dataDir + "/gbrfTauIso_v2.root");
}
double MetLeptonTools::vis(const PFTau *iTau) {
  double lPtTot        = 0.;
  double lChargedPtTot = 0.;
  for(unsigned int i0 = 0; i0 < iTau->NSignalPFCands(); i0++) {
    lPtTot       += iTau->SignalPFCand(i0)->Pt();
    if(iTau->SignalPFCand(i0)->BestTrk() == 0) continue;
    lChargedPtTot += iTau->SignalPFCand(i0)->Pt();
  }
  return lChargedPtTot/lPtTot;
}
Float_t MetLeptonTools::PFIsolation(const ChargedParticle *iLep,const PFCandidateCol *iCands) {
  Double_t lPtSum = 0.;
  for(UInt_t i0 = 0; i0 < iCands->GetEntries(); i0++) {
    const PFCandidate *pCand = iCands->At(i0);
     Double_t pDR = MathUtils::DeltaR(iLep->Mom(), pCand->Mom());
     if(pDR < 0.0001) continue;
     if(pDR < 0.01 && (pCand->PFType() != PFCandidate::eNeutralHadron || pCand->PFType() != PFCandidate::eGamma)) continue;
     if(pCand->Pt() < 0.5 && (pCand->PFType() != PFCandidate::eNeutralHadron || pCand->PFType() != PFCandidate::eGamma)) continue;
     if(pCand->PFType() != PFCandidate::eNeutralHadron && pCand->PFType() != PFCandidate::eGamma && pCand->PFType() != PFCandidate::eHadron ) continue;
     if(pDR         > 0.3)                                 continue;
    lPtSum += pCand->Pt();
  }
  return lPtSum/iLep->Pt();
}
Float_t MetLeptonTools::PFIsolationNoGamma(const ChargedParticle *iLep,const PFCandidateCol *iCands) {
  Double_t lPtSum = 0.;
  for(UInt_t i0 = 0; i0 < iCands->GetEntries(); i0++) {
    const PFCandidate *pCand = iCands->At(i0);
     Double_t pDR = MathUtils::DeltaR(iLep->Mom(), pCand->Mom());
     if(pDR < 0.0001) continue;
     if(pDR < 0.01 && (pCand->PFType() != PFCandidate::eNeutralHadron)) continue;
     if(pCand->Pt() < 0.5 && (pCand->PFType() != PFCandidate::eNeutralHadron)) continue;
     if(pCand->PFType() != PFCandidate::eNeutralHadron && pCand->PFType() != PFCandidate::eHadron ) continue;
     if(pDR         > 0.3)                                 continue;
    lPtSum += pCand->Pt();
  }
  return lPtSum/iLep->Pt();
}
Float_t MetLeptonTools::isoPV(const ChargedParticle *iLep,const PFCandidateCol *iCands,
			      const Vertex *iPV,const VertexCol *iVertices,bool iEle) {
  Float_t lPtSumCharge = 0.;
  for(UInt_t i0 = 0; i0 < iCands->GetEntries(); i0++) {
    const PFCandidate *pCand = iCands->At(i0);
    Double_t pDR = MathUtils::DeltaR(iLep->Mom(), pCand->Mom());
    if(pCand->PFType() != PFCandidate::eHadron) continue;
    if(pDR         > 0.3)                       continue;
    if(pDR         < 0.0001)                    continue;
    if(iEle && pDR         < 0.015 && fabs(iLep->Eta()) > 1.56)   continue;
    if(iEle && pDR         < 0.01  && fabs(iLep->Eta()) < 1.56)   continue;
    if(pCand->HasTrackerTrk() && iPV !=0) {
      if( iPV->HasTrack(pCand->TrackerTrk())) lPtSumCharge += pCand->Pt();
      //if( iPV->HasTrack(pCand->TrackerTrk())) std::cout << "===> Adding ===> " << pCand->Pt() << " --" << pDR << " -- " << pCand->BestTrk()->DzCorrected(*iPV) << std::endl;
      if( iPV->HasTrack(pCand->TrackerTrk())) continue;
    } 
    Double_t pDzMin = 10000;
    Bool_t pVertexFound  = kFALSE;
    const Vertex *pClosestVtx  = 0;
    for(UInt_t i1 = 0; i1 < iVertices->GetEntries(); i1++) {
      const Vertex *pVtx = iVertices->At(i1);
      if(pVtx->HasTrack(pCand->TrackerTrk())) { 
	pClosestVtx  = pVtx; pVertexFound  = kTRUE;  break; 
      }
      Double_t pDz = fabs(pCand->SourceVertex().Z() - pVtx->Z());
      if(pDz < pDzMin) {
	pClosestVtx = pVtx;
	pDzMin = pDz;
      }
    }
    if(pVertexFound  || pClosestVtx != iPV) continue;
    //std::cout << "===> Adding NV ===> " << pCand->Pt() << " --" << pDR << " -- " << pCand->BestTrk()->DzCorrected(*iPV) << std::endl;
    lPtSumCharge += pCand->Pt();
  }
  return lPtSumCharge;// + TMath::Max(lPtSumNeut-lPtSumPU*0.5,0.))/iLep->Pt();
}

// ------ ported from -------------------------------------------------------- // 
// ------ JetMETCorrections-METPUSubtraction/plugins/PFMETProducerMVA.cc ----- // 
Float_t MetLeptonTools::chargedFracInCone(const Photon *iPhoton,const PFCandidateCol *iCands,
                                          const Vertex *iPV) {
  FourVectorM lVis(0,0,0,0);
  for(UInt_t i0 = 0; i0 < iCands->GetEntries(); i0++) {
    const PFCandidate *pCand = iCands->At(i0);
    Double_t pDR = MathUtils::DeltaR(iPhoton->Mom(), pCand->Mom());
    if(pDR         > 0.2)                       continue;
    Double_t pDz = -999; // PH: If no vertex is reconstructed in the event
                         // or PFCandidate has no track, set dZ to -999
    if(iPV != 0 && pCand->HasTrk()) 
      pDz = TMath::Abs(pCand->BestTrk()->DzCorrected(*iPV)); 
    
    if(TMath::Abs(pDz) > 0.1) continue;
    lVis += pCand->Mom();
  }
  
  return lVis.Pt()/iPhoton->Pt();
}
