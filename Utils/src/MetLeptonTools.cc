#include "MitPhysics/Utils/interface/MetLeptonTools.h"
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

bool MetLeptonTools::looseTauId(const PFTau *iTau,const PileupEnergyDensityCol* iPUEnergyDensity) {
  if(iTau->Pt() < 19)                                                 return false;
  if(fabs(iTau->Eta()) > fabs(2.3) )                                  return false;
  if(!iTau->DiscriminationByDecayModeFinding())                       return false;
  if(!iTau->DiscriminationByLooseElectronRejection())                 return false;
  if(!iTau->LooseCombinedIsolationDBSumPtCorr3Hits() < 2.0)           return false;
  if(!iTau->LooseMuonRejection2())                                    return false;
  return true;
}
bool MetLeptonTools::looseEleId(const Electron *iElectron,const PileupEnergyDensityCol* iPUEnergyDensity,
				const PFCandidateCol *iCands,const Vertex *iPV,const VertexCol *iVertices) {
  if(iElectron->SCluster()  == 0)    return false;  
  if(iElectron->Pt()         < 9.5)  return false;
  if(fabs(iElectron->Eta())  > 2.5)  return false;
  if(fabs(iElectron->Eta()) > 1.4442 && fabs(iElectron->Eta()) < 1.566) return false;  
  if(iElectron->GsfTrk() == 0) return false;
  if(iElectron->GsfTrk()) if(iElectron->GsfTrk()->NExpectedHitsInner() > 0) return false;
  if(PFIsolationNoGamma(iElectron,iCands)            > 0.3) return false;
  if(iElectron->TrackIsolationDr03()/iElectron->Et() > 0.2) return false;        
  //Electron Veto Id 
  if(fabs(iElectron->Eta()) < 1.5) { 
    if(fabs(iElectron->DeltaEtaSuperClusterTrackAtVtx()) > 0.007)   return false;
    if(fabs(iElectron->DeltaPhiSuperClusterTrackAtVtx()) > 0.8)     return false;
    if(iElectron->CoviEtaiEta()                          > 0.01)    return false;
    if(iElectron->HadronicOverEm()                       > 0.15)    return false;
    double lE = iElectron->SCluster()->Energy();
    double lP = iElectron->P();
    if(fabs(1./lE-1./lP)                                 > 0.05 )   return false;
  } else { 
    if(fabs(iElectron->DeltaEtaSuperClusterTrackAtVtx()) > 0.009)   return false;
    if(fabs(iElectron->DeltaPhiSuperClusterTrackAtVtx()) > 0.10)    return false;
    if(iElectron->CoviEtaiEta()                          > 0.03)    return false;
    if(iElectron->HadronicOverEm()                       > 0.10)    return false;
    double lE = iElectron->SCluster()->Energy();
    double lP = iElectron->P();
    if(fabs(1./lE-1./lP)                                 > 0.05 )   return false;
  }
  return true;
}
bool MetLeptonTools::looseMuId(const Muon *iMu,const PFCandidateCol *iCands,const Vertex *iPV,const VertexCol *iVertices) {
  if(iMu->TrackerTrk()                         == 0  )          return false;
  if(iMu->Pt()                                  < 9.5)          return false;
  if(fabs(iMu->BestTrk()->Eta())                > 2.5)          return false;
  if(iMu->BestTrk()->D0()                       > 2.0)          return false;
  if(iMu->BestTrk()->RChi2()                    > 10 )          return false;
  if(iMu->TrackerTrk()->NPixelHits()            < 1  )          return false;
  if(iMu->TrackerTrk()->NHits()                 < 6  )          return false;
  if(iMu->NValidHits()                          < 1  )          return false;
  if(iMu->NMatches()                            < 1  )          return false;
  if(PFIsolation(iMu,iCands)                    > 0.3)          return false;
  return true;
}
bool MetLeptonTools::loosePhotonId(const Photon *iPhoton) { 
  if(iPhoton->Pt()        <   20)                                   return false;
  if(fabs(iPhoton->Eta()) > 1.45)                                   return false; 
  if(iPhoton->HadOverEm() > 0.05)                                   return false;
  if(iPhoton->HadOverEm() > 0.05)                                   return false;
  if(iPhoton->HasPixelSeed())                                       return false;
  if(iPhoton->CoviEtaiEta() > 0.01)                                 return false;
  if(iPhoton->HollowConeTrkIsoDr04()    > (2.+0.002*iPhoton->Pt())) return false;
  return true;
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
