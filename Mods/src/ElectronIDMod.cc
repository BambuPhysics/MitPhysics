#include "MitPhysics/Mods/interface/ElectronIDMod.h"
#include "MitAna/DataTree/interface/StableData.h"
#include "MitAna/DataTree/interface/ElectronFwd.h"
#include "MitAna/DataTree/interface/MuonFwd.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/DataTree/interface/TriggerObjectCol.h"
#include "MitAna/DataTree/interface/DecayParticleCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

using namespace mithep;

ClassImp(mithep::ElectronIDMod)

//--------------------------------------------------------------------------------------------------
ElectronIDMod::ElectronIDMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fPrintMVADebugInfo(kFALSE),
  fElectronBranchName(Names::gkElectronBrn),
  fConversionBranchName(Names::gkMvfConversionBrn),
  fGoodElectronsName(ModNames::gkGoodElectronsName),  
  fNonIsolatedMuonsName("random"),  
  fNonIsolatedElectronsName("random"),  
  fVertexName(ModNames::gkGoodVertexesName),
  fBeamSpotName(Names::gkBeamSpotBrn),
  fTrackName(Names::gkTrackBrn),
  fPFCandidatesName(Names::gkPFCandidatesBrn),
  fPFNoPileUpName("PFNoPileUp"),
  fElectronIDType("CustomTight"),
  fElectronIsoType("PFIso"),
  fTrigObjectsName("HLTModTrigObjs"),
  fElectronPtMin(10),
  fElectronEtMin(0.0),  
  fElectronEtaMax(2.5),
  fIDLikelihoodCut(-999.0),
  fTrackIsolationCut(5.0),
  fCaloIsolationCut(5.0),
  fEcalJuraIsoCut(5.0),
  fHcalIsolationCut(5.0),
  fCombIsolationCut(0.1),
  fCombRelativeIsolationCut(0.10),
  fCombRelativeIsolationCut_EE(0.10),
  fPFIsolationCut(-1.0),
  fApplyConvFilterType1(kTRUE),
  fApplyConvFilterType2(kFALSE),
  fNWrongHitsMax(0),
  fNExpectedHitsInnerCut(999),
  fInvertNExpectedHitsInnerCut(kFALSE),
  fCombinedIdCut(kFALSE),
  fApplySpikeRemoval(kTRUE),
  fApplyD0Cut(kTRUE),
  fApplyDZCut(kTRUE),
  fChargeFilter(kTRUE),
  fD0Cut(0.020),
  fDZCut(0.10),
  fWhichVertex(-1),
  fApplyTriggerMatching(kFALSE),
  fApplyEcalSeeded(kFALSE),
  fApplyCombinedIso(kTRUE),
  fApplyEcalFiducial(kFALSE),
  fElectronsFromBranch(kTRUE),
  fElIdType(ElectronTools::kIdUndef),
  fElIsoType(ElectronTools::kIsoUndef),
  fElectrons(0),
  fConversions(0),
  fVertices(0),
  fBeamSpot(0),
  fTracks(0),
  fPFCandidates(0),
  fPFNoPileUpCands(0),
  fIntRadius(0.0),
  fNonIsolatedMuons(0),
  fNonIsolatedElectrons(0),
  fLH(0),
  fPileupEnergyDensityName(Names::gkPileupEnergyDensityBrn),
  fPileupEnergyDensity(0),
  fElectronIDMVA(0),
  fElectronMVAWeights_Subdet0Pt10To20(""),
  fElectronMVAWeights_Subdet1Pt10To20(""),
  fElectronMVAWeights_Subdet2Pt10To20(""),
  fElectronMVAWeights_Subdet0Pt20ToInf(""),
  fElectronMVAWeights_Subdet1Pt20ToInf(""),
  fElectronMVAWeights_Subdet2Pt20ToInf(""),
  fPVName(Names::gkPVBeamSpotBrn),
  fTheRhoType(RhoUtilities::DEFAULT)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
Bool_t ElectronIDMod::PassLikelihoodID(const Electron *ele) const
{

  Double_t LikelihoodValue = ElectronTools::Likelihood(fLH, ele);

  double likCut = fIDLikelihoodCut;
  if (likCut > -900) {
    if (ele->Pt() > 20) {
      if (ele->SCluster()->AbsEta() < 1.479) {
        if (ele->NumberOfClusters() - 1 == 0) likCut = 3.5;
	else                                 likCut = 4.0;
      }
      else  {                                
        if (ele->NumberOfClusters() - 1 == 0) likCut = 4.0;
	else                                 likCut = 4.0;
      }
    }
    else {
      if (ele->SCluster()->AbsEta() < 1.479) {
        if (ele->NumberOfClusters() - 1 == 0) likCut =  4.0;
	else                                 likCut =  4.5;
      }
      else  {                                
        if (ele->NumberOfClusters() - 1 == 0) likCut =  4.0;
	else                                 likCut =  4.0;
      }
    }
  }
  if (LikelihoodValue > likCut) return kTRUE;
  return kFALSE;
}

//--------------------------------------------------------------------------------------------------
Bool_t ElectronIDMod::PassMVAID(const Electron *el, ElectronTools::EElIdType idType, 
                                const Vertex *vertex, const PFCandidateCol *PFCands,
                                const PileupEnergyDensityCol *PileupEnergyDensity) const
{ 
  Bool_t isDebug = kFALSE;
  Double_t MVAValue = 0;
  if     (idType == ElectronTools::kMVAID_BDTG_IDIsoCombined) {
    MVAValue = fElectronIDMVA->MVAValue(el, vertex, PFCands, PileupEnergyDensity, fIntRadius);
  } 
  else if (idType == ElectronTools::kMVAID_BDTG_IDHWW2012TrigV0) {
    ElectronOArr *tempElectrons = new  ElectronOArr;
    MuonOArr     *tempMuons     = new  MuonOArr;
    MVAValue = fElectronIDMVA->MVAValue(el, vertex, PFCands, PileupEnergyDensity, ElectronTools::kEleEANoCorr, 
                                        tempElectrons, tempMuons, isDebug);
    delete tempElectrons;
    delete tempMuons;
  }
  else if (idType == ElectronTools::kMVAID_BDTG_IDIsoCombinedHWW2012TrigV4) {
    MVAValue = fElectronIDMVA->MVAValue(el, vertex, fVertices, PFCands, 
                                        PileupEnergyDensity, ElectronTools::kEleEANoCorr, 
                                        isDebug);
  }
  else {
    MVAValue = fElectronIDMVA->MVAValue(el, vertex);
  }
  
  Double_t eta = el->SCluster()->AbsEta();
  Int_t MVABin = -1;
  if (idType == ElectronTools::kMVAID_BDTG_IDHWW2012TrigV0
     || idType == ElectronTools::kMVAID_BDTG_IDIsoCombinedHWW2012TrigV4
     || idType == ElectronTools::kHggLeptonTagId2012HCP
   ) {
     if (el->Pt() <  20 && eta <  0.800 		    ) MVABin = 0;
     if (el->Pt() <  20 && eta >= 0.800 && fabs(eta) < 1.479) MVABin = 1;
     if (el->Pt() <  20 && eta >= 1.479 		    ) MVABin = 2;
     if (el->Pt() >= 20 && eta <  0.800 		    ) MVABin = 3;
     if (el->Pt() >= 20 && eta >= 0.800 && fabs(eta) < 1.479) MVABin = 4;
     if (el->Pt() >= 20 && eta >= 1.479 		    ) MVABin = 5;
  } else {
     if (el->Pt() <  20 && eta <  1.000 		    ) MVABin = 0;
     if (el->Pt() <  20 && eta >= 1.000 && fabs(eta) < 1.479) MVABin = 1;
     if (el->Pt() <  20 && eta >= 1.479 		    ) MVABin = 2;
     if (el->Pt() >= 20 && eta <  1.000 		    ) MVABin = 3;
     if (el->Pt() >= 20 && eta >= 1.000 && fabs(eta) < 1.479) MVABin = 4;
     if (el->Pt() >= 20 && eta >= 1.479 		    ) MVABin = 5;  
  }
  if (MVABin == -1) assert(0);

  Double_t MVACut = -9999;
  if (idType == ElectronTools::kMVAID_BDTG_NoIPInfo) {
    if      (MVABin == 0) MVACut = 0.133;
    else if (MVABin == 1) MVACut = 0.465;
    else if (MVABin == 2) MVACut = 0.518; 
    else if (MVABin == 3) MVACut = 0.942;
    else if (MVABin == 4) MVACut = 0.947;
    else if (MVABin == 5) MVACut = 0.878 ;
  } else if (idType == ElectronTools::kMVAID_BDTG_WithIPInfo) {
    if      (MVABin == 0) MVACut = 0.139;
    else if (MVABin == 1) MVACut = 0.525;
    else if (MVABin == 2) MVACut = 0.543; 
    else if (MVABin == 3) MVACut = 0.947;
    else if (MVABin == 4) MVACut = 0.950;
    else if (MVABin == 5) MVACut = 0.884;
  } else if (idType == ElectronTools::kMVAID_BDTG_IDIsoCombined) {
    if      (MVABin == 0) MVACut = 0.4202;
    else if (MVABin == 1) MVACut = 0.6206;
    else if (MVABin == 2) MVACut = 0.6190; 
    else if (MVABin == 3) MVACut = 0.9590;
    else if (MVABin == 4) MVACut = 0.9586;
    else if (MVABin == 5) MVACut = 0.9278;
  } else if (idType == ElectronTools::kMVAID_BDTG_IDHWW2012TrigV0) {
    if      (MVABin == 0) MVACut = 0.000;
    else if (MVABin == 1) MVACut = 0.100;
    else if (MVABin == 2) MVACut = 0.620;
    else if (MVABin == 3) MVACut = 0.940;
    else if (MVABin == 4) MVACut = 0.850;
    else if (MVABin == 5) MVACut = 0.920;
  } else if (idType == ElectronTools::kMVAID_BDTG_IDIsoCombinedHWW2012TrigV4) {
    if      (MVABin == 0) MVACut = 0.123;
    else if (MVABin == 1) MVACut = 0.219;
    else if (MVABin == 2) MVACut = 0.509;
    else if (MVABin == 3) MVACut = 0.935;
    else if (MVABin == 4) MVACut = 0.889;
    else if (MVABin == 5) MVACut = 0.871;
  } else if (idType == ElectronTools::kHggLeptonTagId2012HCP) {
    if      (MVABin == 0) MVACut = 0.9;
    else if (MVABin == 1) MVACut = 0.9;
    else if (MVABin == 2) MVACut = 0.9;
    else if (MVABin == 3) MVACut = 0.9;
    else if (MVABin == 4) MVACut = 0.9;
    else if (MVABin == 5) MVACut = 0.9;
  }

  //if (isDebug == kTRUE || true) {
  //  printf("PassElMVAID(%d): %d, pt, eta = %f, %f, rho = %f(%f) : MVA = %f, bin: %d\n",
  //         (MVAValue > MVACut),GetEventHeader()->EvtNum(),el->Pt(), eta,
  //	   fPileupEnergyDensity->At(0)->Rho(),fPileupEnergyDensity->At(0)->RhoKt6PFJets(),MVAValue,MVABin);
  //}

  if (MVAValue > MVACut) return kTRUE;
  return kFALSE;
}

Double_t ElectronIDMod::EvaluateMVAID(const Electron *el, ElectronTools::EElIdType idType, 
                                    const Vertex *vertex, const PFCandidateCol *PFCands,
                                    const PileupEnergyDensityCol *PileupEnergyDensity) const
{ 
  Bool_t isDebug = kFALSE;
  Double_t MVAValue = 0;
  if     (idType == ElectronTools::kMVAID_BDTG_IDIsoCombined) {
    MVAValue = fElectronIDMVA->MVAValue(el, vertex, PFCands, PileupEnergyDensity, fIntRadius);
  } 
  else if (idType == ElectronTools::kMVAID_BDTG_IDHWW2012TrigV0) {
    ElectronOArr *tempElectrons = new  ElectronOArr;
    MuonOArr     *tempMuons     = new  MuonOArr;
    MVAValue = fElectronIDMVA->MVAValue(el, vertex, PFCands, PileupEnergyDensity, ElectronTools::kEleEANoCorr, 
                                        tempElectrons, tempMuons, isDebug);
    delete tempElectrons;
    delete tempMuons;
  }
  else if (idType == ElectronTools::kMVAID_BDTG_IDIsoCombinedHWW2012TrigV4) {

    MVAValue = fElectronIDMVA->MVAValue(el, vertex, fVertices, PFCands, 
                                        PileupEnergyDensity, ElectronTools::kEleEANoCorr, 
                                        isDebug);
  }
  else {
    MVAValue = fElectronIDMVA->MVAValue(el, vertex);
  }

  return MVAValue;
  
}

//--------------------------------------------------------------------------------------------------
Bool_t ElectronIDMod::PassIDCut(const Electron *ele, ElectronTools::EElIdType idType, 
                                const Vertex *vertex) const
{
  Bool_t idcut = kFALSE;
  switch (idType) {
    case ElectronTools::kTight:
      idcut = ele->PassTightID();
      break;
    case ElectronTools::kLoose:
      idcut = ele->PassLooseID();
      break;
    case ElectronTools::kLikelihood:
      idcut = ElectronTools::PassCustomID(ele, ElectronTools::kVBTFWorkingPointFakeableId) &&
              PassLikelihoodID(ele);
      break;
    case ElectronTools::kNoId:
      idcut = kTRUE;
      break;
    case ElectronTools::kCustomIdLoose:
      idcut = ElectronTools::PassCustomID(ele, ElectronTools::kCustomIdLoose);
      break;
    case ElectronTools::kCustomIdTight:
      idcut = ElectronTools::PassCustomID(ele, ElectronTools::kCustomIdTight);
      break;
    case ElectronTools::kVBTFWorkingPointFakeableId:
      idcut = ElectronTools::PassCustomID(ele, ElectronTools::kVBTFWorkingPointFakeableId);
      break;
    case ElectronTools::kVBTFWorkingPoint95Id:
      idcut = ElectronTools::PassCustomID(ele, ElectronTools::kVBTFWorkingPoint95Id);
      break;
    case ElectronTools::kVBTFWorkingPoint90Id:
      idcut = ElectronTools::PassCustomID(ele, ElectronTools::kVBTFWorkingPoint90Id);
      break;
    case ElectronTools::kVBTFWorkingPoint85Id:
      idcut = ElectronTools::PassCustomID(ele, ElectronTools::kVBTFWorkingPoint85Id);
      break;
    case ElectronTools::kVBTFWorkingPoint80Id:
      idcut = ElectronTools::PassCustomID(ele, ElectronTools::kVBTFWorkingPoint80Id);
      break;
    case ElectronTools::kVBTFWorkingPointLowPtId:
      idcut = ElectronTools::PassCustomID(ele, ElectronTools::kVBTFWorkingPointLowPtId);
      break;
    case ElectronTools::kVBTFWorkingPoint70Id:
      idcut = ElectronTools::PassCustomID(ele, ElectronTools::kVBTFWorkingPoint70Id);
      break;
    case ElectronTools::kHggLeptonTagId:
      idcut = ElectronTools::PassHggLeptonTagID(ele);
      break;
    case ElectronTools::kHggLeptonTagId2012:
      idcut = ElectronTools::PassHggLeptonTagID2012(ele);
      break;
    case ElectronTools::kHggLeptonTagId2012HCP:
      {
	idcut = PassMVAID(ele, ElectronTools::kHggLeptonTagId2012HCP, 
			  vertex, fPFCandidates, fPileupEnergyDensity);  
      }
      break;
    case ElectronTools::kMVAID_BDTG_NoIPInfo:
    {
      idcut = ElectronTools::PassCustomID(ele, ElectronTools::kVBTFWorkingPointFakeableId);
      if (idcut) idcut = PassMVAID(ele, ElectronTools::kMVAID_BDTG_NoIPInfo, 
                                   vertex, fPFCandidates, fPileupEnergyDensity);
    }
      break;
    case ElectronTools::kMVAID_BDTG_WithIPInfo:
    {
      idcut = ElectronTools::PassCustomID(ele, ElectronTools::kVBTFWorkingPointFakeableId);
      if (idcut) idcut = PassMVAID(ele, ElectronTools::kMVAID_BDTG_WithIPInfo, 
                                   vertex, fPFCandidates, fPileupEnergyDensity);
    }
      break;
    case ElectronTools::kMVAID_BDTG_IDIsoCombined:
    {
      idcut = ElectronTools::PassCustomID(ele, ElectronTools::kVBTFWorkingPointFakeableId);
      if (idcut) idcut = PassMVAID(ele, ElectronTools::kMVAID_BDTG_IDIsoCombined, 
                                   vertex, fPFCandidates, fPileupEnergyDensity );
    }
      break;
    case ElectronTools::kMVAID_BDTG_IDHWW2012TrigV0:
    {
      idcut = ElectronTools::PassCustomID(ele, ElectronTools::kVBTFWorkingPointFakeableId);
      if (idcut) idcut = PassMVAID(ele, ElectronTools::kMVAID_BDTG_IDHWW2012TrigV0, 
                                   vertex, fPFCandidates, fPileupEnergyDensity );
    }
      break;
    case ElectronTools::kMVAID_BDTG_IDIsoCombinedHWW2012TrigV4:
    {
      idcut = ElectronTools::PassCustomID(ele, ElectronTools::kVBTFWorkingPointFakeableId);
      if (idcut) idcut = PassMVAID(ele, ElectronTools::kMVAID_BDTG_IDIsoCombinedHWW2012TrigV4, 
                                   vertex, fPFCandidates, fPileupEnergyDensity );
    }
      break;
    default:
      break;
  }
  
  return idcut;
}

//--------------------------------------------------------------------------------------------------
Bool_t ElectronIDMod::PassIsolationCut(const Electron *ele, ElectronTools::EElIsoType isoType,
                                       const TrackCol *tracks, const Vertex *vertex, 
				       const Double_t rho, ElectronTools::EElIdType idType) const
{

  Bool_t isocut = kFALSE;
  switch (isoType) {
    case ElectronTools::kTrackCalo:
      isocut = (ele->TrackIsolationDr03() < fTrackIsolationCut) &&
        (ele->CaloIsolation() < fCaloIsolationCut);
      break;
    case ElectronTools::kTrackJura:
      isocut = (ele->TrackIsolationDr03() < ele->Pt()*fTrackIsolationCut) &&
               (ele->EcalRecHitIsoDr03()  < ele->Pt()*fEcalJuraIsoCut) &&
               (ele->HcalTowerSumEtDr03() < ele->Pt()*fHcalIsolationCut);
      break;
    case ElectronTools::kTrackJuraCombined:
      isocut = (ele->TrackIsolationDr03() + ele->EcalRecHitIsoDr03() 
                - 1.5 < fCombIsolationCut);
      break;
    case ElectronTools::kTrackJuraSliding:
    {
      Double_t totalIso = ele->TrackIsolationDr03() + TMath::Max(ele->EcalRecHitIsoDr03() + ele->HcalTowerSumEtDr03() - rho * TMath::Pi() * 0.3 * 0.3, 0.0);
      if (ele->SCluster()->AbsEta() < 1.479) totalIso = ele->TrackIsolationDr03() + TMath::Max(TMath::Max(ele->EcalRecHitIsoDr03() - 1.0, 0.0) + ele->HcalTowerSumEtDr03() - rho * TMath::Pi() * 0.3 * 0.3, 0.0);
      if (totalIso < (ele->Pt()*fCombIsolationCut) )
        isocut = kTRUE;
    }
    break;
    case ElectronTools::kTrackJuraSlidingNoCorrection:
    {
      Double_t totalIso = ele->TrackIsolationDr03() + (ele->EcalRecHitIsoDr03() + ele->HcalTowerSumEtDr03());
      if (ele->SCluster()->AbsEta() < 1.479) totalIso = ele->TrackIsolationDr03() + (TMath::Max(ele->EcalRecHitIsoDr03() - 1.0, 0.0) + ele->HcalTowerSumEtDr03());
      if (totalIso < (ele->Pt()*fCombIsolationCut) )
        isocut = kTRUE;
    }
    break;
    case ElectronTools::kCombinedRelativeConeAreaCorrected:
    {
      Double_t totalIso = ele->TrackIsolationDr03() + TMath::Max(ele->EcalRecHitIsoDr03() + ele->HcalTowerSumEtDr03() - rho * TMath::Pi() * 0.3 * 0.3, 0.0);
      if (ele->SCluster()->AbsEta() < 1.5)  { // Barrel
	if (totalIso < (ele->Pt()*fCombRelativeIsolationCut) )
	  isocut = kTRUE;
      } else {
	if (totalIso < (ele->Pt()*fCombRelativeIsolationCut_EE) )
	  isocut = kTRUE;
      }
    }
    break;
    case ElectronTools::kPFIso:
    {
      Double_t pfIsoCutValue = 9999;
      if (fPFIsolationCut > 0) {
        pfIsoCutValue = fPFIsolationCut;
      } else {
        if (ele->SCluster()->AbsEta() < 1.479) {
          if (ele->Pt() > 20) {
            pfIsoCutValue = 0.13;
          } else {
            pfIsoCutValue = 0.13;
          }
        } else {
          if (ele->Pt() > 20) {
            pfIsoCutValue = 0.09;
          } else {
            pfIsoCutValue = 0.09;
          }
	}
      }
      Double_t totalIso = IsolationTools::PFElectronIsolation(ele, fPFCandidates, vertex, 0.1, 1.0, 0.4, fIntRadius);
      if (totalIso < (ele->Pt()*pfIsoCutValue) )
        isocut = kTRUE;
    }
    break;
    case ElectronTools::kPFIsoNoL:
    {
      Double_t pfIsoCutValue = 9999;
      if (fPFIsolationCut > 0) {
        pfIsoCutValue = fPFIsolationCut;
      } else {
        if (ele->SCluster()->AbsEta() < 1.479) {
          if (ele->Pt() > 20) {
            pfIsoCutValue = 0.13;
          } else {
            pfIsoCutValue = 0.13;
          }
        } else {
          if (ele->Pt() > 20) {
            pfIsoCutValue = 0.09;
          } else {
            pfIsoCutValue = 0.09;
          }
	}
      }
      Double_t totalIso = IsolationTools::PFElectronIsolation(ele, fPFCandidates, fNonIsolatedMuons, fNonIsolatedElectrons, vertex, 0.1, 1.0, 0.4, fIntRadius);
      if (totalIso < (ele->Pt()*pfIsoCutValue) )
        isocut = kTRUE;
    }
    break;
    case ElectronTools::kVBTFWorkingPoint95Iso:
      isocut = ElectronTools::PassCustomIso(ele, ElectronTools::kVBTFWorkingPoint95Iso, fApplyCombinedIso);
      break;
    case ElectronTools::kVBTFWorkingPoint90Iso:
      isocut = ElectronTools::PassCustomIso(ele, ElectronTools::kVBTFWorkingPoint90Iso, fApplyCombinedIso);
      break;
    case ElectronTools::kVBTFWorkingPoint85Iso:
      isocut = ElectronTools::PassCustomIso(ele, ElectronTools::kVBTFWorkingPoint85Iso, fApplyCombinedIso);
      break;
    case ElectronTools::kVBTFWorkingPoint80Iso:
      isocut = ElectronTools::PassCustomIso(ele, ElectronTools::kVBTFWorkingPoint80Iso, fApplyCombinedIso);
      break;
    case ElectronTools::kVBTFWorkingPoint70Iso:
      isocut = ElectronTools::PassCustomIso(ele, ElectronTools::kVBTFWorkingPoint70Iso, fApplyCombinedIso);
      break;
    case ElectronTools::kMVAIso_BDTG_IDIsoCombined:
      isocut = (ele->TrackIsolationDr03() < ele->Pt()*0.2) &&
               (ele->EcalRecHitIsoDr03()  < ele->Pt()*0.2) &&
               (ele->HcalTowerSumEtDr03() < ele->Pt()*0.2);
      break;
    case ElectronTools::kPFIso_HWW2012TrigV0:
    {
      Bool_t isDebug = kFALSE;
      if (isDebug == kTRUE) {
        printf("PFIso_HWW2012TrigV0: %d, pt, eta = %f, %f, rho = %f(%f) : ",
           GetEventHeader()->EvtNum(),ele->Pt(), ele->Eta(),
	   fPileupEnergyDensity->At(0)->Rho(),fPileupEnergyDensity->At(0)->RhoKt6PFJets());
      }
      ElectronOArr *tempIsoElectrons = new  ElectronOArr;
      MuonOArr     *tempIsoMuons     = new  MuonOArr;
      Double_t IsoOverPt = IsolationTools::PFElectronIsolation2012(ele, vertex, fPFNoPileUpCands, 
       fPileupEnergyDensity, ElectronTools::kEleEANoCorr, tempIsoElectrons, tempIsoMuons, 0.4, isDebug);
      delete tempIsoElectrons;
      delete tempIsoMuons;
      Double_t eta = ele->SCluster()->AbsEta();
      Double_t IsoCut = -1;
      if (ele->Pt() <  20 && eta <  0.800		) IsoCut = 0.150;
      if (ele->Pt() <  20 && eta >= 0.800 && eta < 1.479) IsoCut = 0.150;
      if (ele->Pt() <  20 && eta >= 1.479		) IsoCut = 0.150;
      if (ele->Pt() >= 20 && eta <  0.800		) IsoCut = 0.150;
      if (ele->Pt() >= 20 && eta >= 0.800 && eta < 1.479) IsoCut = 0.150;
      if (ele->Pt() >= 20 && eta >= 1.479		) IsoCut = 0.150;
      if (IsoOverPt < IsoCut) isocut = kTRUE;
    }
    break;
  case ElectronTools::kPFIso_HggLeptonTag2012:
    {
      Bool_t isDebug = kFALSE;
      if (isDebug == kTRUE) {
        printf("PFIso_HggLeptonTag2012: %d, pt, eta = %f, %f, rho = %f(%f) : ",
	       GetEventHeader()->EvtNum(),ele->Pt(), ele->Eta(),
	       fPileupEnergyDensity->At(0)->Rho(),fPileupEnergyDensity->At(0)->RhoKt6PFJets());
      }
      ElectronOArr *tempIsoElectrons = new  ElectronOArr;
      MuonOArr     *tempIsoMuons     = new  MuonOArr;
      Double_t IsoOverPt = IsolationTools::PFElectronIsolation2012LepTag(ele, vertex, fPFNoPileUpCands, 
									 fPileupEnergyDensity, ElectronTools::kEleEAData2012, tempIsoElectrons, tempIsoMuons, 0.3, isDebug);
      delete tempIsoElectrons;
      delete tempIsoMuons;
      Double_t eta = ele->SCluster()->AbsEta();
      Double_t IsoCut = -1;
      if (ele->Pt() <  20 && eta <  0.800		) IsoCut = 0.150;
      if (ele->Pt() <  20 && eta >= 0.800 && eta < 1.479) IsoCut = 0.150;
      if (ele->Pt() <  20 && eta >= 1.479		) IsoCut = 0.10;
      if (ele->Pt() >= 20 && eta <  0.800		) IsoCut = 0.150;
      if (ele->Pt() >= 20 && eta >= 0.800 && eta < 1.479) IsoCut = 0.150;
      if (ele->Pt() >= 20 && eta >= 1.479		) IsoCut = 0.150;
      if (IsoOverPt < IsoCut) isocut = kTRUE;
    }
    break;
  case ElectronTools::kPFIso_HggLeptonTag2012HCP:
    {
      
      Bool_t isDebug = kFALSE;
      if (isDebug == kTRUE) {
        printf("PFIso_HggLeptonTag2012HCP: %d, pt, eta, sceta = %f, %f, %f,rho = %f(%f) : ",
	       GetEventHeader()->EvtNum(),ele->Pt(), ele->Eta(),ele->SCluster()->AbsEta(),
	   fPileupEnergyDensity->At(0)->Rho(),fPileupEnergyDensity->At(0)->RhoKt6PFJets());
      }
      ElectronOArr *tempIsoElectrons = new  ElectronOArr;
      MuonOArr     *tempIsoMuons     = new  MuonOArr;
      Double_t IsoOverPt = IsolationTools::PFElectronIsolation2012LepTag(ele, vertex, fPFNoPileUpCands,fPileupEnergyDensity, ElectronTools::kEleEAData2012, tempIsoElectrons, tempIsoMuons, 0.3, isDebug);
      //printf("ming sync check IsoOverPt:%f\n\n",IsoOverPt);
      delete tempIsoElectrons;
      delete tempIsoMuons;
      Double_t eta = ele->SCluster()->AbsEta();
      Double_t IsoCut = -1;
      if (ele->Pt() <  20 && eta <  0.800		) IsoCut = 0.150;
      if (ele->Pt() <  20 && eta >= 0.800 && eta < 1.479) IsoCut = 0.150;
      if (ele->Pt() <  20 && eta >= 1.479		) IsoCut = 0.150;
      if (ele->Pt() >= 20 && eta <  0.800		) IsoCut = 0.150;
      if (ele->Pt() >= 20 && eta >= 0.800 && eta < 1.479) IsoCut = 0.150;
      if (ele->Pt() >= 20 && eta >= 1.479		) IsoCut = 0.150;
      if (IsoOverPt < IsoCut) isocut = kTRUE;
    }
      break;
    case ElectronTools::kMVAIso_BDTG_IDIsoCombinedHWW2012TrigV4:
    {
      if (fabs(ele->SCluster()->Eta()) < 1.479) {
        isocut = ((ele->TrackIsolationDr03() - 1.0) < ele->Pt()*0.2) &&
          (ele->EcalRecHitIsoDr03()  < ele->Pt()*0.2) &&
          (ele->HcalTowerSumEtDr03() < ele->Pt()*0.2);
      } else {
        isocut = (ele->TrackIsolationDr03() < ele->Pt()*0.2) &&
          (ele->EcalRecHitIsoDr03()  < ele->Pt()*0.2) &&
          (ele->HcalTowerSumEtDr03() < ele->Pt()*0.2);
      }
    }
      break;
    case ElectronTools::kNoIso:
      isocut = kTRUE;
      break;
    case ElectronTools::kCustomIso:
      break;
    default:
      break;
  }

  return isocut;
}


//--------------------------------------------------------------------------------------------------
void ElectronIDMod::Process()
{
  // Process entries of the tree. 
  if (fElIdType == ElectronTools::kHggLeptonTagId2012HCP) {
    LoadEventObject(fPVName,fVertices);
  }else{
    fVertices = GetObjThisEvt<VertexOArr>(fVertexName);
  }

  if (fElIsoType != ElectronTools::kPFIsoNoL) {
    LoadEventObject(fElectronBranchName, fElectrons);
  }
  else {
    fElectrons    = GetObjThisEvt<ElectronOArr>(fElectronBranchName);
    fNonIsolatedMuons	  = GetObjThisEvt<MuonCol>(fNonIsolatedMuonsName);
    fNonIsolatedElectrons = GetObjThisEvt<ElectronCol>(fNonIsolatedElectronsName);
  }
  LoadEventObject(fBeamSpotName, fBeamSpot);
  LoadEventObject(fTrackName, fTracks);
  LoadEventObject(fPFCandidatesName, fPFCandidates);

  if (fElIsoType == ElectronTools::kTrackJuraSliding || 
     fElIsoType == ElectronTools::kCombinedRelativeConeAreaCorrected || 
     fElIsoType == ElectronTools::kMVAIso_BDTG_IDIsoCombined || 
     fElIdType  == ElectronTools::kMVAID_BDTG_IDHWW2012TrigV0  || 
     fElIsoType == ElectronTools::kMVAIso_BDTG_IDIsoCombinedHWW2012TrigV4 ||
     fElIsoType == ElectronTools::kPFIso_HWW2012TrigV0      || 
     fElIsoType == ElectronTools::kPFIso_HggLeptonTag2012 ||
     fElIsoType == ElectronTools::kPFIso_HggLeptonTag2012HCP
    ) {
    LoadEventObject(fPileupEnergyDensityName, fPileupEnergyDensity);
    assert(fPileupEnergyDensity);
  }

  if (fElIsoType == ElectronTools::kPFIso_HWW2012TrigV0 || fElIsoType == ElectronTools::kPFIso_HggLeptonTag2012 || fElIsoType == ElectronTools::kPFIso_HggLeptonTag2012HCP) {
    // Name is hardcoded, can be changed if someone feels to do it
    fPFNoPileUpCands = GetObjThisEvt<PFCandidateCol>(fPFNoPileUpName);
  }

  //get trigger object collection if trigger matching is enabled
  const TriggerObjectCol *trigObjs = 0;
  if (fApplyTriggerMatching) {
    trigObjs = GetHLTObjects(fTrigObjectsName);
  }
  
  ElectronOArr *GoodElectrons = new ElectronOArr;
  GoodElectrons->SetName(fGoodElectronsName);

  Double_t MVAValue = -100;
  int NPass = 0;


  std::vector<std::pair<double, int> > elemvaidxs;
  
  //printf("\n");
  //printf("ming sync check number of electrons:%d\n",fElectrons->GetEntries());

  for (UInt_t i=0; i<fElectrons->GetEntries() && fVertices->GetEntries() > 0 ; ++i) {

    //printf("ming sync check electron i:%d\n",i);
   
    const Electron *e = fElectrons->At(i);        

    if (e->SCluster() == 0) 
      continue;
    
    //printf("ming sync check elept:%f  sceta:%f \n",e->Pt(),e->SCluster()->Eta());

    if (e->Pt() < fElectronPtMin) 
      continue;
    
    if (e->SCluster()->Et() < fElectronEtMin)
      continue;    

       
    if (e->AbsEta() > fElectronEtaMax) 
      continue;

    if (fApplyEcalFiducial && ((e->SCluster()->AbsEta()>1.4442 && e->SCluster()->AbsEta()<1.5666) || e->SCluster()->AbsEta()>2.5 )) {
      continue;
    }
    

    if (fApplyEcalSeeded && !e->IsEcalDriven()) {
      continue;
    }    
    
    //apply trigger matching
    Bool_t matchTrigger = fApplyTriggerMatching && ElectronTools::PassTriggerMatching(e,trigObjs);
    if (fApplyTriggerMatching && !matchTrigger)
      continue;
    
    //apply ECAL spike removal    
    Bool_t spikecut = ElectronTools::PassSpikeRemovalFilter(e);
    if (fApplySpikeRemoval && !spikecut)
      continue;

    //apply Isolation Cut
    Double_t Rho = 0.0;
    if (fElIsoType == ElectronTools::kTrackJuraSliding
        || fElIsoType == ElectronTools::kCombinedRelativeConeAreaCorrected 
        || fElIsoType == ElectronTools::kMVAIso_BDTG_IDIsoCombined 
     ) {      
      switch(fTheRhoType) {
      case RhoUtilities::MIT_RHO_VORONOI_HIGH_ETA:
	Rho = fPileupEnergyDensity->At(0)->Rho();
	break;
      case RhoUtilities::MIT_RHO_VORONOI_LOW_ETA:
	Rho = fPileupEnergyDensity->At(0)->RhoLowEta();
	break;
      case RhoUtilities::MIT_RHO_RANDOM_HIGH_ETA:
	Rho = fPileupEnergyDensity->At(0)->RhoRandom();
	break;
      case RhoUtilities::MIT_RHO_RANDOM_LOW_ETA:
	Rho = fPileupEnergyDensity->At(0)->RhoRandomLowEta();
	break;
      case RhoUtilities::CMS_RHO_RHOKT6PFJETS:
	Rho = fPileupEnergyDensity->At(0)->RhoKt6PFJets();
	break;
      default:
	// use the old default
	Rho = fPileupEnergyDensity->At(0)->Rho();
	break;
      }
    }
    
    Bool_t isocut = kFALSE;
    Double_t distVtx = 999.0;
    Int_t closestVtx = 0;
    if (fElIsoType == ElectronTools::kPFIso_HggLeptonTag2012HCP) {
      for (UInt_t nv=0; nv<fVertices->GetEntries(); nv++) {
	double dz = TMath::Abs(e->GsfTrk()->DzCorrected(*fVertices->At(nv)));
	if (dz < distVtx) {
	  distVtx    = dz;
	  closestVtx = nv;
	}
      }
      isocut = PassIsolationCut(e, fElIsoType, fTracks, fVertices->At(closestVtx), Rho, fElIdType);
    }else{
      isocut = PassIsolationCut(e, fElIsoType, fTracks, fVertices->At(0), Rho, fElIdType);
    }

    if (!isocut)
      continue;

    //printf("ming sync check pass isocut\n");

    // apply conversion filters
    Bool_t passConvVetoType1 = kFALSE;
    if (fApplyConvFilterType1) {
      LoadEventObject(fConversionBranchName, fConversions);
      passConvVetoType1 = ElectronTools::PassConversionFilter(e, fConversions, 
                                                         fBeamSpot->At(0), 0, 1e-6, 2.0, kTRUE, kFALSE);      
    }
    else {
      passConvVetoType1 = kTRUE;
    }

    if (passConvVetoType1 == kFALSE) continue;

    Bool_t passConvVetoType2 = kFALSE;
    if (fApplyConvFilterType2) {
      passConvVetoType2 = TMath::Abs(e->ConvPartnerDCotTheta()) >= 0.02 || 
                          TMath::Abs(e->ConvPartnerDist())      >= 0.02;
    }
    else {
      passConvVetoType2 = kTRUE;
    }
    
    if (passConvVetoType2 == kFALSE) continue;
    
    // apply NExpectedHitsInner Cut
    if (fInvertNExpectedHitsInnerCut == kFALSE && fNExpectedHitsInnerCut < 999 && 
       e->CorrectedNExpectedHitsInner() > fNExpectedHitsInnerCut) continue;

    // apply NExpectedHitsInner inverted Cut
    if (fInvertNExpectedHitsInnerCut == kTRUE && fNExpectedHitsInnerCut < 999 && 
       e->CorrectedNExpectedHitsInner() <= fNExpectedHitsInnerCut) continue;

    //printf("ming sync check pass Nexpectedinner cut \n");

    
    // apply d0 cut
    if (fApplyD0Cut) {
      Bool_t passD0cut = kTRUE;
      if (fWhichVertex >= -1) passD0cut = ElectronTools::PassD0Cut(e, fVertices, fD0Cut, fWhichVertex);
      else                   passD0cut = ElectronTools::PassD0Cut(e, fBeamSpot, fD0Cut);
      if (!passD0cut)
        continue;
    }
 
    //printf("ming sync check pass d0 cut \n");

 
    // apply dz cut
    if (fApplyDZCut) {
      Bool_t passDZcut = ElectronTools::PassDZCut(e, fVertices, fDZCut, fWhichVertex);
      if (!passDZcut)
        continue;
    }
 
    //apply id cut
    Bool_t idcut = kFALSE;
    if (fElIdType == ElectronTools::kHggLeptonTagId2012HCP) {
      MVAValue = EvaluateMVAID(e, ElectronTools::kHggLeptonTagId2012HCP, fVertices->At(closestVtx), fPFCandidates, fPileupEnergyDensity);
      idcut = PassIDCut(e, fElIdType, fVertices->At(closestVtx));
      if (idcut) {
        //fill temp vector of indexes and mva values for sorting
        elemvaidxs.push_back(std::make_pair(MVAValue,i));
      }
    } else {
      idcut = PassIDCut(e, fElIdType, fVertices->At(0));
    }

    if (!idcut) 
      continue;
    
    // apply charge filter
    if (fChargeFilter == kTRUE) {
      Bool_t passChargeFilter = ElectronTools::PassChargeFilter(e);
      if (!passChargeFilter) continue;
    }

    // apply full combined id, using Tight cuts
    if (fCombinedIdCut == kTRUE) {
      LoadEventObject(fConversionBranchName, fConversions);
      Int_t result = ElectronTools::PassTightId(e, *&fVertices, fConversions, 2);
      if (result != 15) continue;
    }

    // add good electron
    if (fElIdType!=ElectronTools::kHggLeptonTagId2012HCP) {
      // make sure to mark the selected electron (to be able to pickup by skimmer)
      e->Mark();
      // add it to the good electron collection (re-ordering follows below)
      GoodElectrons->Add(e);
    }

    NPass = NPass+1;
  }

  // sort according to pt
  GoodElectrons->Sort();

  if (fElIdType == ElectronTools::kHggLeptonTagId2012HCP) {
    // sort by mva value (in descending order)
    std::sort(elemvaidxs.begin(), elemvaidxs.end(), std::greater<std::pair<double, int> >());
    
    // fill final list of electrons
    for (std::vector<std::pair<double,int> >::const_iterator it=elemvaidxs.begin(); it!=elemvaidxs.end(); ++it) {
      GoodElectrons->Add(fElectrons->At(it->second));
    }
  }

  // add to event for other modules to use
  AddObjThisEvt(GoodElectrons);
}

//--------------------------------------------------------------------------------------------------
void ElectronIDMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we just request the electron collection branch.

  if (fElectronIDType.CompareTo("Hgg_LeptonTag_2012IdHCP") == 0) {
    ReqEventObject(fPVName,fVertices,true);
  }
 
   // In this case we cannot have a branch
  if (fElectronIsoType.CompareTo("PFIsoNoL") != 0) {
    ReqEventObject(fElectronBranchName, fElectrons,fElectronsFromBranch);
  }
  ReqEventObject(fBeamSpotName, fBeamSpot, kTRUE);
  ReqEventObject(fTrackName, fTracks, kTRUE);
  ReqEventObject(fPFCandidatesName, fPFCandidates, kTRUE);
  if (fElectronIsoType.CompareTo("TrackJuraSliding") == 0 
      || fElectronIsoType.CompareTo("CombinedRelativeConeAreaCorrected") == 0 
      || fElectronIsoType.CompareTo("MVA_BDTG_IDIsoCombined") == 0
      || fElectronIDType.CompareTo("MVA_BDTG_IDHWW2012TrigV0") == 0
      || fElectronIsoType.CompareTo("PFIso_HWW2012TrigV0") == 0
      || fElectronIsoType.CompareTo("PFIso_HggLeptonTag2012") == 0
      || fElectronIDType.CompareTo("MVA_BDTG_IDIsoCombinedHWW2012TrigV4") == 0
      || fElectronIsoType.CompareTo("MVA_BDTG_IDIsoCombinedHWW2012TrigV4") == 0
      || fElectronIsoType.CompareTo("PFIso_HggLeptonTag2012HCP") == 0
     ) {
    ReqEventObject(fPileupEnergyDensityName, fPileupEnergyDensity, kTRUE);
  }

  if (fCombinedIdCut == kTRUE) {
    fElectronIDType  	  = "NoId";
    fElectronIsoType 	  = "NoIso";
    fApplyConvFilterType1 = kFALSE;
    fApplyConvFilterType2 = kFALSE;
    fApplyD0Cut           = kFALSE;
    fApplyDZCut           = kFALSE;
  }

  if (fApplyConvFilterType1 || fCombinedIdCut == kTRUE)
    ReqEventObject(fConversionBranchName, fConversions, kTRUE);

  Setup();

}

//--------------------------------------------------------------------------------------------------
void ElectronIDMod::Setup()
{
  // Set all options properly before execution.

  if (fElectronIDType.CompareTo("Tight") == 0) 
    fElIdType = ElectronTools::kTight;
  else if (fElectronIDType.CompareTo("Loose") == 0) 
    fElIdType = ElectronTools::kLoose;
  else if (fElectronIDType.CompareTo("Likelihood") == 0) {
    if (!fLH) { cout << "Error: Likelihood not initialized.\n"; assert(0); }
    fElIdType = ElectronTools::kLikelihood;
  } else if (fElectronIDType.CompareTo("NoId") == 0) 
    fElIdType = ElectronTools::kNoId;
  else if (fElectronIDType.CompareTo("ZeeId") == 0) 
    fElIdType = ElectronTools::kZeeId;
  else if (fElectronIDType.CompareTo("CustomLoose") == 0) 
    fElIdType = ElectronTools::kCustomIdLoose;
  else if (fElectronIDType.CompareTo("CustomTight") == 0) 
    fElIdType = ElectronTools::kCustomIdTight;
  else if (fElectronIDType.CompareTo("VBTFWorkingPointFakeableId") == 0) 
    fElIdType = ElectronTools::kVBTFWorkingPointFakeableId;
  else if (fElectronIDType.CompareTo("VBTFWorkingPoint95Id") == 0) 
    fElIdType = ElectronTools::kVBTFWorkingPoint95Id;
  else if (fElectronIDType.CompareTo("VBTFWorkingPoint90Id") == 0) 
    fElIdType = ElectronTools::kVBTFWorkingPoint90Id;
  else if (fElectronIDType.CompareTo("VBTFWorkingPoint80Id") == 0) 
    fElIdType = ElectronTools::kVBTFWorkingPoint80Id;
  else if (fElectronIDType.CompareTo("VBTFWorkingPointLowPtId") == 0) 
    fElIdType = ElectronTools::kVBTFWorkingPointLowPtId;
  else if (fElectronIDType.CompareTo("VBTFWorkingPoint85Id") == 0) 
    fElIdType = ElectronTools::kVBTFWorkingPoint85Id;
  else if (fElectronIDType.CompareTo("VBTFWorkingPoint70Id") == 0) 
    fElIdType = ElectronTools::kVBTFWorkingPoint70Id;
  else if (fElectronIDType.CompareTo("MVA_BDTG_NoIPInfo") == 0)
    fElIdType = ElectronTools::kMVAID_BDTG_NoIPInfo; 
  else if (fElectronIDType.CompareTo("MVA_BDTG_WithIPInfo") == 0)
    fElIdType = ElectronTools::kMVAID_BDTG_WithIPInfo; 
  else if (fElectronIDType.CompareTo("MVA_BDTG_IDIsoCombined") == 0)
    fElIdType = ElectronTools::kMVAID_BDTG_IDIsoCombined; 
  else if (fElectronIDType.CompareTo("MVA_BDTG_IDHWW2012TrigV0") == 0)
    fElIdType = ElectronTools::kMVAID_BDTG_IDHWW2012TrigV0; 
  else if (fElectronIDType.CompareTo("Hgg_LeptonTag_WP85Id") == 0)
    fElIdType = ElectronTools::kHggLeptonTagId;
  else if (fElectronIDType.CompareTo("Hgg_LeptonTag_2012Id") == 0)
    fElIdType = ElectronTools::kHggLeptonTagId2012;
  else if (fElectronIDType.CompareTo("Hgg_LeptonTag_2012IdHCP") == 0)
    fElIdType = ElectronTools::kHggLeptonTagId2012HCP;
  else if (fElectronIDType.CompareTo("MVA_BDTG_IDIsoCombinedHWW2012TrigV4") == 0)
    fElIdType = ElectronTools::kMVAID_BDTG_IDIsoCombinedHWW2012TrigV4; 

  else {
    SendError(kAbortAnalysis, "SlaveBegin",
              "The specified electron identification %s is not defined.",
              fElectronIDType.Data());
    return;
  }

  if (fElectronIsoType.CompareTo("TrackCalo") == 0 )
    fElIsoType = ElectronTools::kTrackCalo;
  else if (fElectronIsoType.CompareTo("TrackJura") == 0) 
    fElIsoType = ElectronTools::kTrackJura;
  else if (fElectronIsoType.CompareTo("TrackJuraCombined") == 0)
    fElIsoType = ElectronTools::kTrackJuraCombined;
  else if (fElectronIsoType.CompareTo("TrackJuraSliding") == 0)
    fElIsoType = ElectronTools::kTrackJuraSliding;
  else if (fElectronIsoType.CompareTo("TrackJuraSlidingNoCorrection") == 0)
    fElIsoType = ElectronTools::kTrackJuraSlidingNoCorrection;
  else if (fElectronIsoType.CompareTo("CombinedRelativeConeAreaCorrected") == 0)
    fElIsoType = ElectronTools::kCombinedRelativeConeAreaCorrected;
  else if (fElectronIsoType.CompareTo("PFIso") == 0 )
    fElIsoType = ElectronTools::kPFIso;
  else if (fElectronIsoType.CompareTo("PFIsoNoL") == 0 )
    fElIsoType = ElectronTools::kPFIsoNoL;
  else if (fElectronIsoType.CompareTo("NoIso") == 0 )
    fElIsoType = ElectronTools::kNoIso;
  else if (fElectronIsoType.CompareTo("ZeeIso") == 0 )
    fElIsoType = ElectronTools::kZeeIso;
  else if (fElectronIsoType.CompareTo("VBTFWorkingPoint95Iso") == 0 )
    fElIsoType = ElectronTools::kVBTFWorkingPoint95Iso;
  else if (fElectronIsoType.CompareTo("VBTFWorkingPoint90Iso") == 0 )
    fElIsoType = ElectronTools::kVBTFWorkingPoint90Iso;
  else if (fElectronIsoType.CompareTo("VBTFWorkingPoint85Iso") == 0 )
    fElIsoType = ElectronTools::kVBTFWorkingPoint85Iso;
  else if (fElectronIsoType.CompareTo("VBTFWorkingPoint80Iso") == 0 )
    fElIsoType = ElectronTools::kVBTFWorkingPoint80Iso;
  else if (fElectronIsoType.CompareTo("VBTFWorkingPoint70Iso") == 0 )
    fElIsoType = ElectronTools::kVBTFWorkingPoint70Iso;
  else if (fElectronIsoType.CompareTo("MVA_BDTG_IDIsoCombined") == 0 )
    fElIsoType = ElectronTools::kMVAIso_BDTG_IDIsoCombined;
  else if (fElectronIsoType.CompareTo("PFIso_HWW2012TrigV0") == 0)
    fElIsoType = ElectronTools::kPFIso_HWW2012TrigV0;
  else if (fElectronIsoType.CompareTo("PFIso_HggLeptonTag2012") == 0)
    fElIsoType = ElectronTools::kPFIso_HggLeptonTag2012;
  else if (fElectronIsoType.CompareTo("PFIso_HggLeptonTag2012HCP") == 0)
    fElIsoType = ElectronTools::kPFIso_HggLeptonTag2012HCP;
  else if (fElectronIsoType.CompareTo("MVA_BDTG_IDIsoCombinedHWW2012TrigV4") == 0 )
    fElIsoType = ElectronTools::kMVAIso_BDTG_IDIsoCombinedHWW2012TrigV4;
  else if (fElectronIsoType.CompareTo("Custom") == 0) {
    fElIsoType = ElectronTools::kCustomIso;
    SendError(kWarning, "SlaveBegin",
              "Custom electron isolation is not yet implemented.");
  } else {
    SendError(kAbortAnalysis, "SlaveBegin",
              "The specified electron isolation %s is not defined.",
              fElectronIsoType.Data());
    return;
  }

  // If we use MVA ID, need to load MVA weights
  if (fElIdType == ElectronTools::kMVAID_BDTG_NoIPInfo) {
    fElectronIDMVA = new ElectronIDMVA();
    fElectronIDMVA->Initialize("BDTG method",
                               fElectronMVAWeights_Subdet0Pt10To20,
                               fElectronMVAWeights_Subdet1Pt10To20,
                               fElectronMVAWeights_Subdet2Pt10To20,
                               fElectronMVAWeights_Subdet0Pt20ToInf,
                               fElectronMVAWeights_Subdet1Pt20ToInf,
                               fElectronMVAWeights_Subdet2Pt20ToInf,
                               ElectronIDMVA::kNoIPInfo,
			       fTheRhoType);
  }
  if (fElIdType == ElectronTools::kMVAID_BDTG_WithIPInfo) {
    fElectronIDMVA = new ElectronIDMVA();
    fElectronIDMVA->Initialize("BDTG method",
                               fElectronMVAWeights_Subdet0Pt10To20,
                               fElectronMVAWeights_Subdet1Pt10To20,
                               fElectronMVAWeights_Subdet2Pt10To20,
                               fElectronMVAWeights_Subdet0Pt20ToInf,
                               fElectronMVAWeights_Subdet1Pt20ToInf,
                               fElectronMVAWeights_Subdet2Pt20ToInf,
                               ElectronIDMVA::kWithIPInfo,
			       fTheRhoType);
  }
  if (fElIdType == ElectronTools::kMVAID_BDTG_IDIsoCombined || fElIsoType == ElectronTools::kMVAIso_BDTG_IDIsoCombined) {
    fElectronIDMVA = new ElectronIDMVA();
    fElectronIDMVA->Initialize("BDTG method",
                               fElectronMVAWeights_Subdet0Pt10To20,
                               fElectronMVAWeights_Subdet1Pt10To20,
                               fElectronMVAWeights_Subdet2Pt10To20,
                               fElectronMVAWeights_Subdet0Pt20ToInf,
                               fElectronMVAWeights_Subdet1Pt20ToInf,
                               fElectronMVAWeights_Subdet2Pt20ToInf,
                               ElectronIDMVA::kIDIsoCombined,
			       fTheRhoType);
  }

  if (fElIdType == ElectronTools::kMVAID_BDTG_IDHWW2012TrigV0) {
    fElectronIDMVA = new ElectronIDMVA();
    fElectronIDMVA->Initialize("BDTG method",
                               fElectronMVAWeights_Subdet0Pt10To20,
                               fElectronMVAWeights_Subdet1Pt10To20,
                               fElectronMVAWeights_Subdet2Pt10To20,
                               fElectronMVAWeights_Subdet0Pt20ToInf,
                               fElectronMVAWeights_Subdet1Pt20ToInf,
                               fElectronMVAWeights_Subdet2Pt20ToInf,
                               ElectronIDMVA::kIDHWW2012TrigV0,
			       fTheRhoType);
  }

  if (fElIdType == ElectronTools::kMVAID_BDTG_IDIsoCombinedHWW2012TrigV4
      || fElIsoType == ElectronTools::kMVAIso_BDTG_IDIsoCombinedHWW2012TrigV4) {
    fElectronIDMVA = new ElectronIDMVA();
    fElectronIDMVA->Initialize("BDTG method",
                               fElectronMVAWeights_Subdet0Pt10To20,
                               fElectronMVAWeights_Subdet1Pt10To20,
                               fElectronMVAWeights_Subdet2Pt10To20,
                               fElectronMVAWeights_Subdet0Pt20ToInf,
                               fElectronMVAWeights_Subdet1Pt20ToInf,
                               fElectronMVAWeights_Subdet2Pt20ToInf,
                               ElectronIDMVA::kIDIsoCombinedHWW2012TrigV4,
			       fTheRhoType);
  }

  if (fElIdType == ElectronTools::kHggLeptonTagId2012HCP) {
    fElectronIDMVA = new ElectronIDMVA();
    fElectronIDMVA->Initialize("BDTG method",
                               fElectronMVAWeights_Subdet0Pt10To20,
                               fElectronMVAWeights_Subdet1Pt10To20,
                               fElectronMVAWeights_Subdet2Pt10To20,
                               fElectronMVAWeights_Subdet0Pt20ToInf,
                               fElectronMVAWeights_Subdet1Pt20ToInf,
                               fElectronMVAWeights_Subdet2Pt20ToInf,
                               ElectronIDMVA::kIDEGamma2012NonTrigV1,
			       fTheRhoType);
  }
  
}

//--------------------------------------------------------------------------------------------------
void ElectronIDMod::Terminate()
{
  // Run finishing code on the computer (slave) that did the analysis
  delete fElectronIDMVA;
}
