// $Id: MuonIDMod.cc,v 1.92 2013/10/18 14:09:34 ceballos Exp $

#include "MitPhysics/Mods/interface/MuonIDMod.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataTree/interface/MuonFwd.h"
#include "MitAna/DataTree/interface/ElectronFwd.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitPhysics/Init/interface/ModNames.h"

using namespace mithep;

ClassImp(mithep::MuonIDMod)

//--------------------------------------------------------------------------------------------------
  MuonIDMod::MuonIDMod(const char *name, const char *title) :
  BaseMod(name,title),
  fPrintMVADebugInfo(kFALSE),
  fMuonBranchName(Names::gkMuonBrn),
  fCleanMuonsName(ModNames::gkCleanMuonsName),
  fNonIsolatedMuonsName("random"),
  fNonIsolatedElectronsName("random"),
  fVertexName(ModNames::gkGoodVertexesName),
  fBeamSpotName(Names::gkBeamSpotBrn),
  fTrackName(Names::gkTrackBrn),
  fPFCandidatesName(Names::gkPFCandidatesBrn),
  fPFNoPileUpName("PFNoPileUp"),
  fPFPileUpName("PFPileUp"),
  fMuonIDType("NoId"),
  fMuonIsoType("PFIso"),
  fMuonClassType("GlobalorTracker"),
  fTrackIsolationCut(3.0),
  fCaloIsolationCut(3.0),
  fCombIsolationCut(0.15),
  fCombRelativeIsolationCut(0.15),
  fPFIsolationCut(-999.0),
  fMuonPtMin(10.),
  fApplyD0Cut(kTRUE),
  fApplyDZCut(kTRUE),
  fD0Cut(0.020),
  fDZCut(0.10),
  fWhichVertex(-1),
  fEtaCut(2.4),
  fMuIDType(kIdUndef),
  fMuIsoType(kIsoUndef),
  fMuClassType(kClassUndef),
  fMuons(0),
  fVertices(0),
  fBeamSpot(0),
  fTracks(0),
  fPFCandidates(0),
  fPFNoPileUpCands(0),
  fPFPileUpCands(0),
  fIntRadius(0.0),
  fNonIsolatedMuons(0),
  fNonIsolatedElectrons(0),
  fPileupEnergyDensityName(Names::gkPileupEnergyDensityBrn),
  fPileupEnergyDensity(0),
  fMuonTools(0),
  fMuonIDMVA(0),
  fPVName(Names::gkPVBeamSpotBrn),
  fTheRhoType(RhoUtilities::DEFAULT)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void MuonIDMod::Process()
{
  // Process entries of the tree.

  if (fCleanMuonsName.CompareTo("HggLeptonTagMuons") == 0) {
    LoadEventObject(fPVName,fVertices);
  }
  else{
    fVertices = GetObjThisEvt<VertexOArr>(fVertexName);
  }

  if (fMuIsoType != kPFIsoNoL) {
    LoadEventObject(fMuonBranchName, fMuons);
  }
  else {
    fMuons = GetObjThisEvt<MuonOArr>(fMuonBranchName);
  }
  LoadEventObject(fBeamSpotName, fBeamSpot);
  LoadEventObject(fTrackName, fTracks);
  LoadEventObject(fPFCandidatesName, fPFCandidates);
  if (fMuIsoType == kTrackCaloSliding ||
      fMuIsoType == kCombinedRelativeConeAreaCorrected ||
      fMuIsoType == kPFIsoEffectiveAreaCorrected ||
      fMuIsoType == kMVAIso_BDTG_IDIso ||
      fMuIsoType == kIsoRingsV0_BDTG_Iso ||
      fMuIsoType == kIsoDeltaR
      ) {
    LoadEventObject(fPileupEnergyDensityName, fPileupEnergyDensity);
  }
  if (fMuIsoType == kPFRadialIso ||
      fMuIsoType == kIsoDeltaR   ||
      fMuIsoType == kPFIsoBetaPUCorrected) {
    fPFNoPileUpCands = GetObjThisEvt<PFCandidateCol>(fPFNoPileUpName);
    fPFPileUpCands   = GetObjThisEvt<PFCandidateCol>(fPFPileUpName);
  }

  MuonOArr *CleanMuons = new MuonOArr;
  CleanMuons->SetName(fCleanMuonsName);

  for (UInt_t i=0; i<fMuons->GetEntries() && fVertices->GetEntries() > 0 ; ++i) {
    const Muon *mu = fMuons->At(i);

    Bool_t pass = kFALSE;
    Double_t pt = 0;  // make sure pt is taken from the correct track!
    Double_t eta = 0; // make sure eta is taken from the correct track!
    switch (fMuClassType) {
    case kAll:
      pass = kTRUE;
      if (mu->HasTrk()) {
        pt  = mu->Pt();
        eta = TMath::Abs(mu->Eta());
      }
      break;
    case kGlobal:
      pass = mu->HasGlobalTrk() && mu->IsTrackerMuon();
      if (pass && mu->TrackerTrk()) {
        pt  = mu->TrackerTrk()->Pt();
        eta = TMath::Abs(mu->TrackerTrk()->Eta());
      }
      else {
        pt  = mu->Pt();
        eta = TMath::Abs(mu->Eta());
      }
      break;
    case kGlobalorTracker:
      pass = mu->HasGlobalTrk() || mu->IsTrackerMuon();
      if (pass && mu->TrackerTrk()) {
        pt = mu->TrackerTrk()->Pt();
        eta = TMath::Abs(mu->TrackerTrk()->Eta());
      }
      else{
        pt = mu->Pt();
        eta = TMath::Abs(mu->Eta());
      }
    case kGlobalTracker:
      pass = (mu->HasGlobalTrk() && mu->GlobalTrk()->Chi2()/mu->GlobalTrk()->Ndof() < 10 &&
              (mu->NSegments() > 1 || mu->NMatches() > 1) && mu->NValidHits() > 0) ||
        (mu->IsTrackerMuon() &&
         mu->Quality().Quality(MuonQuality::TMLastStationTight));
      if (pass) {
        pt  = mu->TrackerTrk()->Pt();
        eta = TMath::Abs(mu->TrackerTrk()->Eta());
      }
      else {
        pt  = mu->Pt();
        eta = TMath::Abs(mu->Eta());
      }
      break;
    case kSta:
      pass = mu->HasStandaloneTrk();
      if (pass) {
        pt  = mu->StandaloneTrk()->Pt();
        eta = TMath::Abs(mu->StandaloneTrk()->Eta());
      }
      break;
    case kTrackerMuon:
      pass = mu->HasTrackerTrk() && mu->IsTrackerMuon() &&
        mu->Quality().Quality(MuonQuality::TrackerMuonArbitrated);
      if (pass) {
        pt  = mu->TrackerTrk()->Pt();
        eta = TMath::Abs(mu->TrackerTrk()->Eta());
      }
      break;
    case kCaloMuon:
      pass = mu->HasTrackerTrk() && mu->IsCaloMuon();
      if (pass) {
        pt  = mu->TrackerTrk()->Pt();
        eta = TMath::Abs(mu->TrackerTrk()->Eta());
      }
      break;
    case kTrackerBased:
      pass = mu->HasTrackerTrk();
      if (pass) {
        pt  = mu->TrackerTrk()->Pt();
        eta = TMath::Abs(mu->TrackerTrk()->Eta());
      }
      break;
    case kGlobalOnly:
      pass = mu->HasGlobalTrk();
      if (pass && mu->TrackerTrk()) {
        pt  = mu->TrackerTrk()->Pt();
        eta = TMath::Abs(mu->TrackerTrk()->Eta());
      }
      else {
        pt  = mu->Pt();
        eta = TMath::Abs(mu->Eta());
      }
      break;
    default:
      break;
    }

    if (!pass)
      continue;

    if (pt <= fMuonPtMin)
      continue;

    if (eta >= fEtaCut)
      continue;

    Double_t RChi2 = 0.0;

    if     (mu->HasGlobalTrk())
      RChi2 = mu->GlobalTrk()->Chi2()/mu->GlobalTrk()->Ndof();
    else if (mu->BestTrk() != 0)
      RChi2 = mu->BestTrk()->Chi2()/mu->BestTrk()->Ndof();

    Bool_t idpass = kFALSE;

    switch (fMuIDType) {

    case kWMuId:
      idpass = mu->BestTrk() != 0 &&
                 mu->BestTrk()->NHits() > 10 &&
                 RChi2 < 10.0 &&
                (mu->NSegments() > 1 || mu->NMatches() > 1) &&
                 mu->BestTrk()->NPixelHits() > 0 &&
                 mu->Quality().Quality(MuonQuality::GlobalMuonPromptTight);
      break;
    case kZMuId:
      idpass = mu->BestTrk() != 0 &&
               mu->BestTrk()->NHits() > 10 &&
              (mu->NSegments() > 1 || mu->NMatches() > 1) &&
               mu->BestTrk()->NPixelHits() > 0 &&
               mu->Quality().Quality(MuonQuality::GlobalMuonPromptTight);
      break;
    case kLoose:
      idpass = mu->BestTrk() != 0 &&
               mu->Quality().Quality(MuonQuality::TMOneStationLoose) &&
               mu->Quality().Quality(MuonQuality::TM2DCompatibilityLoose) &&
               mu->BestTrk()->NHits() > 10 &&
               RChi2 < 10.0 &&
               mu->Quality().Quality(MuonQuality::GlobalMuonPromptTight);
      break;
    case kTight:
      idpass = mu->BestTrk() != 0 &&
               mu->NTrkLayersHit() > 5 &&
               mu->IsPFMuon() == kTRUE &&
               mu->BestTrk()->NPixelHits() > 0 &&
               RChi2 < 10.0;
      break;
    case kmuonPOG2012CutBasedIDTight:
      idpass = mu->IsGlobalMuon() &&
               mu->IsPFMuon() &&
               mu->GlobalTrk()->RChi2() < 10 &&
               mu->NValidHits() != 0 &&
               mu->NMatches() > 1    &&
               mu->BestTrk()->NPixelHits() != 0 &&
               mu->NTrkLayersHit() > 5;
     break;
    // 2012 WW analysis for 42x (there is no PFMuon link)
    case kWWMuIdV1:
      idpass = mu->BestTrk() != 0 &&
               mu->NTrkLayersHit() > 5 &&
               mu->BestTrk()->NPixelHits() > 0 &&
               mu->BestTrk()->PtErr()/mu->BestTrk()->Pt() < 0.1 &&
               mu->TrkKink() < 20.0;
      break;
    // 2010 WW analysis
    case kWWMuIdV2:
      idpass = mu->BestTrk() != 0 &&
               mu->BestTrk()->NHits() > 10 &&
               mu->BestTrk()->NPixelHits() > 0 &&
               mu->BestTrk()->PtErr()/mu->BestTrk()->Pt() < 0.1;
      break;
    // 2011 WW analysis
    case kWWMuIdV3:
      idpass = mu->BestTrk() != 0 &&
               mu->BestTrk()->NHits() > 10 &&
               mu->BestTrk()->NPixelHits() > 0 &&
               mu->BestTrk()->PtErr()/mu->BestTrk()->Pt() < 0.1 &&
               mu->TrkKink() < 20.0;
      break;
    // 2012 WW analysis
    case kWWMuIdV4:
      idpass = mu->BestTrk() != 0 &&
               mu->NTrkLayersHit() > 5 &&
               mu->IsPFMuon() == kTRUE &&
               mu->BestTrk()->NPixelHits() > 0 &&
               mu->BestTrk()->PtErr()/mu->BestTrk()->Pt() < 0.1 &&
               mu->TrkKink() < 20.0;
      break;
    case kMVAID_BDTG_IDIso:
      {
        Bool_t passDenominatorM2 = (mu->BestTrk() != 0 &&
                                    mu->BestTrk()->NHits() > 10 &&
                                    mu->BestTrk()->NPixelHits() > 0 &&
                                    mu->BestTrk()->PtErr()/mu->BestTrk()->Pt() < 0.1 &&
                                    MuonTools::PassD0Cut(mu, fVertices, 0.20, 0) &&
                                    MuonTools::PassDZCut(mu, fVertices, 0.10, 0) &&
                                    mu->TrkKink() < 20.0
          );
        idpass =  passDenominatorM2;
        //only evaluate MVA if muon passes M2 denominator to save time
        if (idpass)
	  idpass = PassMuonMVA_BDTG_IdIso(mu, fVertices->At(0), fPileupEnergyDensity);
      }
      break;
    case kNoId:
      {
	idpass = kTRUE;
      }
      break;
    default:
      break;
    }

    if (!idpass)
      continue;

    Double_t Rho = 0.;
    if (fPileupEnergyDensity) {
      const PileupEnergyDensity *rho =  fPileupEnergyDensity->At(0);
      
      switch (fTheRhoType) {
      case RhoUtilities::MIT_RHO_VORONOI_LOW_ETA:
        Rho = rho->RhoLowEta();
        break;
      case RhoUtilities::MIT_RHO_VORONOI_HIGH_ETA:
        Rho = rho->Rho();
        break;
      case RhoUtilities::MIT_RHO_RANDOM_LOW_ETA:
        Rho = rho->RhoRandomLowEta();
        break;
      case RhoUtilities::MIT_RHO_RANDOM_HIGH_ETA:
        Rho = rho->RhoRandom();
        break;
      case RhoUtilities::CMS_RHO_RHOKT6PFJETS:
        Rho = rho->RhoKt6PFJets();
        break;
      default:
        Rho = rho->Rho();
      }
      
      if ((TMath::IsNaN(fPileupEnergyDensity->At(0)->Rho()) ||
	   std::isinf(fPileupEnergyDensity->At(0)->Rho()))    )
	Rho = 0.;
    }
    
    Bool_t isocut = kFALSE;
    switch (fMuIsoType) {
    case kTrackCalo:
      isocut = (mu->IsoR03SumPt() < fTrackIsolationCut) &&
	(mu->IsoR03EmEt() + mu->IsoR03HadEt() < fCaloIsolationCut);
      break;
    case kTrackCaloCombined:
      isocut = (1.0 * mu->IsoR03SumPt() +
		1.0 * mu->IsoR03EmEt()  +
		1.0 * mu->IsoR03HadEt() < fCombIsolationCut);
      break;
    case kTrackCaloSliding:
      {
	Double_t totalIso =  mu->IsoR03SumPt() + TMath::Max(mu->IsoR03EmEt() + mu->IsoR03HadEt()
							    - Rho * TMath::Pi() * 0.3 * 0.3, 0.0);
	// trick to change the signal region cut
	double theIsoCut = fCombIsolationCut;
	if (theIsoCut < 0.20) {
	  if (mu->Pt() >  20.0)
	    theIsoCut = 0.15;
	  else
	    theIsoCut = 0.10;
	}
	if (totalIso < (mu->Pt()*theIsoCut))
	  isocut = kTRUE;
      }
      break;
    case kTrackCaloSlidingNoCorrection:
      {
	Double_t totalIso =  1.0 * mu->IsoR03SumPt() +
	  1.0 * mu->IsoR03EmEt()  +
	  1.0 * mu->IsoR03HadEt();
	// trick to change the signal region cut
	double theIsoCut = fCombIsolationCut;
	if (theIsoCut < 0.20) {
	  if (mu->Pt() > 20.0)
	    theIsoCut = 0.15;
	  else
	    theIsoCut = 0.10;
	}
	if (totalIso < (mu->Pt()*theIsoCut))
	  isocut = kTRUE;
      }
      break;
    case kCombinedRelativeConeAreaCorrected:
      {
        //const PileupEnergyDensity *rho =  fPileupEnergyDensity->At(0); // Fabian: made Rho customable
        Double_t totalIso =  mu->IsoR03SumPt() + TMath::Max(mu->IsoR03EmEt() + mu->IsoR03HadEt()
							    - Rho * TMath::Pi() * 0.3 * 0.3, 0.0);
        double theIsoCut = fCombRelativeIsolationCut;
        if (totalIso < (mu->Pt()*theIsoCut))
	  isocut = kTRUE;
      }
      break;
    case kCombinedRelativeEffectiveAreaCorrected:
      {
        Double_t tmpRho = Rho;   // Fabian: made the Rho type customable.
        //if (!(TMath::IsNaN(fPileupEnergyDensity->At(0)->Rho()) ||
	//      std::isinf(fPileupEnergyDensity->At(0)->Rho()))    )
        //  tmpRho = fPileupEnergyDensity->At(0)->Rho();
	
        isocut = (mu->IsoR03SumPt() + mu->IsoR03EmEt() + mu->IsoR03HadEt()
                  -  tmpRho*MuonTools::MuonEffectiveArea(MuonTools::kMuEMIso03, mu->Eta())
                  -  tmpRho*MuonTools::MuonEffectiveArea(MuonTools::kMuHadIso03, mu->Eta())
		  ) < (mu->Pt()* 0.40);
      }
      break;
    case kPFIso:
      {
        Double_t pfIsoCutValue = 9999;
        if (fPFIsolationCut > 0) {
          pfIsoCutValue = fPFIsolationCut;
        } else {
          if (mu->AbsEta() < 1.479) {
            if (mu->Pt() > 20) {
              pfIsoCutValue = 0.13;
            }
	    else {
              pfIsoCutValue = 0.06;
            }
          }
	  else {
            if (mu->Pt() > 20) {
              pfIsoCutValue = 0.09;
            }
	    else {
              pfIsoCutValue = 0.05;
            }
          }
        }
        Double_t totalIso =  IsolationTools::PFMuonIsolation(mu, fPFCandidates, fVertices->At(0),
							     0.1, 1.0, 0.3, 0.0, fIntRadius);
        if (totalIso < (mu->Pt()*pfIsoCutValue) )
          isocut = kTRUE;
      }
      break;
    case kPFRadialIso:
      {
        Double_t pfIsoCutValue = 9999;
        if (fPFIsolationCut > 0) {
          pfIsoCutValue = fPFIsolationCut;
        } else {
          if (mu->Pt() > 20) {
            pfIsoCutValue = 0.10;
          } else {
              pfIsoCutValue = 0.05;
            }
          }
          Double_t totalIso =  IsolationTools::PFRadialMuonIsolation(mu, fPFNoPileUpCands, 1.0, 0.3);
          if (totalIso < (mu->Pt()*pfIsoCutValue) )
            isocut = kTRUE;
      }
      break;
    case kPFIsoBetaPUCorrected:
      {
        Double_t pfIsoCutValue = 9999;
        if (fPFIsolationCut > 0) {
          pfIsoCutValue = fPFIsolationCut;
        } else {
          if (mu->Pt() > 20) {
            pfIsoCutValue = 0.2;
          } else {
            pfIsoCutValue = 0.2;
          }
        }
        Double_t totalIso =  IsolationTools::BetaMwithPUCorrection(fPFNoPileUpCands, fPFPileUpCands, mu, 0.4);

        if (totalIso < (mu->Pt()*pfIsoCutValue) )
          isocut = kTRUE;
      }
      break;
    case kPFIsoEffectiveAreaCorrected:
      {
        Double_t pfIsoCutValue = 9999;
        if (fPFIsolationCut > 0) {
          pfIsoCutValue = fPFIsolationCut;
        } else {
          pfIsoCutValue = fPFIsolationCut; //leave it like this for now
        }
        Double_t EffectiveAreaCorrectedPFIso =
	  IsolationTools::PFMuonIsolation(mu,fPFCandidates,fVertices->At(0),.1,1.,.3,0.,fIntRadius)
          - Rho * MuonTools::MuonEffectiveArea(MuonTools::kMuNeutralIso03, mu->Eta());
        //- fPileupEnergyDensity->At(0)->Rho() * 
	//  MuonTools::MuonEffectiveArea(MuonTools::kMuNeutralIso03, mu->Eta());
        isocut = EffectiveAreaCorrectedPFIso < (mu->Pt() * pfIsoCutValue);
        break;
      }
    case kPFIsoNoL:
      {
	fNonIsolatedMuons     = GetObjThisEvt<MuonCol>(fNonIsolatedMuonsName);
	fNonIsolatedElectrons = GetObjThisEvt<ElectronCol>(fNonIsolatedElectronsName);
	
	Double_t pfIsoCutValue = 9999;
	if (fPFIsolationCut > 0) {
	  pfIsoCutValue = fPFIsolationCut;
	} else {
	  if (mu->AbsEta() < 1.479) {
	    if (mu->Pt() > 20) {
	      pfIsoCutValue = 0.13;
	    } else {
	      pfIsoCutValue = 0.06;
	    }
	  } else {
	    if (mu->Pt() > 20) {
	      pfIsoCutValue = 0.09;
	    } else {
	      pfIsoCutValue = 0.05;
	    }
	  }
	}
	Double_t totalIso = IsolationTools::PFMuonIsolation(mu, fPFCandidates, fNonIsolatedMuons,
							    fNonIsolatedElectrons,
							    fVertices->At(0), 0.1, 1.0, 0.3, 0.0,
							    fIntRadius);
	if (totalIso < (mu->Pt()*pfIsoCutValue))
	  isocut = kTRUE;
      }
      break;
    case kMVAIso_BDTG_IDIso:
      {
        Double_t totalIso = IsolationTools::PFMuonIsolation(mu, fPFCandidates, fVertices->At(0),
							    0.1, 1.0, 0.3, 0.0, fIntRadius);
        isocut = (totalIso < (mu->Pt()*0.4));
      }
      break;
    case kIsoRingsV0_BDTG_Iso:
      {
        isocut = PassMuonIsoRingsV0_BDTG_Iso(mu, fVertices->At(0), fPileupEnergyDensity);
      }
      break;
    case kIsoDeltaR:
      {
        isocut = PassMuonIsoDeltaR(mu, fVertices->At(0), fPileupEnergyDensity);
      }
      break;
    case kNoIso:
      isocut = kTRUE;
      break;
    case kCustomIso:
    default:
      break;
    }
    
    if (isocut == kFALSE)
      continue;

    // apply d0 cut
    if (fApplyD0Cut) {
      Bool_t passD0cut = kTRUE;
      if (fD0Cut < 0.05) { // trick to change the signal region cut
        if      (mu->Pt() >  20.0)
          fD0Cut = 0.02;
        else if (mu->Pt() <= 20.0)
          fD0Cut = 0.01;
      }
      if (fWhichVertex >= -1)
        passD0cut = MuonTools::PassD0Cut(mu, fVertices, fD0Cut, fWhichVertex);
      else
        passD0cut = MuonTools::PassD0Cut(mu, fBeamSpot, fD0Cut);
      if (!passD0cut)
        continue;
    }

    // apply dz cut
    if (fApplyDZCut) {
      Bool_t passDZcut = MuonTools::PassDZCut(mu, fVertices, fDZCut, fWhichVertex);
      if (!passDZcut)
        continue;
    }
    
    // mark good muon and add it
    mu->Mark();
    CleanMuons->Add(mu);
  }

  // sort according to pt
  CleanMuons->Sort();

  // add objects for other modules to use
  AddObjThisEvt(CleanMuons);
}

//--------------------------------------------------------------------------------------------------
void MuonIDMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here, we just request the
  // muon collection branch.

  using std::string;

  if (fCleanMuonsName.CompareTo("HggLeptonTagMuons") == 0) {
    ReqEventObject(fPVName,fVertices,true);
  }

   // In this case we cannot have a branch
  if (fMuonIsoType.CompareTo("PFIsoNoL") != 0) {
    ReqEventObject(fMuonBranchName, fMuons, kTRUE);
  }
  ReqEventObject(fBeamSpotName, fBeamSpot, kTRUE);
  ReqEventObject(fTrackName, fTracks, kTRUE);
  ReqEventObject(fPFCandidatesName, fPFCandidates, kTRUE);
  if (fMuonIsoType.CompareTo("TrackCaloSliding") == 0
      || fMuonIsoType.CompareTo("CombinedRelativeConeAreaCorrected") == 0
      || fMuonIsoType.CompareTo("CombinedRelativeEffectiveAreaCorrected") == 0
      || fMuonIsoType.CompareTo("PFIsoEffectiveAreaCorrected") == 0
      || fMuonIsoType.CompareTo("MVA_BDTG_IDIso") == 0
      || fMuonIsoType.CompareTo("IsoRingsV0_BDTG_Iso") == 0
      || fMuonIsoType.CompareTo("IsoDeltaR") == 0
      ) {
    ReqEventObject(fPileupEnergyDensityName, fPileupEnergyDensity, kTRUE);
  }


  if (fMuonIDType.CompareTo("WMuId") == 0)
    fMuIDType = kWMuId;
  else if (fMuonIDType.CompareTo("ZMuId") == 0)
    fMuIDType = kZMuId;
  else if (fMuonIDType.CompareTo("Tight") == 0)
    fMuIDType = kTight;
  else if (fMuonIDType.CompareTo("muonPOG2012CutBasedIDTight") == 0)
    fMuIDType = kmuonPOG2012CutBasedIDTight;
  else if (fMuonIDType.CompareTo("Loose") == 0)
    fMuIDType = kLoose;
  else if (fMuonIDType.CompareTo("WWMuIdV1") == 0)
    fMuIDType = kWWMuIdV1;
  else if (fMuonIDType.CompareTo("WWMuIdV2") == 0)
    fMuIDType = kWWMuIdV2;
  else if (fMuonIDType.CompareTo("WWMuIdV3") == 0)
    fMuIDType = kWWMuIdV3;
  else if (fMuonIDType.CompareTo("WWMuIdV4") == 0)
    fMuIDType = kWWMuIdV4;
  else if (fMuonIDType.CompareTo("NoId") == 0)
    fMuIDType = kNoId;
  else if (fMuonIDType.CompareTo("Custom") == 0) {
    fMuIDType = kCustomId;
    SendError(kWarning, "SlaveBegin","Custom muon identification is not yet implemented.");
  }
  else if (fMuonIDType.CompareTo("MVA_BDTG_IDIso") == 0)
    fMuIDType = kMVAID_BDTG_IDIso;
  else {
    SendError(kAbortAnalysis,"SlaveBegin","The specified muon identification %s is not defined.",
              fMuonIDType.Data());
    return;
  }
    
  if (fMuonIsoType.CompareTo("TrackCalo") == 0)
    fMuIsoType = kTrackCalo;
  else if (fMuonIsoType.CompareTo("TrackCaloCombined") == 0)
    fMuIsoType = kTrackCaloCombined;
  else if (fMuonIsoType.CompareTo("TrackCaloSliding") == 0)
    fMuIsoType = kTrackCaloSliding;
  else if (fMuonIsoType.CompareTo("TrackCaloSlidingNoCorrection") == 0)
    fMuIsoType = kTrackCaloSlidingNoCorrection;
  else if (fMuonIsoType.CompareTo("CombinedRelativeConeAreaCorrected") == 0)
    fMuIsoType = kCombinedRelativeConeAreaCorrected;
  else if (fMuonIsoType.CompareTo("CombinedRelativeEffectiveAreaCorrected") == 0)
    fMuIsoType = kCombinedRelativeEffectiveAreaCorrected;
  else if (fMuonIsoType.CompareTo("PFIso") == 0)
    fMuIsoType = kPFIso;
  else if (fMuonIsoType.CompareTo("PFRadialIso") == 0)
    fMuIsoType = kPFRadialIso;
  else if (fMuonIsoType.CompareTo("PFIsoBetaPUCorrected") == 0)
    fMuIsoType = kPFIsoBetaPUCorrected;
  else if (fMuonIsoType.CompareTo("PFIsoEffectiveAreaCorrected") == 0)
    fMuIsoType = kPFIsoEffectiveAreaCorrected;
  else if (fMuonIsoType.CompareTo("PFIsoNoL") == 0)
    fMuIsoType = kPFIsoNoL;
  else if (fMuonIsoType.CompareTo("NoIso") == 0)
    fMuIsoType = kNoIso;
  else if (fMuonIsoType.CompareTo("Custom") == 0) {
    fMuIsoType = kCustomIso;
    SendError(kWarning, "SlaveBegin","Custom muon isolation is not yet implemented.");
  }
  else if (fMuonIsoType.CompareTo("MVA_BDTG_IDIso") == 0) 
    fMuIsoType = kMVAIso_BDTG_IDIso;
  else if (fMuonIsoType.CompareTo("IsoRingsV0_BDTG_Iso") == 0)
    fMuIsoType = kIsoRingsV0_BDTG_Iso;
  else if (fMuonIsoType.CompareTo("IsoDeltaR") == 0)
    fMuIsoType = kIsoDeltaR;
  else {
    SendError(kAbortAnalysis, "SlaveBegin","The specified muon isolation %s is not defined.",
              fMuonIsoType.Data());
    return;
  }

  if (fMuonClassType.CompareTo("All") == 0)
    fMuClassType = kAll;
  else if (fMuonClassType.CompareTo("Global") == 0)
    fMuClassType = kGlobal;
  else if (fMuonClassType.CompareTo("GlobalTracker") == 0)
    fMuClassType = kGlobalTracker;
  else if (fMuonClassType.CompareTo("Standalone") == 0)
    fMuClassType = kSta;
  else if (fMuonClassType.CompareTo("TrackerMuon") == 0)
    fMuClassType = kTrackerMuon;
  else if (fMuonClassType.CompareTo("CaloMuon") == 0)
    fMuClassType = kCaloMuon;
  else if (fMuonClassType.CompareTo("TrackerBased") == 0)
    fMuClassType = kTrackerBased;
  else if (fMuonClassType.CompareTo("GlobalOnly") == 0)
    fMuClassType = kGlobalOnly;
  else if (fMuonClassType.CompareTo("GlobalorTracker") == 0)
    fMuClassType = kGlobalorTracker;
  else {
    SendError(kAbortAnalysis, "SlaveBegin","The specified muon class %s is not defined.",
              fMuonClassType.Data());
    return;
  }

  // If we use MVA ID, need to load MVA weights
  if     (fMuIsoType == kMVAIso_BDTG_IDIso || fMuIDType == kMVAID_BDTG_IDIso) {
    fMuonTools = new MuonTools();
    fMuonIDMVA = new MuonIDMVA();
    fMuonIDMVA->Initialize("BDTG method",
                           string((getenv("MIT_DATA")+string("/MuonMVAWeights/BarrelPtBin0_IDIsoCombined_BDTG.weights.xml"))),
                           string((getenv("MIT_DATA")+string("/MuonMVAWeights/EndcapPtBin0_IDIsoCombined_BDTG.weights.xml"))),
                           string((getenv("MIT_DATA")+string("/MuonMVAWeights/BarrelPtBin1_IDIsoCombined_BDTG.weights.xml"))),
                           string((getenv("MIT_DATA")+string("/MuonMVAWeights/EndcapPtBin1_IDIsoCombined_BDTG.weights.xml"))),
                           string((getenv("MIT_DATA")+string("/MuonMVAWeights/BarrelPtBin2_IDIsoCombined_BDTG.weights.xml"))),
                           string((getenv("MIT_DATA")+string("/MuonMVAWeights/EndcapPtBin2_IDIsoCombined_BDTG.weights.xml"))),
                           MuonIDMVA::kIDIsoCombinedDetIso,
                           fTheRhoType);
  }
  else if (fMuIsoType == kIsoRingsV0_BDTG_Iso) {
    std::vector<std::string> muonidiso_weightfiles;
    muonidiso_weightfiles.push_back(string((getenv("MIT_DATA")+string("/MuonMVAWeights/MuonIsoMVA_BDTG_V0_barrel_lowpt.weights.xml"))));
    muonidiso_weightfiles.push_back(string((getenv("MIT_DATA")+string("/MuonMVAWeights/MuonIsoMVA_BDTG_V0_barrel_highpt.weights.xml"))));
    muonidiso_weightfiles.push_back(string((getenv("MIT_DATA")+string("/MuonMVAWeights/MuonIsoMVA_BDTG_V0_endcap_lowpt.weights.xml"))));
    muonidiso_weightfiles.push_back(string((getenv("MIT_DATA")+string("/MuonMVAWeights/MuonIsoMVA_BDTG_V0_endcap_highpt.weights.xml"))));
    muonidiso_weightfiles.push_back(string((getenv("MIT_DATA")+string("/MuonMVAWeights/MuonIsoMVA_BDTG_V0_tracker.weights.xml"))));
    muonidiso_weightfiles.push_back(string((getenv("MIT_DATA")+string("/MuonMVAWeights/MuonIsoMVA_BDTG_V0_global.weights.xml"))));
    fMuonTools = new MuonTools();
    fMuonIDMVA = new MuonIDMVA();
    fMuonIDMVA->Initialize("MuonIso_BDTG_IsoRings",
			   MuonIDMVA::kIsoRingsV0,
			   kTRUE,
			   muonidiso_weightfiles,
			   fTheRhoType);
  }
  else if (fMuIsoType == kIsoDeltaR) {
    std::vector<std::string> muonidiso_weightfiles;
    muonidiso_weightfiles.push_back(string((getenv("MIT_DATA")+string("/MuonMVAWeights/MuonIsoMVA_santi-V1_LB_BDT.weights.xml"))));
    muonidiso_weightfiles.push_back(string((getenv("MIT_DATA")+string("/MuonMVAWeights/MuonIsoMVA_santi-V1_LE_BDT.weights.xml"))));
    muonidiso_weightfiles.push_back(string((getenv("MIT_DATA")+string("/MuonMVAWeights/MuonIsoMVA_santi-V1_HB_BDT.weights.xml"))));
    muonidiso_weightfiles.push_back(string((getenv("MIT_DATA")+string("/MuonMVAWeights/MuonIsoMVA_santi-V1_HE_BDT.weights.xml"))));
    fMuonTools = new MuonTools();
    fMuonIDMVA = new MuonIDMVA();
    fMuonIDMVA->Initialize("muonHZZ2012IsoDRMVA",
			   MuonIDMVA::kIsoDeltaR,
			   kTRUE,
			   muonidiso_weightfiles,
			   fTheRhoType);
  }
}


//--------------------------------------------------------------------------------------------------
Bool_t MuonIDMod::PassMuonMVA_BDTG_IdIso(const Muon *mu, const Vertex *vertex,
                                         const PileupEnergyDensityCol *PileupEnergyDensity) const
{

  const Track *muTrk=0;
  if      (mu->HasTrackerTrk())
    muTrk = mu->TrackerTrk();
  else if (mu->HasStandaloneTrk())
    muTrk = mu->StandaloneTrk();

  Double_t MVAValue = fMuonIDMVA->MVAValue(mu,vertex,fMuonTools,fPFCandidates,PileupEnergyDensity);

  Int_t subdet = 0;
  if (fabs(muTrk->Eta()) < 1.479)
    subdet = 0;
  else
    subdet = 1;
  Int_t ptBin = 0;
  if (muTrk->Pt() > 14.5)
    ptBin = 1;
  if (muTrk->Pt() > 20.0)
    ptBin = 2;

  Int_t MVABin = -1;
  if      (subdet == 0 && ptBin == 0)
    MVABin = 0;
  else if (subdet == 1 && ptBin == 0)
    MVABin = 1;
  else if (subdet == 0 && ptBin == 1)
    MVABin = 2;
  else if (subdet == 1 && ptBin == 1)
    MVABin = 3;
  else if (subdet == 0 && ptBin == 2)
    MVABin = 4;
  else if (subdet == 1 && ptBin == 2)
    MVABin = 5;
  
  Double_t MVACut = -999;
  if      (MVABin == 0)
    MVACut = -0.5618;
  else if (MVABin == 1)
    MVACut = -0.3002;
  else if (MVABin == 2)
    MVACut = -0.4642;
  else if (MVABin == 3)
    MVACut = -0.2478;
  else if (MVABin == 4)
    MVACut =  0.1706;
  else if (MVABin == 5)
    MVACut =  0.8146;

  if (MVAValue > MVACut)
    return kTRUE;

  return kFALSE;
}

//--------------------------------------------------------------------------------------------------
Bool_t MuonIDMod::PassMuonIsoRingsV0_BDTG_Iso(const Muon *mu, const Vertex *vertex,
                                              const PileupEnergyDensityCol *PileupEnergyDensity)
  const
{

  Bool_t isDebug = kFALSE;
  const Track *muTrk = 0;
  if      (mu->HasTrackerTrk())
    muTrk = mu->TrackerTrk();
  else if (mu->HasStandaloneTrk())
    muTrk = mu->StandaloneTrk();

  ElectronOArr *tempElectrons = new  ElectronOArr;
  MuonOArr     *tempMuons     = new  MuonOArr;
  Double_t MVAValue = fMuonIDMVA->MVAValue(mu,vertex,fMuonTools,fPFCandidates,
                      PileupEnergyDensity,MuonTools::kMuEAFall11MC,tempElectrons,tempMuons,isDebug);
  delete tempElectrons;
  delete tempMuons;

  Int_t MVABin = fMuonIDMVA->GetMVABin(muTrk->Eta(), muTrk->Pt(), mu->IsGlobalMuon(),
				       mu->IsTrackerMuon());

  Double_t MVACut = -1.0;
  Double_t eta = mu->AbsEta();

  if      (mu->Pt() <  20 && eta <  1.479)
    MVACut = 0.86;
  else if (mu->Pt() <  20 && eta >= 1.479)
    MVACut = 0.82;
  else if (mu->Pt() >= 20 && eta <  1.479)
    MVACut = 0.82;
  else if (mu->Pt() >= 20 && eta >= 1.479)
    MVACut = 0.86;

  if (fPFIsolationCut > -1.0)
    MVACut = fPFIsolationCut;

  if (isDebug == kTRUE) {
    printf("PassMuonIsoRingsV0_BDTG_Iso: %d, pt,eta= %f,%f, rho= %f(%f): RingsMVA= %f, bin: %d\n",
           GetEventHeader()->EvtNum(),mu->Pt(), mu->Eta(),
           fPileupEnergyDensity->At(0)->Rho(),fPileupEnergyDensity->At(0)->RhoKt6PFJets(),
	   MVAValue,MVABin);
  }

  if (MVAValue > MVACut)
    return kTRUE;

  return kFALSE;
}

//--------------------------------------------------------------------------------------------------
Bool_t MuonIDMod::PassMuonIsoDeltaR(const Muon *mu, const Vertex *vertex,
                                    const PileupEnergyDensityCol *PileupEnergyDensity) const
{
  const Track *muTrk=0;
  if (mu->HasTrackerTrk())
    muTrk = mu->TrackerTrk();
  else if (mu->HasStandaloneTrk())
    muTrk = mu->StandaloneTrk();

  ElectronOArr *tempElectrons = new  ElectronOArr;
  MuonOArr     *tempMuons     = new  MuonOArr;
  Double_t MVAValue = fMuonIDMVA->MVAValue(mu,vertex,fMuonTools,fPFNoPileUpCands,
                      PileupEnergyDensity,MuonTools::kMuEAFall11MC,tempElectrons,tempMuons,kFALSE);
  delete tempElectrons;
  delete tempMuons;

  Int_t MVABin = fMuonIDMVA->GetMVABin(muTrk->Eta(), muTrk->Pt(), mu->IsGlobalMuon(),
				       mu->IsTrackerMuon());

  Double_t MVACut = -999;
  if      (MVABin == 0)
    MVACut =  0.000;
  else if (MVABin == 1)
    MVACut =  0.000;
  else if (MVABin == 2)
    MVACut =  0.000;
  else if (MVABin == 3)
    MVACut =  0.000;

  if (MVAValue > MVACut)
    return kTRUE;

  return kFALSE;
}

//--------------------------------------------------------------------------------------------------
void MuonIDMod::Terminate()
{
  // Run finishing code on the computer (slave) that did the analysis
  delete fMuonIDMVA;
  delete fMuonTools;
}
