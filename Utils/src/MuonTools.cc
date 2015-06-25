#include "MitPhysics/Utils/interface/MuonTools.h"
#include "MitPhysics/Utils/interface/IsolationTools.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/DataTree/interface/BeamSpotCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityCol.h"
#include <limits>
#include <TFile.h>

ClassImp(mithep::MuonTools)

using namespace mithep;

//--------------------------------------------------------------------------------------------------
MuonTools::MuonTools(const char *mutemp, const char *pitemp) :
  fIsInit(kFALSE),
  fmuon_em_etaEmi(0),
  fmuon_had_etaEmi(0),
  fmuon_had_etaTmi(0),
  fmuon_em_etaB(0),
  fmuon_had_etaB(0),
  fmuon_ho_etaB(0),
  fmuon_had_etaTpl(0),
  fmuon_em_etaEpl(0),
  fmuon_had_etaEpl(0),
  fpion_em_etaEmi(0),
  fpion_had_etaEmi(0),
  fpion_had_etaTmi(0),
  fpion_em_etaB(0),
  fpion_had_etaB(0),
  fpion_ho_etaB(0),
  fpion_had_etaTpl(0),
  fpion_em_etaEpl(0),
  fpion_had_etaEpl(0)
{
  // Constructor.

  if (mutemp && pitemp)
    Init(mutemp, pitemp);
}

//--------------------------------------------------------------------------------------------------
MuonTools::~MuonTools() 
{
  // Destructor.

  DeleteHistos();
}

//--------------------------------------------------------------------------------------------------
void MuonTools::DeleteHistos()
{
  // Delete histograms.

  if (fIsInit) {
    delete fpion_em_etaEmi; 
    delete fpion_had_etaEmi;
    delete fpion_had_etaTmi;
    delete fpion_em_etaB;
    delete fpion_had_etaB;
    delete fpion_ho_etaB;
    delete fpion_had_etaTpl;
    delete fpion_em_etaEpl;
    delete fpion_had_etaEpl;
    delete fmuon_em_etaEmi;
    delete fmuon_had_etaEmi;
    delete fmuon_had_etaTmi;
    delete fmuon_em_etaB;
    delete fmuon_had_etaB;
    delete fmuon_ho_etaB;
    delete fmuon_had_etaTpl;
    delete fmuon_em_etaEpl;
    delete fmuon_had_etaEpl;
    fpion_em_etaEmi  = 0;
    fpion_had_etaEmi = 0;
    fpion_had_etaTmi = 0;
    fpion_em_etaB    = 0;
    fpion_had_etaB   = 0;
    fpion_ho_etaB    = 0;
    fpion_had_etaTpl = 0;
    fpion_em_etaEpl  = 0;
    fpion_had_etaEpl = 0;
    fmuon_em_etaEmi  = 0;
    fmuon_had_etaEmi = 0;
    fmuon_had_etaTmi = 0;
    fmuon_em_etaB    = 0;
    fmuon_had_etaB   = 0;
    fmuon_ho_etaB    = 0;
    fmuon_had_etaTpl = 0;
    fmuon_em_etaEpl  = 0;
    fmuon_had_etaEpl = 0;
    fIsInit = kFALSE;
  }
}

//--------------------------------------------------------------------------------------------------
Double_t MuonTools::GetCaloCompatability(const Muon *iMuon,
                                         Bool_t iEMSpecial, Bool_t iCorrectedHCAL) const
{
  // Get calo compatibility value for given muon based on calorimeter templates.
  // If iEMSpecial is true, then a use different arrangement of ECAL for compatibility.

  Double_t lEta = iMuon->Eta();
  Double_t aEta = TMath::Abs(lEta);
  if (aEta > 2.5) 
    return 0.5; 

  Double_t lP = iMuon->P();
  if (lP >= 2000.) 
    lP = 1999.9;
  if(lP < 0. )           
    return 0.5; 

  Double_t lEM  = -5.;      
  if (!iEMSpecial || iMuon->EmEnergy() != 0.) 
    lEM  = iMuon->EmEnergy();

  Double_t lHad = iMuon->HadEnergy();
  Double_t lHO = iMuon->HoEnergy();;

  TH2D *lTMuonHad = 0;
  TH2D *lTPionHad = 0;
  TH2D *lTMuonHo  = 0;
  TH2D *lTPionHo  = 0;
  TH2D *lTMuonEm  = 0;
  TH2D *lTPionEm  = 0;
    
  if (aEta >= 1.27) {
    if (iCorrectedHCAL) 
      lHad *= 1.8/2.2;
    if (lEta > 0) {
      lTPionHad = fpion_had_etaEpl;
      lTMuonHad = fmuon_had_etaEpl;
    } else {
      lTPionHad = fpion_had_etaEmi;
      lTMuonHad = fmuon_had_etaEmi;
    }
  }

  if (aEta < 1.27 && aEta >= 1.1) {
    if (iCorrectedHCAL)    
      lHad *= (1.8/(-2.2*aEta+5.5));
    if (lEta > 0) {
      lTPionHad  = fpion_had_etaTpl;
      lTMuonHad  = fmuon_had_etaTpl;
    } else {
      lTPionHad  = fpion_had_etaTmi;
      lTMuonHad  = fmuon_had_etaTmi;
    }
  }

  if (aEta < 1.1) {
    if(iCorrectedHCAL)    
      lHad *= TMath::Sin(2*TMath::ATan(TMath::Exp(lEta)));
    lTPionHad  = fpion_had_etaB;
    lTMuonHad  = fmuon_had_etaB;
  }
  if (lEta > 1.479) {
    lTPionEm = fpion_em_etaEpl;
    lTMuonEm = fmuon_em_etaEpl;
  }
  if (aEta <= 1.479) {
    lTPionEm  = fpion_em_etaB;
    lTMuonEm  = fmuon_em_etaB;
  }
  if (lEta < -1.479) {
    lTPionEm  = fpion_em_etaEmi;
    lTMuonEm  = fmuon_em_etaEmi;
  }
  if (aEta < 1.28) {
    lTPionHo  = fpion_ho_etaB;
    lTMuonHo  = fmuon_ho_etaB;
  }
  
  Double_t lPBX = 1.;     
  Double_t lPSX = 1.; 
  Double_t lPBY = 1.;     
  Double_t lPSY = 1.; 
  Double_t lPBZ = 1.;     
  Double_t lPSZ = 1.; 
  if (!Overflow(lTPionEm, lP,lEM))  
    lPBX = lTPionEm ->GetBinContent(lTPionEm ->GetXaxis()->FindBin(lP),
                                    lTPionEm ->GetYaxis()->FindBin(lEM));
  if (!Overflow(lTPionHad,lP,lHad)) 
    lPBY = lTPionHad->GetBinContent(lTPionHad->GetXaxis()->FindBin(lP),
                                    lTPionHad->GetYaxis()->FindBin(lHad));
  if (!Overflow(lTPionHo, lP,lHO))  
    lPBZ = lTPionHo ->GetBinContent(lTPionHo ->GetXaxis()->FindBin(lP),
                                    lTPionHo ->GetYaxis()->FindBin(lHO));
  if (!Overflow(lTMuonEm, lP,lEM )) 
    lPSX = lTMuonEm ->GetBinContent(lTMuonEm ->GetXaxis()->FindBin(lP),
                                    lTMuonEm ->GetYaxis()->FindBin(lEM));
  if (!Overflow(lTMuonHad,lP,lHad)) 
    lPSY = lTMuonHad->GetBinContent(lTMuonHad->GetXaxis()->FindBin(lP),
                                    lTMuonHad->GetYaxis()->FindBin(lHad));
  if (!Overflow(lTMuonHo ,lP,lHO))  
    lPSZ =  lTMuonHo ->GetBinContent(lTMuonHo ->GetXaxis()->FindBin(lP),
                                     lTMuonHo ->GetYaxis()->FindBin(lHO));
  
  if (lPSX == 0. || lPBX == 0. || (lEM <= 0. && !iEMSpecial)) {
    lPSX = 1.; 
    lPBX = 1.;
  } 
  if (lPSY == 0. || lPBY == 0. || lHad == 0.) {
    lPSY = 1.; 
    lPBY = 1.;
  }
  if (lPSZ == 0. || lPBZ == 0. || lHO  == 0.) {
    lPSZ = 1.; 
    lPBZ = 1.;
  }
  if ((lPSX*lPSY*lPSZ+lPBX*lPBY*lPBZ) > 0.) 
    return lPSX*lPSY*lPSZ/(lPSX*lPSY*lPSZ+lPBX*lPBY*lPBZ);

  return 0.5;
}

//--------------------------------------------------------------------------------------------------
Bool_t MuonTools::Init(const char *mutemp, const char *pitemp)
{
  // Read histograms from given files.

  if (fIsInit) {
    DeleteHistos();
  }

  TDirectory::TContext context(0);

  TFile *muon_templates = TFile::Open(mutemp);
  if (!muon_templates) {
    Fatal("Init", "Could not open file %s", mutemp);
    return kFALSE;
  }
  fmuon_em_etaEmi  = LoadHisto("em_etaEmi",  muon_templates);
  fmuon_had_etaEmi = LoadHisto("had_etaEmi", muon_templates);
  fmuon_had_etaTmi = LoadHisto("had_etaTmi", muon_templates);
  fmuon_em_etaB    = LoadHisto("em_etaB",    muon_templates);
  fmuon_had_etaB   = LoadHisto("had_etaB",   muon_templates);
  fmuon_ho_etaB    = LoadHisto("ho_etaB",    muon_templates);
  fmuon_had_etaTpl = LoadHisto("had_etaTpl", muon_templates);
  fmuon_em_etaEpl  = LoadHisto("em_etaEpl",  muon_templates);
  fmuon_had_etaEpl = LoadHisto("had_etaEpl", muon_templates);
  muon_templates->Close();
  delete muon_templates;

  TFile *pion_templates = TFile::Open(pitemp);
  if (!pion_templates) {
    Fatal("Init", "Could not open file %s", pitemp);
    return kFALSE;
  }

  fpion_em_etaEmi  = LoadHisto("em_etaEmi",  pion_templates);
  fpion_had_etaEmi = LoadHisto("had_etaEmi", pion_templates);
  fpion_had_etaTmi = LoadHisto("had_etaTmi", pion_templates);
  fpion_em_etaB    = LoadHisto("em_etaB",    pion_templates);
  fpion_had_etaB   = LoadHisto("had_etaB",   pion_templates);
  fpion_ho_etaB    = LoadHisto("ho_etaB",    pion_templates);
  fpion_had_etaTpl = LoadHisto("had_etaTpl", pion_templates);
  fpion_em_etaEpl  = LoadHisto("em_etaEpl",  pion_templates);
  fpion_had_etaEpl = LoadHisto("had_etaEpl", pion_templates);
  pion_templates->Close();
  delete pion_templates;

  fIsInit = kTRUE;
  return kTRUE;
}

//--------------------------------------------------------------------------------------------------
Bool_t MuonTools::IsGood(const mithep::Muon *iMuon, ESelType iSel) const
{
  // Return true if given muon qualifies given selection criterium.

  Double_t tm2dcut = 0.;

  switch(iSel) {
  case kAllArbitrated:
    if (iMuon->StandaloneTrk() != 0 || iMuon->GlobalTrk()!= 0)  
      return kTRUE;
    if (iMuon->NSegments() > 0) 
      return kTRUE;
    break;
  case kPromptTight: 
    return iMuon->PromptTight(Muon::kAny);
    break;
  case kTMOneStationLoose:
    return iMuon->TMOneStation(999999,999999);
    break;
  case kTMOneStationTight:
    return iMuon->TMOneStation();
    break;
  case kTMLastStationLoose: 
    return iMuon->TMLastStation(999999,999999);
    break;
  case kTMLastStationTight: 
    return iMuon->TMLastStation();
    break;
  case kTM2DCompatibilityLoose: 
    tm2dcut = 0.7;
    break;
  case kTM2DCompatibilityTight: 
    tm2dcut = 1.0;
    break;
  default:
    return kFALSE;
    break;
  }

  Double_t lVal = GetSegmentCompatability(iMuon); 
  if (lVal == 0.5) // exclude this border case
    return kFALSE;

  lVal *= 1.2;
  lVal += 0.8*GetCaloCompatability(iMuon,kTRUE,kTRUE);
  if (lVal > tm2dcut) 
    return kTRUE;

  return kFALSE;

}

//--------------------------------------------------------------------------------------------------
Double_t MuonTools::GetSegmentCompatability(const mithep::Muon *iMuon) const
{
  // Get segment compatability for given muon based on likelihood of well defined 
  // track through chambers.

  Int_t lNStationsCrossed = 0;
  Int_t lNStationsSegment = 0;
  
  Int_t lStSegmentmatch[8];
  Int_t lStCrossed[8];
  Double_t lStBoundary[8];

  Double_t lWeight  = 0.;
  for (Int_t i0 = 0; i0 < 8; ++i0) {
    lStBoundary[i0] = 0.;
    if(iMuon->GetTrackDist(i0) < 999999. ) { 
      lNStationsCrossed++;
      lStCrossed[i0]  = 1;
      if (iMuon->GetTrackDist(i0) > -10. ) 
        lStBoundary[i0] = iMuon->GetTrackDist(i0); 
    } else
      lStCrossed[i0]  = 0;

    if(iMuon->GetDX(i0) < 999999.) {
      lNStationsSegment++;
      lStSegmentmatch[i0] = 1;
    } else
      lStSegmentmatch[i0] = 0;

  }

  if (lNStationsCrossed == 0)
    return 0.5;

  Double_t lStWeight[8];
  Int_t lPCross = -1;
  const Double_t lAtWeight = 0.5; 
  for (Int_t i0 = 0; i0< 8; ++i0) { 
    lStWeight[i0] = 0;
    if (lStCrossed[i0] > 0) { 
      lPCross++;

      switch (lNStationsCrossed) { 
      case 1 : 
        lStWeight[i0] =  1.;
        break;
      case 2 :
        if (lPCross == 0 ) 
          lStWeight[i0] = 0.33;
        else 
          lStWeight[i0] = 0.67;
        break;
      case 3 : 
        if (lPCross == 0) 
          lStWeight[i0] = 0.23;
        else if (lPCross == 1) 
          lStWeight[i0] = 0.33;
        else 
          lStWeight[i0] = 0.44;
        break;
      case 4 : 
        if (lPCross == 0) 
          lStWeight[i0] = 0.10;
        else if (lPCross == 1) 
          lStWeight[i0] = 0.20;
        else if (lPCross == 2) 
          lStWeight[i0] = 0.30;
        else                    
          lStWeight[i0] = 0.40;
        break;
      default : 
        lStWeight[i0] = 1./lNStationsCrossed;
      }
    
      if (lStSegmentmatch[i0] <= 0 && lStBoundary[i0] != 0.)
        lStWeight[i0] *= lAtWeight*0.5*(TMath::Erf(lStBoundary[i0]/6.)+1.); 
      else if (lStSegmentmatch[i0] <= 0 && lStBoundary[i0] == 0) 
        lStWeight[i0] = 0.;
      
      if (lStSegmentmatch[i0] > 0) { 
        Double_t lP2X = TMath::Power(iMuon->GetPullX(i0),2.);
        Double_t lP2Y = TMath::Power(iMuon->GetPullY(i0),2.);
        Double_t lD2X = TMath::Power(iMuon->GetDX(i0),2.);
        Double_t lD2Y = TMath::Power(iMuon->GetDY(i0),2.);
        if (iMuon->GetDY(i0) < 999999 && iMuon->GetDX(i0) < 999999) 
          lStWeight[i0] *= SigWeight(TMath::Sqrt(lD2X+lD2Y),TMath::Sqrt(lP2X+lP2Y));
        else if (iMuon->GetDY(i0) >= 999999 && i0 < 4)
          lStWeight[i0] *= SigWeight(iMuon->GetDX(i0),iMuon->GetPullX(i0));
        else if(i0 < 4)
          lStWeight[i0] *= SigWeight(iMuon->GetDY(i0),iMuon->GetPullY(i0));
      }
    } 
    lWeight += lStWeight[i0];
  }

  return lWeight;
}

//--------------------------------------------------------------------------------------------------
Bool_t
mithep::MuonTools::PassClass(const Muon *mu, EMuClassType classType)
{
  switch (classType) {
  case kAll:
    return true;

  case kGlobal:
    return mu->HasGlobalTrk() && mu->IsTrackerMuon();

  case kGlobalorTracker:
    return mu->HasGlobalTrk() || mu->IsTrackerMuon();

  case kGlobalTracker:
    if (mu->HasGlobalTrk() && mu->GlobalTrk()->Chi2() / mu->GlobalTrk()->Ndof() < 10 &&
        (mu->NSegments() > 1 || mu->NMatches() > 1) && mu->NValidHits() > 0)
      return true;
    if (mu->IsTrackerMuon() && mu->Quality().Quality(MuonQuality::TMLastStationTight))
      return true;

  case kSta:
    return mu->HasStandaloneTrk();

  case kTrackerMuon:
    return mu->HasTrackerTrk() && mu->IsTrackerMuon() &&
      mu->Quality().Quality(MuonQuality::TrackerMuonArbitrated);

  case kCaloMuon:
    return mu->HasTrackerTrk() && mu->IsCaloMuon();

  case kTrackerBased:
    return mu->HasTrackerTrk();

  case kGlobalOnly:
    return mu->HasGlobalTrk();

  case kPFGlobalorTracker:
    return mu->IsPFMuon() && (mu->HasGlobalTrk() || mu->IsTrackerMuon());

  default:
    break;
  }

  return false;
}

//--------------------------------------------------------------------------------------------------
void
mithep::MuonTools::MuonPtEta(const Muon *mu, EMuClassType classType, Double_t& pt, Double_t& absEta)
{
  switch (classType) {
  case kAll:
    pt = mu->Pt();
    absEta = TMath::Abs(mu->Eta());
    return;

  case kGlobal:
    if (mu->TrackerTrk()) {
      pt = mu->TrackerTrk()->Pt();
      absEta = TMath::Abs(mu->TrackerTrk()->Eta());
    }
    else {
      pt = mu->Pt();
      absEta = TMath::Abs(mu->Eta());
    }
    return;

  case kGlobalorTracker:
    if (mu->TrackerTrk()) {
      pt = mu->TrackerTrk()->Pt();
      absEta = TMath::Abs(mu->TrackerTrk()->Eta());
    }
    else {
      pt = mu->Pt();
      absEta = TMath::Abs(mu->Eta());
    }
    return;

  case kGlobalTracker:
    pt = mu->TrackerTrk()->Pt();
    absEta = TMath::Abs(mu->TrackerTrk()->Eta());
    return;

  case kSta:
    pt = mu->StandaloneTrk()->Pt();
    absEta = TMath::Abs(mu->StandaloneTrk()->Eta());
    return;

  case kTrackerMuon:
    pt = mu->TrackerTrk()->Pt();
    absEta = TMath::Abs(mu->TrackerTrk()->Eta());
    return;

  case kCaloMuon:
    pt = mu->TrackerTrk()->Pt();
    absEta = TMath::Abs(mu->TrackerTrk()->Eta());
    return;

  case kTrackerBased:
    pt = mu->TrackerTrk()->Pt();
    absEta = TMath::Abs(mu->TrackerTrk()->Eta());
    return;

  case kGlobalOnly:
    if (mu->TrackerTrk()) {
      pt = mu->TrackerTrk()->Pt();
      absEta = TMath::Abs(mu->TrackerTrk()->Eta());
    }
    else {
      pt = mu->Pt();
      absEta = TMath::Abs(mu->Eta());
    }
    return;

  default:
    break;
  }

  pt = mu->Pt();
  absEta = TMath::Abs(mu->Eta());
}

//--------------------------------------------------------------------------------------------------
Bool_t 
mithep::MuonTools::PassIso(const Muon *mu, EMuIsoType isoType)
{
  switch (isoType) {
  case kTrackCalo:
    return mu->IsoR03SumPt() < 5. && mu->IsoR03EmEt() + mu->IsoR03HadEt() < 5.;

  case kTrackCaloCombined:
    return mu->IsoR03SumPt() + mu->IsoR03EmEt()  + mu->IsoR03HadEt() < 10.;

  case kNoIso:
    return true;

  default:
    break;
  }

  return false;
}

Bool_t
mithep::MuonTools::PassIsoRhoCorr(Muon const* mu, MuonTools::EMuIsoType type, Double_t rho, PFCandidateCol const* pfCandidates/* = 0*/, Vertex const* vertex/* = 0*/)
{
  switch (type) {
  case kTrackCaloSliding:
    return mu->IsoR03SumPt() + TMath::Max(mu->IsoR03EmEt() + mu->IsoR03HadEt() - rho * TMath::Pi() * 0.3 * 0.3, 0.) < mu->Pt() * 0.15;

  case kCombinedRelativeConeAreaCorrected:
    return mu->IsoR03SumPt() + mu->IsoR03EmEt() + mu->IsoR03HadEt()
      - rho * MuonEffectiveArea(kMuEMIso03, mu->Eta())
      - rho * MuonEffectiveArea(kMuHadIso03, mu->Eta())
      < mu->Pt() * 0.40;

  case kPFIsoEffectiveAreaCorrected:
    return IsolationTools::PFMuonIsolation(mu, pfCandidates, vertex, 0.1, 1., 0.3, 0., 0.3)
      - rho * MuonEffectiveArea(kMuNeutralIso03, mu->Eta()) < mu->Pt() * 0.15;

  default:
    break;
  }

  return false;
}

Bool_t
mithep::MuonTools::PassPFIso(Muon const* mu, EMuIsoType type, PFCandidateCol const* pfCandidates,
                             Vertex const* vertex/* = 0*/,
                             PFCandidateCol const* pileupCands/* = 0*/,
                             ElectronCol const* nonisolatedElectrons/* = 0*/,
                             MuonCol const* nonisolatedMuons/* = 0*/)
{
  switch (type) {
  case kPFIso:
    {
      Double_t cutValue = 0.;
      if (mu->AbsEta() < 1.479) {
        if (mu->Pt() > 20)
          cutValue = 0.13;
        else
          cutValue = 0.06;
      }
      else {
        if (mu->Pt() > 20)
          cutValue = 0.09;
        else
          cutValue = 0.05;
      }
      return IsolationTools::PFMuonIsolation(mu, pfCandidates, vertex,
                                             0.1, 1.0, 0.3, 0.0, 0.3) < mu->Pt() * cutValue;
    }

  case kPFRadialIso:
    {
      Double_t cutValue = 0.;
      if (mu->Pt() > 20)
        cutValue = 0.10;
      else
        cutValue = 0.05;

      // pfCandidates here should be NoPileupCandidates
      return IsolationTools::PFRadialMuonIsolation(mu, pfCandidates, 1.0, 0.3) < mu->Pt() * cutValue;
    }

  case kPFIsoBetaPUCorrected:
    // pfCandidates here should be NoPileupCandidates
    return IsolationTools::BetaMwithPUCorrection(pfCandidates, pileupCands, mu, 0.4) < mu->Pt() * 0.2;

  case kPFIsoBetaPUCorrectedTight:
    // pfCandidates here should be NoPileupCandidates
    return IsolationTools::BetaMwithPUCorrection(pfCandidates, pileupCands, mu, 0.4) < mu->Pt() * 0.12;

  case kPFIsoNoL:
    {
      Double_t cutValue = 0.;
      if (mu->AbsEta() < 1.479) {
        if (mu->Pt() > 20)
          cutValue = 0.13;
        else
          cutValue = 0.06;
      }
      else {
        if (mu->Pt() > 20)
          cutValue = 0.09;
        else
          cutValue = 0.05;
      }
      return IsolationTools::PFMuonIsolation(mu, pfCandidates, nonisolatedMuons, nonisolatedElectrons,
                                             vertex, 0.1, 1.0, 0.3, 0., 0.3) < mu->Pt() * cutValue;
    }

  case kMVAIso_BDTG_IDIso:
    return IsolationTools::PFMuonIsolation(mu, pfCandidates, vertex,
                                           0.1, 1.0, 0.3, 0.0, 0.3) < mu->Pt() * 0.4;

  default:
    break;
  }

  return false;
}

//--------------------------------------------------------------------------------------------------
Bool_t
mithep::MuonTools::PassId(const Muon *mu, EMuIdType idType)
{
  Double_t normChi2 = 0.0;
  if (mu->HasGlobalTrk())
    normChi2 = mu->GlobalTrk()->Chi2() / mu->GlobalTrk()->Ndof();
  else if (mu->BestTrk() != 0)
    normChi2 = mu->BestTrk()->Chi2() / mu->BestTrk()->Ndof();

  switch (idType) {
  case kWMuId:
    return mu->BestTrk() != 0 &&
      mu->BestTrk()->NHits() > 10 &&
      normChi2 < 10.0 &&
      (mu->NSegments() > 1 || mu->NMatches() > 1) &&
      mu->BestTrk()->NPixelHits() > 0 &&
      mu->Quality().Quality(MuonQuality::GlobalMuonPromptTight);

  case kZMuId:
    return mu->BestTrk() != 0 &&
      mu->BestTrk()->NHits() > 10 &&
      (mu->NSegments() > 1 || mu->NMatches() > 1) &&
      mu->BestTrk()->NPixelHits() > 0 &&
      mu->Quality().Quality(MuonQuality::GlobalMuonPromptTight);

  case kLoose:
    return mu->BestTrk() != 0 &&
      mu->Quality().Quality(MuonQuality::TMOneStationLoose) &&
      mu->Quality().Quality(MuonQuality::TM2DCompatibilityLoose) &&
      mu->BestTrk()->NHits() > 10 &&
      normChi2 < 10.0 &&
      mu->Quality().Quality(MuonQuality::GlobalMuonPromptTight);

  case kTight:
    return mu->BestTrk() != 0 &&
      mu->NTrkLayersHit() > 5 &&
      mu->IsPFMuon() == kTRUE &&
      mu->BestTrk()->NPixelHits() > 0 &&
      normChi2 < 10.0;

  case kMuonPOG2012CutBasedIdTight:
    return mu->IsGlobalMuon() &&
      mu->IsPFMuon() &&
      mu->GlobalTrk()->RChi2() < 10 &&
      mu->NValidHits() != 0 &&
      mu->NMatches() > 1    &&
      mu->BestTrk()->NPixelHits() != 0 &&
      mu->NTrkLayersHit() > 5;

    // 2012 WW analysis for 42x (there is no PFMuon link)
  case kWWMuIdV1:
    return mu->BestTrk() != 0 &&
      mu->NTrkLayersHit() > 5 &&
      mu->BestTrk()->NPixelHits() > 0 &&
      mu->BestTrk()->PtErr()/mu->BestTrk()->Pt() < 0.1 &&
      mu->TrkKink() < 20.0;

    // 2010 WW analysis
  case kWWMuIdV2:
    return mu->BestTrk() != 0 &&
      mu->BestTrk()->NHits() > 10 &&
      mu->BestTrk()->NPixelHits() > 0 &&
      mu->BestTrk()->PtErr()/mu->BestTrk()->Pt() < 0.1;

    // 2011 WW analysis
  case kWWMuIdV3:
    return mu->BestTrk() != 0 &&
      mu->BestTrk()->NHits() > 10 &&
      mu->BestTrk()->NPixelHits() > 0 &&
      mu->BestTrk()->PtErr()/mu->BestTrk()->Pt() < 0.1 &&
      mu->TrkKink() < 20.0;

    // 2012 WW analysis
  case kWWMuIdV4:
    return mu->BestTrk() != 0 &&
      mu->NTrkLayersHit() > 5 &&
      mu->IsPFMuon() == kTRUE &&
      mu->BestTrk()->NPixelHits() > 0 &&
      mu->BestTrk()->PtErr()/mu->BestTrk()->Pt() < 0.1 &&
      mu->TrkKink() < 20.0;
    break;

  case kNoId:
    return true;

  default:
    break;
  }

  return false;
}

//--------------------------------------------------------------------------------------------------
TH2D *MuonTools::LoadHisto(const char *name, TFile *file) const
{
  // Load histogram with given name from given file and return it.

  TH2D *ret = dynamic_cast<TH2D*>(file->Get(name));
  if (!ret) {
    Fatal("LoadHisto", "Could not load histogram %s from file %s", name, file->GetName());
    return 0;
  }
  ret->SetDirectory(0);
  return ret;
}

Bool_t
MuonTools::PassD0Cut(Muon const* mu, VertexCol const* vertices, EMuIdType idType, Int_t iVertex/* = 0*/)
{
  Track const* mt = mu->BestTrk();
  if (!mt)
    return false;

  if (iVertex >= (int)vertices->GetEntries())
    iVertex = vertices->GetEntries() - 1;

  Double_t d0 = 0.;
  if(iVertex >= 0)
    d0 = TMath::Abs(mt->D0Corrected(*vertices->At(iVertex)));
  else {
    Double_t distVtx = std::numeric_limits<double>::max();
    for(UInt_t nv = 0; nv != vertices->GetEntries(); ++nv){
      double dz = TMath::Abs(mt->DzCorrected(*vertices->At(nv)));
      if(dz < distVtx) {
        distVtx = dz;
        d0 = TMath::Abs(mt->D0Corrected(*vertices->At(nv)));
      }
    }
  }

  return PassD0Cut(mu, d0, idType);
}

Bool_t
MuonTools::PassD0Cut(Muon const* mu, BeamSpotCol const* beamspots, EMuIdType idType)
{
  Track const* mt = mu->BestTrk();
  if (!mt)
    return false;

  // d0 cut
  Double_t d0 = std::numeric_limits<double>::max();
  for (UInt_t i0 = 0; i0 < beamspots->GetEntries(); ++i0) {
    Double_t pD0 = mt->D0Corrected(*beamspots->At(i0));
    if(TMath::Abs(pD0) < d0)
      d0 = TMath::Abs(pD0);
  }

  return PassD0Cut(mu, d0, idType);
}

Bool_t
MuonTools::PassD0Cut(Muon const*, Double_t d0, EMuIdType idType)
{
  switch (idType) {
  case kMuonPOG2012CutBasedIdTight:
  case kMVAID_BDTG_IDIso:
    return d0 < 0.2;
    break;
  default:
    return d0 < 0.3;
    break;
  }
}

Bool_t
MuonTools::PassDZCut(Muon const* mu, VertexCol const* vertices, EMuIdType idType, Int_t iVertex/* = 0*/)
{
  const Track *mt = mu->BestTrk();
  if (!mt)
    return false;

  if (iVertex >= (int) vertices->GetEntries())
    iVertex = vertices->GetEntries()-1;

  Double_t dz;

  if (iVertex >= 0)
    dz = TMath::Abs(mt->DzCorrected(*vertices->At(iVertex)));
  else {
    dz = std::numeric_limits<double>::max();
    for(UInt_t nv = 0; nv != vertices->GetEntries(); ++nv) {
      double test = TMath::Abs(mt->DzCorrected(*vertices->At(nv)));
      if(test < dz)
        dz = test;
    }
  }

  return PassDZCut(mu, dz, idType);
}

Bool_t
MuonTools::PassDZCut(Muon const*, Double_t dz, EMuIdType idType)
{
  switch (idType) {
  case kMuonPOG2012CutBasedIdTight:
    return dz < 0.5;

  case kMVAID_BDTG_IDIso:
    return dz < 0.1;

  default:
    return dz < 20.;
  }
}

//--------------------------------------------------------------------------------------------------
Bool_t MuonTools::PassSoftMuonCut(const Muon *mu, const VertexCol *vertices, const Double_t,
                                  const Bool_t applyIso) 
{
  if(mu->Pt() <= 3.0) return kFALSE;

  if(!mu->IsTrackerMuon()) return kFALSE;
  
  if(!mu->Quality().Quality(MuonQuality::TMLastStationAngTight)) return kFALSE;

  if(mu->NTrkLayersHit() <= 5) return kFALSE;

  if(!PassD0Cut(mu, vertices, kMuonPOG2012CutBasedIdTight, 0)) return kFALSE;

  if(!PassDZCut(mu, vertices, kMuonPOG2012CutBasedIdTight, 0)) return kFALSE;

  if(applyIso == kTRUE){
    Double_t totalIso = 1.0 * mu->IsoR03SumPt() + 
      1.0 * mu->IsoR03EmEt() + 
      1.0 * mu->IsoR03HadEt();
    if (totalIso < (mu->Pt()*0.10) && mu->Pt() > 20.0) return kFALSE;
  }

  return kTRUE;
}

Double_t MuonTools::MuonEffectiveArea(EMuonEffectiveAreaType type, Double_t Eta, 
                                      EMuonEffectiveAreaTarget EffectiveAreaTarget) {

  Double_t EffectiveArea = 0;
  if (fabs(Eta) < 1.0) {
    if (type == kMuChargedIso03) EffectiveArea = 0.000;
    if (type == kMuNeutralIso03) EffectiveArea = 0.080;
    if (type == kMuChargedIso04) EffectiveArea = 0.000;
    if (type == kMuNeutralIso04) EffectiveArea = 0.163;
    if (type == kMuHadEnergy)    EffectiveArea = 0.000;
    if (type == kMuHoEnergy)     EffectiveArea = 0.000;
    if (type == kMuEmEnergy)     EffectiveArea = 0.000;
    if (type == kMuHadS9Energy)  EffectiveArea = 0.016;
    if (type == kMuHoS9Energy)   EffectiveArea = 0.000;
    if (type == kMuEmS9Energy)   EffectiveArea = 0.000;
    if (type == kMuTrkIso03)     EffectiveArea = 0.000;
    if (type == kMuEMIso03)      EffectiveArea = 0.080;
    if (type == kMuHadIso03)     EffectiveArea = 0.025;
    if (type == kMuTrkIso05)     EffectiveArea = 0.000;
    if (type == kMuEMIso05)      EffectiveArea = 0.290;
    if (type == kMuHadIso05)     EffectiveArea = 0.091;
  } else if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) {
    if (type == kMuChargedIso03) EffectiveArea = 0.000;
    if (type == kMuNeutralIso03) EffectiveArea = 0.083;
    if (type == kMuChargedIso04) EffectiveArea = 0.000;
    if (type == kMuNeutralIso04) EffectiveArea = 0.168;
    if (type == kMuHadEnergy)    EffectiveArea = 0.005;
    if (type == kMuHoEnergy)     EffectiveArea = 0.000;
    if (type == kMuEmEnergy)     EffectiveArea = 0.000;
    if (type == kMuHadS9Energy)  EffectiveArea = 0.041;
    if (type == kMuHoS9Energy)   EffectiveArea = 0.000;
    if (type == kMuEmS9Energy)   EffectiveArea = 0.000;
    if (type == kMuTrkIso03)     EffectiveArea = 0.000;
    if (type == kMuEMIso03)      EffectiveArea = 0.043;
    if (type == kMuHadIso03)     EffectiveArea = 0.028;
    if (type == kMuTrkIso05)     EffectiveArea = 0.000;
    if (type == kMuEMIso05)      EffectiveArea = 0.184;
    if (type == kMuHadIso05)     EffectiveArea = 0.106;
  } else if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) {
    if (type == kMuChargedIso03) EffectiveArea = 0.000;
    if (type == kMuNeutralIso03) EffectiveArea = 0.060;
    if (type == kMuChargedIso04) EffectiveArea = 0.000;
    if (type == kMuNeutralIso04) EffectiveArea = 0.131;
    if (type == kMuHadEnergy)    EffectiveArea = 0.020;
    if (type == kMuHoEnergy)     EffectiveArea = 0.000;
    if (type == kMuEmEnergy)     EffectiveArea = 0.000;
    if (type == kMuHadS9Energy)  EffectiveArea = 0.072;
    if (type == kMuHoS9Energy)   EffectiveArea = 0.000;
    if (type == kMuEmS9Energy)   EffectiveArea = 0.000;
    if (type == kMuTrkIso03)     EffectiveArea = 0.000;
    if (type == kMuEMIso03)      EffectiveArea = 0.025;
    if (type == kMuHadIso03)     EffectiveArea = 0.036;
    if (type == kMuTrkIso05)     EffectiveArea = 0.000;
    if (type == kMuEMIso05)      EffectiveArea = 0.124;
    if (type == kMuHadIso05)     EffectiveArea = 0.140;
  } else if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.25 ) {
    if (type == kMuChargedIso03) EffectiveArea = 0.000;
    if (type == kMuNeutralIso03) EffectiveArea = 0.066;
    if (type == kMuChargedIso04) EffectiveArea = 0.000;
    if (type == kMuNeutralIso04) EffectiveArea = 0.149;
    if (type == kMuHadEnergy)    EffectiveArea = 0.056;
    if (type == kMuHoEnergy)     EffectiveArea = 0.000;
    if (type == kMuEmEnergy)     EffectiveArea = 0.000;
    if (type == kMuHadS9Energy)  EffectiveArea = 0.148;
    if (type == kMuHoS9Energy)   EffectiveArea = 0.000;
    if (type == kMuEmS9Energy)   EffectiveArea = 0.000;
    if (type == kMuTrkIso03)     EffectiveArea = 0.000;
    if (type == kMuEMIso03)      EffectiveArea = 0.025;
    if (type == kMuHadIso03)     EffectiveArea = 0.050;
    if (type == kMuTrkIso05)     EffectiveArea = 0.000;
    if (type == kMuEMIso05)      EffectiveArea = 0.120;
    if (type == kMuHadIso05)     EffectiveArea = 0.186;
  } else if (fabs(Eta) >= 2.25 && fabs(Eta) < 2.4 ) {
    if (type == kMuChargedIso03) EffectiveArea = 0.000;
    if (type == kMuNeutralIso03) EffectiveArea = 0.098;
    if (type == kMuChargedIso04) EffectiveArea = 0.000;
    if (type == kMuNeutralIso04) EffectiveArea = 0.200;
    if (type == kMuHadEnergy)    EffectiveArea = 0.093;
    if (type == kMuHoEnergy)     EffectiveArea = 0.000;
    if (type == kMuEmEnergy)     EffectiveArea = 0.000;
    if (type == kMuHadS9Energy)  EffectiveArea = 0.260;
    if (type == kMuHoS9Energy)   EffectiveArea = 0.000;
    if (type == kMuEmS9Energy)   EffectiveArea = 0.000;
    if (type == kMuTrkIso03)     EffectiveArea = 0.000;
    if (type == kMuEMIso03)      EffectiveArea = 0.027;
    if (type == kMuHadIso03)     EffectiveArea = 0.060;
    if (type == kMuTrkIso05)     EffectiveArea = 0.000;
    if (type == kMuEMIso05)      EffectiveArea = 0.139;
    if (type == kMuHadIso05)     EffectiveArea = 0.228;
  }

  if (EffectiveAreaTarget == kMuEANoCorr) {
    return 0.0;
  }
  
  //2012 Data Effective Areas
  else if (EffectiveAreaTarget == kMuEAData2012) {
    if (type == kMuGammaIsoDR0p0To0p1) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.005;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.002;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.005;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.023;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.000;
      if (fabs(Eta) >= 2.3 ) EffectiveArea = 0.000;
    }
    if (type == kMuGammaIsoDR0p1To0p2) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.013;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.007;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.006;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.010;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.010;
      if (fabs(Eta) >= 2.3 ) EffectiveArea = 0.016;
    }
    if (type == kMuGammaIsoDR0p2To0p3) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.027;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.021;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.012;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.018;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.019;
      if (fabs(Eta) >= 2.3  ) EffectiveArea = 0.025;
    }
    if (type == kMuGammaIsoDR0p3To0p4) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.044;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.031;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.019;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.024;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.021;
      if (fabs(Eta) >= 2.3  ) EffectiveArea = 0.040;
    }
    if (type == kMuGammaIsoDR0p4To0p5) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.060;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.045;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.027;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.037;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.043;
      if (fabs(Eta) >= 2.3  ) EffectiveArea = 0.052;
    }
    if (type == kMuNeutralHadronIsoDR0p0To0p1) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.000;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.003;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.006;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.000;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.000;
      if (fabs(Eta) >= 2.3  ) EffectiveArea = 0.003;
    }
    if (type == kMuNeutralHadronIsoDR0p1To0p2) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.005;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.005;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.003;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.004;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.002;
      if (fabs(Eta) >= 2.3  ) EffectiveArea = 0.007;
    }
    if (type == kMuNeutralHadronIsoDR0p2To0p3) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.006;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.009;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.009;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.008;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.013;
      if (fabs(Eta) >= 2.3  ) EffectiveArea = 0.015;
    }
    if (type == kMuNeutralHadronIsoDR0p3To0p4) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.008;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.013;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.013;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.013;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.013;
      if (fabs(Eta) >= 2.3  ) EffectiveArea = 0.013;
    }
    if (type == kMuNeutralHadronIsoDR0p4To0p5) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.014;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.016;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.015;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.024;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.023;
      if (fabs(Eta) >= 2.3  ) EffectiveArea = 0.077;
    }

        
        
    if (type == kMuGammaIso04){
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 )   EffectiveArea = 0.50419;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.30582;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.19765;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 )   EffectiveArea = 0.28723;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 )   EffectiveArea = 0.52529;
      if (fabs(Eta) >= 2.3 )                  EffectiveArea = 0.48818;
    }
    if (type == kMuNeutralHadronIso04){
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 )   EffectiveArea = 0.16580;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.25904;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.24695;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 )   EffectiveArea = 0.22021;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 )   EffectiveArea = 0.34045;
      if (fabs(Eta) >= 2.3 )                  EffectiveArea = 0.21592;
    }
    if (type == kMuGammaAndNeutralHadronIso04){
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 )   EffectiveArea = 0.674;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.565;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.442;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 )   EffectiveArea = 0.515;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 )   EffectiveArea = 0.821;
      if (fabs(Eta) >= 2.3 )                  EffectiveArea = 0.660;
    }
    if (type == kMuGammaAndNeutralHadronIso03){
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 )   EffectiveArea = 0.382;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.317;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.242;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 )   EffectiveArea = 0.326;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 )   EffectiveArea = 0.462;
      if (fabs(Eta) >= 2.3 )                  EffectiveArea = 0.372;
    }
    if (type == kMuGammaAndNeutralHadronIso04Tight){
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 )   EffectiveArea = 0.340;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.310;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.315;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 )   EffectiveArea = 0.415;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 )   EffectiveArea = 0.658;
      if (fabs(Eta) >= 2.3 )                  EffectiveArea = 0.405;
    }
    if (type == kMuGammaAndNeutralHadronIso03Tight){
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 )   EffectiveArea = 0.207;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.183;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.177;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 )   EffectiveArea = 0.271;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 )   EffectiveArea = 0.348;
      if (fabs(Eta) >= 2.3 )                  EffectiveArea = 0.246;
    }
  }

  //2011 Data Effective Areas
  else if (EffectiveAreaTarget == kMuEAData2011) {
    
    if (type == kMuGammaIsoDR0p0To0p1) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.004;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.002;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.002;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.000;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.000;
      if (fabs(Eta) >= 2.3 ) EffectiveArea = 0.005;
    }
    if (type == kMuGammaIsoDR0p1To0p2) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.011;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.008;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.005;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.008;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.008;
      if (fabs(Eta) >= 2.3 ) EffectiveArea = 0.011;
    }
    if (type == kMuGammaIsoDR0p2To0p3) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.023;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.016;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.010;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.014;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.017;
      if (fabs(Eta) >= 2.3  ) EffectiveArea = 0.021;
    }
    if (type == kMuGammaIsoDR0p3To0p4) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.036;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.026;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.017;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.023;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.028;
      if (fabs(Eta) >= 2.3  ) EffectiveArea = 0.032;
    }
    if (type == kMuGammaIsoDR0p4To0p5) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.051;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.037;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.028;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.033;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.042;
      if (fabs(Eta) >= 2.3  ) EffectiveArea = 0.052;
    }
    if (type == kMuNeutralHadronIsoDR0p0To0p1) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.002;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.001;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.001;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.001;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.005;
      if (fabs(Eta) >= 2.3  ) EffectiveArea = 0.007;
    }
    if (type == kMuNeutralHadronIsoDR0p1To0p2) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.005;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.008;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.009;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.009;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.010;
      if (fabs(Eta) >= 2.3  ) EffectiveArea = 0.014;
    }
    if (type == kMuNeutralHadronIsoDR0p2To0p3) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.010;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.015;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.017;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.017;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.019;
      if (fabs(Eta) >= 2.3  ) EffectiveArea = 0.024;
    }
    if (type == kMuNeutralHadronIsoDR0p3To0p4) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.015;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.021;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.024;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.032;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.038;
      if (fabs(Eta) >= 2.3  ) EffectiveArea = 0.038;
    }
    if (type == kMuNeutralHadronIsoDR0p4To0p5) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.020;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.026;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.033;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.045;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.051;
      if (fabs(Eta) >= 2.3  ) EffectiveArea = 0.114;
    }
    /// BEGIN FROM SLIDE 11 OF  https://indico.cern.ch/getFile.py/access?contribId=1&resId=0&materialId=slides&confId=188494
    /// NOTE: to be used with the rho from ALL pf candidates within |eta|<2.5
    if (type == kMuGammaIso03){
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.049;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.5 ) EffectiveArea = 0.030;
      if (fabs(Eta) >= 1.5 && fabs(Eta) < 2.0 ) EffectiveArea = 0.022;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.034;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.041;
      if (fabs(Eta) >= 2.3 )                EffectiveArea = 0.048;
    }
    if (type == kMuGammaIso04){
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.085;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.5 ) EffectiveArea = 0.052;
      if (fabs(Eta) >= 1.5 && fabs(Eta) < 2.0 ) EffectiveArea = 0.038;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.055;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.070;
      if (fabs(Eta) >= 2.3 )                EffectiveArea = 0.081;
    }
    if (type == kMuNeutralHadronIso03){
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.027;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.5 ) EffectiveArea = 0.039;
      if (fabs(Eta) >= 1.5 && fabs(Eta) < 2.0 ) EffectiveArea = 0.044;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.047;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.055;
      if (fabs(Eta) >= 2.3 )                  EffectiveArea = 0.065;
    }
    if (type == kMuNeutralHadronIso04){
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.046;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.5 ) EffectiveArea = 0.067;
      if (fabs(Eta) >= 1.5 && fabs(Eta) < 2.0 ) EffectiveArea = 0.074;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.083;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.095;
      if (fabs(Eta) >= 2.3 )                  EffectiveArea = 0.105;
    }
    if (type == kMuGammaAndNeutralHadronIso03){
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.076;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.5 ) EffectiveArea = 0.070;
      if (fabs(Eta) >= 1.5 && fabs(Eta) < 2.0 ) EffectiveArea = 0.067;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.082;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.097;
      if (fabs(Eta) >= 2.3 )                  EffectiveArea = 0.115;
    }
    if (type == kMuGammaAndNeutralHadronIso04){
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.132;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.5 ) EffectiveArea = 0.120;
      if (fabs(Eta) >= 1.5 && fabs(Eta) < 2.0 ) EffectiveArea = 0.114;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.139;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.168;
      if (fabs(Eta) >= 2.3 )                  EffectiveArea = 0.189;
    }
    /// END FROM SLIDE 11 OF  https://indico.cern.ch/getFile.py/access?contribId=1&resId=0&materialId=slides&confId=188494

    if (type == kMuGammaIso05){
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.05317;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.03502;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.03689;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.05221;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.06668;
      if (fabs(Eta) >= 2.3  ) EffectiveArea = 0.0744;
    }
    if (type == kMuNeutralIso05) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.06408;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.07557;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.08864;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.11492;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.13784;
      if (fabs(Eta) >= 2.3  ) EffectiveArea = 0.18745;
    }
  } 
  
  //Summer11 MC Effective Areas
  else if (EffectiveAreaTarget == kMuEASummer11MC) {
    if (type == kMuGammaIsoDR0p0To0p1) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.000;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.000;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.000;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.000;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.000;
      if (fabs(Eta) >= 2.3  ) EffectiveArea = 0.006;
    }
    if (type == kMuGammaIsoDR0p1To0p2) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.012;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.007;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.006;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.008;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.019;
      if (fabs(Eta) >= 2.3  ) EffectiveArea = 0.015;
    }
    if (type == kMuGammaIsoDR0p2To0p3) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.023;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.018;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.013;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.016;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.024;
      if (fabs(Eta) >= 2.3  ) EffectiveArea = 0.036;
    }
    if (type == kMuGammaIsoDR0p3To0p4) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.038;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.027;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.019;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.033;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.041;
      if (fabs(Eta) >= 2.3  ) EffectiveArea = 0.062;
    }
    if (type == kMuGammaIsoDR0p4To0p5) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.055;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.038;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.032;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.052;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.066;
      if (fabs(Eta) >= 2.3  ) EffectiveArea = 0.093;
    }
    if (type == kMuNeutralHadronIsoDR0p0To0p1) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.002;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.005;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.000;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.000;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.000;
      if (fabs(Eta) >= 2.3  ) EffectiveArea = 0.003;
    }
    if (type == kMuNeutralHadronIsoDR0p1To0p2) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.005;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.006;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.009;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.008;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.009;
      if (fabs(Eta) >= 2.3  ) EffectiveArea = 0.013;
    }
    if (type == kMuNeutralHadronIsoDR0p2To0p3) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.009;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.013;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.015;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.016;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.020;
      if (fabs(Eta) >= 2.3  ) EffectiveArea = 0.024;
    }
    if (type == kMuNeutralHadronIsoDR0p3To0p4) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.012;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.019;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.021;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.025;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.030;
      if (fabs(Eta) >= 2.3  ) EffectiveArea = 0.044;
    }
    if (type == kMuNeutralHadronIsoDR0p4To0p5) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.016;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.026;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.030;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.038;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.048;
      if (fabs(Eta) >= 2.3  ) EffectiveArea = 0.118;
    }
  } 
  
  //Fall11 MC Effective Areas
  else if (EffectiveAreaTarget == kMuEAFall11MC) {
    if (type == kMuGammaIsoDR0p0To0p1) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.004;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.002;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.003;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.009;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.003;
      if (fabs(Eta) >= 2.3  ) EffectiveArea = 0.011;
    }
    if (type == kMuGammaIsoDR0p1To0p2) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.012;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.008;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.006;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.012;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.019;
      if (fabs(Eta) >= 2.3  ) EffectiveArea = 0.024;
    }
    if (type == kMuGammaIsoDR0p2To0p3) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.026;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.020;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.012;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.022;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.027;
      if (fabs(Eta) >= 2.3  ) EffectiveArea = 0.034;
    }
    if (type == kMuGammaIsoDR0p3To0p4) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.042;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.033;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.022;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.036;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.059;
      if (fabs(Eta) >= 2.3  ) EffectiveArea = 0.068;
    }
    if (type == kMuGammaIsoDR0p4To0p5) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.060;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.043;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.036;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.055;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.092;
      if (fabs(Eta) >= 2.3  ) EffectiveArea = 0.115;
    }
    if (type == kMuNeutralHadronIsoDR0p0To0p1) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.002;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.004;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.004;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.004;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.010;
      if (fabs(Eta) >= 2.3  ) EffectiveArea = 0.014;
    }
    if (type == kMuNeutralHadronIsoDR0p1To0p2) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 )   EffectiveArea = 0.005;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.007;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.009;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 )   EffectiveArea = 0.009;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 )   EffectiveArea = 0.015;
      if (fabs(Eta) >= 2.3  )                 EffectiveArea = 0.017;
    }
    if (type == kMuNeutralHadronIsoDR0p2To0p3) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.009;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.015;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.016;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.018;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.022;
      if (fabs(Eta) >= 2.3  ) EffectiveArea = 0.026;
    }
    if (type == kMuNeutralHadronIsoDR0p3To0p4) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.013;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.021;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.026;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.032;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.037;
      if (fabs(Eta) >= 2.3  ) EffectiveArea = 0.042;
    }
    if (type == kMuNeutralHadronIsoDR0p4To0p5) {
      if (fabs(Eta) >= 0.0 && fabs(Eta) < 1.0 ) EffectiveArea = 0.017;
      if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) EffectiveArea = 0.026;
      if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) EffectiveArea = 0.035;
      if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.2 ) EffectiveArea = 0.046;
      if (fabs(Eta) >= 2.2 && fabs(Eta) < 2.3 ) EffectiveArea = 0.063;
      if (fabs(Eta) >= 2.3  ) EffectiveArea = 0.135;
    }
  }

  return EffectiveArea;  

}

