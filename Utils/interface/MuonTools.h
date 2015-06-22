//--------------------------------------------------------------------------------------------------
// $Id: MuonTools.h,v 1.22 2012/04/15 12:07:00 sixie Exp $
//
// MuonTools
//
// This class allows you to classify a given muon according to defined criteria. For this purpose
// is loads histograms from two ROOT files, specified in the constructor. The main function then
// is "IsGood(*muon, selection)" which returns true if the given muon fulfills the selection
// criteria. 
// 
// Logically, the code has been put together by Phil who took most of the ideas from
//
//  http://cmslxr.fnal.gov/lxr/source/RecoMuon/MuonIdentification/
//  http://cmslxr.fnal.gov/lxr/source/DataFormats/MuonReco/
//
// Authors: P.Harris, C.Loizides
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_UTILS_MUONTOOLS_H
#define MITPHYSICS_UTILS_MUONTOOLS_H

#include "MitAna/DataTree/interface/Muon.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/DataTree/interface/BeamSpotCol.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "TH2D.h"

namespace mithep {
  class MuonTools {
    public:
      MuonTools(const char *mutemp="$MIT_DATA/MuonCaloTemplate.root", 
                const char *pitemp="$MIT_DATA/PionCaloTemplate.root");
      virtual ~MuonTools();

      enum EMuIdType {
        kIdUndef = 0,       //not defined
        kWMuId,             //"WMuId"
        kZMuId,             //"ZMuId"
        kTight,             //"Tight"
        kmuonPOG2012CutBasedIDTight,             //"muonPOG2012CutBasedIDTight"
        kLoose,             //"Loose"
        kWWMuIdV1,          //"WWMuIdV1"
        kWWMuIdV2,          //"WWMuIdV2"
        kWWMuIdV3,          //"WWMuIdV3"
        kWWMuIdV4,          //"WWMuIdV4"
        kNoId,              //"NoId"        
        kCustomId,          //"Custom"
        kMVAID_BDTG_IDIso   //"BDTG ID + Iso03, Iso04 Combined"
      };
      enum EMuIsoType {
        kIsoUndef = 0,                      //"not defined"
        kTrackCalo,                         //"TrackCalo"
        kTrackCaloCombined,                 //"TrackCaloCombined"
        kTrackCaloSliding,                  //"TrackCaloSliding"
        kTrackCaloSlidingNoCorrection,      //"TrackCaloSlidingNoCorrection"
        kCombinedRelativeConeAreaCorrected, //"CombinedRelativeConeAreaCorrected"        
        kCombinedRelativeEffectiveAreaCorrected,
        kCustomIso,                         //"Custom"
        kPFIso,                             //"PFIso"
        kPFRadialIso,                       //"PFRadialIso"
        kPFIsoBetaPUCorrected,              //"PFISo with PUcorrection using delta Beta
        kPFIsoEffectiveAreaCorrected,       //"PFIso with EffectiveArea Pileup Correction"
        kPFIsoNoL,                          //"PFIsoNoL"
        kNoIso,                             //"NoIso"
        kMVAIso_BDTG_IDIso,                 //"BDTG ID + Iso03, Iso04 Combined"
        kIsoRingsV0_BDTG_Iso,               //"BDTG Iso Rings"
        kIsoDeltaR                          //"BGDT Iso dR"              
      };
      enum EMuClassType {
        kClassUndef = 0,    //not defined
        kAll,               //"All"
        kGlobal,            //"Global"
        kGlobalorTracker,   //"Global or Tracker Muon"
        kGlobalTracker,     //"GlobalTracker"
        kSta,               //"Standalone"
        kTrackerMuon,       //"TrackerMuon"
        kCaloMuon,          //"CaloMuon"
        kTrackerBased,      //"TrackerMuon or CaloMuon"
        kGlobalOnly         //"GlobalOnly"
      };

      enum ESelType { 
        kAllArbitrated,          //All arbitration (DT/CSC/RPC Hits) put on at least one 
                                 //  segments given a global Muon
        kPromptTight,            //Standard global muon identification
        kTMLastStationLoose,     //Loose matching requirements on lastmost muon station of reco 
        kTMLastStationTight,     //Tight matching requirements on lastmost muon station of reco 
        kTMOneStationLoose,      //Loose matching requirements on at least one muon station of reco 
        kTMOneStationTight,      //Tight matching requirements on at least one muon station of reco 
        kTM2DCompatibilityLoose, //Loose requirement on sum of compatabiliity variables 
                                 //  ===> 1.2 Segment compatability + 0.8 calo compatability > 0.8
        kTM2DCompatibilityTight  //Tight requirement on sum of compatabiliity variables 
                                 //  ===> 1.2 Segment compatability + 0.8 calo compatability > 1.2
      };

      enum EMuonEffectiveAreaType {
	kMuTrkIso03, 
	kMuEcalIso03, 
	kMuHcalIso03, 
	kMuTrkIso05, 
	kMuEcalIso05, 
	kMuHcalIso05, 
	kMuChargedIso03, 
	kMuGammaIso03, 
	kMuNeutralHadronIso03, 
	kMuGammaAndNeutralHadronIso03,
	kMuGammaIso03Tight, 
	kMuNeutralHadronIso03Tight, 
	kMuGammaAndNeutralHadronIso03Tight,
	kMuChargedIso04, 
	kMuGammaIso04, 
	kMuNeutralHadronIso04, 
	kMuGammaAndNeutralHadronIso04,
	kMuGammaIso04Tight, 
	kMuNeutralHadronIso04Tight, 
	kMuGammaAndNeutralHadronIso04Tight,
	kMuGammaIsoDR0p0To0p1,
	kMuGammaIsoDR0p1To0p2,
	kMuGammaIsoDR0p2To0p3,
	kMuGammaIsoDR0p3To0p4,
	kMuGammaIsoDR0p4To0p5,
	kMuNeutralHadronIsoDR0p0To0p1,
	kMuNeutralHadronIsoDR0p1To0p2,
	kMuNeutralHadronIsoDR0p2To0p3,
	kMuNeutralHadronIsoDR0p3To0p4,
	kMuNeutralHadronIsoDR0p4To0p5,
	kMuGammaIso05,
	kMuNeutralIso05,
        kMuNeutralIso03, 
        kMuNeutralIso04, 
        kMuHadEnergy, 
        kMuHoEnergy, 
        kMuEmEnergy, 
        kMuHadS9Energy, 
        kMuHoS9Energy, 
        kMuEmS9Energy,
        kMuEMIso03, 
        kMuHadIso03, 
        kMuEMIso05, 
        kMuHadIso05
      };

      enum EMuonEffectiveAreaTarget {
	kMuEANoCorr,
	kMuEAData2011,
	kMuEASummer11MC,
	kMuEAFall11MC,
	kMuEAData2012
      };

      Bool_t          Init(const char *mutemp, const char *pitemp);
      Bool_t          IsGood(const mithep::Muon *iMuon, ESelType iSel) const;
      Double_t        GetCaloCompatability(const mithep::Muon *iMuon,
                                         Bool_t iEMSpecial, Bool_t iCorrectedHCAL) const; 
      Double_t        GetSegmentCompatability(const mithep::Muon *iMuon)             const;
      static Bool_t   PassD0Cut(const Muon *mu, const VertexCol *vertices, Double_t fD0Cut, Int_t nVertex = 0); 
      static Bool_t   PassD0Cut(const Muon *mu, const BeamSpotCol *beamspots, Double_t fD0Cut);
      static Bool_t   PassDZCut(const Muon *mu, const VertexCol *vertices, Double_t fDZCut, Int_t nVertex = 0);
      static Bool_t   PassSoftMuonCut(const Muon *mu, const VertexCol *vertices, const Double_t fDZCut = 0.2,
                                    const Bool_t applyIso = kTRUE);
      static Double_t MuonEffectiveArea(EMuonEffectiveAreaType type, Double_t Eta, 
                                        EMuonEffectiveAreaTarget EffectiveAreaTarget = kMuEAData2011);

    protected:
      void        DeleteHistos();
      Bool_t      Overflow(const TH2D *iHist, Double_t lVal0, Double_t lVal1)    const; 
      Double_t    SigWeight(Double_t iVal0, Double_t iVal1)                      const; 
   
    private:
      Bool_t      fIsInit;              //!=true if histograms are loaded
      TH2D       *fmuon_em_etaEmi;      //!Neg Endcap EM       Calo Deposit Template for Muons
      TH2D       *fmuon_had_etaEmi;     //!Neg Endcap Hadronic Calo Deposit Template for Muons
      TH2D       *fmuon_had_etaTmi;     //!Neg Transition Hadronic Calo Deposit Template for Muons
      TH2D       *fmuon_em_etaB;        //!Barrel EM       Calo Deposit Template for Muons
      TH2D       *fmuon_had_etaB;       //!Barrel Hadronic Calo Deposit Template for Muons
      TH2D       *fmuon_ho_etaB;        //!Barrel HO       Calo Deposit Template for Muons
      TH2D       *fmuon_had_etaTpl;     //!Plus Transition Hadronic Calo Deposit Template for Muons
      TH2D       *fmuon_em_etaEpl;      //!Plus Endcap EM       Calo Deposit Template for Muons
      TH2D       *fmuon_had_etaEpl;     //!Plus Endcap Hadronic Calo Deposit Template for Muons
      TH2D       *fpion_em_etaEmi;      //!Neg  Endcap EM       Calo Deposit Template for Pions
      TH2D       *fpion_had_etaEmi;     //!Neg  Endcap Hadronic Calo Deposit Template for Pions
      TH2D       *fpion_had_etaTmi;     //!Neg Transition Hadronic Calo Deposit Template for Pions
      TH2D       *fpion_em_etaB;        //!Barrel EM       Calo Deposit Template for Pions
      TH2D       *fpion_had_etaB;       //!Barrel Hadronic Calo Deposit Template for Pions
      TH2D       *fpion_ho_etaB;        //!Barrel HO       Calo Deposit Template for Pions
      TH2D       *fpion_had_etaTpl;     //!Plus Transition Hadronic Calo Deposit Template for Pions
      TH2D       *fpion_em_etaEpl;      //!Plus Endcap EM       Calo Deposit Template for Pions
      TH2D       *fpion_had_etaEpl;     //!Plus Endcap Hadronic Calo Deposit Template for Pions

      TH2D       *LoadHisto(const char *fname, TFile *file)                      const;

    ClassDef(MuonTools, 1) // Muon tools
  };
}

//--------------------------------------------------------------------------------------------------
inline Double_t mithep::MuonTools::SigWeight(Double_t iVal0, Double_t iVal1) const
{
  // Returns weighted uncertainty given segment matching uncertainty (iVal0) and
  // segment matching pull (iVal1).

  if (iVal1 < 1.)  //if pull defined and within range
    return 1.;
  if (iVal0 < 3. && iVal1 > 3.) {       //if pull not well defined and uncertainty defined
    Double_t lVal = TMath::Max(iVal0,1.);
    return 1./TMath::Power(lVal,0.25);
  }

  Double_t lVal = TMath::Max(iVal1,1.); //if pull > 1 and pull < 3 return 1/pull^4
  return 1./TMath::Power(lVal,0.25);
}
//--------------------------------------------------------------------------------------------------
inline Bool_t mithep::MuonTools::Overflow(const TH2D *iHist, Double_t lVal0, Double_t lVal1) const
{
  // Check if values are in overflow bins of given histogram.

  if(iHist == 0)
    return kTRUE;

  if (iHist ->GetXaxis()->FindBin(lVal0) == 0                  || 
      iHist ->GetXaxis()->FindBin(lVal0) >  iHist->GetNbinsX() || 
      iHist ->GetYaxis()->FindBin(lVal1) == 0                  || 
      iHist ->GetYaxis()->FindBin(lVal1) >  iHist->GetNbinsY()) {
    return kTRUE;
  }
  return kFALSE;
}
#endif
