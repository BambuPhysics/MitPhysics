//--------------------------------------------------------------------------------------------------
// MuonIDMVA
//
// Helper Class for Muon Identification MVA
//
// Authors: S.Xie
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_UTILS_MuonIDMVA_H
#define MITPHYSICS_UTILS_MuonIDMVA_H

#include "MitAna/DataTree/interface/MuonFwd.h"
#include "MitAna/DataTree/interface/VertexFwd.h"
#include "MitAna/DataTree/interface/TrackFwd.h"
#include "MitAna/DataTree/interface/Muon.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityCol.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitPhysics/Utils/interface/MuonTools.h"

class TRandom3;
namespace TMVA {
  class Reader;
}

namespace mithep {
  class MuonIDMVA {
  public:
    MuonIDMVA();
    virtual ~MuonIDMVA(); 

    enum MVAType {
      kUninitialized,
      kV2,
      kV3,
      kV8,
      kIDIsoCombinedDetIso,
      kIsoRingsV0,
      kIDIsoCombinedIsoRingsV0,
      kIDV0,
      kIsoDeltaR
    };


    void     Initialize(std::string const& methodName,
                        std::string const& weightsFile,
                        MuonIDMVA::MVAType,
                        mithep::PileupEnergyDensity::Algo = mithep::PileupEnergyDensity::kHighEta);
    void     Initialize(std::string const& methodName,
                        MuonIDMVA::MVAType,
                        Bool_t useBinnedVersion,
                        std::vector<std::string> const& weightsfiles,
                        mithep::PileupEnergyDensity::Algo = mithep::PileupEnergyDensity::kHighEta);
    void     Initialize(TString const& methodName,
                        TString const& Subdet0Pt10To14p5Weights , 
                        TString const& Subdet1Pt10To14p5Weights , 
                        TString const& Subdet0Pt14p5To20Weights,
                        TString const& Subdet1Pt14p5To20Weights, 
                        TString const& Subdet0Pt20ToInfWeights, 
                        TString const& Subdet1Pt20ToInfWeights,
                        MuonIDMVA::MVAType,
                        mithep::PileupEnergyDensity::Algo = mithep::PileupEnergyDensity::kHighEta);

    Bool_t   IsInitialized() const { return fIsInitialized; }
    UInt_t   GetMVABin(double eta,double pt,
                       Bool_t isGlobal = kTRUE, Bool_t isTrackerMuon = kTRUE) const;
    Double_t MVAValue(const Muon *mu, const Vertex *vertex, MuonTools *fMuonTools,
                      const PFCandidateCol *PFCands, 
                      const PileupEnergyDensityCol *PileupEnergyDensity, 
                      Bool_t printDebug = kFALSE);
    Double_t MVAValue(const Muon *mu, const Vertex *vertex, MuonTools *fMuonTools,
                      const PFCandidateCol *PFCands, 
                      const PileupEnergyDensityCol *PileupEnergyDensity, 
                      MuonTools::EMuonEffectiveAreaTarget EffectiveAreaTarget,
                      const ElectronCol *goodElectrons,
                      const MuonCol *goodMuons,            
                      Bool_t printDebug = kFALSE);
    Double_t MVAValue( Double_t MuPt , Double_t MuEta,
                       Double_t                   MuTkNchi2, 
                       Double_t                   MuGlobalNchi2, 
                       Double_t                   MuNValidHits, 
                       Double_t                   MuNTrackerHits, 
                       Double_t                   MuNPixelHits, 
                       Double_t                   MuNMatches, 
                       Double_t                   MuD0, 
                       Double_t                   MuIP3d, 
                       Double_t                   MuIP3dSig, 
                       Double_t                   MuTrkKink, 
                       Double_t                   MuSegmentCompatibility, 
                       Double_t                   MuCaloCompatibility, 
                       Double_t                   MuHadEnergyOverPt, 
                       Double_t                   MuHoEnergyOverPt, 
                       Double_t                   MuEmEnergyOverPt, 
                       Double_t                   MuHadS9EnergyOverPt, 
                       Double_t                   MuHoS9EnergyOverPt, 
                       Double_t                   MuEmS9EnergyOverPt, 
                       Double_t                   MuChargedIso03OverPt,
                       Double_t                   MuNeutralIso03OverPt,
                       Double_t                   MuChargedIso04OverPt,
                       Double_t                   MuNeutralIso04OverPt,
                       Bool_t                     printDebug = kFALSE
                       );
    Double_t MVAValue( Double_t MuPt , Double_t MuEta,
                       Double_t                   MuTkNchi2, 
                       Double_t                   MuGlobalNchi2, 
                       Double_t                   MuNValidHits, 
                       Double_t                   MuNTrackerHits, 
                       Double_t                   MuNPixelHits, 
                       Double_t                   MuNMatches, 
                       Double_t                   MuD0, 
                       Double_t                   MuIP3d, 
                       Double_t                   MuIP3dSig, 
                       Double_t                   MuTrkKink, 
                       Double_t                   MuSegmentCompatibility, 
                       Double_t                   MuCaloCompatibility, 
                       Double_t                   MuHadEnergyOverPt, 
                       Double_t                   MuHoEnergyOverPt, 
                       Double_t                   MuEmEnergyOverPt, 
                       Double_t                   MuHadS9EnergyOverPt, 
                       Double_t                   MuHoS9EnergyOverPt, 
                       Double_t                   MuEmS9EnergyOverPt, 
                       Double_t                   MuTrkIso03OverPt,
                       Double_t                   MuEMIso03OverPt,
                       Double_t                   MuHadIso03OverPt,
                       Double_t                   MuTrkIso05OverPt,
                       Double_t                   MuEMIso05OverPt,
                       Double_t                   MuHadIso05OverPt,
                       Bool_t                     printDebug = kFALSE
                       );
    Double_t MVAValue_IsoRings( Double_t MuPt,
                                Double_t MuEta,
                                Bool_t MuIsGlobal,
                                Bool_t MuIsTracker,
                                Double_t ChargedIso_DR0p0To0p1,
                                Double_t ChargedIso_DR0p1To0p2,
                                Double_t ChargedIso_DR0p2To0p3,
                                Double_t ChargedIso_DR0p3To0p4,
                                Double_t ChargedIso_DR0p4To0p5,
                                Double_t GammaIso_DR0p0To0p1,
                                Double_t GammaIso_DR0p1To0p2,
                                Double_t GammaIso_DR0p2To0p3,
                                Double_t GammaIso_DR0p3To0p4,
                                Double_t GammaIso_DR0p4To0p5,
                                Double_t NeutralHadronIso_DR0p0To0p1,
                                Double_t NeutralHadronIso_DR0p1To0p2,
                                Double_t NeutralHadronIso_DR0p2To0p3,
                                Double_t NeutralHadronIso_DR0p3To0p4,
                                Double_t NeutralHadronIso_DR0p4To0p5,
                                Bool_t printDebug = kFALSE);
    Double_t MVAValue_ID( Double_t MuPt, 
                          Double_t MuEta,                             
                          Bool_t MuIsGlobal,
                          Bool_t MuIsTracker,
                          Double_t MuTkNchi2, 
                          Double_t MuGlobalNchi2, 
                          Double_t MuNValidHits, 
                          Double_t MuNTrackerHits, 
                          Double_t MuNPixelHits, 
                          Double_t MuNMatches, 
                          Double_t MuTrkKink, 
                          Double_t MuSegmentCompatibility, 
                          Double_t MuCaloCompatibility, 
                          Double_t MuHadEnergy, 
                          Double_t MuEmEnergy, 
                          Double_t MuHadS9Energy, 
                          Double_t MuEmS9Energy, 
                          Bool_t printDebug = kFALSE);

  protected:      
    std::vector<TMVA::Reader*> fTMVAReader;
    TString                   fMethodname;
    Bool_t                    fIsInitialized;
    MVAType                   fMVAType;
    Bool_t                    fUseBinnedVersion;
    UInt_t                    fNMVABins;
    mithep::PileupEnergyDensity::Algo fTheRhoAlgo;

    Float_t                   fMVAVar_MuPt; 
    Float_t                   fMVAVar_MuEta; 
    Float_t                   fMVAVar_MuTkNchi2; 
    Float_t                   fMVAVar_MuGlobalNchi2; 
    Float_t                   fMVAVar_MuNValidHits; 
    Float_t                   fMVAVar_MuNTrackerHits; 
    Float_t                   fMVAVar_MuNPixelHits; 
    Float_t                   fMVAVar_MuNMatches; 
    Float_t                   fMVAVar_MuD0; 
    Float_t                   fMVAVar_MuIP3d; 
    Float_t                   fMVAVar_MuIP3dSig; 
    Float_t                   fMVAVar_MuTrkKink; 
    Float_t                   fMVAVar_MuSegmentCompatibility; 
    Float_t                   fMVAVar_MuCaloCompatibility; 
    Float_t                   fMVAVar_MuHadEnergyOverPt; 
    Float_t                   fMVAVar_MuHoEnergyOverPt; 
    Float_t                   fMVAVar_MuEmEnergyOverPt; 
    Float_t                   fMVAVar_MuHadS9EnergyOverPt; 
    Float_t                   fMVAVar_MuHoS9EnergyOverPt; 
    Float_t                   fMVAVar_MuEmS9EnergyOverPt; 
    Float_t                   fMVAVar_MuHadEnergy; 
    Float_t                   fMVAVar_MuEmEnergy; 
    Float_t                   fMVAVar_MuHadS9Energy; 
    Float_t                   fMVAVar_MuEmS9Energy; 
    Float_t                   fMVAVar_MuChargedIso03OverPt;
    Float_t                   fMVAVar_MuNeutralIso03OverPt;
    Float_t                   fMVAVar_MuChargedIso04OverPt;
    Float_t                   fMVAVar_MuNeutralIso04OverPt;
    Float_t                   fMVAVar_MuTrkIso03OverPt;
    Float_t                   fMVAVar_MuEMIso03OverPt;
    Float_t                   fMVAVar_MuHadIso03OverPt;
    Float_t                   fMVAVar_MuTrkIso05OverPt;
    Float_t                   fMVAVar_MuEMIso05OverPt;
    Float_t                   fMVAVar_MuHadIso05OverPt;

    Float_t                   fMVAVar_ChargedIso_DR0p0To0p1;
    Float_t                   fMVAVar_ChargedIso_DR0p1To0p2;
    Float_t                   fMVAVar_ChargedIso_DR0p2To0p3;
    Float_t                   fMVAVar_ChargedIso_DR0p3To0p4;
    Float_t                   fMVAVar_ChargedIso_DR0p4To0p5;
    Float_t                   fMVAVar_ChargedIso_DR0p5To0p7;
    Float_t                   fMVAVar_GammaIso_DR0p0To0p1;
    Float_t                   fMVAVar_GammaIso_DR0p1To0p2;
    Float_t                   fMVAVar_GammaIso_DR0p2To0p3;
    Float_t                   fMVAVar_GammaIso_DR0p3To0p4;
    Float_t                   fMVAVar_GammaIso_DR0p4To0p5;
    Float_t                   fMVAVar_GammaIso_DR0p5To0p7;
    Float_t                   fMVAVar_NeutralHadronIso_DR0p0To0p1;
    Float_t                   fMVAVar_NeutralHadronIso_DR0p1To0p2;
    Float_t                   fMVAVar_NeutralHadronIso_DR0p2To0p3;
    Float_t                   fMVAVar_NeutralHadronIso_DR0p3To0p4;
    Float_t                   fMVAVar_NeutralHadronIso_DR0p4To0p5;
    Float_t                   fMVAVar_NeutralHadronIso_DR0p5To0p7;

    // isolation variables II
    Float_t                   fMVAVar_MuRelIsoPFCharged;
    Float_t                   fMVAVar_MuRelIsoPFNeutral;
    Float_t                   fMVAVar_MuRelIsoPFPhotons;
    Float_t                   fMVAVar_MuDeltaRMean;
    Float_t                   fMVAVar_MuDeltaRSum;
    Float_t                   fMVAVar_MuDensity;


    ClassDef(MuonIDMVA, 0) // Muon MVA
  };
}

#endif
