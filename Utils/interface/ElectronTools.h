//--------------------------------------------------------------------------------------------------
// ElectronTools
//
// Helper Class for electron Identification decisions.
//
// Authors: S.Xie
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_UTILS_ELECTRONTOOLS_H
#define MITPHYSICS_UTILS_ELECTRONTOOLS_H

#include "MitAna/DataTree/interface/Electron.h"
#include "MitAna/DataTree/interface/DecayParticleCol.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/DataTree/interface/BeamSpotCol.h"
#include "MitAna/DataTree/interface/TriggerObjectCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitPhysics/ElectronLikelihood/interface/ElectronLikelihood.h"
#include <utility>

namespace mithep {

  class ElectronTools {
  public:
    ElectronTools();
    virtual ~ElectronTools() {}
  
    enum EElIdType {
      kIdUndef = 0,       //not defined
      kTight,             //"Tight"
      kLoose,             //"Loose"
      kLikelihood,        //"Likelihood"
      kNoId,              //"NoId"
      kZeeId,             //"ZeeId"
      kCustomIdLoose,     //"CustomLoose"
      kCustomIdTight,     //"CustomTight"
      kVBTFWorkingPointFakeableId,
      kVBTFWorkingPoint95Id,
      kVBTFWorkingPoint90Id,
      kVBTFWorkingPoint85Id,
      kVBTFWorkingPoint80Id,
      kVBTFWorkingPointLowPtId,
      kVBTFWorkingPoint70Id,
      kMVAID_BDTG_NoIPInfo,
      kMVAID_BDTG_WithIPInfo,
      kMVAID_BDTG_IDIsoCombined,
      kHggLeptonTagId,
      kHggLeptonTagId2012,
      kHggLeptonTagId2012HCP,
      kMVAID_BDTG_IDHWW2012TrigV0,
      kMVAID_BDTG_IDIsoCombinedHWW2012TrigV4,
      kPhys14Veto,
      kPhys14Loose,
      kPhys14Medium,
      kPhys14Tight
    };

    enum EElIsoType {
      kIsoUndef = 0,      	       //"not defined"
      kTrackCalo,         	       //"TrackCalo"
      kTrackJura,         	       //"TrackJura"
      kTrackJuraCombined, 	       //"TrackJuraCombined"
      kTrackJuraSliding,  	       //"TrackJuraSliding"
      kTrackJuraSlidingNoCorrection, //"TrackJuraSlidingNoCorrection"
      kCombinedRelativeConeAreaCorrected, //"CombinedRelativeConeAreaCorrected"
      kNoIso,             	       //"NoIso"
      kPFIso,             	       //"PFIso"
      kPFIsoRhoCorrected,
      kPFIsoNoL,          	       //"PFIsoNoL"
      kZeeIso,            	       //"ZeeIso"
      kCustomIso,         	       //"Custom"
      kVBTFWorkingPoint95IndividualIso,
      kVBTFWorkingPoint90IndividualIso,
      kVBTFWorkingPoint85IndividualIso,
      kVBTFWorkingPoint80IndividualIso,
      kVBTFWorkingPoint70IndividualIso,
      kVBTFWorkingPoint95CombinedIso,
      kVBTFWorkingPoint90CombinedIso,
      kVBTFWorkingPoint85CombinedIso,
      kVBTFWorkingPoint80CombinedIso,
      kVBTFWorkingPoint70CombinedIso,
      kMVAIso_BDTG_IDIsoCombined,
      kPFIso_HWW2012TrigV0,
      kPFIso_HggLeptonTag2012,
      kPFIso_HggLeptonTag2012HCP,
      kMVAIso_BDTG_IDIsoCombinedHWW2012TrigV4,
      kPhys14VetoIso,
      kPhys14LooseIso,
      kPhys14MediumIso,
      kPhys14TightIso
    };

    enum EElectronIsoVar {
      kTrackIsolation,
      kCaloIsolation,
      kEcalJuraIsolation,
      kHcalIsolation,
      kCombinedIsolationEB,
      kCombinedIsolationEE,
      nEElectronIsoVars
    };

    enum EElectronEffectiveAreaType {
      kEleChargedIso03, 
      kEleNeutralHadronIso03, 
      kEleGammaAndNeutralHadronIso03,
      kEleGammaIso03, 
      kEleGammaIsoVetoEtaStrip03, 
      kEleNeutralIso03,
      kEleChargedIso04, 
      kEleNeutralHadronIso04, 
      kEleGammaAndNeutralHadronIso04,
      kEleGammaIso04, 
      kEleGammaIsoVetoEtaStrip04, 
      kEleNeutralHadronIso007, 
      kEleNeutralIso04,
      kEleHoverE, 
      kEleHcalDepth1OverEcal, 
      kEleHcalDepth2OverEcal,
      kEleGammaIsoDR0p0To0p1,
      kEleGammaIsoDR0p1To0p2,
      kEleGammaIsoDR0p2To0p3,
      kEleGammaIsoDR0p3To0p4,
      kEleGammaIsoDR0p4To0p5,
      kEleNeutralHadronIsoDR0p0To0p1,
      kEleNeutralHadronIsoDR0p1To0p2,
      kEleNeutralHadronIsoDR0p2To0p3,
      kEleNeutralHadronIsoDR0p3To0p4,
      kEleNeutralHadronIsoDR0p4To0p5
    };
      
    enum EElectronEffectiveAreaTarget {
      kEleEANoCorr,
      kEleEAData2011,
      kEleEAData2012,
      kEleEASummer11MC,
      kEleEAFall11MC,
      kEleEAPhys14
    };

    static Bool_t       PassChargeFilter(const Electron *el);
    static Bool_t       PassConversionFilter(const Electron *el, const DecayParticleCol *conversions, 
                                             const BaseVertex *vtx, UInt_t nWrongHitsMax=0, Double_t probMin=1e-6,
                                             Double_t lxyMin = 2.0, Bool_t matchCkf = kTRUE, Bool_t requireArbitratedMerged = kFALSE, Double_t trkptMin = -99.);
    static Bool_t       PassConversionFilterPFAOD(const Electron *el, const DecayParticleCol *conversions, 
                                                  const BaseVertex *vtx, UInt_t nWrongHitsMax=0, Double_t probMin=1e-6,
                                                  Double_t lxyMin = 2.0, Bool_t matchCkf = kTRUE, Bool_t requireArbitratedMerged = kFALSE, Double_t trkptMin = -99.);
    static Bool_t       PassNExpectedHits(Electron const*, EElIdType, Bool_t invert = false);
    static Bool_t       PassCustomID(const Electron *el, EElIdType idType);
    static Bool_t       PassCustomIso(const Electron *el, EElIsoType isoType);
    static Bool_t       PassID(Electron const*, EElIdType); // new implementation to get rid of PassCustomID
    static Bool_t       PassIso(Electron const*, EElIsoType); // new implementation to get rid of PassCustomIso
    static Bool_t       PassPFIso(Electron const*, EElIsoType, PFCandidateCol const*, Vertex const*, MuonCol const* = 0, ElectronCol const* = 0);
    static Bool_t       PassIsoRhoCorr(Electron const*, EElIsoType, Double_t rho, PFCandidateCol const* = 0, Vertex const* = 0); // new implementation to get rid of PassCustomIso
    static Bool_t       PassD0Cut(Electron const*, VertexCol const*, EElIdType, Int_t iVertex = 0);
    static Bool_t       PassD0Cut(Electron const*, BeamSpotCol const*, EElIdType);
    static Bool_t       PassDZCut(Electron const*, VertexCol const*, EElIdType, Int_t iVertex = 0);
    static Bool_t       PassSpikeRemovalFilter(Electron const*);
    static Bool_t       PassTriggerMatching(Electron const*, TriggerObjectCol const*);
    static Int_t        Classify(const Electron *ele);
    static Int_t        PassTightId(const Electron *ele, const VertexCol *vertices, 
                                    const DecayParticleCol *conversions, const Int_t typeCuts,
                                    Double_t beta = 1.0);
    static bool         compute_cut(double x, double et, double cut_min, double cut_max, bool gtn=false);
    static Double_t     Likelihood(ElectronLikelihood *LH, const Electron *ele);
    static Double_t     ElectronEffectiveArea(EElectronEffectiveAreaType type, Double_t Eta, 
                                              EElectronEffectiveAreaTarget EffectiveAreaTarget = kEleEAData2011);
    static std::pair<Double_t,Double_t> ComputeEPCombination( const Electron * ele, const float regression_energy, 
                                                              const float regression_energy_error);

    static Bool_t       PassHggLeptonTagID(const Electron *el);
    static Bool_t       PassHggLeptonTagID2012(const Electron *el);

  private:
    static Bool_t       PassD0Cut(const Electron *el, Double_t d0, EElIdType);
    static Bool_t       PassDZCut(const Electron *el, Double_t dz, EElIdType);
    
    ClassDef(ElectronTools, 1) // Muon tools
  };
}

#endif
