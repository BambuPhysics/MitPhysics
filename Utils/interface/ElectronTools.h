//--------------------------------------------------------------------------------------------------
// $Id $
//
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
        kMVAID_BDTG_IDIsoCombinedHWW2012TrigV4
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
        kPFIsoNoL,          	       //"PFIsoNoL"
        kZeeIso,            	       //"ZeeIso"
        kCustomIso,         	       //"Custom"
        kVBTFWorkingPoint95Iso,
        kVBTFWorkingPoint90Iso,
        kVBTFWorkingPoint85Iso,
        kVBTFWorkingPoint80Iso,
        kVBTFWorkingPoint70Iso,
        kMVAIso_BDTG_IDIsoCombined,
        kPFIso_HWW2012TrigV0,
	kPFIso_HggLeptonTag2012,
	kPFIso_HggLeptonTag2012HCP,
        kMVAIso_BDTG_IDIsoCombinedHWW2012TrigV4
      };

      enum EElectronEffectiveAreaType {
        kEleChargedIso03, 
        kEleNeutralHadronIso03, 
	kEleGammaAndNeutralHadronIso03,
        kEleGammaIso03, 
        kEleGammaIsoVetoEtaStrip03, 
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
        kEleEAFall11MC
      };

      static Bool_t       PassChargeFilter(const Electron *el);
      static Bool_t       PassConversionFilter(const Electron *el, const DecayParticleCol *conversions, 
                                               const BaseVertex *vtx, UInt_t nWrongHitsMax=0, Double_t probMin=1e-6,
                                               Double_t lxyMin = 2.0, Bool_t matchCkf = kTRUE, Bool_t requireArbitratedMerged = kFALSE, Double_t trkptMin = -99.);
      static Bool_t       PassConversionFilterPFAOD(const Electron *el, const DecayParticleCol *conversions, 
						    const BaseVertex *vtx, UInt_t nWrongHitsMax=0, Double_t probMin=1e-6,
						    Double_t lxyMin = 2.0, Bool_t matchCkf = kTRUE, Bool_t requireArbitratedMerged = kFALSE, Double_t trkptMin = -99.);
      static Bool_t       PassCustomID(const Electron *el, EElIdType idType);
      static Bool_t       PassCustomIso(const Electron *el, EElIsoType isoType,
                                        Bool_t useCombineIso = kTRUE);
      static Bool_t       PassD0Cut(const Electron *el, const VertexCol *vertices, Double_t fD0Cut, Int_t nVertex = 0);
      static Bool_t       PassD0Cut(const Electron *el, const BeamSpotCol *beamspots, Double_t fD0Cut);
      static Bool_t       PassDZCut(const Electron *el, const VertexCol *vertices, Double_t fDZCut, Int_t nVertex = 0);
      static Bool_t       PassSpikeRemovalFilter(const Electron *ele);
      static Bool_t       PassTriggerMatching(const Electron *ele, const TriggerObjectCol *trigobjs);
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
    
    ClassDef(ElectronTools, 1) // Muon tools
  };
}

#endif
