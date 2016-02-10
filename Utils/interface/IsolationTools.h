//--------------------------------------------------------------------------------------------------
// IsolationTools
//
// Isolation functions to compute various kinds of isolation.
//
// Authors: S.Xie
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_UTILS_ISOLATIONTOOLS_H
#define MITPHYSICS_UTILS_ISOLATIONTOOLS_H

#include <TMath.h>
#include "MitAna/DataTree/interface/TrackFwd.h"
#include "MitAna/DataTree/interface/PhotonFwd.h"
#include "MitAna/DataTree/interface/BasicClusterFwd.h"
#include "MitAna/DataTree/interface/SuperClusterFwd.h"
#include "MitAna/DataTree/interface/CaloTowerFwd.h"
#include "MitAna/DataTree/interface/VertexFwd.h"
#include "MitAna/DataTree/interface/MuonFwd.h"
#include "MitAna/DataTree/interface/ElectronFwd.h"
#include "MitAna/DataTree/interface/PFCandidateFwd.h"
#include "MitAna/DataTree/interface/TrackFwd.h"
#include "MitAna/DataTree/interface/DecayParticleFwd.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityFwd.h"
#include "MitPhysics/Utils/interface/ElectronTools.h"
#include "MitPhysics/Utils/interface/MuonTools.h"
#include "MitPhysics/Utils/interface/PhotonTools.h"

namespace mithep
{
  class IsolationTools {
  public:
    static Double_t TrackIsolation(const mithep::Track *p, Double_t extRadius, 
                                   Double_t intRadius, Double_t ptLow, Double_t maxVtxZDist, 
                                   const mithep::TrackCol *tracks); 
    static Double_t EcalIsolation(const SuperCluster *sc, Double_t coneSize, Double_t etLow, 
                                  const mithep::BasicClusterCol *basicClusters);
    static Double_t CaloTowerHadIsolation(const ThreeVector *p,  Double_t extRadius, 
                                          Double_t intRadius, Double_t etLow, 
                                          const mithep::CaloTowerCol
                                          *caloTowers);
    static Double_t CaloTowerEmIsolation(const ThreeVector *p, Double_t extRadius, 
                                         Double_t intRadius, Double_t etLow, 
                                         const mithep::CaloTowerCol *caloTowers);
    static Double_t PFRadialMuonIsolation(const Muon *p, const PFCandidateCol *PFCands, 
                                          Double_t ptMin = 1.0, Double_t extRadius = 0.3);
    static Double_t PFMuonIsolation(const Muon *p, const PFCandidateCol *PFCands, const Vertex *vertex, 
                                    Double_t  delta_z = 0.1, Double_t ptMin = 1.0,
                                    Double_t extRadius = 0.4, Double_t intRadiusGamma = 0.07, Double_t intRadius = 0.0);
    static Double_t PFMuonIsolation(const Muon *p, const PFCandidateCol *PFCands,
                                    const MuonCol *goodMuons, const ElectronCol *goodElectrons, 
                                    const Vertex *vertex, Double_t  delta_z, Double_t ptMin,
                                    Double_t extRadius, Double_t intRadiusGamma, Double_t intRadius);
    static Double_t PFElectronIsolation(const Electron *p, const PFCandidateCol *PFCands, 
                                        const Vertex *vertex, Double_t  delta_z, Double_t ptMin,
                                        Double_t extRadius, Double_t intRadius, Int_t PFCandidateType = -1);
    static Double_t PFElectronIsolation(const Electron *p, const PFCandidateCol *PFCands, 
                                        const MuonCol *goodMuons, const ElectronCol *goodElectrons,
                                        const Vertex *vertex, Double_t  delta_z, Double_t ptMin,
                                        Double_t extRadius, Double_t intRadius, Int_t PFCandidateType = -1);
    static Double_t PFElectronIsolation2012(const Electron *ele, const Vertex *vertex, 
                                            const PFCandidateCol *PFCands, 
                                            Double_t rho,
                                            ElectronTools::EElectronEffectiveAreaTarget,
                                            const ElectronCol *goodElectrons = 0,
                                            const MuonCol *goodMuons = 0,
                                            Double_t dRMax = 0.4, Bool_t isDebug = kFALSE);
    static Double_t PFElectronIsolation2012LepTag(const Electron *ele, const Vertex *vertex, 
                                                  const PFCandidateCol *PFCands, 
                                                  Double_t rho,
                                                  ElectronTools::EElectronEffectiveAreaTarget,
                                                  const ElectronCol *goodElectrons = 0,
                                                  const MuonCol *goodMuons = 0,
                                                  Double_t dRMax = 0.3, Bool_t isDebug=kFALSE);

    static Double_t PFMuonIsolationRhoCorr(Muon const*, Double_t rho, MuonTools::EMuonEffectiveAreaTarget);

    static void PFEGIsoFootprintRemoved(Particle const*, Vertex const*, PFCandidateCol const*, Double_t dR, Double_t& chIso, Double_t& nhIso, Double_t& phIso);

    static Double_t BetaM(const TrackCol *tracks, const Muon *p, const Vertex *vertex, 
                          Double_t ptMin, Double_t  delta_z, Double_t extRadius,
                          Double_t intRadius);
    static Double_t BetaMwithPUCorrection(const PFCandidateCol *PFNoPileUP, const PFCandidateCol *PFPileUP, const Muon *p, Double_t extRadius, Double_t* isoArr = 0);
    static Double_t BetaE(const TrackCol *tracks, const Electron *p, const Vertex *vertex, 
                          Double_t ptMin, Double_t  delta_z, Double_t extRadius,
                          Double_t intRadius);

    // method added by F.Stoeckli: computes the track isolation with NO constraint on the OV-track compatibility
    static Double_t TrackIsolationNoPV(const mithep::Particle*, const BaseVertex*, 
                                       Double_t extRadius, 
                                       Double_t intRadius, 
                                       Double_t ptLow, 
                                       Double_t etaStrip,
                                       Double_t maxD0,
                                       mithep::TrackQuality::EQuality,
                                       const mithep::Collection<mithep::Track> *tracks,
                                       UInt_t maxNExpectedHitsInner = 999,
                                       const mithep::DecayParticleCol *conversions = 0);

    // methods for Hgg BaseLien Selection. These isolation are stupid, but what can we do.... ;(
    static Double_t CiCTrackIsolation(const mithep::Photon*, 
                                      const BaseVertex*, 
                                      Double_t extRadius, 
                                      Double_t intRadius, 
                                      Double_t ptLow, 
                                      Double_t etaStrip,
                                      Double_t maxD0,
                                      Double_t maxDZ,
                                      const mithep::Collection<mithep::Track> *tracks,
                                      unsigned int* worstVtxIdx = NULL,
                                      const mithep::Collection<mithep::Vertex> *vtxs = NULL,
                                      const mithep::Collection<mithep::Electron> *eles = NULL,
                                      bool print=false,
                                      double* ptmax=NULL, 
                                      double* dRmax=NULL);
      
    static Double_t PFChargedIsolation(const mithep::Photon *p, 
                                       const BaseVertex *theVtx, 
                                       Double_t extRadius,
                                       Double_t intRadius,
                                       const PFCandidateCol *PFCands,
                                       unsigned int* worstVtxIndex = NULL,
                                       const mithep::Collection<mithep::Vertex> *vtxs = NULL,
                                       bool print = false);

    static Double_t PFGammaIsolation(const mithep::Photon *p, 
                                     Double_t extRadius,
                                     Double_t intRadius,
                                     const PFCandidateCol *PFCands);                                         

    static Double_t PFNeutralHadronIsolation(const mithep::Photon *p, 
                                             Double_t extRadius,
                                             Double_t intRadius,
                                             const PFCandidateCol *PFCands);                                                                                  
                                         
    static Float_t PFChargedCount(const mithep::Photon*, 
                                  const BaseVertex*, 
                                  Double_t extRadius, 
                                  Double_t intRadius, 
                                  Double_t ptLow, 
                                  Double_t etaStrip,
                                  Double_t maxD0,
                                  Double_t maxDZ,
                                  const PFCandidateCol *PFCands,
                                  unsigned int* worstVtxIndex = NULL,
                                  const mithep::Collection<mithep::Vertex> *vtxs = NULL,
                                  const mithep::Collection<mithep::Electron> *eles = NULL,
                                  bool print = false,
                                  double* ptmax = NULL, 
                                  double* dRmax = NULL);
  };
}
#endif
