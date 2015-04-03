//--------------------------------------------------------------------------------------------------
// $Id: 
//
// MetTools
//
// Authors: M. Zanetti
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_UTILS_METTOOLS_H
#define MITPHYSICS_UTILS_METTOOLS_H

#include "MitAna/DataTree/interface/ParticleCol.h"
#include "MitAna/DataTree/interface/MetCol.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitAna/DataTree/interface/GenericParticleCol.h"
#include "MitCommon/MathTools/interface/MathUtils.h"

namespace mithep {

  class MetTools {

  public:

    MetTools(const MuonCol *fMuons, const PFCandidateCol *fPFCandidates, 
             const Vertex *fVertex, float deltaZCut = 0.1, float ptCut = 8.0, float etaCut = 5.0);
    MetTools(const ElectronCol *fElectrons, const PFCandidateCol *fPFCandidates, 
             const Vertex *fVertex, float deltaZCut = 0.1, float ptCut = 8.0, float etaCut = 5.0);

    MetTools(const MuonCol *fMuons, const PFCandidateCol *fPFCandidates, const PFJetCol *fPFJets,
             const Vertex *fVertex, float deltaZCut = 0.1, float ptCut = 8.0, float etaCut = 5.0);
    MetTools(const ElectronCol *fElectrons, const PFCandidateCol *fPFCandidates, const PFJetCol *fPFJets, 
             const Vertex *fVertex, float deltaZCut = 0.1, float ptCut = 8.0, float etaCut = 5.0);

    MetTools(const MuonCol *fMuons, const ElectronCol *fElectrons, const PFCandidateCol *fPFCandidates, 
             const Vertex *fVertex, float deltaZCut = 0.1, float ptCut = 8.0, float etaCut = 5.0, float intRadius = 0.0,
	     const GenericParticle *genP = NULL);

    MetTools(const MuonCol *fMuons, const ElectronCol *fElectrons, const PFCandidateCol *fPFCandidates, const PFJetCol *fPFJets, 
             const Vertex *fVertex, float deltaZCut = 0.1, float ptCut = 8.0, float etaCut = 5.0, float intRadius = 0.0,
	     const GenericParticle *genP = NULL);


    virtual ~MetTools() {}
    
    void AddToCorrectedTrackMet                   ( const Particle *p, bool debug = false);
    void AddToCorrectedMet                        ( const Particle *p );
    void AddToRecoil                              ( const Particle *p );
    void RemoveParticleInIsoConeFromTrackMet      ( const Particle *p, 
                                                    const PFCandidateCol *fPFCandidates, 
                                                    const Vertex *fVertex, float deltaZCut, 
                                                    float deltaR, bool debug = false );
    void RemoveParticleInIsoConeFromCorrectedMet  ( const Particle *p, 
                                                    const PFCandidateCol *fPFCandidates, 
                                                    const Vertex *fVertex, 
                                                    float deltaZCut, float ptCut, float etaCut, 
                                                    float deltaR);
    void RemoveParticleInIsoConeFromRecoil        ( const Particle *p, 
                                                    const PFCandidateCol *fPFCandidates, 
                                                    const Vertex *fVertex, 
                                                    float deltaZCut, float ptCut, float etaCut, 
                                                    float deltaR);
    Met  GetMinimumMet                            (const Met *UncorrectedMet);
    Met  GetMinimumTrackMet                       (const Met *UncorrectedMet);
    Met  GetCorrectedMet()         { return fCorrectedMet; }
    Met  GetCorrectedTrackMet()    { return fCorrectedTrackMet; }
    Met  GetCHSMet()               { return fCHSMet; }
    Met  GetNHSMet()               { return fNHSMet; }
    FourVectorM  Recoil()          { return fRecoil; }
    FourVectorM  ChargedRecoil()   { return fChargedRecoil; }

    template<class V>
    double GetProjectedMet(const V *fV, const Met *UncorrectedMet);
    template<class V>
    double GetProjectedMet(const V *fV);
    template<class V>
    double GetProjectedTrackMet(const V *fV);

  private:
    Met fCorrectedMet;
    Met fCorrectedTrackMet;
    Met fCHSMet;
    Met fNHSMet;
    FourVectorM fRecoil;
    FourVectorM fChargedRecoil;
    
    ClassDef(MetTools, 0) // Met tools
      };

  template<class V>
  double MetTools::GetProjectedMet(const V *fV, const Met *UncorrectedMet) {
    double projectedMet = UncorrectedMet->Pt();
    double minDPhi = 999;
    for (UInt_t m = 0; m < fV->GetEntries(); ++m) {
      if (fabs(MathUtils::DeltaPhi(UncorrectedMet->Phi(), fV->At(m)->Phi())) < minDPhi) {
	minDPhi = fabs(MathUtils::DeltaPhi(UncorrectedMet->Phi(), fV->At(m)->Phi()));
      }
    }
    if (minDPhi < TMath::Pi()/2.) return projectedMet = projectedMet * sin(minDPhi);
    return projectedMet;
  }

  template<class V>
  double MetTools::GetProjectedMet(const V *fV) {
    double projectedMet = fCorrectedMet.Pt();
    double minDPhi = 999;
    for (UInt_t m = 0; m < fV->GetEntries(); ++m) {
      if (fabs(MathUtils::DeltaPhi(fCorrectedMet.Phi(), fV->At(m)->Phi())) < minDPhi) {
	minDPhi = fabs(MathUtils::DeltaPhi(fCorrectedMet.Phi(), fV->At(m)->Phi()));
      }
    }
    if (minDPhi < TMath::Pi()/2.) return projectedMet = projectedMet * sin(minDPhi);
    return projectedMet;
  }

  template<class V>
  double MetTools::GetProjectedTrackMet(const V *fV) {
    double projectedMet = fCorrectedTrackMet.Pt();
    double minDPhi = 999;
    for (UInt_t m = 0; m < fV->GetEntries(); ++m) {
      if (fabs(MathUtils::DeltaPhi(fCorrectedTrackMet.Phi(), fV->At(m)->Phi())) < minDPhi) {
	minDPhi = fabs(MathUtils::DeltaPhi(fCorrectedTrackMet.Phi(), fV->At(m)->Phi()));
      }
    }
    if (minDPhi < TMath::Pi()/2.) return projectedMet = projectedMet * sin(minDPhi);
    return projectedMet;
  }

}

#endif
