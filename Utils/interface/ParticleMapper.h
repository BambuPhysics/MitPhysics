//--------------------------------------------------------------------------------------------------
// Particle Mapper
//
// Tool to avoid looping over entire sets of PF Candidates in an event
// Initialize with a collection of PFCandidates. These fill a grid in eta-phi space with either 
// a default resolution or one specified by the user. Then the object can be used to return the 
// indices of particles in the same bin or all nearby bins. This second option guarantees returning
// all particles within a distance of the resolution (plus some extras). This is much faster than
// looping over all other particles to determine which ones are close.
//
// Authors: D.Abercrombie
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_UTILS_PARTICLEMAPPER_H
#define MITPHYSICS_UTILS_PARTICLEMAPPER_H

#include "MitAna/DataTree/interface/PFCandidateCol.h"

namespace mithep {

  class ParticleMapper {
  public:

    ParticleMapper();
    virtual ~ParticleMapper();

    void Initialize(const PFCandidateCol &Particles,Double_t DeltaEta = 0.3,Double_t DeltaPhi = 0.3,Double_t EtaMax = 5.0);
    
    std::vector<Int_t> GetSurrounding( Int_t index );                // Returns a vector of particle indices in particle bin and adjacent bins
    std::vector<Int_t> GetNearEtaPhi( Double_t eta, Double_t phi);   // Returns a vector of particle indices in area of eta-phi spot given

  private:
    
    Double_t  fDeltaEta;
    Double_t  fDeltaPhi;
    Int_t     fNumParticles;
    Int_t     fNumEtaBins;
    Int_t     fNumPhiBins;
    Int_t     fNumTotBins;

    Int_t *fParticleLocation;
    std::vector<Int_t> *fBinContents;

    std::vector<Int_t> ReturnNear( Int_t etaBin, Int_t phiBin );     // This is the function that the "Get" functions ultimately call

    ClassDef(ParticleMapper, 0)
  };
}

#endif
