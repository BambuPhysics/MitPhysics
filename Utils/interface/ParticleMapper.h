//--------------------------------------------------------------------------------------------------
// Particle Mapper
//
// Tool to avoid looping over entire sets of PF Candidates in an event
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
    
    std::vector<Int_t> GetInBin( Int_t index );               // Returns a vector of particle indices in particle bin
    std::vector<Int_t> GetSurrounding( Int_t index );         // Returns a vector of particle indices in particle bin and adjacent bins

  private:
    
    Int_t fNumParticles;
    Int_t fNumEtaBins;
    Int_t fNumPhiBins;
    Int_t fNumTotBins;

    Int_t *fParticleLocation;
    std::vector<Int_t> *fBinContents;

    ClassDef(ParticleMapper, 0)
  };
}

#endif
