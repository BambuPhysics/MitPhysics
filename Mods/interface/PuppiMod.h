//--------------------------------------------------------------------------------------------------
// PuppiMod
//
// This mod makes Puppi weights
//
// Authors: D.Abercrombie
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_PUPPIMOD_H
#define MITPHYSICS_MODS_PUPPIMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/PFMetCol.h"
#include "MitAna/DataCont/interface/Types.h"
#include "MitPhysics/Utils/interface/ParticleMapper.h"

#include "TFormula.h"

namespace mithep 
{
  class PuppiMod : public BaseMod
  {
    public:
      PuppiMod(const char *name="PuppiMod", 
               const char *title="Puppi module");
     ~PuppiMod();

      const char   *GetVertexesName()              const     { return fVertexesName;       }
      const char   *GetInputName()                 const     { return fPFCandidatesName;   }   
      const char   *GetOutputName()                const     { return fPuppiParticlesName; }
      void SetVertexesName( const char *name )               { fVertexesName = name;       }
      void SetInputName( const char *name )                  { fPFCandidatesName = name;   }
      void SetOutputName( const char *name )                 { fPuppiParticlesName = name; }

      void SetRMin( Double_t RMin )                          { fRMin = RMin;               }
      void SetR0( Double_t R0 )                              { fR0 = R0;                   }
      void SetBeta( Double_t Beta )                          { fBeta = Beta;               }
      
      void SetD0Cut( Double_t cut )                          { fD0Cut = cut;               }
      void SetDZCut( Double_t cut )                          { fDZCut = cut;               }

    protected:
      void                  SlaveBegin();
      void                  SlaveTerminate();
      void                  Process();
      
      TString               fVertexesName;        // Name of vertices collection used for PV
      TString               fPFCandidatesName;    // Name of PFCandidate collection (input)
      TString               fPuppiParticlesName;  // Name of Puppi Particle collection (output)

      const VertexCol *fVertexes;                 // Vertex branch
      const PFCandidateCol *fPFCandidates;        // Particle flow branch

      Double_t fRMin;                             // Minimum dR cut for summing up surrounding particles
      Double_t fR0;                               // Maximum dR cut for summing up surrounding particles
      Double_t fBeta;                             // Parameter for weighting dR of surrounding particles

      Double_t fD0Cut;                            // D0 cut for charged particle vertex matching
      Double_t fDZCut;                            // DZ cut for charged particle vertex matching

      ClassDef(PuppiMod, 1)                       // met correction module
  };
}
#endif
