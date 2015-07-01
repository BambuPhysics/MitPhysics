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

namespace mithep 
{
  class PuppiMod : public BaseMod
  {
    public:
      PuppiMod(const char *name="PuppiMod", 
               const char *title="Puppi module");
     ~PuppiMod();

      Int_t GetParticleType( const PFCandidate *cand );
      Int_t GetEtaBin( const PFCandidate *cand );
      Double_t Chi2fromDZ( Double_t dz );

      const char   *GetVertexesName()              const     { return fVertexesName;       }
      const char   *GetInputName()                 const     { return fPFCandidatesName;   }   
      const char   *GetOutputName()                const     { return fPuppiParticlesName; }
      void SetEtaConfigName( const char *name )              { fEtaConfigName = name;      }
      void SetVertexesName( const char *name )               { fVertexesName = name;       }
      void SetInputName( const char *name )                  { fPFCandidatesName = name;   }
      void SetOutputName( const char *name )                 { fPuppiParticlesName = name; }

      void SetRMin( Double_t RMin )                          { fRMin = RMin;               }
      void SetR0( Double_t R0 )                              { fR0 = R0;                   }
      void SetAlpha( Double_t Alpha )                        { fAlpha = Alpha;             }
      void SetBeta( Double_t Beta )                          { fBeta = Beta;               }
      
      void SetD0Cut( Double_t cut )                          { fD0Cut = cut;               }
      void SetDZCut( Double_t cut )                          { fDZCut = cut;               }

      void SetMinWeightCut( Double_t cut )                   { fMinWeightCut = cut;        }

      void SetRMSScaleFactor( Double_t fact )                { fRMSScaleFactor = fact;     }
      void SetTrackUncertainty( Double_t sig )               { fTrackUncertainty = sig;    }

      void SetKeepPileup( Bool_t keep )                      { fKeepPileup = keep;         }
      void SetInvert( Bool_t invert )                        { fInvert = invert;           }
      void SetApplyCHS( Bool_t apply )                       { fApplyCHS = apply;          }

    protected:
      void                  SlaveBegin();
      void                  SlaveTerminate();
      void                  Process();
      
      TString               fEtaConfigName;       // Name of the configuration file with eta tables
      TString               fVertexesName;        // Name of vertices collection used for PV
      TString               fPFCandidatesName;    // Name of PFCandidate collection (input)
      TString               fPuppiParticlesName;  // Name of Puppi Particle collection (output)

      const VertexCol *fVertexes;                 // Vertex branch
      const PFCandidateCol *fPFCandidates;        // Particle flow branch
      PFCandidateArr *fPuppiParticles;            // The output collection

      Double_t fRMin;                             // Minimum dR cut for summing up surrounding particles
      Double_t fR0;                               // Maximum dR cut for summing up surrounding particles
      Double_t fAlpha;                            // Parameter for weighting pt of surrounding particles
      Double_t fBeta;                             // Parameter for weighting dR of surrounding particles

      Double_t fD0Cut;                            // D0 cut for charged particle vertex matching
      Double_t fDZCut;                            // DZ cut for charged particle vertex matching

      Double_t fMinWeightCut;                     // Minimum weight to drop it to zero

      Double_t fRMSScaleFactor;                   // A scale factor for RMS
      Double_t fTrackUncertainty;                 // The experimental uncertainty in the track fit to vertex distance

      Bool_t fKeepPileup;                         // Keep pileup with zero weight (for debugging)
      Bool_t fInvert;                             // Option to invert weights
      Bool_t fApplyCHS;                           // This will force weights to 0 or 1 for tracked particles

      // These are parameters that are functions of Eta hopefully we can be more clever some day
      Int_t fNumEtaBins;                          // This is the number of eta regions we are dividing into
      std::vector<Double_t> fMaxEtas;             // These are the maximum etas for each region of the table
      std::vector<Double_t> fMinPts;              // Various Pt cuts
      std::vector<Double_t> fMinNeutralPts;       // Minimum Pt cut on neutral particles (after weighting)
      std::vector<Double_t> fMinNeutralPtSlopes;  // Predicted slope of neutral particles as function of PU
      std::vector<Double_t> fRMSEtaSFs;           // Scale factor for RMS as function of Eta
      std::vector<Double_t> fMedEtaSFs;           // Scale factor for median as a function of Eta
      std::vector<Double_t> fEtaMaxExtraps;       // I think this is the maximum eta to calculate median alphas?

      ClassDef(PuppiMod, 1)                       // Puppi module
  };
}
#endif
