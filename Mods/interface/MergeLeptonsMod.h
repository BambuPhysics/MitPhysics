//--------------------------------------------------------------------------------------------------
// MergeLeptonsMod
//
// This module merges muon and electron collections. (Note if need be this can easily be
// generalized).
//
// Authors: C.Loizides
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_MERGELEPTONSMOD_H
#define MITPHYSICS_MODS_MERGELEPTONSMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/ElectronFwd.h"
#include "MitAna/DataTree/interface/MuonFwd.h"
#include "MitAna/DataTree/interface/ParticleFwd.h"
#include "TH1D.h"

namespace mithep 
{
  class MergeLeptonsMod : public BaseMod 
  {
    public:
      MergeLeptonsMod(const char *name="MergeLeptonsMod", 
                      const char *title="Merging leptons module");

      const char              *GetElectronsName()           const { return fElName;         }
      const char              *GetMergedName()              const { return fMergedName;     }
      const char              *GetMuonsName()               const { return fMuName;         }
      const char              *GetOutputName()              const { return GetMergedName(); }
      void                     SetElectronsName(const char *n)    { fElName=n;              }
      void                     SetMergedName(const char *n)       { fMergedName=n;          }
      void                     SetMuonsName(const char *n)        { fMuName=n;              }
      void                     SetOutputName(const char *n)       { SetMergedName(n);       }

    protected:
      void                     BeginRun();
      void                     Process();
      void                     SlaveBegin();
      void                     SlaveTerminate();

      TString                  fElName;        //name of electrons collection (input)
      TString                  fMuName;        //name of muons collection (input)
      TString                  fMergedName;    //name of merged collection (output)
      ParticleOArr            *fColOut;        //!pointer to merged collection
      TH1D*		       fRecoWMuons;
      TH1D*		       fRecoWElectrons;

    ClassDef(MergeLeptonsMod, 1) // Merging leptons module
  };
}
#endif
