//--------------------------------------------------------------------------------------------------
// RemoveMuonsMod
//
// This mod removes trigger muons from an event's PFCandidates. This is useful for MET studies, 
// especially for PUPPI, where simply subtracting the muons from the recoil does not give the 
// correct result.
//
// Authors: D.Abercrombie
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_REMOVEMUONSMOD_H
#define MITPHYSICS_MODS_REMOVEMUONSMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataCont/interface/Types.h"
#include "MitAna/DataTree/interface/TriggerObjectFwd.h"

namespace mithep
{
  class RemoveMuonsMod : public BaseMod
  {
  public:
    RemoveMuonsMod( const char *name = "RemoveMuonsMod",
                   const char *title = "Remove Muon Module" );
   ~RemoveMuonsMod();

    const char    *GetTriggerMatchName()     const   { return fTriggerMatchName;   }
    const char    *GetTriggerObjectsName()   const   { return fTriggerObjectsName; }
    const char    *GetInputName()            const   { return fInputName;          }
    const char    *GetOutputName()           const   { return fOutputName;         }

    void SetTriggerMatchName( const char *name )     { fTriggerMatchName = name;   }
    void SetTriggerObjectsName( const char *name )   { fTriggerObjectsName = name; }
    void SetInputName( const char *name )            { fInputName = name;          }
    void SetOutputName( const char *name )           { fOutputName = name;         }

    void SetInputFromBranch( Bool_t from )           { fInputFromBranch = from;    }
    void SetDeltaR( Double_t dr )                    { fDeltaR = dr;               }

  protected:
    void    SlaveBegin();
    void    SlaveTerminate();
    void    Process();

    const TriggerObjectCol *fTriggerObjects;
    const PFCandidateCol   *fPFCandidates;

    TString fTriggerMatchName;
    TString fTriggerObjectsName;
    TString fInputName;
    TString fOutputName;

    Bool_t  fInputFromBranch;
    Bool_t  fDeltaR;

    ClassDef(RemoveMuonsMod, 1)
  };
}
#endif
