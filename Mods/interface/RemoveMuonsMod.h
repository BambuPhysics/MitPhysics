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
#include "MitAna/DataTree/interface/MuonFwd.h"
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

    const char    *GetTriggerMatchName()     const   { return fTriggerMatchName;      }
    const char    *GetTriggerObjectsName()   const   { return fTriggerObjectsName;    }
    const char    *GetPFCandidateName()      const   { return fPFCandidateName;       }
    const char    *GetMuonName()             const   { return fMuonName;              }
    const char    *GetOutputName()           const   { return fOutputName;            }

    void SetTriggerMatchName( const char *name )     { fTriggerMatchName = name;      }
    void SetTriggerObjectsName( const char *name )   { fTriggerObjectsName = name;    }
    void SetPFCandidateName( const char *name )      { fPFCandidateName = name;       }
    void SetOutputName( const char *name )           { fOutputName = name;            }

    void SetPFCandidateFromBranch( Bool_t from )     { fPFCandidateFromBranch = from; }
    void SetMuonFromBranch( Bool_t from )            { fMuonFromBranch = from;        }
    void SetDeltaR( Double_t dr )                    { fDeltaR = dr;                  }
    void GettingTriggerMatch( Bool_t get )           { fGettingTriggerMatch = get;    }

  protected:
    void    SlaveBegin();
    void    SlaveTerminate();
    void    Process();

    const TriggerObjectCol *fTriggerObjects;
    const PFCandidateCol   *fPFCandidates;
    const MuonOArr         *fMuons;

    TString fTriggerMatchName;
    TString fTriggerObjectsName;
    TString fPFCandidateName;
    TString fMuonName;
    TString fOutputName;

    Bool_t   fPFCandidateFromBranch;
    Bool_t   fMuonFromBranch;
    Double_t fDeltaR;
    Bool_t   fGettingTriggerMatch;

    ClassDef(RemoveMuonsMod, 1)
  };
}
#endif
