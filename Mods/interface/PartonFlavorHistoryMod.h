//--------------------------------------------------------------------------------------------------
// PartonFlavorHistoryMod
//
// This module looks at the generator information and determines the flavor history of partons.
// It is based on the FlavorHistory modules from CMSSW.
//
// Authors: S.Xie
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_PARTONFLAVORHISTORYMOD_H
#define MITPHYSICS_MODS_PARTONFLAVORHISTORYMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 

class TH1D;
class TH2D;

namespace mithep 
{
  class PartonFlavorHistoryMod : public BaseMod
  {
    public:
      PartonFlavorHistoryMod(const char *name="PartonFlavorHistoryMod", 
                             const char *title="FlavorHistory information module");

      enum MCType {
        kMCTypeUndef = 0,    //not defined
        kMCTypeVLightJets,   //W+Light Jets
        kMCTypeWc,           //W+c
        kMCTypeVQQ           //V+2 heavy flavor jets
      };

      const char          *GetMCPartName()               const { return fMCPartName;              }
      const char          *GetMCSampleType()             const { return fMCSampleType;            }
      const Bool_t         GetApplyPartonFlavorFilter()  const { return fApplyPartonFlavorFilter; }

      void                 SetMCPartName(const char *s)	        { fMCPartName              = s;   }
      void                 SetMCSampleType(const char *s)       { fMCSampleType            = s;   }
      void                 SetApplyPartonFlavorFilter(Bool_t b) { fApplyPartonFlavorFilter = b;   }
       
    protected:
      void                 Process();
      void                 SlaveBegin();

      TString              fMCPartName;               //name of MCParticle branch
      TString              fMCSampleType;             //name of MCSampleType
      Bool_t               fApplyPartonFlavorFilter;  //=true then we apply the filter
      MCType               fMCType;                   //!type of MC 
      TH1D                *fFlavorClassification;     //!histos for flavor history

    ClassDef(PartonFlavorHistoryMod, 1) // Module to gather flavor history information
  };
}
#endif
