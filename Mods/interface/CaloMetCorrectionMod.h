//--------------------------------------------------------------------------------------------------
// CaloMetCorrectionMod
//
// This module applies jet energy corrections to CaloMet on the fly at analysis level, reading
// the jet corrections directly out of the bambu jet objects
//
// Authors: J.Bendavid
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_CALOMETCORRECTIONMOD_H
#define MITPHYSICS_MODS_CALOMETCORRECTIONMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 

namespace mithep 
{
  class CaloMetCorrectionMod : public BaseMod
  {
    public:
      CaloMetCorrectionMod(const char *name="CaloMetCorrectionMod", 
               const char *title="Met correction module");

      const char       *GetInputName()                 const      { return fMetName;               }   
      const char       *GetCorrectedName()             const      { return GetCorrectedMetName();  }     
      const char       *GetCorrectedJetsName()         const      { return fCorrectedJetsName;     }  
      const char       *GetCorrectedMetName()          const      { return fCorrectedMetName;      }    
      const char       *GetOutputName()                const      { return GetCorrectedMetName();  }    
      void              SetCorrectedJetsName(const char *name)    { fCorrectedJetsName = name;     }
      void              SetCorrectedMetName(const char *name)     { fCorrectedMetName = name;      }
      void              SetCorrectedName(const char *name)        { SetCorrectedMetName(name);     }    
      void              SetInputName(const char *name)            { fMetName = name;               }
      void              SetMaxJetEMF(Double_t x)                  { fMaxJetEMF = x;                }
      void              SetMinJetPt(Double_t x)                   { fMinJetPt = x;                 }
      void              SetOutputName(const char *name)           { SetCorrectedMetName(name);     }          

    protected:
      void              Process();

      TString           fMetName;               //name of CaloMet collection (input)
      TString           fCorrectedMetName;      //name of corrected met collection (output)
      TString           fCorrectedJetsName;     //corrected CaloJets used for met correction
      Double_t          fMinJetPt;              //min uncorrected pt for input jets
      Double_t          fMaxJetEMF;             //max emf for input jets

      ClassDef(CaloMetCorrectionMod, 1) // Met correction module
  };
}
#endif
