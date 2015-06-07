//--------------------------------------------------------------------------------------------------
// JetIDMod
//
// This module applies jet identification criteria and exports a pointer to a collection
// of "good jet" according to the specified identification scheme.
//
// Authors: S.Xie
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_JETIDMOD_H
#define MITPHYSICS_MODS_JETIDMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitPhysics/Utils/interface/JetIDMVA.h"

namespace mithep 
{
  class JetIDMod : public BaseMod
  {
    public:
      JetIDMod(const char *name="JetIDMod", 
               const char *title="Jet identification module");

      const char       *GetInputName()                 const { return fJetsName;            }
      const char       *GetGoodName()                  const { return GetGoodJetsName();    }
      const char       *GetGoodJetsName()              const { return fGoodJetsName;        }
      const char       *GetOutputName()                const { return GetGoodJetsName();    }
      Double_t          GetPtCut()                     const { return fJetPtCut;            }
      Bool_t            GetUseCorrection()             const { return fUseJetCorrection;    }
      Double_t          GetEtaMaxCut()                 const { return fJetEtaMaxCut;        }
      Double_t          GetJetEEMFractionMinCut()      const { return fJetEEMFractionMinCut;}
      Bool_t            GetApplyBetaCut()              const { return fApplyBetaCut;        }
      Bool_t            GetApplyMVACut()               const { return fApplyMVACut;         }
      void              SetGoodJetsName(const char *name)    { fGoodJetsName = name;       }
      void              SetGoodName(const char *name)        { SetGoodJetsName(name);      }
      void              SetInputName(const char *name)       { fJetsName = name;           }
      void              SetOutputName(const char *name)      { SetGoodJetsName(name);      }
      void              SetPtCut(Double_t cut)               { fJetPtCut = cut;            }
      void              SetUseCorrection(Bool_t b)           { fUseJetCorrection = b;      }
      void              SetEtaMaxCut(Double_t cut)           { fJetEtaMaxCut = cut;        }
      void              SetJetEEMFractionMinCut(Double_t cut){ fJetEEMFractionMinCut = cut;}
      void              SetApplyBetaCut(Bool_t b)            { fApplyBetaCut = b;          }
      void              SetApplyMVACut(Bool_t b)             { fApplyMVACut = b;           }
      void              SetApplyMVACHS(Bool_t b)             { fApplyMVACHS = b;           }

    protected:
      void              Process();
      void              SlaveBegin();
      void              SlaveTerminate();

      TString           fJetsName;              //name of jet collection (input)
      TString           fGoodJetsName;          //name of good jets collection (output)
      TString           fVertexName;	        //name of vertex collection
      Bool_t            fUseJetCorrection;      //=true then use corrected energy
      Double_t          fJetPtCut;              //jet pt cut
      Double_t          fJetEtaMaxCut;          //jet eta max cut
      Double_t          fJetEEMFractionMinCut;  //jet Eem fraction min cut
      Bool_t            fApplyBetaCut;          //=true then apply beta cut
      Bool_t            fApplyMVACut;           //=true then apply MVA cut
      Bool_t            fApplyMVACHS;           //=true then apply MVA for CHS
      const VertexCol  *fVertices;	        //Vertices branches
      JetIDMVA         *fJetIDMVA;
      ClassDef(JetIDMod, 1) // Jet identification module
  };
}
#endif
