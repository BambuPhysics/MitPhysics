//--------------------------------------------------------------------------------------------------
// MetCorrectionMod
//
// This module applies Type0/1 and XY shift MET corrections on the fly at analysis level.
// The methods are synchronized with JetMET POG 2012 studies, documented in this twiki
// https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMetAnalysis#7_7_6_MET_Corrections
//
// Authors: L.Di Matteo
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_METCORRECTIONMOD_H
#define MITPHYSICS_MODS_METCORRECTIONMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/PFMetCol.h"
#include "MitAna/DataCont/interface/Types.h"

#include "TFormula.h"

namespace mithep 
{
  class MetCorrectionMod : public BaseMod
  {
    public:
      MetCorrectionMod(const char *name="MetCorrectionMod", 
		       const char *title="Met correction module");
     ~MetCorrectionMod();

      const char   *GetInputName()                 const     { return fMetName;          }   
      const char   *GetCorrectedName()             const     { return fCorrectedMetName; }     
      const char   *GetJetsName()                  const     { return fJetsName;         }   
      const char   *GetCorrectedJetsName()         const     { return fCorrectedJetsName;}     
      const char   *GetOutputName()                const     { return GetCorrectedName();}
      Double_t     GetMinDz()                      const     { return fMinDz;            }
      const char   *GetExprType0()                 const     { return fExprType0;        }
      const char   *GetExprShiftDataPx()           const     { return fExprShiftDataPx;  }
      const char   *GetExprShiftDataPy()           const     { return fExprShiftDataPy;  }
      const char   *GetExprShiftMCPx()             const     { return fExprShiftMCPx;    }
      const char   *GetExprShiftMCPy()             const     { return fExprShiftMCPy;    }
      void         SetInputName(const char *name)            { fMetName = name;          }
      void         SetCorrectedName(const char *name)        { fCorrectedMetName = name; }    
      void         SetOutputName(const char *name)           { fCorrectedMetName = name; }    
      void         SetJetsName(const char *name)             { fJetsName = name;         }     
      void         SetCorrectedJetsName(const char *name)    { fCorrectedJetsName = name;}     
      void         SetMinDz(Double_t d)                      { fMinDz = d;               }     
      void         ApplyType0(bool b)                        { fApplyType0 = b;          }     
      void         ApplyType1(bool b)                        { fApplyType1 = b;          }
      void         ApplyShift(bool b)                        { fApplyShift = b;          }
      void         SetExprType0(const char *expr)            { fExprType0 = expr;        }
      void         SetExprShiftDataPx(const char *expr)      { fExprShiftDataPx = expr;  }
      void         SetExprShiftDataPy(const char *expr)      { fExprShiftDataPy = expr;  }
      void         SetExprShiftMCPx(const char *expr)        { fExprShiftMCPx = expr;    }
      void         SetExprShiftMCPy(const char *expr)        { fExprShiftMCPy = expr;    }   
      void         IsData(bool b)                            { fIsData = b;              }
      void         SetPrint(bool b)                          { fPrint = b;               }

    protected:
      void                  SlaveBegin();
      void                  SlaveTerminate();
      void                  Process();
      
      TString               fMetName;            //name of met collection (input)
      TString               fCorrectedMetName;   //name of corrected met collection (output)
      TString               fJetsName;           //name of uncorrected jet collection (input)
      TString               fCorrectedJetsName;  //name of corrected jet collection (input)
      TString               fPFCandidatesName;   //name of PF candidates collection (input)
      TString               fVertexName;         //name of vertex collection (input)
      Double_t              fMinDz;              //delta Z for separating PU verteces from PV
                           
      bool                  fApplyType0;         //switch on type 0 correction
      bool                  fApplyType1;         //switch on type 1 correction
      bool                  fApplyShift;         //switch on XY shift correction
                                    
      TString               fExprType0;          //expr for type0 correction
      TFormula             *fFormulaType0;       //formula for type0 correction
      TString               fExprShiftDataPx;    //expr for shift correction (data)
      TFormula             *fFormulaShiftDataPx; //formula for type0 correction
      TString               fExprShiftDataPy;    //expr for shift correction (data)
      TFormula             *fFormulaShiftDataPy; //formula for type0 correction
      TString               fExprShiftMCPx;      //expr for shift correction (MC)
      TFormula             *fFormulaShiftMCPx;   //formula for type0 correction
      TString               fExprShiftMCPy;      //expr for shift correction (MC)
      TFormula             *fFormulaShiftMCPy;   //formula for type0 correction
                                    
      bool                  fIsData;             //flag for data/MC distinction
      bool                  fPrint;              //flag for debug print out
                           
      const PFMetCol       *fPFMet;              //met branch
      const PFCandidateCol *fPFCandidates;       //particle flow branch

      ClassDef(MetCorrectionMod, 1)              // met correction module
  };
}
#endif
