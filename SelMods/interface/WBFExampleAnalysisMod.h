//--------------------------------------------------------------------------------------------------
// $Id $
//
// WBFExampleAnalysisMod
//
// A Module to select qqH events
// and produces some distributions.
//
//
// Authors: Si Xie
//--------------------------------------------------------------------------------------------------

#ifndef ANA_SELMODS_WBFExampleAnalysisMod_H
#define ANA_SELMODS_WBFExampleAnalysisMod_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/MetCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/SuperClusterCol.h"
#include "MitAna/DataTree/interface/DecayParticleCol.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitPhysics/Utils/interface/MuonTools.h"

class TH1D;
class TH2D;

namespace mithep 
{
  class WBFExampleAnalysisMod : public BaseMod
  {
    public:
    WBFExampleAnalysisMod(const char *name="WBFExampleAnalysisMod", 
		 const char *title="Example analysis module with all branches");
      ~WBFExampleAnalysisMod() {}

      const char *GetCleanJetsName()                 const { return fCleanJetsName;                }

      void  SetCleanJetsName(const char *name)             { fCleanJetsName                = name; }
      void  SetJetPtMax(double x)                          { fJetPtMax                        = x; }
      void  SetJetPtMin(double x)                          { fJetPtMin                        = x; }
      void  SetDeltaEtaMin(double x)                       { fDeltaEtaMin                     = x; }
      void  SetDiJetMassMin(double x)                      { fDiJetMassMin                    = x; }

    protected:
      TString                  fCleanJetsName;           //name of clean central jets collection
      double       	       fJetPtMax;
      double       	       fJetPtMin;
      double       	       fDeltaEtaMin;
      double       	       fDiJetMassMin;

      TH1D                    *fWBFSelection;
      TH1D                    *fWBFDeltaPhi;

      TH1D                    *fWBFPtJetMax_NMinusOne;
      TH1D                    *fWBFPtJetMin_NMinusOne;
      TH1D                    *fWBFEta12_NMinusOne;
      TH1D                    *fWBFdeltaEta_NMinusOne;
      TH1D                    *fWBFdijetMass_NMinusOne;
      TH1D                    *fWBFZVar_NMinusOne;

      TH1D                    *fWBFSSMass_afterCuts;   
      TH1D                    *fWBFSSDeltaPhi_afterCuts;    
      TH1D                    *fWBFOSMass_afterCuts;   
      TH1D                    *fWBFOSDeltaPhi_afterCuts;    

      void         Begin();
      void         Process();
      void         SlaveBegin();
      void         SlaveTerminate();
      void         Terminate();      

      ClassDef(WBFExampleAnalysisMod,1) // TAM example analysis module
  };
}
#endif
