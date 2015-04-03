//--------------------------------------------------------------------------------------------------
// $Id: MVAMetMod.h,v 1.4 2012/04/07 16:45:04 ceballos Exp $
//
// MVAMetMod
//
// Example on how to call regressed MET
//
// Authors: P.Harris
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_MVAMETMOD_H
#define MITPHYSICS_MODS_MVAMETMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/DataTree/interface/PFMetCol.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityCol.h"
#include "MitPhysics/Utils/interface/MVAMet.h"

namespace mithep 
{
  class MVAMetMod : public BaseMod
  {
    public:
      MVAMetMod(const char *name="MVAMetMod", 
		const char *title="MVAMet example");
      void   SetJetsName(TString s) { fJetsName = s;}   

    protected:
      void                Process();
      void                SlaveBegin();
      void                SlaveTerminate();

      TString             fMVAMetName;
      TString             fJetsName  ;
      TString             fPFCandName;
      TString             fVertexName;
      TString             fPFMetName ;
      TString             fRhoName   ;

      const JetCol                 *fJets;
      const PFCandidateCol         *fCands;
      const VertexCol              *fVertices;
      const PFMetCol               *fPFMet;
      const PileupEnergyDensityCol *fRhoCol;

      MVAMet         *fMVAMet;
      ClassDef(MVAMetMod, 1) // Jet identification module
  };
}
#endif
