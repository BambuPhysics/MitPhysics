//--------------------------------------------------------------------------------------------------
// $Id: $
//
// VTagMod
//
// This module analyzes input jets and tries to identify them as jets corresponding to Vector
// bosons. So called fat jets are identified and groomed to identify substructure within each jet
// to ultimately lead to a so called VTag (Vector boson tag), where the vector boson decays
// hadronically.
//
// Documentation here: 
//
// * Experiment - CMS
//   https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsJME13006
//   https://indico.cern.ch/getFile.py/access?contribId=6&resId=0&materialId=slides&confId=266244
//
// * Theory
//   http://arxiv.org/abs/1011.2268
//
// Authors: C.Paus
//--------------------------------------------------------------------------------------------------
#ifndef MITPHYSICS_MODS_VTAGMOD_H
#define MITPHYSICS_MODS_VTAGMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"

namespace fastjet {
  class Pruner;
  class JetDefinition;
  class GhostedAreaSpec;
  class AreaDefinition;
  class PseudoJet;
}

namespace mithep
{
  class VTagMod : public BaseMod
  {
    public:
      VTagMod(const char *name  = "VTagMod",
              const char *title = "Vector Boson Tagging module");

      const char                   *GetGoodVTagsName()              const { return fGoodVTagsName;     }
      const char                   *GetOutputName()                 const { return GetGoodVTagsName(); }

      void                          SetGoodVTagsName(const char *n)       { fGoodVTagsName = n;        }
			            
    protected:		            
      void                          Process();
      void                          SlaveBegin();
      void                          Terminate();

    private:
      float                         GetTau(fastjet::PseudoJet &iJet,int iN, float iKappa);
			            
      TString                       fJetsBranchName;        //(i) name of jet collection to analyze
      TString                       fPFCandidatesName;      //(i) name of PF candidates colleciont
      TString                       fGoodVTagsName;         //(o) name of VTags collection

      // Bambu inputs needed
      const PFCandidateCol         *fPFCandidates;          //Particle Flow candidates

      // Objects/configs from fastjet we want to use
      double                        fConeSize;
      fastjet::Pruner              *fPruner;
      fastjet::JetDefinition       *fCAJetDef;
      fastjet::GhostedAreaSpec     *fActiveArea;
      fastjet::AreaDefinition      *fAreaDefinition;

    ClassDef(VTagMod, 1) // Vector boson tagging module
  };
}
#endif
