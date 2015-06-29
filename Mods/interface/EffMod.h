//--------------------------------------------------------------------------------------------------
// EffMod
//
// This module calculates reconstruction efficiency (and fake rate) between reconstructed
// objects and the MC truth by simple matching in DeltaR.
//
// Authors: C.Loizides
//--------------------------------------------------------------------------------------------------

#ifndef MITMODS_MODS_EFFMOD_H
#define MITMODS_MODS_EFFMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/MCParticle.h"

class TH1D;

namespace mithep 
{
  class EffMod : public BaseMod
  {
    public:
      EffMod(const char *name="EffMod", 
             const char *title="Efficiency analysis module");

      void                     SetCol1Name(const char *n)       { fCol1Name = n;  }
      void                     SetCol2Name(const char *n)       { fCol2Name = n;  }
      void                     SetMinPt(Double_t pt)            { fMinPt    = pt; }
      void                     SetMaxEta(Double_t e)            { fMaxEta   = e;  }
      void                     SetRadius(Double_t r)            { fRadius = r; }
      void                     SetType(MCParticle::EPartType t) { fPartType = t; }

    protected:
      void                     Process();
      void                     SlaveBegin();
      void                     Terminate();

      TString                  fCol1Name;           //first  collection name (input)
      TString                  fCol2Name;           //second collection name (input)
      Double_t                 fMinPt;              //minimum pt for col1
      Double_t                 fMaxEta;             //maximum absolute eta for col1
      Double_t                 fRadius;             //radius used for matching
      MCParticle::EPartType    fPartType;           //particle type (default = kNone)
      TH1D                    *fCol1Pt;             //!pt efficiency denominator
      TH1D                    *fCol1Eta;            //!eta efficiency denominator
      TH1D                    *fCol2Pt;             //!pt efficiency numerator
      TH1D                    *fCol2Eta;            //!eta efficiency numerator
      TH1D                    *fGoodPt;             //!pt fake denominator
      TH1D                    *fGoodEta;            //!eta fake denominator
      TH1D                    *fFakePt;             //!pt fake numerator
      TH1D                    *fFakeEta;            //!eta fake numerator

    ClassDef(EffMod, 1) // Efficiency analysis module
  };
}
#endif
