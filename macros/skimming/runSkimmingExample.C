// $Id: runSkimmingExample.C,v 1.5 2009/07/17 20:48:23 loizides Exp $

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <TSystem.h>
#include <TRandom.h>
#include <TParameter.h>
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/TreeMod/interface/Analysis.h"
#include "MitAna/TreeMod/interface/OutputMod.h"
#endif

namespace mithep 
{
  class Sel : public BaseMod
  {
    public:
      Sel(Int_t r=10) : fRand(r), fP("RValue",0) {}

    protected:
      void               Process() { Double_t r = gRandom->Rndm(); 
                                     fP.SetVal(r);
                                     if (r>1./fRand) SkipEvent(); }
      void               SlaveBegin() { PublishObj(&fP);          }

      Int_t              fRand; //inverse fraction of to be skipped events
      TParameter<double> fP;    //parameter holding the random number

    ClassDef(Sel, 1)
  };
}

//--------------------------------------------------------------------------------------------------
void runSkimmingExample(const char *files, 
                        const Double_t fkeep=0.1, 
                        const char *prefix="skimtest", 
                        UInt_t nev=0)
{
  using namespace mithep;
  gDebugMask  = Debug::kAnalysis;
  gDebugLevel = 1;

  Sel *smod = new Sel((Int_t)(1./fkeep));

  OutputMod *omod = new OutputMod;
  omod->SetFileName(prefix);
  //omod->SetCheckBrDep(kFALSE);
  omod->AddNewBranch("RValue");

  omod->Keep("*");
  if (0) { // more complex case
    omod->Drop("*");
    omod->Keep("*HLT*");
    omod->Keep("Muons");
    omod->Keep("MCParticles");
    omod->Keep("Electrons");
    omod->Keep("*Tracks*");
    omod->Keep("*Clusters*");
  }
  smod->Add(omod);

  // set up analysis
  Analysis *ana = new Analysis;
  if (0) {
    ana->SetHLTTreeName("HLT_E29");
    ana->SetHLTObjsName("HLTObjects_E29");
  }

  ana->SetSuperModule(smod);
  if (nev)
    ana->SetProcessNEvents(nev);
  ana->AddFile(files);
  
  // run the analysis after successful initialisation
  ana->Run(0);
}
