// $Id: runMCParticlesVal.C,v 1.1 2009/03/23 09:09:05 loizides Exp $

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <TROOT.h>
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/TreeMod/interface/Analysis.h"
#include "MitPhysics/Validation/interface/MCParticlesValMod.h"
#endif

//--------------------------------------------------------------------------------------------------
void runMCParticlesVal(const char *files, UInt_t nev=0)
{
  using namespace mithep;
  gDebugMask  = Debug::kAnalysis;
  gDebugLevel = 1;

  MCParticlesValMod *mod = new MCParticlesValMod;

  // set up analysis
  Analysis *ana = new Analysis;
  ana->SetSuperModule(mod);
  ana->AddFile(files);
  ana->SetUseHLT(0); 
  if (nev>0) 
    ana->SetProcessNEvents(nev);
  if (gROOT->IsBatch()) 
    ana->SetOutputName("mcparticlesval.root");
  
  // run the analysis after successful initialisation
  ana->Run(!gROOT->IsBatch());
}
