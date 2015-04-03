// $Id: runHgg.C,v 1.5 2012/03/22 15:54:08 bendavid Exp $
#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <TSystem.h>
#include <TProfile.h>
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/Catalog/interface/Catalog.h"
#include "MitAna/TreeMod/interface/Analysis.h"
#include "MitAna/TreeMod/interface/HLTMod.h"
#include "MitAna/PhysicsMod/interface/RunLumiSelectionMod.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Mods/interface/GoodPVFilterMod.h"
#include "MitPhysics/Mods/interface/MuonIDMod.h"
#include "MitPhysics/Mods/interface/ElectronIDMod.h"
#include "MitPhysics/Mods/interface/ElectronCleaningMod.h"
#include "MitPhysics/Mods/interface/PhotonIDMod.h"
#include "MitPhysics/Mods/interface/PhotonCleaningMod.h"
#include "MitPhysics/Mods/interface/MergeLeptonsMod.h"
#include "MitPhysics/Mods/interface/JetCorrectionMod.h"
#include "MitPhysics/SelMods/interface/TwoPhotonSelMod.h"
#include "MitPhysics/SelMods/interface/MCParticleSelMod.h"
#include "MitAna/PhysicsMod/interface/SkimMod.h"
#include "MitAna/TreeMod/interface/OutputMod.h"

#endif

//--------------------------------------------------------------------------------------------------
void runH2gSkim(const char *fileset    = "",
		const char *skim       = "noskim",
		const char *dataset    = "f11--h120gg-gf-v14b-pu",
		const char *book       = "cern/filefi/025",
		const char *catalogDir = "/home/mitprod/catalog",
		const char *outputName = "h2g",
		int         nEvents    = -1)
{
  //------------------------------------------------------------------------------------------------
  // some parameters get passed through the environment
  //------------------------------------------------------------------------------------------------
  char json[1024], overlap[1024];
  float overlapCut = -1;

  if (gSystem->Getenv("MIT_PROD_JSON"))
    sprintf(json,   "%s",gSystem->Getenv("MIT_PROD_JSON"));
  else {
    sprintf(json, "%s", "~");
    //printf(" JSON file was not properly defined. EXIT!\n");
    //return;
  } 

  TString jsonFile = TString("/home/cmsprod/cms/json/") + TString(json);
  Bool_t  isData   = ( (jsonFile.CompareTo("/home/cmsprod/cms/json/~") != 0) );
  
  if (gSystem->Getenv("MIT_PROD_OVERLAP")) {
    sprintf(overlap,"%s",gSystem->Getenv("MIT_PROD_OVERLAP"));
    if (EOF == sscanf(overlap,"%f",&overlapCut)) {
      printf(" Overlap was not properly defined. EXIT!\n");
      return;
    }
  }
  else {
     sprintf(overlap,"%s", "-1.0");
    //printf(" OVERLAP file was not properly defined. EXIT!\n");
    //return;
  } 

  printf("\n Initialization worked. \n\n");

  //------------------------------------------------------------------------------------------------
  // some global setups
  //------------------------------------------------------------------------------------------------
  using namespace mithep;
  gDebugMask  = Debug::kGeneral;
  gDebugLevel = 3;

  //------------------------------------------------------------------------------------------------
  // set up information
  //------------------------------------------------------------------------------------------------
  RunLumiSelectionMod *runLumiSel = new RunLumiSelectionMod;
  runLumiSel->SetAcceptMC(kTRUE);                          // Monte Carlo events are always accepted

  
  // only select on run- and lumisection numbers when valid json file present
  if ((jsonFile.CompareTo("/home/cmsprod/cms/json/~") != 0) &&
      (jsonFile.CompareTo("/home/cmsprod/cms/json/-") != 0)   ) {
    runLumiSel->AddJSONFile(jsonFile.Data());
  }
  if ((jsonFile.CompareTo("/home/cmsprod/cms/json/-") == 0)   ) {
    printf("\n WARNING -- Looking at data without JSON file: always accept.\n\n");
    runLumiSel->SetAbortIfNotAccepted(kFALSE);   // accept all events if there is no valid JSON file
  }


  TwoPhotonSelMod *twophotonsel = new TwoPhotonSelMod;
  twophotonsel->SetPtMin(20.);
  twophotonsel->SetMassMin(55.);

  MCParticleSelMod *mcselmod  = new MCParticleSelMod;
  mcselmod->AddAbsPdgId(11);
  mcselmod->AddAbsPdgId(12);
  mcselmod->AddAbsPdgId(13);
  mcselmod->AddAbsPdgId(14);
  mcselmod->AddAbsPdgId(15);
  mcselmod->AddAbsPdgId(16);
  mcselmod->AddAbsPdgId(23);
  mcselmod->AddAbsPdgId(24);
  mcselmod->AddAbsPdgId(25);


  SkimMod<MCParticle> *skimMod = 0;
  if (! isData) {
    skimMod = new SkimMod<MCParticle>;
    skimMod->SetBranchName("MCParticles");
  }

  OutputMod *outMod = new OutputMod;
  outMod->Keep("*");
  outMod->Drop("CaloTowers");
  outMod->Drop("MCParticles");
  outMod->Drop("*Jets*");
  outMod->Keep("AKT5GenJets");
  outMod->Keep("AKt5PFJets");

  if (! isData)
    outMod->AddNewBranch("SkmMCParticles");


  TString skmRootFile = TString(outputName);
  skmRootFile += TString("_") + TString(dataset) + TString("_") + TString(skim);
  if (TString(fileset) != TString(""))
    skmRootFile += TString("_") + TString(fileset);
  outMod->SetFileName(skmRootFile);
  outMod->SetPathName(".");


  //------------------------------------------------------------------------------------------------
  // making analysis chain
  //------------------------------------------------------------------------------------------------
  // this is how it always starts
  runLumiSel      ->Add(twophotonsel);

  if (isData) { 
    twophotonsel  ->Add(outMod);
  }
  else {
    twophotonsel  ->Add(mcselmod);
    mcselmod      ->Add(skimMod);
    skimMod       ->Add(outMod);  
  }
  
  //------------------------------------------------------------------------------------------------
  // setup analysis
  //------------------------------------------------------------------------------------------------
  Analysis *ana = new Analysis;
  ana->SetUseHLT       (kTRUE);
  ana->SetKeepHierarchy(kTRUE);
  ana->SetSuperModule  (runLumiSel);
  ana->SetPrintScale   (100);
  if (nEvents >= 0)
    ana->SetProcessNEvents(nEvents);

  //------------------------------------------------------------------------------------------------
  // organize input
  //------------------------------------------------------------------------------------------------
  Catalog *c = new Catalog(catalogDir);
  TString skimdataset = TString(dataset)+TString("/") +TString(skim);
  Dataset *d = NULL;
  TString bookstr = book;
  if (TString(skim).CompareTo("noskim") == 0)
    d = c->FindDataset(bookstr,dataset,fileset);
  else 
    d = c->FindDataset(bookstr,skimdataset.Data(),fileset);
  ana->AddDataset(d);
  
  //------------------------------------------------------------------------------------------------
  // organize output
  //------------------------------------------------------------------------------------------------
  TString rootFile = TString(outputName);
  rootFile += TString("_") + TString(dataset) + TString("_") + TString(skim);
  if (TString(fileset) != TString(""))
    rootFile += TString("_") + TString(fileset);
  rootFile += TString(".root");
  //ana->SetOutputName(rootFile.Data());
  //ana->SetCacheSize(64*1024*1024);
  //ana->SetCacheSize(0);

  ana->SetOutputName("test.root");
  
  //------------------------------------------------------------------------------------------------
  // Say what we are doing
  //------------------------------------------------------------------------------------------------
  printf("\n==== PARAMETER SUMMARY FOR THIS JOB ====\n");
  printf("\n JSON file: %s\n  and overlap cut: %f (%s)\n",jsonFile.Data(),overlapCut,overlap);
  printf("\n Rely on Catalog: %s\n",catalogDir);
  printf("  -> Book: %s  Dataset: %s  Skim: %s  Fileset: %s <-\n",book,dataset,skim,fileset);
  printf("\n Root output: %s\n\n",rootFile.Data());  
  printf("\n========================================\n");

  //------------------------------------------------------------------------------------------------
  // run the analysis after successful initialisation
  //------------------------------------------------------------------------------------------------
  ana->Run(!gROOT->IsBatch());

  return;
}
