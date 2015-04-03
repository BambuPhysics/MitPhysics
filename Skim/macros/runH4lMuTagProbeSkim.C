// $Id: runH4lSkim.C,v 1.1 2012/06/02 20:46:24 paus Exp $
#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <TSystem.h>
#include <TProfile.h>
#include <TFile.h>
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/Catalog/interface/Catalog.h"
#include "MitAna/TreeMod/interface/Analysis.h"
#include "MitAna/TreeMod/interface/HLTMod.h"
#include "MitAna/TreeMod/interface/OutputMod.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitAna/PhysicsMod/interface/PublisherMod.h"
#include "MitAna/PhysicsMod/interface/RunLumiSelectionMod.h"
#include "MitAna/PhysicsMod/interface/MCProcessSelectionMod.h"
#include "MitPhysics/Mods/interface/GoodPVFilterMod.h"
#include "MitPhysics/Mods/interface/PhotonPairSelector.h"
#include "MitPhysics/Mods/interface/PhotonTreeWriter.h"
#include "MitPhysics/Skim/interface/H4lMuTagProbeSkim.h"
#endif

using namespace mithep;

int  decodeEnv(char* json, char* overlap, float overlapCut, char* path);

//--------------------------------------------------------------------------------------------------
void runH4lMuTagProbeSkim(const char *fileset    = "0000",
			  const char *skim       = "noskim",
			  const char *dataset    = "r11b-pho-pr-v1",
			  const char *book       = "local/filefi/025",
			  const char *catalogDir = "/home/cmsprod/catalog",
			  const char *outputName = "mu_tagprobe",
			  int         nEvents    = 1000)
{
  //------------------------------------------------------------------------------------------------
  // some parameters get passed through the environment
  //------------------------------------------------------------------------------------------------
  char json[1024], overlap[1024], path[1024];
  float overlapCut = -1;

  if (decodeEnv(json,overlap,overlapCut,path) != 0)
    return;

  TString jsonFile = TString("/home/cmsprod/cms/json/") + TString(json);
  Bool_t  isData   = (jsonFile.CompareTo("/home/cmsprod/cms/json/~") != 0);

  printf("\n Initialization worked: \n\n");
  printf("   JSON   : %s (file: %s)\n",  json,jsonFile.Data());
  printf("   OVERLAP: %s\n\n",overlap);
  printf("   PATH   : %s\n",  path);
  printf("   isData : %d\n\n",isData);

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
    printf("\n Jason file added: %s \n\n",jsonFile.Data());
    runLumiSel->AddJSONFile(jsonFile.Data());
  }
  if ((jsonFile.CompareTo("/home/cmsprod/cms/json/-") == 0)   ) {
    printf("\n WARNING -- Looking at data without JSON file: always accept.\n\n");
    runLumiSel->SetAbortIfNotAccepted(kFALSE);   // accept all events if there is no valid JSON file
  }
  printf("\n Run lumi worked. \n\n");

  //------------------------------------------------------------------------------------------------
  // the skimmer
  //------------------------------------------------------------------------------------------------
  H4lMuTagProbeSkim *skmMod = new H4lMuTagProbeSkim();

  //------------------------------------------------------------------------------------------------
  // making analysis chain
  //------------------------------------------------------------------------------------------------
  runLumiSel->Add(skmMod);

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
  TString skimDataset = TString(dataset)+TString("/") +TString(skim);
  Dataset *d = NULL;
  if (TString(skim).CompareTo("noskim") == 0)
    d = c->FindDataset(book,dataset,fileset);
  else
    d = c->FindDataset(book,skimDataset.Data(),fileset);
  ana->AddDataset(d);

  //------------------------------------------------------------------------------------------------
  // organize hist/ntuple output
  //------------------------------------------------------------------------------------------------
  ana->SetOutputName("test.root");
  ana->SetCacheSize(64*1024*1024);

  //------------------------------------------------------------------------------------------------
  // organize root output
  //------------------------------------------------------------------------------------------------
  OutputMod *outMod = new OutputMod;
  outMod->Drop("*");
  outMod->Keep("EventHeader");
  outMod->Keep("PileupInfo");
  outMod->Keep("HLTObjects");
  outMod->Keep("HLTObjectsRelation");
  outMod->Keep("HLTBits");
  outMod->Keep("Muons");
  outMod->Keep("Electrons");
  outMod->Keep("PrimaryVertexes");
  outMod->Keep("BeamSpot");
  outMod->Keep("AKt5PFJets");
  outMod->Keep("HLTBits");
  outMod->Keep("PFMet");
  outMod->Keep("Photons");
  outMod->Keep("MergedConversions");
  outMod->Keep("Rho");
  outMod->Keep("PFCandidates");
  outMod->Keep("Tracks");
  outMod->Keep("GlobalMuonTracks");
  outMod->Keep("StandaloneMuonTracksWVtxConstraint");
  outMod->Keep("GsfTracks");
  outMod->Keep("BarrelSuperClusters");
  outMod->Keep("EndcapSuperClusters");
  outMod->Keep("BarrelBasicClusters");
  outMod->Keep("MergedConversions_StableDatas");
  outMod->Keep("ElectronsStable");
  outMod->Keep("GsfElectronsStable");
  outMod->Keep("EndcapBasicClusters");
  outMod->Keep("PFSuperClusters");
  outMod->Keep("PFBasicClusters");
  outMod->Keep("ConversionOutInElectronsStable");
  outMod->Keep("ConversionOutInTracks");
  outMod->Keep("StandaloneMuonTracks");
  outMod->Keep("ConversionInOutElectronsStable");
  outMod->Keep("ConversionInOutTracks");
  if(!isData) {
    outMod->Keep("MCEventInfo");
    outMod->Keep("MCParticles");
  }

  TString rootFile = TString(outputName);
  rootFile += TString("_") + TString(dataset) + TString("_") + TString(skim);
  if (TString(fileset) != TString(""))
    rootFile += TString("_") + TString(fileset);
  outMod->SetFileName(rootFile);
  outMod->SetPathName(".");

  // Last step is the output module
  skmMod->Add(outMod);

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

int decodeEnv(char* json, char* overlap, float overlapCut, char* path)
{
  if (gSystem->Getenv("MIT_PROD_JSON"))
    sprintf(json,   "%s",gSystem->Getenv("MIT_PROD_JSON"));
  else {
    printf(" JSON file was not properly defined. EXIT!\n");
    return -1;
  }
  if (gSystem->Getenv("MIT_PROD_OVERLAP")) {
    sprintf(overlap,"%s",gSystem->Getenv("MIT_PROD_OVERLAP"));
    if (EOF == sscanf(overlap,"%f",&overlapCut)) {
      printf(" Overlap was not properly defined. EXIT!\n");
      return -1;
    }
  }
  else {
    printf(" OVERLAP file was not properly defined. EXIT!\n");
    return -1;
  }
  if (gSystem->Getenv("CMSSW_BASE"))
    sprintf(path,   "%s/src/MitPhysics/data/",gSystem->Getenv("CMSSW_BASE"));
  else {
    printf(" PATH was not properly defined. EXIT!\n");
    return -1;
  }

  return 0;
}
