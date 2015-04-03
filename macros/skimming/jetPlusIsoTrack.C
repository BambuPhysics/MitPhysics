#if !defined(__CLING__) || defined(__ROOTCLING__)
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/DataTree/interface/Names.h"
#include "MitAna/Catalog/interface/Catalog.h"
#include "MitAna/TreeMod/interface/Analysis.h"
#include "MitAna/TreeMod/interface/OutputMod.h"
#include "MitAna/PhysicsMod/interface/PublisherMod.h"
#include "MitAna/PhysicsMod/interface/PlotKineMod.h"
#include "MitPhysics/Mods/interface/MuonIDMod.h"
#include "MitPhysics/Mods/interface/MuonIDMod.h"
#include "MitPhysics/Mods/interface/ElectronIDMod.h"
#include "MitPhysics/Mods/interface/ElectronCleaningMod.h"
#include "MitPhysics/Mods/interface/JetIDMod.h"
#include "MitPhysics/Mods/interface/JetCleaningMod.h"
#include "MitPhysics/SelMods/interface/JetPlusIsoTrackSelMod.h"
#endif

//--------------------------------------------------------------------------------------------------
void jetPlusIsoTrack(const char *fileset    = "0000",
                     const char *dataset    = "s8-wm-id9",
                     const char *book       = "mit/filler/006",
                     const char *catalogDir = "/home/mitprod/catalog",
                     Int_t       nEvents    = -1)
{
  TString skimName("jetPlusIsoTrack");
  using namespace mithep;
  gDebugMask  = Debug::kAnalysis;
  gDebugLevel = 1;

  //------------------------------------------------------------------------------------------------
  // organize selection
  //------------------------------------------------------------------------------------------------
  //const char     *jetInput   = Names::gkSC5JetBrn;
  //const char     *jetInput   = Names::gkPFJetBrn;
  const char     *jetInput   = "AKt5PFJetsCHS";
  const char     *gsfTracks  = "GsfTracks";
  const Double_t  jetPtMin   = 30;
  const Double_t  trackPtMin = 10;


  PublisherMod<PFJet,Jet> *pubJet = new PublisherMod<PFJet,Jet>("JetPub");
  //pubJet->SetInputName(Names::gkPFJetBrn);
  //pubJet->SetInputName("AKt5PFJetsCHS");
  pubJet->SetInputName("AKt5PFJetsCHS");
  pubJet->SetOutputName("PubAKt5PFJetsCHS");

  JetIDMod *jetId = new JetIDMod;  
  jetId->SetInputName(pubJet->GetOutputName());
  jetId->SetUseCorrection(kFALSE);
  jetId->SetPtCut(jetPtMin); 


  JetPlusIsoTrackSelMod *selMod = new JetPlusIsoTrackSelMod;
  selMod->SetTrackPtMin(trackPtMin);
  selMod->SetJetColName(jetId->GetOutputName());
  selMod->SetTrackerTrackColName(Names::gkTrackBrn);
  selMod->SetGsfTrackColName(gsfTracks);

  //------------------------------------------------------------------------------------------------
  // link modules together
  //------------------------------------------------------------------------------------------------
  pubJet->Add(jetId);

  //------------------------------------------------------------------------------------------------
  // organize output
  //------------------------------------------------------------------------------------------------
  OutputMod *outMod = new OutputMod;
  outMod->Keep("*");
  //selMod->Add(outMod);
  TString rootFile = skimName;
  if (TString(fileset) != TString(""))
    rootFile += TString("_") + TString(fileset);
  printf("\nRoot output: %s\n\n",rootFile.Data());  
  outMod->SetFileName(rootFile);
  outMod->SetPathName(".");

  //------------------------------------------------------------------------------------------------
  // set up analysis
  //------------------------------------------------------------------------------------------------
  Analysis *ana = new Analysis;
  ana->AddSuperModule(pubJet);
  //ana->AddSuperModule(jetId);
  //ana->AddSuperModule(selMod);
  if (nEvents>0)
    ana->SetProcessNEvents(nEvents);

  //------------------------------------------------------------------------------------------------
  // organize input
  //------------------------------------------------------------------------------------------------
  printf("\nRely on Catalog: %s\n",catalogDir);
  printf("  -> Book: %s  Dataset: %s  Fileset: %s <-\n\n",book,dataset,fileset);
  //Catalog *c = new Catalog(catalogDir);
  //Dataset *d = c->FindDataset(book,dataset,fileset);
  //ana->AddDataset(d);
  ana->AddFile("/mnt/hadoop/cmsprod/XX-MITDATASET-XX_000.root");
  //ana->AddFile("root://xrootd.cmsaf.mit.edu//store/user/paus/filefi/032/s12-ttj-v1-v7a/047D197D-0EE2-E111-9C87-0030487E4ED3.root");

  //------------------------------------------------------------------------------------------------
  // run the analysis after successful initialisation
  //------------------------------------------------------------------------------------------------
  ana->Run(kFALSE);
}
