// $Id: leptonPlusIsoTrack.C,v 1.3 2009/03/22 09:04:13 loizides Exp $

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/DataTree/interface/Names.h"
#include "MitAna/Catalog/interface/Catalog.h"
#include "MitAna/TreeMod/interface/Analysis.h"
#include "MitAna/TreeMod/interface/OutputMod.h"
#include "MitAna/PhysicsMod//interface/PlotKineMod.h"
#include "MitPhysics/Mods/interface/MuonIDMod.h"
#include "MitPhysics/Mods/interface/MuonIDMod.h"
#include "MitPhysics/Mods/interface/ElectronIDMod.h"
#include "MitPhysics/Mods/interface/ElectronCleaningMod.h"
#include "MitPhysics/Mods/interface/MergeLeptonsMod.h"
#include "MitPhysics/SelMods/interface/LeptonPlusIsoTrackSelMod.h"
#endif

//--------------------------------------------------------------------------------------------------
void leptonPlusIsoTrack(const char *fileset    = "",
                        const char *dataset    = "s8-wm-id9",
                        const char *book       = "mit/filler/006",
                        const char *catalogDir = "/home/mitprod/catalog",
                        Int_t       nEvents    = -1)
{
  TString skimName("leptonPlusIsoTrack");
  using namespace mithep;
  gDebugMask  = Debug::kAnalysis;
  gDebugLevel = 1;

  //------------------------------------------------------------------------------------------------
  // organize selection
  //------------------------------------------------------------------------------------------------
  const char     *muInput   = Names::gkMuonBrn;
  const char     *elInput   = Names::gkElectronBrn;
  const char     *gsfTracks = "GsfTracks";
  const Double_t  ptMin     = 10;

  MuonIDMod *muId = new MuonIDMod;  
  muId->SetInputName  (muInput);
  muId->SetIDType     ("Loose");
  muId->SetIsoType    ("TrackCalo");
  muId->SetTrackIsoCut( 10.0);
  muId->SetCaloIsoCut (100.0);
  muId->SetPtMin      (ptMin);

  ElectronIDMod *elId = new ElectronIDMod;
  elId->SetInputName  (elInput);
  elId->SetIDType     ("NoId");
  elId->SetIsoType    ("TrackCalo");
  elId->SetTrackIsoCut( 10.0);
  elId->SetCaloIsoCut (100.0);
  elId->SetPtMin      (ptMin);

  ElectronCleaningMod *elCl = new ElectronCleaningMod;
  elCl->SetGoodElectronsName(elId->GetOutputName());
  elCl->SetCleanMuonsName   (muId->GetOutputName());
  elId->Add(elCl);

  MergeLeptonsMod *merger = new MergeLeptonsMod;
  merger->SetMuonsName    (muId->GetOutputName());
  merger->SetElectronsName(elCl->GetOutputName());

  LeptonPlusIsoTrackSelMod *selMod = new LeptonPlusIsoTrackSelMod;
  selMod->SetLeptonPtMin(ptMin);
  selMod->SetTrackPtMin(ptMin);
  selMod->SetLeptonColName(merger->GetOutputName());
  selMod->SetTrackerTrackColName(Names::gkTrackBrn);
  selMod->SetGsfTrackColName(gsfTracks);
  merger->Add(selMod);

  //------------------------------------------------------------------------------------------------
  // organize output
  //------------------------------------------------------------------------------------------------
  OutputMod *outMod = new OutputMod;
  outMod->Keep("*");
  selMod->Add(outMod);
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
  ana->AddSuperModule(muId);
  ana->AddSuperModule(elId);
  ana->AddSuperModule(merger);
  if (nEvents>0)
    ana->SetProcessNEvents(nEvents);

  //------------------------------------------------------------------------------------------------
  // organize input
  //------------------------------------------------------------------------------------------------
  printf("\nRely on Catalog: %s\n",catalogDir);
  printf("  -> Book: %s  Dataset: %s  Fileset: %s <-\n\n",book,dataset,fileset);
  Catalog *c = new Catalog(catalogDir);
  Dataset *d = c->FindDataset(book,dataset,fileset);
  ana->AddDataset(d);

  //------------------------------------------------------------------------------------------------
  // run the analysis after successful initialisation
  //------------------------------------------------------------------------------------------------
  ana->Run(kFALSE);
}
