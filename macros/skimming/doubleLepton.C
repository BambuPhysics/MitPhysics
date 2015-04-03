// $Id: doubleLepton.C,v 1.11 2012/08/10 19:09:37 ceballos Exp $

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/DataTree/interface/Names.h"
#include "MitAna/Catalog/interface/Catalog.h"
#include "MitAna/TreeMod/interface/Analysis.h"
#include "MitAna/TreeMod/interface/OutputMod.h"
#include "MitPhysics/Mods/interface/MuonIDMod.h"
#include "MitPhysics/Mods/interface/ElectronIDMod.h"
#include "MitPhysics/Mods/interface/ElectronCleaningMod.h"
#include "MitPhysics/Mods/interface/MergeLeptonsMod.h"
#include "MitPhysics/Mods/interface/GoodPVFilterMod.h"
#include "MitPhysics/SelMods/interface/DilepSelMod.h"
#endif

//--------------------------------------------------------------------------------------------------
void doubleLepton(const char *catalogDir   = "~/scratch0/catalog",
                  const char *book	   = "cern/filefi/028",
                  const char *dataset	   = "r12a-mueg-pr-v1",
                  const char *fileset	   = "0003",
                  int NEvents = 1000)
{
  TString skimName("doubleLepton");
  using namespace mithep;
  gDebugMask  = Debug::kAnalysis;
  gDebugLevel = 1;

  //------------------------------------------------------------------------------------------------
  // organize selection
  //------------------------------------------------------------------------------------------------
  GoodPVFilterMod *goodPVFilterMod = new GoodPVFilterMod;
  goodPVFilterMod->SetMinVertexNTracks(0);
  goodPVFilterMod->SetMinNDof(4);
  goodPVFilterMod->SetMaxAbsZ(24.0);
  goodPVFilterMod->SetMaxRho(2.0);
  goodPVFilterMod->SetVertexesName("PrimaryVertexes");

  const char     *muInput  = Names::gkMuonBrn;
  const char     *elInput  = Names::gkElectronBrn;

  MuonIDMod *muId = new MuonIDMod;  
  muId->SetInputName  (muInput);
  muId->SetPtMin      (10.0);
  muId->SetApplyD0Cut (kTRUE);
  muId->SetD0Cut      (0.20);
  muId->SetApplyDZCut (kTRUE);
  muId->SetDZCut      (0.50);
  muId->SetWhichVertex(0);
  muId->SetClassType  ("All");
  muId->SetIDType     ("NoId");
  muId->SetIsoType    ("NoIso");

  ElectronIDMod *elId = new ElectronIDMod;
  elId->SetInputName                 (elInput);
  elId->SetPtMin                     (10.0);
  elId->SetApplyConversionFilterType1(kFALSE);
  elId->SetApplyConversionFilterType2(kFALSE);
  elId->SetChargeFilter              (kFALSE);
  elId->SetApplySpikeRemoval         (kFALSE);
  elId->SetApplyD0Cut                (kTRUE);
  elId->SetD0Cut                     (0.20);
  elId->SetApplyDZCut                (kTRUE);
  elId->SetDZCut                     (0.50);
  elId->SetWhichVertex               (0);
  elId->SetNExpectedHitsInnerCut     (999);
  elId->SetIDType                    ("VBTFWorkingPointFakeableId");
  elId->SetIsoType                   ("NoIso");

  ElectronCleaningMod *elCl = new ElectronCleaningMod;
  elCl->SetGoodElectronsName(elId->GetOutputName());
  elCl->SetCleanMuonsName   (muId->GetOutputName());
  elId->Add(elCl);

  MergeLeptonsMod *merger = new MergeLeptonsMod;
  merger->SetMuonsName    (muId->GetOutputName());
  merger->SetElectronsName(elCl->GetOutputName());

  DilepSelMod *selMod = new DilepSelMod;
  selMod->SetCleanLeptonsName(merger->GetOutputName());
  selMod->SetMinPt(0.0);
  selMod->SetMinDilMass(0.1);
  selMod->SetFillHist(kFALSE);
  merger->Add(selMod);

  //------------------------------------------------------------------------------------------------
  // organize output
  //------------------------------------------------------------------------------------------------
  OutputMod *outMod = new OutputMod;
  outMod->Keep("*");
  outMod->Drop("CorMuonMet");
  selMod->Add(outMod);
  TString rootFile = "";
  rootFile += TString("skims/") + dataset + TString("/");
  rootFile += skimName;
  if (TString(fileset) != TString(""))
    rootFile += TString("_") + TString(fileset);
  printf("\nRoot output: %s\n\n",rootFile.Data());  
  outMod->SetFileName(rootFile);
  outMod->SetPathName(".");

  //------------------------------------------------------------------------------------------------
  // set up analysis
  //------------------------------------------------------------------------------------------------
  Analysis *ana = new Analysis;
  ana->AddSuperModule(goodPVFilterMod);
  goodPVFilterMod->Add(muId);
  muId->Add(elId);
  elId->Add(merger);
  if (NEvents>0)
    ana->SetProcessNEvents(NEvents);

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
