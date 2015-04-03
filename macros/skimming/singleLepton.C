// $Id: singleLepton.C,v 1.3 2009/04/30 12:13:37 ceballos Exp $

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/DataTree/interface/Names.h"
#include "MitAna/Catalog/interface/Catalog.h"
#include "MitAna/TreeMod/interface/Analysis.h"
#include "MitAna/TreeMod/interface/OutputMod.h"
#include "MitAna/PhysicsMod/interface/PlotKineMod.h"
#include "MitPhysics/Mods/interface/MuonIDMod.h"
#include "MitPhysics/Mods/interface/MuonIDMod.h"
#include "MitPhysics/Mods/interface/ElectronIDMod.h"
#include "MitPhysics/Mods/interface/ElectronCleaningMod.h"
#include "MitPhysics/Mods/interface/MergeLeptonsMod.h"
#include "MitPhysics/SelMods/interface/GenericSelMod.h"
#endif

//--------------------------------------------------------------------------------------------------
void singleLepton(const char *fileset    = "0000",
		  const char *dataset    = "p10-wjets-mg-v26",
		  const char *book       = "cern/filler/014a",
		  const char *catalogDir = "/home/ceballos/catalog",
		  Int_t       nEvents    = 1000)
{
  TString skimName("singlelepton");
  using namespace mithep;
  gDebugMask  = Debug::kAnalysis;
  gDebugLevel = 1;

  //------------------------------------------------------------------------------------------------
  // organize selection
  //------------------------------------------------------------------------------------------------
  MuonIDMod *muId = new MuonIDMod;  
  muId->SetPtMin     (17.0);
  muId->SetApplyD0Cut(kFALSE);
  muId->SetClassType ("Global");
  muId->SetIDType    ("NoId");
  muId->SetIsoType   ("NoIso");

  ElectronIDMod       *electronIDWP95       = new ElectronIDMod;
  electronIDWP95->SetPtMin(0.0);
  electronIDWP95->SetEtaMax(999.);
  electronIDWP95->SetEtMin(17.0);
  electronIDWP95->SetIDType(TString("VBTFWorkingPoint95Id"));
  electronIDWP95->SetIsoType(TString("VBTFWorkingPoint95Iso"));
  electronIDWP95->SetNExpectedHitsInnerCut(1);
  electronIDWP95->SetApplyConversionFilterType1(kFALSE);
  electronIDWP95->SetApplyConversionFilterType2(kFALSE);
  electronIDWP95->SetApplyD0Cut(kFALSE);
  electronIDWP95->SetChargeFilter(kFALSE);
  electronIDWP95->SetApplyTriggerMatching(kFALSE);
  electronIDWP95->SetApplyEcalFiducial(kFALSE);
  electronIDWP95->SetApplyEcalSeeded(kFALSE);
  electronIDWP95->SetApplyCombinedIso(kFALSE);
  electronIDWP95->SetApplySpikeRemoval(kFALSE);
  electronIDWP95->SetOutputName("GoodElectronsWP95");

  ElectronIDMod       *electronIDWP80Id       = new ElectronIDMod;
  electronIDWP80Id->SetPtMin(0.0);
  electronIDWP80Id->SetEtaMax(999.);
  electronIDWP80Id->SetEtMin(17.0);
  electronIDWP80Id->SetIDType(TString("VBTFWorkingPoint80Id"));
  electronIDWP80Id->SetIsoType(TString("NoIso"));
  electronIDWP80Id->SetNExpectedHitsInnerCut(1);
  electronIDWP80Id->SetApplyConversionFilterType1(kFALSE);
  electronIDWP80Id->SetApplyConversionFilterType2(kTRUE);
  electronIDWP80Id->SetApplyD0Cut(kFALSE);
  electronIDWP80Id->SetChargeFilter(kFALSE);
  electronIDWP80Id->SetApplyTriggerMatching(kFALSE);
  electronIDWP80Id->SetApplyEcalFiducial(kFALSE);
  electronIDWP80Id->SetApplyEcalSeeded(kFALSE);
  electronIDWP80Id->SetApplyCombinedIso(kFALSE);
  electronIDWP80Id->SetApplySpikeRemoval(kFALSE);
  electronIDWP80Id->SetOutputName("GoodElectronsWP80Id");

  ElectronIDMod       *electronIDWP80Iso       = new ElectronIDMod;
  electronIDWP80Iso->SetPtMin(0.0);
  electronIDWP80Iso->SetEtaMax(999.);
  electronIDWP80Iso->SetEtMin(17.0);
  electronIDWP80Iso->SetIDType(TString("NoId"));
  electronIDWP80Iso->SetIsoType(TString("VBTFWorkingPoint80Iso"));
  electronIDWP80Iso->SetNExpectedHitsInnerCut(999);
  electronIDWP80Iso->SetApplyConversionFilterType1(kFALSE);
  electronIDWP80Iso->SetApplyConversionFilterType2(kFALSE);
  electronIDWP80Iso->SetApplyD0Cut(kFALSE);
  electronIDWP80Iso->SetChargeFilter(kFALSE);
  electronIDWP80Iso->SetApplyTriggerMatching(kFALSE);
  electronIDWP80Iso->SetApplyEcalFiducial(kFALSE);
  electronIDWP80Iso->SetApplyEcalSeeded(kFALSE);
  electronIDWP80Iso->SetApplyCombinedIso(kFALSE);
  electronIDWP80Iso->SetApplySpikeRemoval(kFALSE);
  electronIDWP80Iso->SetOutputName("GoodElectronsWP80Iso");

  GenericSelMod<mithep::Muon> *selMod0 = new GenericSelMod<mithep::Muon>;
  selMod0->SetPtMin(0.0);
  selMod0->SetMinCounts(1);
  selMod0->SetColName(muId->GetOutputName());

  GenericSelMod<mithep::Electron> *selMod1 = new GenericSelMod<mithep::Electron>;
  selMod1->SetPtMin(0.0);
  selMod1->SetMinCounts(1);
  selMod1->SetColName(electronIDWP95->GetOutputName());

  GenericSelMod<mithep::Electron> *selMod2 = new GenericSelMod<mithep::Electron>;
  selMod2->SetPtMin(0.0);
  selMod2->SetMinCounts(1);
  selMod2->SetColName(electronIDWP80Id->GetOutputName());

  GenericSelMod<mithep::Electron> *selMod3 = new GenericSelMod<mithep::Electron>;
  selMod3->SetPtMin(0.0);
  selMod3->SetMinCounts(1);
  selMod3->SetColName(electronIDWP80Iso->GetOutputName());

  //------------------------------------------------------------------------------------------------
  // organize output
  //------------------------------------------------------------------------------------------------
  OutputMod *outMod = new OutputMod;
  outMod->Keep("*");

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
  ana->AddSuperModule(electronIDWP95);
  ana->AddSuperModule(electronIDWP80Id);
  ana->AddSuperModule(electronIDWP80Iso);
  muId->Add(selMod0);
  electronIDWP95->Add(selMod1);
  electronIDWP80Id->Add(selMod2);
  electronIDWP80Iso->Add(selMod3);
  selMod0->Add(outMod);
  selMod1->Add(outMod);
  selMod2->Add(outMod);
  selMod3->Add(outMod);
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
