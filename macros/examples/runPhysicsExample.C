//root -l -q -b $CMSSW_BASE/src/MitHiggs/macros/runMacros/runHwwExampleAnalysis.C+\(\"0000\",\"noskim\",\"s8-h190ww2l-gf-mc3\",\"mit/filler/011\",\"/home/mitprod/catalog\",\"HwwExampleAnalysis\",1000,1\)

// $Id: runPhysicsExample.C,v 1.16 2011/04/29 16:18:16 ceballos Exp $

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <TROOT.h>
#include <exception>
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/Catalog/interface/Catalog.h"
#include "MitAna/TreeMod/interface/Analysis.h"
#include "MitAna/TreeMod/interface/HLTMod.h"
#include "MitAna/PhysicsMod/interface/PublisherMod.h"
#include "MitAna/PhysicsMod/interface/RunLumiSelectionMod.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Mods/interface/GoodPVFilterMod.h"
#include "MitPhysics/Mods/interface/GeneratorMod.h"
#include "MitPhysics/Mods/interface/HKFactorProducer.h"
#include "MitPhysics/Mods/interface/JetCorrectionMod.h"
#include "MitPhysics/Mods/interface/CaloMetCorrectionMod.h"
#include "MitPhysics/Mods/interface/MuonIDMod.h"
#include "MitPhysics/Mods/interface/ElectronIDMod.h"
#include "MitPhysics/Mods/interface/PhotonIDMod.h"
#include "MitPhysics/Mods/interface/TauIDMod.h"
#include "MitPhysics/Mods/interface/JetIDMod.h"
#include "MitPhysics/Mods/interface/ElectronCleaningMod.h"
#include "MitPhysics/Mods/interface/PhotonCleaningMod.h"
#include "MitPhysics/Mods/interface/TauCleaningMod.h"
#include "MitPhysics/Mods/interface/JetCleaningMod.h"
#include "MitPhysics/Mods/interface/MergeLeptonsMod.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitAna/DataTree/interface/MetCol.h" 
#include "MitAna/DataTree/interface/PFMetCol.h" 
#include "MitAna/DataTree/interface/CaloMetCol.h"
#include "MitPhysics/SelMods/interface/HwwExampleAnalysisMod.h"
#include "MitPhysics/SelMods/interface/WBFExampleAnalysisMod.h"
#endif

//--------------------------------------------------------------------------------------------------
void runPhysicsExample(const char *catalogDir = "/home/mitprod/catalog",
		       const char *book	      = "cern/filefi/020",
                       const char *dataset    = "p11-ww2l-v1g1-pu",
                       const char *fileset    = "0000",
                       const char *skim       = "noskim",
                       const char *outputName = "histo",
                       int   sampleID	      = -1,
                       int   nEvents	      = 10000)
{
  //------------------------------------------------------------------------------------------------
  // some global setups
  //------------------------------------------------------------------------------------------------
  using namespace mithep;
  gDebugMask  = Debug::kAnalysis;
  gDebugLevel = 3;

  //------------------------------------------------------------------------------------------------
  // set up information
  //------------------------------------------------------------------------------------------------
  Bool_t applyISRFilter     = kFALSE;
  Bool_t applyMllGenCut     = kFALSE;
  Bool_t isData             = kFALSE;
  Bool_t isDataMuonElectron = kFALSE;
  Bool_t isDataDMuon        = kFALSE;
  Bool_t isDataSMuon        = kFALSE;
  Bool_t isDataDElectron    = kFALSE;
  Bool_t isDataSElectron    = kFALSE;
  int processId         = -999999999; // use 999 for MCatNLO MC sample, 102 for H->WW
  TString fInputFilenameKF = "/home/ceballos/releases/CMSSW_4_2_2/src/MitPhysics/data/HWW_KFactors_160_10TeV.dat";

  if(sampleID >= 1000) isData             = kTRUE;

  if(sampleID >= 1000) isDataMuonElectron = kTRUE;
  if(sampleID >= 2000) isDataDMuon	  = kTRUE;
  if(sampleID >= 3000) isDataSMuon	  = kTRUE;
  if(sampleID >= 4000) isDataDElectron	  = kTRUE;
  if(sampleID >= 5000) isDataSElectron	  = kTRUE;

  //------------------------------------------------------------------------------------------------
  // generator information
  //------------------------------------------------------------------------------------------------
  GeneratorMod *generatorMod = new GeneratorMod;
  generatorMod->SetPrintDebug(kFALSE);
  generatorMod->SetPtLeptonMin(0.0);
  generatorMod->SetEtaLeptonMax(2.7);
  generatorMod->SetPtPhotonMin(15.0);
  generatorMod->SetEtaPhotonMax(2.7);
  generatorMod->SetPtRadPhotonMin(10.0);
  generatorMod->SetEtaRadPhotonMax(2.7);
  generatorMod->SetIsData(isData);
  generatorMod->SetFillHist(!isData);
  if(applyMllGenCut == kTRUE){
    generatorMod->SetPdgIdCut(23);
    generatorMod->SetMassMinCut( 0.);
    generatorMod->SetMassMaxCut(50.);
  }
  generatorMod->SetApplyISRFilter(applyISRFilter);

  HKFactorProducer *hKFactorProducer = new HKFactorProducer;
  hKFactorProducer->SetProcessID(processId);
  hKFactorProducer->SetInputFilename(fInputFilenameKF);
  hKFactorProducer->SetIsData(isData);
  hKFactorProducer->SetFillHist(!isData);

  //------------------------------------------------------------------------------------------------
  // Run RunLumiSelectionMod
  //------------------------------------------------------------------------------------------------
  RunLumiSelectionMod *runLumiSelectionMod = new RunLumiSelectionMod;
  runLumiSelectionMod->SetAcceptMC(!isData);    
  runLumiSelectionMod->AddJSONFile("/home/ceballos/releases/CMSSW_4_2_2/src/json/Cert_160404-166502_7TeV_PromptReco_Collisions11_JSON.txt");

  //------------------------------------------------------------------------------------------------
  // PV filter selection
  //------------------------------------------------------------------------------------------------
  GoodPVFilterMod *goodPVFilterMod = new GoodPVFilterMod;
  goodPVFilterMod->SetMinVertexNTracks(0);
  goodPVFilterMod->SetMinNDof(4);
  goodPVFilterMod->SetMaxAbsZ(24.0);
  goodPVFilterMod->SetMaxRho(2.0);
  goodPVFilterMod->SetVertexesName("DAPrimaryVertexes");

  //------------------------------------------------------------------------------------------------
  // HLT information
  //------------------------------------------------------------------------------------------------
  HLTMod *hltmod = new HLTMod;

  if     (isData == kFALSE){
    hltmod->AddTrigger("HLT_Mu15_v1");
    hltmod->AddTrigger("!HLT_Mu15_v1");
  }
  else if(isData == true && isDataMuonElectron == true) {
    hltmod->AddTrigger("HLT_Mu17_Ele8_CaloIdL_v1",150000,161176);
    hltmod->AddTrigger("HLT_Mu8_Ele17_CaloIdL_v1",150000,161176);
    hltmod->AddTrigger("HLT_Mu17_Ele8_CaloIdL_v2",161179,163261);
    hltmod->AddTrigger("HLT_Mu8_Ele17_CaloIdL_v2",161179,163261);
    hltmod->AddTrigger("HLT_Mu17_Ele8_CaloIdL_v3",163262,999999);
    hltmod->AddTrigger("HLT_Mu8_Ele17_CaloIdL_v3",163262,999999);
  }
  else if(isData == true && isDataDMuon == true) {
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v1&!HLT_Mu17_Ele8_CaloIdL_v1&HLT_DoubleMu7_v1",150000,161176);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v2&!HLT_Mu17_Ele8_CaloIdL_v2&HLT_DoubleMu7_v1",161179,163261);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v3&!HLT_Mu17_Ele8_CaloIdL_v3&HLT_DoubleMu7_v2",163262,999999);
  }
  else if(isData == true && isDataSMuon == true) {
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v1&!HLT_Mu17_Ele8_CaloIdL_v1&!HLT_DoubleMu7_v1&HLT_Mu15_v2",150000,161176);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v2&!HLT_Mu17_Ele8_CaloIdL_v2&!HLT_DoubleMu7_v1&HLT_Mu15_v2",161179,163261);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v3&!HLT_Mu17_Ele8_CaloIdL_v3&!HLT_DoubleMu7_v2&HLT_Mu24_v2",163262,999999);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v3&!HLT_Mu17_Ele8_CaloIdL_v3&!HLT_DoubleMu7_v2&HLT_IsoMu17_v6",163262,999999);
  }
  else if(isData == true && isDataDElectron == true) {
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v1&!HLT_Mu17_Ele8_CaloIdL_v1&!HLT_DoubleMu7_v1&!HLT_Mu15_v2&HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v1",150000,161176);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v2&!HLT_Mu17_Ele8_CaloIdL_v2&!HLT_DoubleMu7_v1&!HLT_Mu15_v2&HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2",161179,163261);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v3&!HLT_Mu17_Ele8_CaloIdL_v3&!HLT_DoubleMu7_v2&!HLT_Mu24_v2&!HLT_IsoMu17_v6&HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v3",163262,999999);
  }
  else if(isData == true && isDataSElectron == true) {
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v1&!HLT_Mu17_Ele8_CaloIdL_v1&!HLT_DoubleMu7_v1&!HLT_Mu15_v2&!HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v1&HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1",150000,161176);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v2&!HLT_Mu17_Ele8_CaloIdL_v2&!HLT_DoubleMu7_v1&!HLT_Mu15_v2&!HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2&HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2",161179,163261);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v3&!HLT_Mu17_Ele8_CaloIdL_v3&!HLT_DoubleMu7_v2&!HLT_Mu24_v2&!HLT_IsoMu17_v6&!HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v3&HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3",163262,999999);
  }
  hltmod->SetTrigObjsName("myhltobjs");

  //------------------------------------------------------------------------------------------------
  // publisher Mod
  //------------------------------------------------------------------------------------------------
  PublisherMod<PFJet,Jet> *pubJet = new PublisherMod<PFJet,Jet>("JetPub");
  pubJet->SetInputName("AKt5PFJets");
  pubJet->SetOutputName("PubAKt5PFJets");

  PublisherMod<PFMet,Met> *pubMet = new PublisherMod<PFMet,Met>("MetPub");
  pubMet->SetInputName("PFMet");
  pubMet->SetOutputName("PubPFMet");

  //------------------------------------------------------------------------------------------------
  // Apply Jet Corrections
  //------------------------------------------------------------------------------------------------
  JetCorrectionMod *jetCorr = new JetCorrectionMod;
  jetCorr->AddCorrectionFromFile("/home/ceballos/releases/CMSSW_4_2_2/src/MitPhysics/data/START41_V0_AK5PF_L1FastJet.txt"); 
  jetCorr->AddCorrectionFromFile("/home/ceballos/releases/CMSSW_4_2_2/src/MitPhysics/data/START41_V0_AK5PF_L2Relative.txt"); 
  jetCorr->AddCorrectionFromFile("/home/ceballos/releases/CMSSW_4_2_2/src/MitPhysics/data/START41_V0_AK5PF_L3Absolute.txt");
  if(isData == true){
    jetCorr->AddCorrectionFromFile("/home/ceballos/releases/CMSSW_4_2_2/src/MitPhysics/data/START41_V0_AK5PF_L2L3Residual.txt");  
  }
  jetCorr->SetInputName(pubJet1->GetOutputName());
  jetCorr->SetCorrectedName("CorrectedJets");

  //------------------------------------------------------------------------------------------------
  // object id and cleaning sequence
  //------------------------------------------------------------------------------------------------
  MuonIDMod *muonID = new MuonIDMod;
  muonID->SetClassType("GlobalTracker");
  muonID->SetIDType("WWMuIdV2");
  muonID->SetIsoType("PFIso");
  muonID->SetApplyD0Cut(kTRUE);
  muonID->SetApplyDZCut(kTRUE);
  muonID->SetWhichVertex(0);

  ElectronIDMod *electronID = new ElectronIDMod;
  electronID->SetIDType("VBTFWorkingPointLowPtId");
  electronID->SetIsoType("PFIso");
  electronID->SetApplyConversionFilterType1(kTRUE);
  electronID->SetApplyConversionFilterType2(kFALSE);
  electronID->SetChargeFilter(kFALSE);
  electronID->SetApplyD0Cut(kTRUE);
  electronID->SetApplyDZCut(kTRUE);
  electronID->SetWhichVertex(0);
  electronID->SetNExpectedHitsInnerCut(0);

  // Object ID and Cleaning Sequence
  PhotonIDMod         *photonID      = new PhotonIDMod;
  TauIDMod            *tauID         = new TauIDMod;
  JetIDMod            *jetID         = new JetIDMod;
  jetID->SetInputName(jetCorr->GetOutputName());
  jetID->SetPtCut(30.0);
  jetID->SetEtaMaxCut(5.0);
  jetID->SetJetEEMFractionMinCut(0.0);
  jetID->SetOutputName("GoodJets");
  jetID->SetApplyBetaCut(kFALSE);

  ElectronCleaningMod *electronCleaning = new ElectronCleaningMod;
  PhotonCleaningMod   *photonCleaning   = new PhotonCleaningMod;
  TauCleaningMod      *tauCleaning      = new TauCleaningMod;
  JetCleaningMod      *jetCleaning      = new JetCleaningMod;
  jetCleaning->SetGoodJetsName("GoodJets");
  jetCleaning->SetCleanJetsName("CleanJets");

  JetIDMod            *jetIDNoPtCut     = new JetIDMod;
  jetIDNoPtCut->SetInputName(jetCorr->GetOutputName());
  jetIDNoPtCut->SetPtCut(0.0);
  jetIDNoPtCut->SetEtaMaxCut(5.0);
  jetIDNoPtCut->SetJetEEMFractionMinCut(0.0);
  jetIDNoPtCut->SetOutputName("GoodJetsNoPtCut");
  jetIDNoPtCut->SetApplyBetaCut(kFALSE);

  JetCleaningMod      *jetCleaningNoPtCut = new JetCleaningMod;
  jetCleaningNoPtCut->SetGoodJetsName("GoodJetsNoPtCut");
  jetCleaningNoPtCut->SetCleanJetsName("CleanJetsNoPtCut");

  //------------------------------------------------------------------------------------------------
  // merge modules
  //------------------------------------------------------------------------------------------------
  MergeLeptonsMod *mergeLeptonsMod = new MergeLeptonsMod;
  mergeLeptonsMod->SetMuonsName(muonID->GetOutputName());
  mergeLeptonsMod->SetElectronsName(electronCleaning->GetOutputName());

  //------------------------------------------------------------------------------------------------
  // analyses modules
  //------------------------------------------------------------------------------------------------
  HwwExampleAnalysisMod *WWanalysisMod = new HwwExampleAnalysisMod;
  WWanalysisMod->SetMetName(pubMet->GetOutputName());
  WWanalysisMod->SetCleanJetsName(jetCleaning->GetOutputName());
  WWanalysisMod->SetCleanJetsNoPtCutName(jetCleaningNoPtCut->GetOutputName());

  WBFExampleAnalysisMod *WBFanalysisMod = new WBFExampleAnalysisMod;
  WBFanalysisMod->SetMetName(pubMet->GetOutputName());
  WBFanalysisMod->SetCleanJetsName(jetCleaning->GetOutputName());

  //------------------------------------------------------------------------------------------------
  // making analysis chain
  //------------------------------------------------------------------------------------------------
  generatorMod->Add(hKFactorProducer);
  hKFactorProducer->Add(runLumiSelectionMod);
  runLumiSelectionMod->Add(goodPVFilterMod);
  goodPVFilterMod->Add(hltmod);
  hltmod->Add(muonID);
  muonID->Add(electronID);
  electronID->Add(photonID);
  photonID->Add(tauID);
  tauID->Add(pubJet);
  pubJet->Add(pubMet); 
  pubMet->Add(jetCorr);
  jetCorr->Add(jetID);
  jetID->Add(electronCleaning);
  electronCleaning->Add(photonCleaning);
  photonCleaning->Add(tauCleaning);
  tauCleaning->Add(jetCleaning);
  jetCleaning->Add(jetIDNoPtCut);
  jetIDNoPtCut->Add(jetCleaningNoPtCut);
  jetCleaningNoPtCut->Add(mergeLeptonsMod);
  mergeLeptonsMod->Add(WWanalysisMod);
  WWanalysisMod->Add(WBFanalysisMod);

  //------------------------------------------------------------------------------------------------
  // setup analysis
  //------------------------------------------------------------------------------------------------
  Analysis *ana = new Analysis;
  ana->SetUseHLT(kTRUE);
  ana->SetKeepHierarchy(kFALSE);
  if (nEvents >= 0)
    ana->SetProcessNEvents(nEvents);
  ana->SetSuperModule(generatorMod);
  ana->SetPrintScale(100);

  //------------------------------------------------------------------------------------------------
  // organize input
  //------------------------------------------------------------------------------------------------
  printf("\nRely on Catalog: %s\n",catalogDir);
  printf("  -> Book: %s  Dataset: %s  Skim: %s  Fileset: %s <-\n\n",book,dataset,skim,fileset);
  Catalog *c = new Catalog(catalogDir);
  TString skimdataset = TString(dataset)+TString("/") +TString(skim);
  Dataset *d = NULL;
  if (TString(skim).CompareTo("noskim") == 0)
    d = c->FindDataset(book,dataset,fileset);
  else 
    d = c->FindDataset(book,skimdataset.Data(),fileset);
  ana->AddDataset(d);

  //------------------------------------------------------------------------------------------------
  // organize output
  //------------------------------------------------------------------------------------------------
  TString rootFile = TString(outputName);
  rootFile += TString("_") + TString(dataset) + TString("_") + TString(skim);
  if (TString(fileset) != TString(""))
    rootFile += TString("_") + TString(fileset);
  rootFile += TString(".root");
  printf("\nRoot output: %s\n\n",rootFile.Data());  
  ana->SetOutputName(rootFile.Data());
  ana->SetCacheSize(64*1024*1024);

  //------------------------------------------------------------------------------------------------
  // run the analysis after successful initialisation
  //------------------------------------------------------------------------------------------------
  ana->Run(!gROOT->IsBatch());  

  return;
}
