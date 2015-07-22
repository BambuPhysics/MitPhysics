#if !defined(__CINT__) || defined(__MAKECINT__)
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/Catalog/interface/Catalog.h"
#include "MitAna/TreeMod/interface/Analysis.h"
#include "MitAna/TreeMod/interface/OutputMod.h"
#include "MitAna/TreeMod/interface/HLTMod.h"
#include "MitAna/PhysicsMod/interface/RunLumiSelectionMod.h"
#include "MitPhysics/Mods/interface/GeneratorMod.h"
#include "MitPhysics/Mods/interface/GoodPVFilterMod.h"
#include "MitPhysics/Mods/interface/SeparatePileUpMod.h"
#include "MitPhysics/Mods/interface/MuonIdMod.h"
#include "MitPhysics/Utils/interface/MuonTools.h"
#include "MitPhysics/Mods/interface/ElectronIdMod.h"
#include "MitPhysics/Utils/interface/ElectronTools.h"
#include "MitPhysics/Mods/interface/PhotonIdMod.h"
#include "MitPhysics/Mods/interface/PFTauIdMod.h"
#include "MitPhysics/Mods/interface/JetIdMod.h"
#include "MitPhysics/Mods/interface/JetCorrectionMod.h"
#include "MitPhysics/Mods/interface/MetCorrectionMod.h"
#include "MitPhysics/Skim/interface/MonoJetAnalysisMod.h"

#include "TSystem.h"

#include <iostream>
#endif

using namespace mithep;

//--------------------------------------------------------------------------------------------------
void runMonoJetSkim(const char *fileset    = "0000",
		    const char *skim       = "noskim",
		    const char *dataset    = "WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM",
		    const char *book       = "t2mit/filefi/041",
		    const char *catalogDir = "/home/cmsprod/catalog",
		    const char *outputLabel = "monojet",
		    int         nEvents    = 1000)
{
  //------------------------------------------------------------------------------------------------
  // json parameters get passed through the environment
  // for MC, the value must be "~"
  //------------------------------------------------------------------------------------------------
  TString json(gSystem->Getenv("MIT_PROD_JSON"));
  if (json.Length() == 0) {
    printf(" JSON file was not properly defined. EXIT!\n");
    return;
  }

  TString jsonFile = TString("/home/cmsprod/cms/json/") + json;
  Bool_t  isData   = (json != "~");

  TString MitData(gSystem->Getenv("MIT_DATA"));
  if(MitData.Length() == 0){
    MitData = gSystem->Getenv("CMSSW_BASE");
    MitData += "/src/MitPhysics/data";
  }

  printf("\n Initialization worked: \n\n");
  printf("   JSON   : %s (file: %s)\n",  json.Data(), jsonFile.Data());
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
  std::vector<mithep::BaseMod*> modules;

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

  modules.push_back(runLumiSel);

  // Generator info
  GeneratorMod *generatorMod = 0;
  if (!isData) {
    generatorMod = new GeneratorMod;
    generatorMod->SetPrintDebug(kFALSE);
    generatorMod->SetPtLeptonMin(0.0);
    generatorMod->SetEtaLeptonMax(2.7);
    generatorMod->SetPtPhotonMin(0.0);
    generatorMod->SetEtaPhotonMax(2.7);
    generatorMod->SetPtRadPhotonMin(0.0);
    generatorMod->SetEtaRadPhotonMax(2.7);
    generatorMod->SetIsData(isData);
    generatorMod->SetFillHist(! isData);
    generatorMod->SetApplyISRFilter(kFALSE);
    generatorMod->SetApplyVVFilter(kFALSE);
    generatorMod->SetApplyVGFilter(kFALSE);
    generatorMod->SetFilterBTEvents(kFALSE);

    modules.push_back(generatorMod);
  }

  //-----------------------------------------------------------------------------------------------------------
  // HLT information : trigger not applied (neither for data nor for MC, store info to apply selection offline
  //-----------------------------------------------------------------------------------------------------------
  HLTMod *hltMod = new HLTMod();

  // monojet triggers
  TString triggers[] = {"HLT_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight_v*",
                        "HLT_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight_v*"};
  for (auto& trigger : triggers)
    hltMod->AddTrigger(trigger);

  hltMod->SetBitsName("HLTBits");
  hltMod->SetTrigObjsName("MonoJetTriggerObjects");
  hltMod->SetAbortIfNotAccepted(isData);
  hltMod->SetAbortIfNoData(kFALSE);

  modules.push_back(hltMod);

  JetIdMod* leadJetId = new JetIdMod("LeadJetId");
  leadJetId->SetMinChargedHadronFraction(0.2);
  leadJetId->SetMaxNeutralHadronFraction(0.7);
  leadJetId->SetMaxNeutralEMFraction(0.7);
  leadJetId->SetApplyPFLooseId(kTRUE);
  leadJetId->SetPtMin(90.);
  leadJetId->SetMinOutput(1);

  modules.push_back(leadJetId);

  //------------------------------------------------------------------------------------------------
  // select events with a good primary vertex
  //------------------------------------------------------------------------------------------------
  GoodPVFilterMod *goodPvMod = new GoodPVFilterMod;
  goodPvMod->SetMinVertexNTracks(0);
  goodPvMod->SetMinNDof(4.0);
  goodPvMod->SetMaxAbsZ(24.0);
  goodPvMod->SetMaxRho(2.0);
  goodPvMod->SetIsMC(!isData);
  goodPvMod->SetVertexesName("PrimaryVertexes");

  modules.push_back(goodPvMod);

  //------------------------------------------------------------------------------------------------
  // split pfcandidates to PFPU and PFnoPU
  //------------------------------------------------------------------------------------------------
  SeparatePileUpMod* SepPUMod = new SeparatePileUpMod;
  SepPUMod->SetPFNoPileUpName("pfnopileupcands");
  SepPUMod->SetPFPileUpName("pfpileupcands");
  SepPUMod->SetCheckClosestZVertex(kFALSE);

  modules.push_back(SepPUMod);
 
  //-----------------------------------
  // Lepton Selection 
  //-----------------------------------
  MuonIdMod *vetoMuonIdMod = new MuonIdMod("VetoMuonId");
  vetoMuonIdMod->SetMuonClassType(mithep::MuonTools::kPFGlobalorTracker);
  vetoMuonIdMod->SetIdType(mithep::MuonTools::kNoId);
  vetoMuonIdMod->SetPFNoPileupCandidatesName(SepPUMod->GetPFNoPileUpName());
  vetoMuonIdMod->SetPFPileupCandidatesName(SepPUMod->GetPFPileUpName());
  vetoMuonIdMod->SetIsoType(mithep::MuonTools::kPFIsoBetaPUCorrected);
  vetoMuonIdMod->SetApplyD0Cut(kTRUE);
  vetoMuonIdMod->SetApplyDZCut(kTRUE);
  vetoMuonIdMod->SetWhichVertex(0);
  vetoMuonIdMod->SetPtMin(10.);
  vetoMuonIdMod->SetEtaMax(2.4);
  vetoMuonIdMod->SetOutputName("VetoMuons");

  modules.push_back(vetoMuonIdMod);

  MuonIdMod *muonIdMod = new MuonIdMod("GoodMuonId");
  muonIdMod->SetInputName(vetoMuonIdMod->GetOutputName());
  muonIdMod->SetMuonClassType(mithep::MuonTools::kPFGlobalorTracker);
  muonIdMod->SetIdType(mithep::MuonTools::kMuonPOG2012CutBasedIdTight);
  muonIdMod->SetPFNoPileupCandidatesName(SepPUMod->GetPFNoPileUpName());
  muonIdMod->SetPFPileupCandidatesName(SepPUMod->GetPFPileUpName());
  muonIdMod->SetIsoType(mithep::MuonTools::kPFIsoBetaPUCorrectedTight);
  muonIdMod->SetApplyD0Cut(kTRUE);
  muonIdMod->SetApplyDZCut(kTRUE);
  muonIdMod->SetWhichVertex(0);
  muonIdMod->SetPtMin(20.);
  muonIdMod->SetEtaMax(2.4);
  muonIdMod->SetIsFilterMode(kFALSE);
  muonIdMod->SetOutputName("GoodMuons");

  modules.push_back(muonIdMod);

  ElectronIdMod* vetoEleIdMod = new ElectronIdMod("VetoElectronId");
  vetoEleIdMod->SetPtMin(10.);  
  vetoEleIdMod->SetEtaMax(2.4);
  vetoEleIdMod->SetApplyEcalFiducial(true);
  vetoEleIdMod->SetIdType(mithep::ElectronTools::kPhys14Veto);
  vetoEleIdMod->SetIsoType(mithep::ElectronTools::kPhys14VetoIso);
  vetoEleIdMod->SetConversionsName("Conversions");
  vetoEleIdMod->SetApplyConversionFilterType1(kTRUE);
  vetoEleIdMod->SetApplyConversionFilterType2(kFALSE);
  vetoEleIdMod->SetApplyD0Cut(kTRUE);
  vetoEleIdMod->SetApplyDZCut(kTRUE);
  vetoEleIdMod->SetWhichVertex(0);
  vetoEleIdMod->SetOutputName("VetoElectrons");

  modules.push_back(vetoEleIdMod);

  ElectronIdMod* eleIdMod = new ElectronIdMod("GoodElectronId");
  eleIdMod->SetInputName(vetoEleIdMod->GetOutputName());
  eleIdMod->SetPtMin(20.);
  eleIdMod->SetEtaMax(2.4);
  eleIdMod->SetApplyEcalFiducial(true);
  eleIdMod->SetIdType(mithep::ElectronTools::kPhys14Medium);
  eleIdMod->SetIsoType(mithep::ElectronTools::kPhys14MediumIso);
  eleIdMod->SetConversionsName("Conversions");
  eleIdMod->SetApplyConversionFilterType1(kTRUE);
  eleIdMod->SetApplyConversionFilterType2(kFALSE);
  eleIdMod->SetApplyD0Cut(kTRUE);
  eleIdMod->SetApplyDZCut(kTRUE);
  eleIdMod->SetWhichVertex(0);
  eleIdMod->SetIsFilterMode(kFALSE);
  eleIdMod->SetOutputName("GoodElectrons");

  modules.push_back(eleIdMod);

  //-----------------------------------
  // Photon Id 
  //-----------------------------------

  PhotonIdMod *vetoPhotonIdMod = new PhotonIdMod("VetoPhotonId");
  vetoPhotonIdMod->SetPtMin(15.);
  vetoPhotonIdMod->SetOutputName("VetoPhotons");
  vetoPhotonIdMod->SetIdType(mithep::PhotonTools::kPhys14Loose);
  vetoPhotonIdMod->SetIsoType(mithep::PhotonTools::kPhys14LooseIso);
  vetoPhotonIdMod->SetApplyElectronVeto(kTRUE);

  modules.push_back(vetoPhotonIdMod);

  PhotonIdMod *photonIdMod = new PhotonIdMod("GoodPhotonId");
  photonIdMod->SetInputName(vetoPhotonIdMod->GetOutputName());
  photonIdMod->SetPtMin(90.);
  photonIdMod->SetOutputName("GoodPhotons");
  photonIdMod->SetIdType(mithep::PhotonTools::kPhys14Medium);
  photonIdMod->SetIsoType(mithep::PhotonTools::kPhys14MediumIso);
  photonIdMod->SetApplyElectronVeto(kTRUE);
  photonIdMod->SetIsFilterMode(kFALSE);

  modules.push_back(photonIdMod);

  PFTauIdMod *pfTauIdMod = new PFTauIdMod;
  pfTauIdMod->SetPtMin(18.);
  pfTauIdMod->SetEtaMax(2.3);
  pfTauIdMod->SetInputName("HPSTaus");
  pfTauIdMod->AddDiscriminator(mithep::PFTau::kDiscriminationByDecayModeFindingNewDMs);
  pfTauIdMod->AddCutDiscriminator(mithep::PFTau::kDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits,5);
  pfTauIdMod->SetOutputName("GoodTaus");

  modules.push_back(pfTauIdMod);

  JetCorrectionMod *jetCorr = new JetCorrectionMod;
  if (isData){ 
    jetCorr->AddCorrectionFromFile((MitData+TString("/Summer13_V1_DATA_L1FastJet_AK5PF.txt")).Data()); 
    jetCorr->AddCorrectionFromFile((MitData+TString("/Summer13_V1_DATA_L2Relative_AK5PF.txt")).Data()); 
    jetCorr->AddCorrectionFromFile((MitData+TString("/Summer13_V1_DATA_L3Absolute_AK5PF.txt")).Data()); 
    jetCorr->AddCorrectionFromFile((MitData+TString("/Summer13_V1_DATA_L2L3Residual_AK5PF.txt")).Data());
  }                                                                                      
  else {                                                                                 
    jetCorr->AddCorrectionFromFile((MitData+TString("/MCRUN2_74_V9_L1FastJet_AK4PFchs.txt")).Data()); 
    jetCorr->AddCorrectionFromFile((MitData+TString("/MCRUN2_74_V9_L2Relative_AK4PFchs.txt")).Data()); 
    jetCorr->AddCorrectionFromFile((MitData+TString("/MCRUN2_74_V9_L3Absolute_AK4PFchs.txt")).Data()); 
  }
  jetCorr->SetInputName("AKt4PFJetsCHS");
  jetCorr->SetCorrectedName("CorrectedJets");    

  modules.push_back(jetCorr);

  JetIdMod *jetId = new JetIdMod;
  jetId->SetInputName(jetCorr->GetOutputName());
  jetId->SetOutputName("GoodJets");
  jetId->SetPtMin(30.0);
  jetId->SetEtaMax(2.5);
  jetId->SetApplyPFLooseId(kTRUE);
  jetId->SetMVATrainingSet(JetIDMVA::k53BDTCHSFullPlusRMS);
  jetId->SetMVACutWP(JetIDMVA::kLoose);
  jetId->SetMVACutsFile(MitData + "/jetIDCuts_121221.dat");
  jetId->SetMVAWeightsFile(MitData + "/TMVAClassification_5x_BDT_chsFullPlusRMS.weights.xml");

  modules.push_back(jetId);

  MetCorrectionMod *type1MetCorr = new MetCorrectionMod;
  type1MetCorr->SetInputName("PFMet");
  type1MetCorr->SetOutputName("PFType1CorrectedMet");
  type1MetCorr->SetJetsName("AKt4PFJets");
  type1MetCorr->SetRhoAlgo(PileupEnergyDensity::kFixedGridFastjetAll);
  type1MetCorr->SetMaxEMFraction(0.9);
  type1MetCorr->SetSkipMuons(kTRUE);
  type1MetCorr->ApplyType0(kFALSE);
  type1MetCorr->ApplyType1(kTRUE);
  type1MetCorr->ApplyShift(kFALSE);
  type1MetCorr->AddJetCorrectionFromFile(MitData + "/MCRUN2_74_V9_L1FastJet_AK4PF.txt");
  type1MetCorr->AddJetCorrectionFromFile(MitData + "/MCRUN2_74_V9_L2Relative_AK4PF.txt");
  type1MetCorr->AddJetCorrectionFromFile(MitData + "/MCRUN2_74_V9_L3Absolute_AK4PF.txt");
  type1MetCorr->IsData(isData);

  modules.push_back(type1MetCorr);

  //------------------------------------------------------------------------------------------------
  // select events
  //------------------------------------------------------------------------------------------------

  float maxJetEta       = 4.7;
  float minMet          = 150.;
  float minLeadJetPt    = 90.;

  MonoJetAnalysisMod *monojetSel = new MonoJetAnalysisMod("MonoJetSelector");
  monojetSel->SetMetName(type1MetCorr->GetOutputName());
  monojetSel->SetJetsName(jetId->GetOutputName());
  monojetSel->SetVetoElectronsName(vetoEleIdMod->GetOutputName());
  monojetSel->SetElectronMaskName(eleIdMod->GetOutputName());
  monojetSel->SetVetoMuonsName(vetoMuonIdMod->GetOutputName());
  monojetSel->SetMuonMaskName(muonIdMod->GetOutputName());
  monojetSel->SetVetoTausName(pfTauIdMod->GetOutputName());
  monojetSel->SetVetoPhotonsName(vetoPhotonIdMod->GetOutputName());
  monojetSel->SetPhotonMaskName(photonIdMod->GetOutputName());
  monojetSel->SetCategoryFlagsName("MonoJetEventCategories");

  // Using uniform setup for all categories
  for (unsigned iCat = 0; iCat != MonoJetAnalysisMod::nMonoJetCategories; ++iCat) {
    monojetSel->SetCategoryActive(iCat, kTRUE);
    monojetSel->SetMaxNumJets(iCat, 3);
    monojetSel->SetMinLeadJetPt(iCat, minLeadJetPt);
    monojetSel->SetMaxJetEta(iCat, maxJetEta);
    monojetSel->SetMinMetPt(iCat, minMet);
    monojetSel->SetMinChargedHadronFrac(iCat, 0.2); 
    monojetSel->SetMaxNeutralHadronFrac(iCat, 0.7);
    monojetSel->SetMaxNeutralEmFrac(iCat, 0.7);
  }

  modules.push_back(monojetSel);

  //------------------------------------------------------------------------------------------------
  // skim output
  //------------------------------------------------------------------------------------------------
  TString outputName = TString(outputLabel);
  outputName += TString("_") + TString(dataset) + TString("_") + TString(skim);
  if (TString(fileset) != TString(""))
    outputName += TString("_") + TString(fileset);
  
  OutputMod *skimOutput = new OutputMod;
  skimOutput->Drop("*");
  skimOutput->Keep("HLT*");

  skimOutput->Keep("MC*");
  skimOutput->Keep("PileupInfo");
  skimOutput->Keep("Rho");
  skimOutput->Keep("EvtSelData");
  skimOutput->Keep("BeamSpot");
  skimOutput->Keep("PrimaryVertexes");
  skimOutput->Keep("PFMet");
  skimOutput->Keep("AKt4PFJetsCHS");
  skimOutput->Keep("AKt8PFJetsCHS");
  skimOutput->Keep("Electrons");
  skimOutput->Keep("Conversions");
  skimOutput->Keep("*Stable*");
  skimOutput->Keep("Muons");
  skimOutput->Keep("HPSTaus");
  skimOutput->Keep("Photons");
  skimOutput->Keep("AKT4GenJets");
  skimOutput->AddNewBranch(monojetSel->GetCategoryFlagsName());

  skimOutput->SetMaxFileSize(10 * 1024); // 10 GB - should never exceed
  skimOutput->SetFileName(outputName);
  skimOutput->SetPathName(".");

  modules.push_back(skimOutput);

  //------------------------------------------------------------------------------------------------
  // making the analysis chain
  //------------------------------------------------------------------------------------------------
  auto mItr(modules.begin());
  while (true) {
    auto* mod = *(mItr++);
    if (mItr == modules.end())
      break;
    mod->Add(*mItr);
  }

  //------------------------------------------------------------------------------------------------
  // setup analysis
  //------------------------------------------------------------------------------------------------
  Analysis *ana = new Analysis;
  ana->SetUseCacher(1);
  ana->SetUseHLT(kTRUE);
  ana->SetAllowNoHLTTree(kTRUE); // for private MC with no HLT info
  ana->SetKeepHierarchy(kFALSE);
  ana->SetPrintScale(100);
  ana->SetOutputName(outputName + "_hist.root");
  if (nEvents >= 0)
    ana->SetProcessNEvents(nEvents);

  ana->AddSuperModule(modules.front());

  //------------------------------------------------------------------------------------------------
  // organize input
  //------------------------------------------------------------------------------------------------
  Catalog *c = new Catalog(catalogDir);
  TString skimDataset = TString(dataset)+TString("/") +TString(skim);
  Dataset *d = NULL;
  if (TString(skim).CompareTo("noskim") == 0)
    d = c->FindDataset(book,dataset,fileset, 1); // 1 to use smartcache
  else
    d = c->FindDataset(book,skimDataset.Data(),fileset, 1);
  ana->AddDataset(d);
  ana->SetCacheSize(0);

  //------------------------------------------------------------------------------------------------
  // Say what we are doing
  //------------------------------------------------------------------------------------------------
  printf("\n==== PARAMETER SUMMARY FOR THIS JOB ====\n");
  printf("\n JSON file: %s\n",jsonFile.Data());
  printf("\n Rely on Catalog: %s\n",catalogDir);
  printf("  -> Book: %s  Dataset: %s  Skim: %s  Fileset: %s <-\n",book,dataset,skim,fileset);
  printf("\n========================================\n");

  std::cout << std::endl;
  std::cout << "+++++ ANALYSIS FLOW +++++" << std::endl;
  ana->PrintModuleTree();
  std::cout << std::endl;
  std::cout << "+++++++++++++++++++++++++" << std::endl;

  //------------------------------------------------------------------------------------------------
  // run the analysis after successful initialisation
  //------------------------------------------------------------------------------------------------
  ana->Run(false);

  delete ana; // all modules deleted recursively

  // rename the output file so that condor can see it
  gSystem->Rename("./" + outputName + "_000.root", "./" + outputName + ".root");

  return;
}

