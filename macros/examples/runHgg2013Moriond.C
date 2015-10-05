// $Id: runHgg2012HCP.C,v 1.2 2012/10/24 14:47:27 mingyang Exp $
#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <TSystem.h>
#include <TProfile.h>
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/Catalog/interface/Catalog.h"
#include "MitAna/TreeMod/interface/Analysis.h"
#include "MitAna/TreeMod/interface/HLTMod.h"
#include "MitAna/PhysicsMod/interface/RunLumiSelectionMod.h"
#include "MitAna/PhysicsMod/interface/MCProcessSelectionMod.h"
#include "MitAna/PhysicsMod/interface/PublisherMod.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitAna/DataTree/interface/Names.h"
#include "MitPhysics/Mods/interface/GoodPVFilterMod.h"
#include "MitPhysics/Mods/interface/MuonIDMod.h"
#include "MitPhysics/Mods/interface/ElectronIDMod.h"
#include "MitPhysics/Mods/interface/ElectronCleaningMod.h"
#include "MitPhysics/Mods/interface/PhotonIDMod.h"
#include "MitPhysics/Mods/interface/PhotonPairSelector.h"
#include "MitPhysics/Mods/interface/PhotonTreeWriter.h"
#include "MitPhysics/Mods/interface/PhotonCleaningMod.h"
#include "MitPhysics/Mods/interface/MergeLeptonsMod.h"
#include "MitPhysics/Mods/interface/JetCorrectionMod.h"
#include "MitPhysics/Mods/interface/PhotonMvaMod.h"
#include "MitPhysics/Mods/interface/MVASystematicsMod.h"
#include "MitPhysics/Mods/interface/SeparatePileUpMod.h"

#endif

//--------------------------------------------------------------------------------------------------
void runHgg2013Moriond(const char *fileset    = "0000",
		       const char *skim       = "noskim",
		       //const char *dataset = "f11--h120gg-gf-v14b-pu",
		       //const char *dataset = "r11a-pho-j16-v1",   
		       //const char *dataset = "meridiani2012-diphoj-v9",
		       //const char *dataset = "s12-h125gg-gf-v9",
		       //const char *dataset = "s12-pj40-2em-v9",
		       //const char *dataset = "s12-diphoj-3-v9",
		       //const char *dataset = "s12-zllm50-v9",
		       //const char *dataset = "f11--h121gg-gf-v14b-pu",
		       //const char *dataset = "r12a-pho-pr-v1",
		       //const char *dataset = "s12-pj40-2em-v9",
		       //const char *dataset = "s12-h125gg-vbf-v7a", 
		       //const char *dataset = "s12-h125gg-vh-v7a", 
		       const char *dataset = "s12-h125gg-gf-v7a", 
		       //const char *dataset = "r12a-pho-j13-v1", 
		       //const char *dataset = "r12b-dph-j13-v1", 
		       //const char *dataset = "s12-zllm50-v7a", 
		       //const char *book       = "t2mit/filefi/029",
		       const char *book       = "local/filefi/029",
		       const char *catalogDir = "/home/cmsprod/catalog",
		       const char *outputName = "hgg",
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
  
  TString jsonFile = TString("/home/mingyang/cms/json/") + TString(json);
  //TString jsonFile = TString("/home/mingyang/cms/json/") + TString("Cert_136033-149442_7TeV_Dec22ReReco_Collisions10_JSON_v4.txt");
  Bool_t  isData   = ( (jsonFile.CompareTo("/home/mingyang/cms/json/~") != 0) );
  
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

  MCProcessSelectionMod *mcselmod = new MCProcessSelectionMod;
  
  MVASystematicsMod *sysMod = new MVASystematicsMod;
  sysMod->SetMCR9Scale(1.0035, 1.0035);  
  sysMod->SetIsData(isData);
  
  // only select on run- and lumisection numbers when valid json file present
  if ((jsonFile.CompareTo("/home/mingyang/cms/json/~") != 0) &&
      (jsonFile.CompareTo("/home/mingyang/cms/json/-") != 0)   ) {
    runLumiSel->AddJSONFile(jsonFile.Data());
  }
  if ((jsonFile.CompareTo("/home/mingyang/cms/json/-") == 0)   ) {
    printf("\n WARNING -- Looking at data without JSON file: always accept.\n\n");
    runLumiSel->SetAbortIfNotAccepted(kFALSE);   // accept all events if there is no valid JSON file
  }

  printf("\n Run lumi worked. \n\n");

  //------------------------------------------------------------------------------------------------
  // HLT information
  //------------------------------------------------------------------------------------------------
  HLTMod *hltModM = new HLTMod("HLTModM");
  hltModM->AddTrigger("HLT_Mu9");
  hltModM->AddTrigger("HLT_Mu11");
  hltModM->AddTrigger("HLT_Mu15_v1");
  hltModM->SetTrigObjsName("MyHltMuonObjs");
  hltModM->SetAbortIfNotAccepted(kFALSE);

  HLTMod *hltModE = new HLTMod("HLTModE");
  hltModE->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v1",150000,161176);
  hltModE->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2",161179,163261);
  hltModE->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v3",163262,164237);
  hltModE->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v4",165085,165888);
  hltModE->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v5",165900,166967);
  hltModE->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v6",166968,170053);
  hltModE->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6",170054,170759);
  hltModE->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7",170760,173198);
  hltModE->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8",173199,178380);
  hltModE->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9",178381,179889);
  hltModE->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v10",179890,189999);
  hltModE->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",190000,999999);
  hltModE->SetTrigObjsName("MyHltElecObjs");
  hltModE->SetAbortIfNotAccepted(isData);

  HLTMod *hltModES = new HLTMod("HLTModES");
  hltModES->AddTrigger("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v*",160000,189999);
  hltModES->AddTrigger("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2",160000,189999);
  hltModES->AddTrigger("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3",160000,189999);
  hltModES->AddTrigger("HLT_Ele25_WP80_PFMT40_v*",160000,189999);
  hltModES->AddTrigger("HLT_Ele27_WP80_PFMT50_v*",160000,189999);
  hltModES->AddTrigger("HLT_Ele27_WP80_v*",190000,999999);
  hltModES->SetTrigObjsName("MyHltElecSObjs");
  hltModES->SetAbortIfNotAccepted(isData);

  HLTMod *hltModP = new HLTMod("HLTModP");
  
  hltModP->AddTrigger("HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v*",160000,161176);
  
  hltModP->AddTrigger("HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v*",161216,165633);
  hltModP->AddTrigger("HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v*",161216,165633);
  hltModP->AddTrigger("HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v*",161216,165633);
  hltModP->AddTrigger("HLT_Photon20_R9Id_Photon18_R9Id_v*",161216,165633);

  hltModP->AddTrigger("HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v*",165970,173198);
  hltModP->AddTrigger("HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v*",165970,173198);
  hltModP->AddTrigger("HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v*",165970,173198);
  hltModP->AddTrigger("HLT_Photon26_R9Id_Photon18_R9Id_v*",165970,173198);
  hltModP->AddTrigger("HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL_v*",165970,173198);
  hltModP->AddTrigger("HLT_Photon36_R9Id_Photon22_CaloIdL_IsoVL_v*",165970,173198);
  hltModP->AddTrigger("HLT_Photon36_R9Id_Photon22_R9Id_v*",165970,173198);
  
  hltModP->AddTrigger("HLT_Photon36_CaloId_IsoVL_Photon22_R9Id_v*",165970,166967);
  hltModP->AddTrigger("HLT_Photon36_CaloIdL_IsoVL_Photon22_R9Id_v*",167039,173198);
  
  hltModP->AddTrigger("HLT_Photon26_CaloIdXL_IsoXL_Photon18_CaloIdXL_IsoXL_v*",173236,178380);
  hltModP->AddTrigger("HLT_Photon26_CaloIdXL_IsoXL_Photon18_R9Id_v*",173236,178380);
  hltModP->AddTrigger("HLT_Photon26_R9Id_Photon18_CaloIdXL_IsoXL_v*",173236,178380);
  hltModP->AddTrigger("HLT_Photon26_R9Id_Photon18_R9Id_v*",173236,178380);
  hltModP->AddTrigger("HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL_v*",173236,178380);
  hltModP->AddTrigger("HLT_Photon36_CaloIdL_IsoVL_Photon22_R9Id_v*",173236,178380);
  hltModP->AddTrigger("HLT_Photon36_R9Id_Photon22_CaloIdL_IsoVL_v*",173236,178380);
  hltModP->AddTrigger("HLT_Photon36_R9Id_Photon22_R9Id_v*",173236,178380);
  
  hltModP->AddTrigger("HLT_Photon26_CaloIdXL_IsoXL_Photon18_CaloIdXL_IsoXL_Mass60_v*",178420,189999);
  hltModP->AddTrigger("HLT_Photon26_CaloIdXL_IsoXL_Photon18_R9IdT_Mass60_v*",178420,189999);
  hltModP->AddTrigger("HLT_Photon26_R9IdT_Photon18_CaloIdXL_IsoXL_Mass60_v*",178420,189999);
  hltModP->AddTrigger("HLT_Photon26_R9IdT_Photon18_R9IdT_Mass60_v*",178420,189999);
  hltModP->AddTrigger("HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL_v*",178420,189999);
  hltModP->AddTrigger("HLT_Photon36_CaloIdL_IsoVL_Photon22_R9Id_v*",178420,189999);
  hltModP->AddTrigger("HLT_Photon36_R9Id_Photon22_CaloIdL_IsoVL_v*",178420,189999);
  hltModP->AddTrigger("HLT_Photon36_R9Id_Photon22_R9Id_v*",178420,189999);

  //hltModP->AddTrigger("HLT_Photon26_CaloId10_Iso50_Photon18_CaloId10_Iso50_Mass60_v*",190000,999999); 
  //hltModP->AddTrigger("HLT_Photon26_CaloId10_Iso50_Photon18_R9Id85_Mass60_v*",190000,999999); 
  //hltModP->AddTrigger("HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass60_v*",190000,999999); 
  //hltModP->AddTrigger("HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass70_v*",190000,999999); 
  //hltModP->AddTrigger("HLT_Photon26_R9Id85_Photon18_CaloId10_Iso50_Mass60_v*",190000,999999); 
  //hltModP->AddTrigger("HLT_Photon26_R9Id85_Photon18_R9Id85_Mass60_v*",190000,999999); 
  //hltModP->AddTrigger("HLT_Photon36_CaloId10_Iso50_Photon22_CaloId10_Iso50_v*",190000,999999); 
  //hltModP->AddTrigger("HLT_Photon36_CaloId10_Iso50_Photon22_R9Id85_v*",190000,999999); 
  //hltModP->AddTrigger("HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_R9Id85_OR_CaloId10_Iso50_v*",190000,999999); 
  //hltModP->AddTrigger("HLT_Photon36_R9Id85_Photon22_CaloId10_Iso50_v*",190000,999999); 
  //hltModP->AddTrigger("HLT_Photon36_R9Id85_Photon22_R9Id85_v*",190000,999999); 

  hltModP->AddTrigger("HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass*",190000,999999); 
  hltModP->AddTrigger("HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_R9Id85_OR_CaloId10_Iso50_v*",190000,999999); 
    
  hltModP->SetTrigObjsName("MyHltPhotObjs");
  hltModP->SetAbortIfNotAccepted(isData);
  //------------------------------------------------------------------------------------------------
  // split pfcandidates to PFPU and PFnoPU
  //------------------------------------------------------------------------------------------------
  SeparatePileUpMod* SepPUMod = new SeparatePileUpMod;
  //  SepPUMod->SetUseAllVerteces(kFALSE);
  // SepPUMod->SetVertexName("OutVtxCiC");
  SepPUMod->SetPFNoPileUpName("pfnopileupcands");
  SepPUMod->SetPFPileUpName("pfpileupcands");
  SepPUMod->SetCheckClosestZVertex(kTRUE);
   
  //------------------------------------------------------------------------------------------------
  // select events with a good primary vertex
  //------------------------------------------------------------------------------------------------
  GoodPVFilterMod *goodPVFilterMod = new GoodPVFilterMod;
  goodPVFilterMod->SetMinVertexNTracks(0);
  goodPVFilterMod->SetMinNDof         (4.0);
  goodPVFilterMod->SetMaxAbsZ         (24.0);
  goodPVFilterMod->SetMaxRho          (2.0);
  goodPVFilterMod->SetAbortIfNotAccepted(kFALSE);
  goodPVFilterMod->SetIsMC(!isData);

  GoodPVFilterMod *goodPVFilterModE = new GoodPVFilterMod("GoodPVFilterModE");
  goodPVFilterModE->SetOutputName("GoodVertexesE");
  goodPVFilterModE->SetMinVertexNTracks(0);
  goodPVFilterModE->SetMinNDof         (4.0);
  goodPVFilterModE->SetMaxAbsZ         (24.0);
  goodPVFilterModE->SetMaxRho          (2.0);  
  goodPVFilterModE->SetIsMC(!isData);
  
  GoodPVFilterMod *goodPVFilterModES = new GoodPVFilterMod("GoodPVFilterModES");
  goodPVFilterModES->SetOutputName("GoodVertexesES");
  goodPVFilterModES->SetMinVertexNTracks(0);
  goodPVFilterModES->SetMinNDof         (4.0);
  goodPVFilterModES->SetMaxAbsZ         (24.0);
  goodPVFilterModES->SetMaxRho          (2.0);    
  

  //------------------------------------------------------------------------------------------------
  // object id and cleaning sequence
  //------------------------------------------------------------------------------------------------
  //-----------------------------------
  // Lepton Selection 
  //-----------------------------------
  ElectronIDMod* eleIdMod = new ElectronIDMod;
  eleIdMod -> SetPtMin(20);  
  eleIdMod -> SetEtaMax(2.5);
  eleIdMod -> SetApplyEcalFiducial(true);
  eleIdMod -> SetIDType("Hgg_LeptonTag_2012IdHCP");  
  eleIdMod -> SetElectronMVAWeightsSubdet0Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat1.weights.xml"))));
  eleIdMod -> SetElectronMVAWeightsSubdet1Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat2.weights.xml"))));  
  eleIdMod -> SetElectronMVAWeightsSubdet2Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat3.weights.xml"))));  
  eleIdMod -> SetElectronMVAWeightsSubdet0Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat4.weights.xml")))); 
  eleIdMod -> SetElectronMVAWeightsSubdet1Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat5.weights.xml")))); 
  eleIdMod -> SetElectronMVAWeightsSubdet2Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat6.weights.xml")))); 
  eleIdMod -> SetWhichVertex(-1);
  eleIdMod -> SetD0Cut(0.02);
  eleIdMod -> SetDZCut(0.2); //h  
  eleIdMod -> SetIsoType("PFIso_HggLeptonTag2012HCP"); //h
  eleIdMod -> SetOutputName("HggLeptonTagElectrons");
  eleIdMod -> SetRhoType(RhoUtilities::CMS_RHO_RHOKT6PFJETS);
  eleIdMod -> SetPFNoPileUpName("pfnopileupcands");
  eleIdMod -> SetInvertNExpectedHitsInnerCut(kFALSE);
  eleIdMod -> SetNExpectedHitsInnerCut(1);   
  eleIdMod -> SetApplyConversionFilterType1(kTRUE);
  eleIdMod -> SetPVName(Names::gkPVBeamSpotBrn);   
   
  MuonIDMod* muonIdMod = new MuonIDMod;
  // base kinematics
  muonIdMod -> SetPtMin(20.);
  muonIdMod -> SetEtaCut(2.4);
  // base ID
  //muonIdMod -> SetIDType("Tight");
  muonIdMod -> SetIDType("muonPOG2012CutBasedIDTight");
  muonIdMod -> SetWhichVertex(-1); // this is a 'hack'.. but hopefully good enough...
  muonIdMod -> SetD0Cut(0.2);
  muonIdMod -> SetDZCut(0.5);
  muonIdMod -> SetIsoType("PFIsoBetaPUCorrectedTight"); //h
  muonIdMod -> SetPFIsoCut(0.2); //h
  muonIdMod -> SetOutputName("HggLeptonTagMuons");
  muonIdMod -> SetPFNoPileUpName("pfnopileupcands");
  muonIdMod -> SetPFPileUpName("pfpileupcands");
  muonIdMod -> SetPVName(Names::gkPVBeamSpotBrn); 
  
  PublisherMod<PFJet,Jet> *pubJet = new PublisherMod<PFJet,Jet>("JetPub");
  pubJet->SetInputName("AKt5PFJets");
  pubJet->SetOutputName("PubAKt5PFJets");
  
  PublisherMod<PFJet,Jet> *pubJetOpen = new PublisherMod<PFJet,Jet>("JetPubOpen");
  pubJetOpen->SetInputName("AKt5PFJets");
  pubJetOpen->SetOutputName("PubAKt5PFJetsOpen");  
  
  JetCorrectionMod *jetCorr = new JetCorrectionMod;
  /*if(isData){ 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/GR_P_V42_AN3_L1FastJet_AK5PF.txt")).Data())); 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/GR_P_V42_AN3_L2Relative_AK5PF.txt")).Data())); 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/GR_P_V42_AN3_L3Absolute_AK5PF.txt")).Data())); 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/GR_P_V42_AN3_L2L3Residual_AK5PF.txt")).Data())); 
  }
  else {
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/START53_V15_L1FastJet_AK5PF.txt")).Data())); 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/START53_V15_L2Relative_AK5PF.txt")).Data())); 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/START53_V15_L3Absolute_AK5PF.txt")).Data())); 
    }*/
  if(isData){ 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer12_V7_DATA_L1FastJet_AK5PF.txt")).Data())); 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer12_V7_DATA_L2Relative_AK5PF.txt")).Data())); 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer12_V7_DATA_L3Absolute_AK5PF.txt")).Data())); 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer12_V7_DATA_L2L3Residual_AK5PF.txt")).Data())); 
  }
  else {
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer12_V7_MC_L1FastJet_AK5PF.txt")).Data())); 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer12_V7_MC_L2Relative_AK5PF.txt")).Data())); 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer12_V7_MC_L3Absolute_AK5PF.txt")).Data())); 
  }
  jetCorr->SetInputName(pubJet->GetOutputName());
  jetCorr->SetCorrectedName("CorrectedJets");    
  
  
  Bool_t excludedoubleprompt = kFALSE;
  if (TString(dataset).Contains("-pj")) {
    mcselmod->ExcludeProcess(18);
    mcselmod->ExcludeProcess(114);
    excludedoubleprompt = kTRUE;
  }
  
  if (TString(dataset).Contains("-qcd2em") || TString(dataset).Contains("-qcd-2em")) {
    excludedoubleprompt = kTRUE;
  }
  
  Bool_t is25 = kFALSE;
  if (TString(book).Contains("025")) is25 = kTRUE;
  
  PhotonMvaMod *photreg = new PhotonMvaMod;
  photreg->SetRegressionVersion(3);
  photreg->SetRegressionWeights(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/gbrv3ph_52x.root")).Data()));
  photreg->SetOutputName("GoodPhotonsRegr");
  photreg->SetApplyShowerRescaling(kTRUE);
  photreg->SetIsData(isData);

  
  PhotonPairSelector         *photcic = new PhotonPairSelector("PhotonPairSelectorCiC");
  photcic->SetOutputName("GoodPhotonsCIC");
  photcic->SetOutputVtxName("OutVtxCiC");        
  photcic->SetPhotonSelType("CiCPFSelection");
  photcic->SetVertexSelType("CiCMVA2012Selection");
  photcic->SetUseSingleLegConversions(kFALSE);
  photcic->DoMCSmear(kTRUE);
  photcic->DoDataEneCorr(kTRUE);
  photcic->SetPhotonsFromBranch(kFALSE);
  photcic->SetInputPhotonsName(photreg->GetOutputName());
  
  //------------------------------------------2012 HCP--------------------------------------------------------------
  photcic->SetMCSmearFactors2012HCP(0.0111,0.0111,0.0107,0.0107,0.0155,0.0194,0.0295,0.0276,0.037,0.0371);
  photcic->SetMCSmearFactors2012HCPMVA(0.01,0.011,0.0106,0.0106,0.0186,0.0196,0.0283,0.0267,0.0343,0.0345);
  photcic->UseSpecialSmearForDPMVA(true);
  photcic->AddEnCorrPerRun2012HCP(190645,190781,0.9964,0.9964,1.0020,1.0020,0.9893,1.0028,0.9871,0.9937,0.9839,0.9958);
  photcic->AddEnCorrPerRun2012HCP(190782,191042,1.0024,1.0024,1.0079,1.0079,0.9923,1.0058,0.9911,0.9977,0.9886,1.0005);
  photcic->AddEnCorrPerRun2012HCP(191043,193555,0.9935,0.9935,0.9991,0.9991,0.9861,0.9997,0.9894,0.9960,0.9864,0.9982);
  photcic->AddEnCorrPerRun2012HCP(193556,194150,0.9920,0.9920,0.9976,0.9976,0.9814,0.9951,0.9896,0.9962,0.9872,0.9990);
  photcic->AddEnCorrPerRun2012HCP(194151,194532,0.9925,0.9925,0.9981,0.9981,0.9826,0.9963,0.9914,0.9980,0.9874,0.9993);
  photcic->AddEnCorrPerRun2012HCP(194533,195113,0.9927,0.9927,0.9983,0.9983,0.9844,0.9981,0.9934,0.9999,0.9878,0.9996);
  photcic->AddEnCorrPerRun2012HCP(195114,195915,0.9929,0.9929,0.9984,0.9984,0.9838,0.9974,0.9942,1.0007,0.9878,0.9997);
  photcic->AddEnCorrPerRun2012HCP(195916,198115,0.9919,0.9919,0.9975,0.9975,0.9827,0.9964,0.9952,1.0017,0.9869,0.9987);
  photcic->AddEnCorrPerRun2012HCP(198116,199803,0.9955,0.9955,1.0011,1.0011,0.9859,0.9995,0.9893,0.9959,0.9923,1.0041);
  photcic->AddEnCorrPerRun2012HCP(199804,200048,0.9967,0.9967,1.0023,1.0023,0.9870,1.0006,0.9893,0.9959,0.9937,1.0055);
  photcic->AddEnCorrPerRun2012HCP(200049,200151,0.9980,0.9980,1.0036,1.0036,0.9877,1.0012,0.9910,0.9976,0.9980,1.0097);
  photcic->AddEnCorrPerRun2012HCP(200152,200490,0.9958,0.9958,1.0013,1.0013,0.9868,1.0004,0.9922,0.9988,0.9948,1.0065);
  photcic->AddEnCorrPerRun2012HCP(200491,200531,0.9979,0.9979,1.0035,1.0035,0.9876,1.0012,0.9915,0.9981,0.9979,1.0096);
  photcic->AddEnCorrPerRun2012HCP(200532,201656,0.9961,0.9961,1.0017,1.0017,0.9860,0.9996,0.9904,0.9970,0.9945,1.0063);
  photcic->AddEnCorrPerRun2012HCP(201657,202305,0.9969,0.9969,1.0025,1.0025,0.9866,1.0002,0.9914,0.9980,0.9999,1.0116);
  photcic->AddEnCorrPerRun2012HCP(202305,203002,0.9982,0.9982,1.0038,1.0038,0.9872,1.0008,0.9934,1.0000,1.0018,1.0135);
  photcic->AddEnCorrPerRun2012HCP(203003,203984,1.0006,1.0006,1.0061,1.0061,0.9880,1.0017,0.9919,0.9988,0.9992,1.0104);	
  photcic->AddEnCorrPerRun2012HCP(203985,205085,0.9993,0.9993,1.0048,1.0048,0.9903,1.0040,0.9928,0.9997,0.9987,1.0099);	
  photcic->AddEnCorrPerRun2012HCP(205086,205310,1.0004,1.0004,1.0059,1.0059,0.9901,1.0037,0.9987,1.0055,1.0091,1.0202);	
  photcic->AddEnCorrPerRun2012HCP(205311,206207,1.0000,1.0000,1.0055,1.0055,0.9891,1.0028,0.9948,1.0017,1.0032,1.0144);	
  photcic->AddEnCorrPerRun2012HCP(206208,206483,1.0003,1.0003,1.0058,1.0058,0.9895,1.0032,0.9921,0.9989,1.0056,1.0167);	
  photcic->AddEnCorrPerRun2012HCP(206484,206597,1.0005,1.0005,1.0060,1.0060,0.9895,1.0032,0.9968,1.0036,1.0046,1.0158);	
  photcic->AddEnCorrPerRun2012HCP(206598,206896,1.0006,1.0006,1.0061,1.0061,0.9881,1.0017,0.9913,0.9982,1.0050,1.0162);	
  photcic->AddEnCorrPerRun2012HCP(206897,207220,1.0006,1.0006,1.0061,1.0061,0.9884,1.0021,0.9909,0.9978,1.0053,1.0165);	
  photcic->AddEnCorrPerRun2012HCP(207221,208686,1.0006,1.0006,1.0061,1.0061,0.9894,1.0030,0.9951,1.0020,1.0060,1.0172);	

  //-----------------------------------------------------------------------------------------------------------------
  
  //  photcic->SetDoMCR9Scaling(kTRUE);
  //  photcic->SetMCR9Scale(1.0035, 1.0035);
  photcic->SetDoShowerShapeScaling(kTRUE); 
  photcic->SetShowerShapeType("2012ShowerShape");
  //photcic->SetDoMCErrScaling(kTRUE);
  //photcic->SetMCErrScale(1.07, 1.045); 
  //photcic->SetMCErrScale(1, 1); //ming:scale(sigE/E)
  photcic->SetJetsName(jetCorr->GetOutputName());    
  //photcic->SetRescaledBeamspotWidth(5.0);
  photcic->SetIsData(isData);
  photcic->SetApplyLeptonTag(kTRUE);
  photcic->SetLeptonTagElectronsName("HggLeptonTagElectrons");
  photcic->SetLeptonTagMuonsName("HggLeptonTagMuons");  
  photcic->Set2012HCP(kTRUE);
  photcic->SetIdMVAType("2012IdMVA_globe");
  
  PhotonPairSelector         *photcicnoeleveto = new PhotonPairSelector("PhotonPairSelectorCiCInvertEleVeto");
  photcicnoeleveto->SetOutputName("GoodPhotonsCICNoEleVeto");
  photcicnoeleveto->SetOutputVtxName("OutVtxCiCInvertEleVeto");      
  photcicnoeleveto->SetPhotonSelType("CiCPFSelection");
  photcicnoeleveto->SetVertexSelType("CiCMVA2012Selection");
  photcicnoeleveto->SetUseSingleLegConversions(kFALSE);
  photcicnoeleveto->DoMCSmear(kTRUE);
  photcicnoeleveto->DoDataEneCorr(kTRUE);
  photcicnoeleveto->SetPhotonsFromBranch(kFALSE);
  photcicnoeleveto->SetInputPhotonsName(photreg->GetOutputName());
  
  //------------------------------------------2012 HCP--------------------------------------------------------------
  photcicnoeleveto->SetMCSmearFactors2012HCP(0.0111,0.0111,0.0107,0.0107,0.0155,0.0194,0.0295,0.0276,0.037,0.0371);
  photcicnoeleveto->SetMCSmearFactors2012HCPMVA(0.01,0.011,0.0106,0.0106,0.0186,0.0196,0.0283,0.0267,0.0343,0.0345);
  photcicnoeleveto->UseSpecialSmearForDPMVA(true);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(190645,190781,0.9964,0.9964,1.0020,1.0020,0.9893,1.0028,0.9871,0.9937,0.9839,0.9958);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(190782,191042,1.0024,1.0024,1.0079,1.0079,0.9923,1.0058,0.9911,0.9977,0.9886,1.0005);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(191043,193555,0.9935,0.9935,0.9991,0.9991,0.9861,0.9997,0.9894,0.9960,0.9864,0.9982);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(193556,194150,0.9920,0.9920,0.9976,0.9976,0.9814,0.9951,0.9896,0.9962,0.9872,0.9990);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(194151,194532,0.9925,0.9925,0.9981,0.9981,0.9826,0.9963,0.9914,0.9980,0.9874,0.9993);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(194533,195113,0.9927,0.9927,0.9983,0.9983,0.9844,0.9981,0.9934,0.9999,0.9878,0.9996);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(195114,195915,0.9929,0.9929,0.9984,0.9984,0.9838,0.9974,0.9942,1.0007,0.9878,0.9997);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(195916,198115,0.9919,0.9919,0.9975,0.9975,0.9827,0.9964,0.9952,1.0017,0.9869,0.9987);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(198116,199803,0.9955,0.9955,1.0011,1.0011,0.9859,0.9995,0.9893,0.9959,0.9923,1.0041);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(199804,200048,0.9967,0.9967,1.0023,1.0023,0.9870,1.0006,0.9893,0.9959,0.9937,1.0055);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(200049,200151,0.9980,0.9980,1.0036,1.0036,0.9877,1.0012,0.9910,0.9976,0.9980,1.0097);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(200152,200490,0.9958,0.9958,1.0013,1.0013,0.9868,1.0004,0.9922,0.9988,0.9948,1.0065);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(200491,200531,0.9979,0.9979,1.0035,1.0035,0.9876,1.0012,0.9915,0.9981,0.9979,1.0096);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(200532,201656,0.9961,0.9961,1.0017,1.0017,0.9860,0.9996,0.9904,0.9970,0.9945,1.0063);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(201657,202305,0.9969,0.9969,1.0025,1.0025,0.9866,1.0002,0.9914,0.9980,0.9999,1.0116);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(202305,203002,0.9982,0.9982,1.0038,1.0038,0.9872,1.0008,0.9934,1.0000,1.0018,1.0135);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(203003,203984,1.0006,1.0006,1.0061,1.0061,0.9880,1.0017,0.9919,0.9988,0.9992,1.0104);	
  photcicnoeleveto->AddEnCorrPerRun2012HCP(203985,205085,0.9993,0.9993,1.0048,1.0048,0.9903,1.0040,0.9928,0.9997,0.9987,1.0099);	
  photcicnoeleveto->AddEnCorrPerRun2012HCP(205086,205310,1.0004,1.0004,1.0059,1.0059,0.9901,1.0037,0.9987,1.0055,1.0091,1.0202);	
  photcicnoeleveto->AddEnCorrPerRun2012HCP(205311,206207,1.0000,1.0000,1.0055,1.0055,0.9891,1.0028,0.9948,1.0017,1.0032,1.0144);	
  photcicnoeleveto->AddEnCorrPerRun2012HCP(206208,206483,1.0003,1.0003,1.0058,1.0058,0.9895,1.0032,0.9921,0.9989,1.0056,1.0167);	
  photcicnoeleveto->AddEnCorrPerRun2012HCP(206484,206597,1.0005,1.0005,1.0060,1.0060,0.9895,1.0032,0.9968,1.0036,1.0046,1.0158);	
  photcicnoeleveto->AddEnCorrPerRun2012HCP(206598,206896,1.0006,1.0006,1.0061,1.0061,0.9881,1.0017,0.9913,0.9982,1.0050,1.0162);	
  photcicnoeleveto->AddEnCorrPerRun2012HCP(206897,207220,1.0006,1.0006,1.0061,1.0061,0.9884,1.0021,0.9909,0.9978,1.0053,1.0165);	
  photcicnoeleveto->AddEnCorrPerRun2012HCP(207221,208686,1.0006,1.0006,1.0061,1.0061,0.9894,1.0030,0.9951,1.0020,1.0060,1.0172);	

  //-----------------------------------------------------------------------------------------------------------------
  
  //photcicnoeleveto->SetDoMCR9Scaling(kTRUE);
  //photcicnoeleveto->SetMCR9Scale(1.0035, 1.0035);
  photcicnoeleveto->SetDoShowerShapeScaling(kTRUE);
  photcicnoeleveto->SetShowerShapeType("2012ShowerShape");
  //photcicnoeleveto->SetDoMCErrScaling(kTRUE);
  //photcicnoeleveto->SetMCErrScale(1.07, 1.045);    
  //photcicnoeleveto->SetMCErrScale(1, 1);    
  photcicnoeleveto->SetApplyEleVeto(kFALSE);
  photcicnoeleveto->SetInvertElectronVeto(kTRUE);
  photcicnoeleveto->SetJetsName(jetCorr->GetOutputName());  
  //photcicnoeleveto->SetRescaledBeamspotWidth(5.0);
  photcicnoeleveto->SetIsData(isData);
  photcicnoeleveto->SetApplyLeptonTag(kTRUE);
  photcicnoeleveto->SetLeptonTagElectronsName("HggLeptonTagElectrons");
  photcicnoeleveto->SetLeptonTagMuonsName("HggLeptonTagMuons");  
  photcicnoeleveto->Set2012HCP(kTRUE);  
  photcicnoeleveto->SetIdMVAType("2012IdMVA_globe");
 
  PhotonPairSelector         *photpresel = new PhotonPairSelector("PhotonPairSelectorPresel");
  photpresel->SetOutputName("GoodPhotonsPresel");
  photpresel->SetPhotonSelType("MITPFSelectionNoEcal");
  photpresel->SetVertexSelType("CiCMVA2012Selection");
  photpresel->SetUseSingleLegConversions(kFALSE);
  photpresel->SetIdMVAType("2012IdMVA_globe");
  photpresel->DoMCSmear(kTRUE);
  photpresel->DoDataEneCorr(kTRUE);
  photpresel->SetPhotonsFromBranch(kFALSE);
  photpresel->SetInputPhotonsName(photreg->GetOutputName());
  //------------------------------------------2012 HCP--------------------------------------------------------------
  photpresel->SetMCSmearFactors2012HCP(0.0111,0.0111,0.0107,0.0107,0.0155,0.0194,0.0295,0.0276,0.037,0.0371);
  photpresel->SetMCSmearFactors2012HCPMVA(0.01,0.011,0.0106,0.0106,0.0186,0.0196,0.0283,0.0267,0.0343,0.0345);
  photpresel->UseSpecialSmearForDPMVA(true);
  photpresel->AddEnCorrPerRun2012HCP(190645,190781,0.9964,0.9964,1.0020,1.0020,0.9893,1.0028,0.9871,0.9937,0.9839,0.9958);
  photpresel->AddEnCorrPerRun2012HCP(190782,191042,1.0024,1.0024,1.0079,1.0079,0.9923,1.0058,0.9911,0.9977,0.9886,1.0005);
  photpresel->AddEnCorrPerRun2012HCP(191043,193555,0.9935,0.9935,0.9991,0.9991,0.9861,0.9997,0.9894,0.9960,0.9864,0.9982);
  photpresel->AddEnCorrPerRun2012HCP(193556,194150,0.9920,0.9920,0.9976,0.9976,0.9814,0.9951,0.9896,0.9962,0.9872,0.9990);
  photpresel->AddEnCorrPerRun2012HCP(194151,194532,0.9925,0.9925,0.9981,0.9981,0.9826,0.9963,0.9914,0.9980,0.9874,0.9993);
  photpresel->AddEnCorrPerRun2012HCP(194533,195113,0.9927,0.9927,0.9983,0.9983,0.9844,0.9981,0.9934,0.9999,0.9878,0.9996);
  photpresel->AddEnCorrPerRun2012HCP(195114,195915,0.9929,0.9929,0.9984,0.9984,0.9838,0.9974,0.9942,1.0007,0.9878,0.9997);
  photpresel->AddEnCorrPerRun2012HCP(195916,198115,0.9919,0.9919,0.9975,0.9975,0.9827,0.9964,0.9952,1.0017,0.9869,0.9987);
  photpresel->AddEnCorrPerRun2012HCP(198116,199803,0.9955,0.9955,1.0011,1.0011,0.9859,0.9995,0.9893,0.9959,0.9923,1.0041);
  photpresel->AddEnCorrPerRun2012HCP(199804,200048,0.9967,0.9967,1.0023,1.0023,0.9870,1.0006,0.9893,0.9959,0.9937,1.0055);
  photpresel->AddEnCorrPerRun2012HCP(200049,200151,0.9980,0.9980,1.0036,1.0036,0.9877,1.0012,0.9910,0.9976,0.9980,1.0097);
  photpresel->AddEnCorrPerRun2012HCP(200152,200490,0.9958,0.9958,1.0013,1.0013,0.9868,1.0004,0.9922,0.9988,0.9948,1.0065);
  photpresel->AddEnCorrPerRun2012HCP(200491,200531,0.9979,0.9979,1.0035,1.0035,0.9876,1.0012,0.9915,0.9981,0.9979,1.0096);
  photpresel->AddEnCorrPerRun2012HCP(200532,201656,0.9961,0.9961,1.0017,1.0017,0.9860,0.9996,0.9904,0.9970,0.9945,1.0063);
  photpresel->AddEnCorrPerRun2012HCP(201657,202305,0.9969,0.9969,1.0025,1.0025,0.9866,1.0002,0.9914,0.9980,0.9999,1.0116);
  photpresel->AddEnCorrPerRun2012HCP(202305,203002,0.9982,0.9982,1.0038,1.0038,0.9872,1.0008,0.9934,1.0000,1.0018,1.0135);
  photpresel->AddEnCorrPerRun2012HCP(203003,203984,1.0006,1.0006,1.0061,1.0061,0.9880,1.0017,0.9919,0.9988,0.9992,1.0104);	
  photpresel->AddEnCorrPerRun2012HCP(203985,205085,0.9993,0.9993,1.0048,1.0048,0.9903,1.0040,0.9928,0.9997,0.9987,1.0099);	
  photpresel->AddEnCorrPerRun2012HCP(205086,205310,1.0004,1.0004,1.0059,1.0059,0.9901,1.0037,0.9987,1.0055,1.0091,1.0202);	
  photpresel->AddEnCorrPerRun2012HCP(205311,206207,1.0000,1.0000,1.0055,1.0055,0.9891,1.0028,0.9948,1.0017,1.0032,1.0144);	
  photpresel->AddEnCorrPerRun2012HCP(206208,206483,1.0003,1.0003,1.0058,1.0058,0.9895,1.0032,0.9921,0.9989,1.0056,1.0167);	
  photpresel->AddEnCorrPerRun2012HCP(206484,206597,1.0005,1.0005,1.0060,1.0060,0.9895,1.0032,0.9968,1.0036,1.0046,1.0158);	
  photpresel->AddEnCorrPerRun2012HCP(206598,206896,1.0006,1.0006,1.0061,1.0061,0.9881,1.0017,0.9913,0.9982,1.0050,1.0162);	
  photpresel->AddEnCorrPerRun2012HCP(206897,207220,1.0006,1.0006,1.0061,1.0061,0.9884,1.0021,0.9909,0.9978,1.0053,1.0165);	
  photpresel->AddEnCorrPerRun2012HCP(207221,208686,1.0006,1.0006,1.0061,1.0061,0.9894,1.0030,0.9951,1.0020,1.0060,1.0172);	

  //-----------------------------------------------------------------------------------------------------------------
  //photpresel->SetDoMCR9Scaling(kTRUE);
  photpresel->SetDoShowerShapeScaling(kTRUE);
  photpresel->SetShowerShapeType("2012ShowerShape");
  //photpresel->SetDoMCErrScaling(kTRUE);
  //photpresel->SetMCErrScale(1.07, 1.045); 
  //photpresel->SetMCErrScale(1, 1); 
  photpresel->SetJetsName(jetCorr->GetOutputName()); 
  //photpresel->SetRescaledBeamspotWidth(5.0);
  photpresel->SetIsData(isData);
  photpresel->SetApplyLeptonTag(kTRUE);
  photpresel->SetLeptonTagElectronsName("HggLeptonTagElectrons");
  photpresel->SetLeptonTagMuonsName("HggLeptonTagMuons");  
  photpresel->Set2012HCP(kTRUE); 
 
  PhotonPairSelector         *photpreselinverteleveto = new PhotonPairSelector("PhotonPairSelectorPreselInvertEleVeto");
  photpreselinverteleveto->SetOutputName("GoodPhotonsPreselInvertEleVeto");
  photpreselinverteleveto->SetOutputVtxName("OutVtxPreselInvertEleVeto");    
  photpreselinverteleveto->SetPhotonSelType("MITPFSelectionNoEcal");
  photpreselinverteleveto->SetIdMVAType("2012IdMVA_globe");
  //photpreselinverteleveto->SetVertexSelType("CiCMVA2012Selection");//ming change
  photpreselinverteleveto->SetUseSingleLegConversions(kFALSE);
  photpreselinverteleveto->DoMCSmear(kTRUE);
  photpreselinverteleveto->DoDataEneCorr(kTRUE);
  photpreselinverteleveto->SetPhotonsFromBranch(kFALSE);
  photpreselinverteleveto->SetInputPhotonsName(photreg->GetOutputName());
  //------------------------------------------2012 HCP--------------------------------------------------------------
  photpreselinverteleveto->SetMCSmearFactors2012HCP(0.0111,0.0111,0.0107,0.0107,0.0155,0.0194,0.0295,0.0276,0.037,0.0371);
  photpreselinverteleveto->SetMCSmearFactors2012HCPMVA(0.01,0.011,0.0106,0.0106,0.0186,0.0196,0.0283,0.0267,0.0343,0.0345);
  photpreselinverteleveto->UseSpecialSmearForDPMVA(true);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(190645,190781,0.9964,0.9964,1.0020,1.0020,0.9893,1.0028,0.9871,0.9937,0.9839,0.9958);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(190782,191042,1.0024,1.0024,1.0079,1.0079,0.9923,1.0058,0.9911,0.9977,0.9886,1.0005);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(191043,193555,0.9935,0.9935,0.9991,0.9991,0.9861,0.9997,0.9894,0.9960,0.9864,0.9982);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(193556,194150,0.9920,0.9920,0.9976,0.9976,0.9814,0.9951,0.9896,0.9962,0.9872,0.9990);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(194151,194532,0.9925,0.9925,0.9981,0.9981,0.9826,0.9963,0.9914,0.9980,0.9874,0.9993);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(194533,195113,0.9927,0.9927,0.9983,0.9983,0.9844,0.9981,0.9934,0.9999,0.9878,0.9996);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(195114,195915,0.9929,0.9929,0.9984,0.9984,0.9838,0.9974,0.9942,1.0007,0.9878,0.9997);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(195916,198115,0.9919,0.9919,0.9975,0.9975,0.9827,0.9964,0.9952,1.0017,0.9869,0.9987);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(198116,199803,0.9955,0.9955,1.0011,1.0011,0.9859,0.9995,0.9893,0.9959,0.9923,1.0041);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(199804,200048,0.9967,0.9967,1.0023,1.0023,0.9870,1.0006,0.9893,0.9959,0.9937,1.0055);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(200049,200151,0.9980,0.9980,1.0036,1.0036,0.9877,1.0012,0.9910,0.9976,0.9980,1.0097);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(200152,200490,0.9958,0.9958,1.0013,1.0013,0.9868,1.0004,0.9922,0.9988,0.9948,1.0065);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(200491,200531,0.9979,0.9979,1.0035,1.0035,0.9876,1.0012,0.9915,0.9981,0.9979,1.0096);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(200532,201656,0.9961,0.9961,1.0017,1.0017,0.9860,0.9996,0.9904,0.9970,0.9945,1.0063);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(201657,202305,0.9969,0.9969,1.0025,1.0025,0.9866,1.0002,0.9914,0.9980,0.9999,1.0116);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(202305,203002,0.9982,0.9982,1.0038,1.0038,0.9872,1.0008,0.9934,1.0000,1.0018,1.0135);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(203003,203984,1.0006,1.0006,1.0061,1.0061,0.9880,1.0017,0.9919,0.9988,0.9992,1.0104);	
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(203985,205085,0.9993,0.9993,1.0048,1.0048,0.9903,1.0040,0.9928,0.9997,0.9987,1.0099);	
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(205086,205310,1.0004,1.0004,1.0059,1.0059,0.9901,1.0037,0.9987,1.0055,1.0091,1.0202);	
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(205311,206207,1.0000,1.0000,1.0055,1.0055,0.9891,1.0028,0.9948,1.0017,1.0032,1.0144);	
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(206208,206483,1.0003,1.0003,1.0058,1.0058,0.9895,1.0032,0.9921,0.9989,1.0056,1.0167);	
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(206484,206597,1.0005,1.0005,1.0060,1.0060,0.9895,1.0032,0.9968,1.0036,1.0046,1.0158);	
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(206598,206896,1.0006,1.0006,1.0061,1.0061,0.9881,1.0017,0.9913,0.9982,1.0050,1.0162);	
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(206897,207220,1.0006,1.0006,1.0061,1.0061,0.9884,1.0021,0.9909,0.9978,1.0053,1.0165);	
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(207221,208686,1.0006,1.0006,1.0061,1.0061,0.9894,1.0030,0.9951,1.0020,1.0060,1.0172);

  //-----------------------------------------------------------------------------------------------------------------
  photpreselinverteleveto->SetShowerShapeType("2012ShowerShape");
  photpreselinverteleveto->SetDoShowerShapeScaling(kTRUE);  
  //photpreselinverteleveto->SetDoMCErrScaling(kTRUE);
  //photpreselinverteleveto->SetMCErrScale(1.07, 1.045);  
  //photpreselinverteleveto->SetMCErrScale(1, 1);  
  photpreselinverteleveto->SetApplyEleVeto(kFALSE);
  photpreselinverteleveto->SetInvertElectronVeto(kTRUE);
  photpreselinverteleveto->SetJetsName(jetCorr->GetOutputName());  
  //photpreselinverteleveto->SetRescaledBeamspotWidth(5.0);   
  photpreselinverteleveto->SetIsData(isData);
  photpreselinverteleveto->SetApplyLeptonTag(kTRUE);
  photpreselinverteleveto->SetLeptonTagElectronsName("HggLeptonTagElectrons");
  photpreselinverteleveto->SetLeptonTagMuonsName("HggLeptonTagMuons");    
  photpreselinverteleveto->Set2012HCP(kTRUE);  
  photpreselinverteleveto->SetLeadingPtMin(20.);
  photpreselinverteleveto->SetTrailingPtMin(20.);

   PhotonPairSelector         *photpreselnosmear = new PhotonPairSelector("PhotonPairSelectorPreselNoSmear");
   photpreselnosmear->SetOutputName("GoodPhotonsPreselNoSmear");
   photpreselnosmear->SetPhotonSelType("MITPFSelectionNoEcal");
   photpreselnosmear->SetVertexSelType("CiCMVA2012Selection");
   photpreselnosmear->SetUseSingleLegConversions(kFALSE);  
   photpreselnosmear->SetIdMVAType("2012IdMVA_globe");
   photpreselnosmear->SetShowerShapeType("2012ShowerShape");
   photpreselnosmear->SetDoShowerShapeScaling(kTRUE);  
   photpreselnosmear->SetPhotonsFromBranch(kFALSE);
   photpreselnosmear->SetInputPhotonsName(photreg->GetOutputName());
   photpreselnosmear->SetJetsName(jetCorr->GetOutputName());  
   photpreselnosmear->SetOutputVtxName("OutVtxNoSmear");  
   photpreselnosmear->SetLeadingPtMin(30.);
   photpreselnosmear->SetTrailingPtMin(22.);
   //photpreselnosmear->SetRescaledBeamspotWidth(5.0);      
   photpreselnosmear->SetIsData(isData);  
   photpreselnosmear->SetApplyLeptonTag(kTRUE);
   photpreselnosmear->SetLeptonTagElectronsName("HggLeptonTagElectrons");
   photpreselnosmear->SetLeptonTagMuonsName("HggLeptonTagMuons");
   photpreselnosmear->Set2012HCP(kTRUE); 
   photpreselnosmear->DoMCSmear(kFALSE);
   photpreselnosmear->DoDataEneCorr(kFALSE);

   PhotonPairSelector         *photcicnosmear = new PhotonPairSelector("PhotonPairSelectorCiCNoSmear");
   photcicnosmear->SetOutputName("GoodPhotonsCICNoSmear");
   photcicnosmear->SetOutputVtxName("OutVtxCiCNoSmear");        
   photcicnosmear->SetPhotonSelType("CiCPFSelection");
   photcicnosmear->SetVertexSelType("CiCMVA2012Selection");
   photcicnosmear->SetUseSingleLegConversions(kFALSE);
   photcicnosmear->DoMCSmear(kFALSE);
   photcicnosmear->DoDataEneCorr(kFALSE);
   photcicnosmear->SetPhotonsFromBranch(kFALSE);
   photcicnosmear->SetInputPhotonsName(photreg->GetOutputName());
   //  photcicnosmear->SetDoMCR9Scaling(kTRUE);
   //  photcicnosmear->SetMCR9Scale(1.0035, 1.0035);
   photcicnosmear->SetDoShowerShapeScaling(kTRUE); 
   photcicnosmear->SetShowerShapeType("2012ShowerShape");
   //photcicnosmear->SetDoMCErrScaling(kTRUE);
   //photcicnosmear->SetMCErrScale(1.07, 1.045); 
   //photcicnosmear->SetMCErrScale(1, 1); //ming:scale(sigE/E)
   photcicnosmear->SetJetsName(jetCorr->GetOutputName());    
   //photcicnosmear->SetRescaledBeamspotWidth(5.0);
   photcicnosmear->SetIsData(isData);
   photcicnosmear->SetApplyLeptonTag(kTRUE);
   photcicnosmear->SetLeptonTagElectronsName("HggLeptonTagElectrons");
   photcicnosmear->SetLeptonTagMuonsName("HggLeptonTagMuons");  
   photcicnosmear->Set2012HCP(kTRUE);
   
   PhotonPairSelector         *photpreselinvertelevetonosmear = new PhotonPairSelector("PhotonPairSelectorPreselInvertEleVetoNoSmear");
   photpreselinvertelevetonosmear->SetOutputName("GoodPhotonsPreselInvertEleVetoNoSmear");
   photpreselinvertelevetonosmear->SetPhotonSelType("MITPFSelectionNoEcal");
   //photpreselinvertelevetonosmear->SetVertexSelType("CiCMVA2012Selection");//ming change
   photpreselinvertelevetonosmear->SetUseSingleLegConversions(kFALSE);  
   photpreselinvertelevetonosmear->SetIdMVAType("2012IdMVA_globe");
   photpreselinvertelevetonosmear->SetShowerShapeType("2012ShowerShape");
   photpreselinvertelevetonosmear->SetDoShowerShapeScaling(kTRUE);  
   photpreselinvertelevetonosmear->SetPhotonsFromBranch(kFALSE);
   photpreselinvertelevetonosmear->SetInputPhotonsName(photreg->GetOutputName());
   photpreselinvertelevetonosmear->SetJetsName(jetCorr->GetOutputName());  
   photpreselinvertelevetonosmear->SetOutputVtxName("OutVtxPreselInvertEleVetoNoSmear");  
   //photpreselinvertelevetonosmear->SetLeadingPtMin(30.);
   photpreselinvertelevetonosmear->SetLeadingPtMin(20.);//ming
   //photpreselinvertelevetonosmear->SetTrailingPtMin(22.);
   photpreselinvertelevetonosmear->SetTrailingPtMin(20.);//ming
   //photpreselinvertelevetonosmear->SetRescaledBeamspotWidth(5.0);      
   photpreselinvertelevetonosmear->SetIsData(isData);  
   photpreselinvertelevetonosmear->SetApplyLeptonTag(kTRUE);
   photpreselinvertelevetonosmear->SetLeptonTagElectronsName("HggLeptonTagElectrons");
   photpreselinvertelevetonosmear->SetLeptonTagMuonsName("HggLeptonTagMuons");
   photpreselinvertelevetonosmear->Set2012HCP(kTRUE); 
   photpreselinvertelevetonosmear->DoMCSmear(kFALSE);
   photpreselinvertelevetonosmear->DoDataEneCorr(kFALSE);
   photpreselinvertelevetonosmear->SetApplyEleVeto(kFALSE);
   photpreselinvertelevetonosmear->SetInvertElectronVeto(kTRUE);
   
   PhotonTreeWriter *phottreecic = new PhotonTreeWriter("PhotonTreeWriterCiC");
   phottreecic->SetPhotonsFromBranch(kFALSE);
   phottreecic->SetInputPhotonsName(photcic->GetOutputName());
   phottreecic->SetEnableJets(kTRUE);
   phottreecic->SetApplyJetId(kTRUE);
   phottreecic->SetPFJetsFromBranch(kFALSE);
   phottreecic->SetPFJetName(jetCorr->GetOutputName());
   phottreecic->SetExcludeDoublePrompt(excludedoubleprompt);
   phottreecic->SetIsData(isData);
   if (is25) phottreecic->SetEnablePFPhotons(kFALSE);
   phottreecic->SetApplyLeptonTag(kTRUE);
   phottreecic->SetLeptonTagElectronsName("HggLeptonTagElectrons");
   phottreecic->SetLeptonTagMuonsName("HggLeptonTagMuons");
   phottreecic->SetApplyVBFTag(kTRUE);
   phottreecic->SetApplyPFMetCorr(kTRUE);
   phottreecic->SetBeamspotWidth(5.0);
   
   PhotonTreeWriter *phottreecicnoeleveto = new PhotonTreeWriter("PhotonTreeWriterCiCInvertEleVeto");
   phottreecicnoeleveto->SetPhotonsFromBranch(kFALSE);
   phottreecicnoeleveto->SetInputPhotonsName(photcicnoeleveto->GetOutputName());
   phottreecicnoeleveto->SetEnableJets(kTRUE);
   phottreecicnoeleveto->SetApplyJetId(kTRUE);  
   phottreecicnoeleveto->SetPFJetsFromBranch(kFALSE);
   phottreecicnoeleveto->SetPFJetName(jetCorr->GetOutputName());
   phottreecicnoeleveto->SetApplyElectronVeto(kFALSE);
   phottreecicnoeleveto->SetExcludeDoublePrompt(excludedoubleprompt);
   phottreecicnoeleveto->SetIsData(isData);  
   if (is25) phottreecicnoeleveto->SetEnablePFPhotons(kFALSE);
   phottreecicnoeleveto->SetApplyLeptonTag(kTRUE);
   phottreecicnoeleveto->SetLeptonTagElectronsName("HggLeptonTagElectrons");
   phottreecicnoeleveto->SetLeptonTagMuonsName("HggLeptonTagMuons");
   phottreecicnoeleveto->SetApplyVBFTag(kTRUE);
   phottreecicnoeleveto->SetApplyPFMetCorr(kTRUE);
   phottreecicnoeleveto->SetBeamspotWidth(5.0);
   
   PhotonTreeWriter *phottreepresel = new PhotonTreeWriter("PhotonTreeWriterPresel");
   phottreepresel->SetPhotonsFromBranch(kFALSE);
   phottreepresel->SetInputPhotonsName(photpresel->GetOutputName());
   phottreepresel->SetEnableJets(kTRUE);
   phottreepresel->SetApplyJetId(kTRUE);    
   phottreepresel->SetPFJetsFromBranch(kFALSE);
   phottreepresel->SetPFJetName(jetCorr->GetOutputName());  
   phottreepresel->SetExcludeDoublePrompt(excludedoubleprompt);  
   phottreepresel->SetIsData(isData);  
   if (is25) phottreepresel->SetEnablePFPhotons(kFALSE);
   phottreepresel->SetApplyLeptonTag(kTRUE);
   phottreepresel->SetLeptonTagElectronsName("HggLeptonTagElectrons");
   phottreepresel->SetLeptonTagMuonsName("HggLeptonTagMuons");
   phottreepresel->SetApplyVBFTag(kTRUE);
   phottreepresel->SetApplyPFMetCorr(kTRUE);
   phottreepresel->SetBeamspotWidth(5.0);
   phottreepresel->SetElectronMVAWeightsSubdet0Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat1.weights.xml"))));
   phottreepresel->SetElectronMVAWeightsSubdet1Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat2.weights.xml"))));  
   phottreepresel->SetElectronMVAWeightsSubdet2Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat3.weights.xml"))));  
   phottreepresel->SetElectronMVAWeightsSubdet0Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat4.weights.xml")))); 
   phottreepresel->SetElectronMVAWeightsSubdet1Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat5.weights.xml")))); 
   phottreepresel->SetElectronMVAWeightsSubdet2Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat6.weights.xml"))));
   phottreepresel->SetDoSynching(kTRUE);
    
   PhotonTreeWriter *phottreepreselinverteleveto = new PhotonTreeWriter("PhotonTreeWriterPreselInvertEleVeto");
   phottreepreselinverteleveto->SetPhotonsFromBranch(kFALSE);
   phottreepreselinverteleveto->SetInputPhotonsName(photpreselinverteleveto->GetOutputName());
   phottreepreselinverteleveto->SetEnableJets(kTRUE);
   phottreepreselinverteleveto->SetApplyJetId(kTRUE);    
   phottreepreselinverteleveto->SetPFJetsFromBranch(kFALSE);
   phottreepreselinverteleveto->SetPFJetName(jetCorr->GetOutputName());  
   phottreepreselinverteleveto->SetApplyElectronVeto(kFALSE);  
   phottreepreselinverteleveto->SetExcludeDoublePrompt(excludedoubleprompt);    
   phottreepreselinverteleveto->SetIsData(isData); 
   if (is25) phottreepreselinverteleveto->SetEnablePFPhotons(kFALSE);  
   phottreepreselinverteleveto->SetApplyLeptonTag(kTRUE);
   phottreepreselinverteleveto->SetLeptonTagElectronsName("HggLeptonTagElectrons");
   phottreepreselinverteleveto->SetLeptonTagMuonsName("HggLeptonTagMuons");
   phottreepreselinverteleveto->SetApplyVBFTag(kTRUE);
   phottreepreselinverteleveto->SetApplyPFMetCorr(kTRUE);
   phottreepreselinverteleveto->SetBeamspotWidth(5.0);
  
   PhotonTreeWriter *phottreepreselnosmear = new PhotonTreeWriter("PhotonTreeWriterPreselNoSmear");
   phottreepreselnosmear->SetPhotonsFromBranch(kFALSE);
   phottreepreselnosmear->SetInputPhotonsName(photpreselnosmear->GetOutputName());
   phottreepreselnosmear->SetEnableJets(kTRUE);
   phottreepreselnosmear->SetApplyJetId(kTRUE);  
   phottreepreselnosmear->SetPFJetsFromBranch(kFALSE);
   phottreepreselnosmear->SetPFJetName(jetCorr->GetOutputName());  
   phottreepreselnosmear->SetExcludeDoublePrompt(excludedoubleprompt);  
   phottreepreselnosmear->SetIsData(isData);    
   if (is25) phottreepreselnosmear->SetEnablePFPhotons(kFALSE); 
   phottreepreselnosmear->SetApplyLeptonTag(kTRUE);
   phottreepreselnosmear->SetLeptonTagElectronsName("HggLeptonTagElectrons");
   phottreepreselnosmear->SetLeptonTagMuonsName("HggLeptonTagMuons");
   phottreepreselnosmear->SetApplyVBFTag(kTRUE);
   phottreepreselnosmear->SetApplyPFMetCorr(kTRUE);
   phottreepreselnosmear->SetBeamspotWidth(5.0); 

   PhotonTreeWriter *phottreecicnosmear = new PhotonTreeWriter("PhotonTreeWriterCiCNoSmear");
   phottreecicnosmear->SetPhotonsFromBranch(kFALSE);
   phottreecicnosmear->SetInputPhotonsName(photcicnosmear->GetOutputName());
   phottreecicnosmear->SetEnableJets(kTRUE);
   phottreecicnosmear->SetApplyJetId(kTRUE);
   phottreecicnosmear->SetPFJetsFromBranch(kFALSE);
   phottreecicnosmear->SetPFJetName(jetCorr->GetOutputName());
   phottreecicnosmear->SetExcludeDoublePrompt(excludedoubleprompt);
   phottreecicnosmear->SetIsData(isData);
   if (is25) phottreecicnosmear->SetEnablePFPhotons(kFALSE);
   phottreecicnosmear->SetApplyLeptonTag(kTRUE);
   phottreecicnosmear->SetLeptonTagElectronsName("HggLeptonTagElectrons");
   phottreecicnosmear->SetLeptonTagMuonsName("HggLeptonTagMuons");
   phottreecicnosmear->SetApplyVBFTag(kTRUE);
   phottreecicnosmear->SetApplyPFMetCorr(kTRUE);
   phottreecicnosmear->SetBeamspotWidth(5.0);

   PhotonTreeWriter *phottreepreselinvertelevetonosmear = new PhotonTreeWriter("PhotonTreeWriterPreselInvertEleVetoNoSmear");
   phottreepreselinvertelevetonosmear->SetPhotonsFromBranch(kFALSE);
   phottreepreselinvertelevetonosmear->SetInputPhotonsName(photpreselinvertelevetonosmear->GetOutputName());
   phottreepreselinvertelevetonosmear->SetEnableJets(kTRUE);
   phottreepreselinvertelevetonosmear->SetApplyJetId(kTRUE);  
   phottreepreselinvertelevetonosmear->SetPFJetsFromBranch(kFALSE);
   phottreepreselinvertelevetonosmear->SetPFJetName(jetCorr->GetOutputName());  
   phottreepreselinvertelevetonosmear->SetExcludeDoublePrompt(excludedoubleprompt);  
   phottreepreselinvertelevetonosmear->SetIsData(isData);    
   if (is25) phottreepreselinvertelevetonosmear->SetEnablePFPhotons(kFALSE); 
   phottreepreselinvertelevetonosmear->SetApplyLeptonTag(kTRUE);
   phottreepreselinvertelevetonosmear->SetLeptonTagElectronsName("HggLeptonTagElectrons");
   phottreepreselinvertelevetonosmear->SetLeptonTagMuonsName("HggLeptonTagMuons");
   phottreepreselinvertelevetonosmear->SetApplyVBFTag(kTRUE);
   phottreepreselinvertelevetonosmear->SetApplyPFMetCorr(kTRUE);
   phottreepreselinvertelevetonosmear->SetBeamspotWidth(5.0); 
   phottreepreselinvertelevetonosmear->SetApplyElectronVeto(kFALSE); 
   
   PhotonIDMod         *photidpresel = new PhotonIDMod("PhotonIDModPresel");
   photidpresel->SetPtMin(25.0);
   photidpresel->SetOutputName("GoodPhotonsPreselid");
   photidpresel->SetIDType("MITPFSelection");
   photidpresel->SetApplyElectronVeto(kTRUE);
   photidpresel->SetIsData(isData);
   
  PhotonIDMod         *photidpreselinvert = new PhotonIDMod("PhotonIDModPreselInvert");
  photidpreselinvert->SetPtMin(25.0);
  photidpreselinvert->SetOutputName("GoodPhotonsPreselidInvert");
  photidpreselinvert->SetIDType("MITPFSelection");
  photidpreselinvert->SetApplyElectronVeto(kFALSE);
  photidpreselinvert->SetInvertElectronVeto(kTRUE);
  photidpreselinvert->SetIsData(isData);  
  
  PhotonTreeWriter *phottreesingle = new PhotonTreeWriter("PhotonTreeWriterSingle");
  phottreesingle->SetWriteDiphotonTree(kFALSE);
  phottreesingle->SetPhotonsFromBranch(kFALSE);
  phottreesingle->SetInputPhotonsName(photidpresel->GetOutputName());  
  phottreesingle->SetEnableJets(kTRUE);
  phottreesingle->SetPFJetsFromBranch(kFALSE);
  phottreesingle->SetPFJetName(jetCorr->GetOutputName());      
  phottreesingle->SetBeamspotWidth(5.0);
  phottreesingle->SetIsData(isData);
  if (is25) phottreesingle->SetEnablePFPhotons(kFALSE);  
  
  PhotonTreeWriter *phottreesingleinvert = new PhotonTreeWriter("PhotonTreeWriterSingleInvert");
  phottreesingleinvert->SetWriteDiphotonTree(kTRUE);
  phottreesingleinvert->SetPhotonsFromBranch(kFALSE);
  phottreesingleinvert->SetInputPhotonsName(photidpreselinvert->GetOutputName());  
  phottreesingleinvert->SetEnableJets(kTRUE);
  phottreesingleinvert->SetPFJetsFromBranch(kFALSE);
  phottreesingleinvert->SetPFJetName(pubJetOpen->GetOutputName());
  phottreesingleinvert->SetApplyElectronVeto(kFALSE);
  phottreesingleinvert->SetBeamspotWidth(5.0);
  phottreesingleinvert->SetIsData(isData);  
  if (is25) phottreesingleinvert->SetEnablePFPhotons(kFALSE);    
  

  //------------------------------------------------------------------------------------------------
  // making analysis chain
  //------------------------------------------------------------------------------------------------
  // this is how it always starts
  runLumiSel      ->Add(mcselmod);
    
  if (TString(dataset).Contains("-h")) {    
    mcselmod        ->Add(sysMod);
  }
  
  // high level trigger is always first
  //mcselmod         ->Add(hltModES);
  //mcselmod         ->Add(hltModE);
  
  if (!TString(dataset).Contains("meridiani")) {      
    mcselmod         ->Add(hltModP);
    hltModP         ->Add(goodPVFilterMod);
  }
  else {
    mcselmod->Add(goodPVFilterMod);
  }
  //hltModE         ->Add(goodPVFilterModE);
  //hltModES        ->Add(goodPVFilterModES);
  
  //goodPVFilterMod ->Add(muonId);
  goodPVFilterMod->Add(photreg);
  photreg->Add(pubJet);
  pubJet->Add(jetCorr);
  
  // simple object id modules
  //goodPVFilterModE -> Add(pubJetOpen);
  //pubJetOpen       -> Add(photidpreselinvert);
  //goodPVFilterModES ->Add(elecIdS);
  jetCorr          ->Add(SepPUMod); 
  SepPUMod         ->Add(eleIdMod);
  eleIdMod         ->Add(muonIdMod);
   
  muonIdMod          ->Add(photcic);
  muonIdMod          ->Add(photcicnoeleveto);  
  muonIdMod          ->Add(photpresel);  
  muonIdMod          ->Add(photpreselinverteleveto);  
  muonIdMod          ->Add(photpreselnosmear);  
  muonIdMod          ->Add(photpreselinvertelevetonosmear);  
  muonIdMod          ->Add(photcicnosmear);  

  photcic         ->Add(phottreecic);
  photcicnoeleveto       ->Add(phottreecicnoeleveto);
  photpresel    ->Add(phottreepresel);
  photpreselinverteleveto    ->Add(phottreepreselinverteleveto);
  photpreselnosmear    ->Add(phottreepreselnosmear);
  photpreselinvertelevetonosmear    ->Add(phottreepreselinvertelevetonosmear);
  photcicnosmear    ->Add(phottreecicnosmear);

  //jetCorr          ->Add(photidpresel);
  //photidpresel       ->Add(phottreesingle);
  
  //photidpreselinvert       ->Add(phottreesingleinvert);  
  //elecIdS->Add(phottreeES);
  

  //TFile::SetOpenTimeout(0);
  //TFile::SetCacheFileDir("./rootfilecache",kTRUE,kTRUE);
  TFile::SetReadaheadSize(128*1024*1024);
  
  //------------------------------------------------------------------------------------------------
  // setup analysis
  //------------------------------------------------------------------------------------------------
  Analysis *ana = new Analysis;
  ana->SetUseHLT(kTRUE);
  ana->SetKeepHierarchy(kTRUE);
  ana->SetSuperModule(runLumiSel);
  ana->SetPrintScale(100);
  if (nEvents >= 0)
    ana->SetProcessNEvents(nEvents);

  if (TString(dataset).Contains("meridiani")) {    
    ana->SetUseHLT(kFALSE);
  }
  
  //------------------------------------------------------------------------------------------------
  // organize input
  //------------------------------------------------------------------------------------------------
  Catalog *c = new Catalog(catalogDir);
  TString skimdataset = TString(dataset)+TString("/") +TString(skim);
  Dataset *d = NULL;
  TString bookstr = book;
  //if (TString(dataset).Contains("s11-h")) bookstr.ReplaceAll("local","t2mit");
  if (TString(skim).CompareTo("noskim") == 0)
    d = c->FindDataset(bookstr,dataset,fileset);
  else 
    d = c->FindDataset(bookstr,skimdataset.Data(),fileset);
 
  ana->AddDataset(d);
  //ana->AddFile("/mnt/hadoop/cmsprod/filefi/029/r12a-pho-j13-v1/385B77DE-58D0-E111-B925-001E67396928.root");//79737729
  //ana->AddFile("/mnt/hadoop/cmsprod/filefi/029/r12b-dph-j13-v1/FE308E0A-23D2-E111-8B2D-00266CFAE7D0.root");//871378986
  //ana->AddFile("/mnt/hadoop/cmsprod/filefi/029/r12d-dph-pr-v1/54E239DC-2725-E211-A52C-003048D373F6.root");//528937923
  //ana->AddFile("/mnt/hadoop/cmsprod/filefi/029/r12d-dph-pr-v1/EC80B121-7E30-E211-BFFF-5404A640A63D.root");//261543921
  //ana->AddFile("/mnt/hadoop/cmsprod/filefi/029/r12d-dph-pr-v1/66CE9443-472C-E211-B08F-001D09F2437B.root");//483561562
  //ana->AddFile("/mnt/hadoop/cmsprod/filefi/029/r12d-dph-pr-v1/F03B81B2-0B29-E211-B28F-BCAEC518FF67.root");
  //ana->AddFile("/mnt/hadoop/cmsprod/filefi/029/r12d-dph-pr-v1/A6EC91DD-CB39-E211-8897-BCAEC53296F8.root");
  //ana->AddFile("/mnt/hadoop/cmsprod/filefi/029/r12c-dph-pr-v2/64F7EAB0-CBEA-E111-B8AD-003048F117B6.root");
  //ana->AddFile("/mnt/hadoop/cmsprod/filefi/029/r12c-dph-pr-v2/0615B16F-ECDA-E111-A266-003048F117B6.root"); 
  //ana->AddFile("/mnt/hadoop/cmsprod/filefi/029/r12a-pho-j13-v1/385B77DE-58D0-E111-B925-001E67396928.root"); 
  //ana->AddFile("/mnt/hadoop/cmsprod/filefi/029/r12c-dph-pr-v2/068B257D-81D0-E111-A835-5404A63886B4.root");//369441614 
  //ana->AddFile("/mnt/hadoop/cmsprod/filefi/029/r12d-dph-pr-v1/5443DBCB-562C-E211-8042-00237DDBE0E2.root");//333643114
  //ana->AddFile("/mnt/hadoop/cmsprod/filefi/029/r12d-dph-pr-v1/B2698ACF-843A-E211-8D9E-001D09F244DE.root");//89022540
  //ana->AddFile("/mnt/hadoop/cmsprod/filefi/029/r12b-dph-j13-v1/1E050D78-6DD4-E111-81D2-00266CFAE8D0.root");//223740859
  //ana->AddFile("/mnt/hadoop/cmsprod/filefi/029/r12b-dph-j13-v1/2CD64E47-79D3-E111-90C5-00A0D1EE8E60.root");//541603559
  //ana->AddFile("/mnt/hadoop/cmsprod/filefi/029/r12b-dph-j13-v1/147EF338-B2D3-E111-B1A6-00266CF256CC.root");//8983064 
  //ana->AddFile("/mnt/hadoop/cmsprod/filefi/029/r12d-dph-pr-v1/8843BD9F-DC3E-E211-9E2F-003048678110.root");//876316897 
  //ana->AddFile("/mnt/hadoop/cmsprod/filefi/029/r12c-dph-pr-v2/CAA12077-F9F6-E111-8DF2-001D09F28D54.root");
 
  //------------------------------------------------------------------------------------------------
  // organize output
  //------------------------------------------------------------------------------------------------
  TString rootFile = TString(outputName);
  rootFile += TString("_") + TString(dataset) + TString("_") + TString(skim);
  if (TString(fileset) != TString(""))
    rootFile += TString("_") + TString(fileset);
  rootFile += TString(".root");
  ana->SetOutputName(rootFile.Data());
  //ana->SetCacheSize(64*1024*1024);
  ana->SetCacheSize(0);
  
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
