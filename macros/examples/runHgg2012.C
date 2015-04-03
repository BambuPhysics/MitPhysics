// $Id: runHgg2012.C,v 1.14 2012/06/26 23:24:57 bendavid Exp $
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

#endif

//--------------------------------------------------------------------------------------------------
void runHgg2012(const char *fileset    = "0000",
            const char *skim       = "noskim",
	    //const char *dataset = "f11--h120gg-gf-v14b-pu",
            //const char *dataset = "r11a-pho-j16-v1",   
            //const char *dataset = "meridiani2012-diphoj-v9",
	    //const char *dataset = "s12-h125gg-gf-v9",
            //const char *dataset = "s12-pj40-2em-v9",
            const char *dataset = "s12-diphoj-3-v9",
            //const char *dataset = "s12-zllm50-v9",
	    //const char *dataset = "f11--h121gg-gf-v14b-pu",
	    //const char *dataset = "r12a-pho-pr-v1",
	    //const char *dataset = "s12-pj40-2em-v9", 
            const char *book       = "t2mit/filefi/028",
            const char *catalogDir = "/home/cmsprod/catalog",
            const char *outputName = "hgg",
            int         nEvents    = 500)
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

  TString jsonFile = TString("/home/bendavid/cms/json/") + TString(json);
  //TString jsonFile = TString("/home/bendavid/cms/json/") + TString("Cert_136033-149442_7TeV_Dec22ReReco_Collisions10_JSON_v4.txt");
  Bool_t  isData   = ( (jsonFile.CompareTo("/home/bendavid/cms/json/~") != 0) );
  
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
  if ((jsonFile.CompareTo("/home/bendavid/cms/json/~") != 0) &&
      (jsonFile.CompareTo("/home/bendavid/cms/json/-") != 0)   ) {
    runLumiSel->AddJSONFile(jsonFile.Data());
  }
  if ((jsonFile.CompareTo("/home/bendavid/cms/json/-") == 0)   ) {
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

  hltModP->AddTrigger("HLT_Photon26_CaloId10_Iso50_Photon18_CaloId10_Iso50_Mass60_v*",190000,999999); 
  hltModP->AddTrigger("HLT_Photon26_CaloId10_Iso50_Photon18_R9Id85_Mass60_v*",190000,999999); 
  hltModP->AddTrigger("HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass60_v*",190000,999999); 
  hltModP->AddTrigger("HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass70_v*",190000,999999); 
  hltModP->AddTrigger("HLT_Photon26_R9Id85_Photon18_CaloId10_Iso50_Mass60_v*",190000,999999); 
  hltModP->AddTrigger("HLT_Photon26_R9Id85_Photon18_R9Id85_Mass60_v*",190000,999999); 
  hltModP->AddTrigger("HLT_Photon36_CaloId10_Iso50_Photon22_CaloId10_Iso50_v*",190000,999999); 
  hltModP->AddTrigger("HLT_Photon36_CaloId10_Iso50_Photon22_R9Id85_v*",190000,999999); 
  hltModP->AddTrigger("HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_R9Id85_OR_CaloId10_Iso50_v*",190000,999999); 
  hltModP->AddTrigger("HLT_Photon36_R9Id85_Photon22_CaloId10_Iso50_v*",190000,999999); 
  hltModP->AddTrigger("HLT_Photon36_R9Id85_Photon22_R9Id85_v*",190000,999999); 
    
  hltModP->SetTrigObjsName("MyHltPhotObjs");
  hltModP->SetAbortIfNotAccepted(isData);
 
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
  MuonIDMod *muonId = new MuonIDMod;  
  muonId->SetClassType ("Global");
  muonId->SetIDType    ("ZMuId");
  muonId->SetIsoType   ("TrackCaloSliding");
  muonId->SetApplyD0Cut(kTRUE);

  
  PublisherMod<PFJet,Jet> *pubJet = new PublisherMod<PFJet,Jet>("JetPub");
  pubJet->SetInputName("AKt5PFJets");
  pubJet->SetOutputName("PubAKt5PFJets");
  
  PublisherMod<PFJet,Jet> *pubJetOpen = new PublisherMod<PFJet,Jet>("JetPubOpen");
  pubJetOpen->SetInputName("AKt5PFJets");
  pubJetOpen->SetOutputName("PubAKt5PFJetsOpen");  
  
  JetCorrectionMod *jetCorr = new JetCorrectionMod;
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
  //photcic->SetUseSingleLegConversions(kTRUE);
  photcic->DoMCSmear(kTRUE);
  photcic->DoDataEneCorr(kTRUE);
  photcic->SetPhotonsFromBranch(kFALSE);
  photcic->SetInputPhotonsName(photreg->GetOutputName());
  //2012 ICHEP promptreco 5.3/fb
  photcic->SetMCSmearFactors(0.00825,0.01035,0.0116,0.0222,0.0138,0.0292,0.0275,0.0322,0.0310);
  photcic->AddEnCorrPerRun(190450, 190781, 0.9944, 0.9944, 1.0006, 0.9948, 1.0083, 0.9881, 0.9956, 0.9878, 0.9991);
  photcic->AddEnCorrPerRun(190782, 190949, 1.0085, 1.0085, 1.0146, 0.9858, 0.9993, 1.0047, 1.0120, 0.9962, 1.0074);
  photcic->AddEnCorrPerRun(190950, 191833, 0.9939, 0.9939, 1.0001, 0.9879, 1.0014, 0.9906, 0.9981, 0.9874, 0.9986);
  photcic->AddEnCorrPerRun(191834, 193686, 0.9916, 0.9916, 0.9978, 0.9855, 0.9990, 0.9907, 0.9982, 0.9893, 1.0005);
  photcic->AddEnCorrPerRun(193746, 194210, 0.9921, 0.9921, 0.9983, 0.9834, 0.9970, 0.9921, 0.9996, 0.9883, 0.9995);
  photcic->AddEnCorrPerRun(194211, 194479, 0.9926, 0.9926, 0.9988, 0.9845, 0.9980, 0.9928, 1.0003, 0.9895, 1.0007);
  photcic->AddEnCorrPerRun(194480, 195147, 0.9927, 0.9927, 0.9989, 0.9862, 0.9998, 0.9954, 1.0028, 0.9891, 1.0003);
  photcic->AddEnCorrPerRun(195148, 195350, 0.9928, 0.9928, 0.9990, 0.9852, 0.9988, 1.0007, 1.0081, 0.9943, 1.0054);
  photcic->AddEnCorrPerRun(195351, 195395, 0.9931, 0.9931, 0.9992, 0.9861, 0.9996, 0.9933, 1.0008, 0.9900, 1.0012);
  photcic->AddEnCorrPerRun(195396, 195530, 0.9932, 0.9932, 0.9994, 0.9861, 0.9996, 0.9951, 1.0025, 0.9869, 0.9982);
  photcic->AddEnCorrPerRun(195531, 999999, 0.9924, 0.9924, 0.9986, 0.9857, 0.9992, 0.9962, 1.0036, 0.9885, 0.9997);
  
  
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
  
  PhotonPairSelector         *photcicnoeleveto = new PhotonPairSelector("PhotonPairSelectorCiCInvertEleVeto");
  photcicnoeleveto->SetOutputName("GoodPhotonsCICNoEleVeto");
  photcicnoeleveto->SetOutputVtxName("OutVtxCiCInvertEleVeto");      
  photcicnoeleveto->SetPhotonSelType("CiCPFSelection");
  photcicnoeleveto->SetVertexSelType("CiCMVA2012Selection");
  //photcicnoeleveto->SetUseSingleLegConversions(kTRUE);
  photcicnoeleveto->DoMCSmear(kTRUE);
  photcicnoeleveto->DoDataEneCorr(kTRUE);
  photcicnoeleveto->SetPhotonsFromBranch(kFALSE);
  photcicnoeleveto->SetInputPhotonsName(photreg->GetOutputName());
  photcicnoeleveto->SetMCSmearFactors(0.00825,0.01035,0.0116,0.0222,0.0138,0.0292,0.0275,0.0322,0.0310);
  photcicnoeleveto->AddEnCorrPerRun(190450, 190781, 0.9944, 0.9944, 1.0006, 0.9948, 1.0083, 0.9881, 0.9956, 0.9878, 0.9991);
  photcicnoeleveto->AddEnCorrPerRun(190782, 190949, 1.0085, 1.0085, 1.0146, 0.9858, 0.9993, 1.0047, 1.0120, 0.9962, 1.0074);
  photcicnoeleveto->AddEnCorrPerRun(190950, 191833, 0.9939, 0.9939, 1.0001, 0.9879, 1.0014, 0.9906, 0.9981, 0.9874, 0.9986);
  photcicnoeleveto->AddEnCorrPerRun(191834, 193686, 0.9916, 0.9916, 0.9978, 0.9855, 0.9990, 0.9907, 0.9982, 0.9893, 1.0005);
  photcicnoeleveto->AddEnCorrPerRun(193746, 194210, 0.9921, 0.9921, 0.9983, 0.9834, 0.9970, 0.9921, 0.9996, 0.9883, 0.9995);
  photcicnoeleveto->AddEnCorrPerRun(194211, 194479, 0.9926, 0.9926, 0.9988, 0.9845, 0.9980, 0.9928, 1.0003, 0.9895, 1.0007);
  photcicnoeleveto->AddEnCorrPerRun(194480, 195147, 0.9927, 0.9927, 0.9989, 0.9862, 0.9998, 0.9954, 1.0028, 0.9891, 1.0003);
  photcicnoeleveto->AddEnCorrPerRun(195148, 195350, 0.9928, 0.9928, 0.9990, 0.9852, 0.9988, 1.0007, 1.0081, 0.9943, 1.0054);
  photcicnoeleveto->AddEnCorrPerRun(195351, 195395, 0.9931, 0.9931, 0.9992, 0.9861, 0.9996, 0.9933, 1.0008, 0.9900, 1.0012);
  photcicnoeleveto->AddEnCorrPerRun(195396, 195530, 0.9932, 0.9932, 0.9994, 0.9861, 0.9996, 0.9951, 1.0025, 0.9869, 0.9982);
  photcicnoeleveto->AddEnCorrPerRun(195531, 999999, 0.9924, 0.9924, 0.9986, 0.9857, 0.9992, 0.9962, 1.0036, 0.9885, 0.9997);
  //photcicnoeleveto->SetDoMCR9Scaling(kTRUE);
  //photcicnoeleveto->SetMCR9Scale(1.0035, 1.0035);
  photcicnoeleveto->SetDoShowerShapeScaling(kTRUE);
  photcicnoeleveto->SetShowerShapeType("2012ShowerShape");
  //photcicnoeleveto->SetDoMCErrScaling(kTRUE);
  //photcicnoeleveto->SetMCErrScale(1.07, 1.045);    
  //photcicnoeleveto->SetMCErrScale(1, 1);    
  photcicnoeleveto->SetApplyEleVeto(kFALSE);
  photcicnoeleveto->SetInvertElectronVeto(kFALSE);
  photcicnoeleveto->SetJetsName(jetCorr->GetOutputName());  
  //photcicnoeleveto->SetRescaledBeamspotWidth(5.0);
  photcicnoeleveto->SetIsData(isData);
  
  
  PhotonPairSelector         *photpresel = new PhotonPairSelector("PhotonPairSelectorPresel");
  photpresel->SetOutputName("GoodPhotonsPresel");
  photpresel->SetPhotonSelType("MITPFSelection");
  photpresel->SetVertexSelType("CiCMVA2012Selection");
  //photpresel->SetUseSingleLegConversions(kTRUE);
  photpresel->SetIdMVAType("2012IdMVA_globe");
  //photpresel->SetVertexSelType("MetSigSelection");
  photpresel->DoMCSmear(kTRUE);
  photpresel->DoDataEneCorr(kTRUE);
  photpresel->SetPhotonsFromBranch(kFALSE);
  photpresel->SetInputPhotonsName(photreg->GetOutputName());
  photpresel->SetMCSmearFactors(0.00825,0.01035,0.0116,0.0222,0.0138,0.0292,0.0275,0.0322,0.0310);
  photpresel->AddEnCorrPerRun(190450, 190781, 0.9944, 0.9944, 1.0006, 0.9948, 1.0083, 0.9881, 0.9956, 0.9878, 0.9991);
  photpresel->AddEnCorrPerRun(190782, 190949, 1.0085, 1.0085, 1.0146, 0.9858, 0.9993, 1.0047, 1.0120, 0.9962, 1.0074);
  photpresel->AddEnCorrPerRun(190950, 191833, 0.9939, 0.9939, 1.0001, 0.9879, 1.0014, 0.9906, 0.9981, 0.9874, 0.9986);
  photpresel->AddEnCorrPerRun(191834, 193686, 0.9916, 0.9916, 0.9978, 0.9855, 0.9990, 0.9907, 0.9982, 0.9893, 1.0005);
  photpresel->AddEnCorrPerRun(193746, 194210, 0.9921, 0.9921, 0.9983, 0.9834, 0.9970, 0.9921, 0.9996, 0.9883, 0.9995);
  photpresel->AddEnCorrPerRun(194211, 194479, 0.9926, 0.9926, 0.9988, 0.9845, 0.9980, 0.9928, 1.0003, 0.9895, 1.0007);
  photpresel->AddEnCorrPerRun(194480, 195147, 0.9927, 0.9927, 0.9989, 0.9862, 0.9998, 0.9954, 1.0028, 0.9891, 1.0003);
  photpresel->AddEnCorrPerRun(195148, 195350, 0.9928, 0.9928, 0.9990, 0.9852, 0.9988, 1.0007, 1.0081, 0.9943, 1.0054);
  photpresel->AddEnCorrPerRun(195351, 195395, 0.9931, 0.9931, 0.9992, 0.9861, 0.9996, 0.9933, 1.0008, 0.9900, 1.0012);
  photpresel->AddEnCorrPerRun(195396, 195530, 0.9932, 0.9932, 0.9994, 0.9861, 0.9996, 0.9951, 1.0025, 0.9869, 0.9982);
  photpresel->AddEnCorrPerRun(195531, 999999, 0.9924, 0.9924, 0.9986, 0.9857, 0.9992, 0.9962, 1.0036, 0.9885, 0.9997);
  
  //photcicnoeleveto->SetDoMCR9Scaling(kTRUE);
  photpresel->SetDoShowerShapeScaling(kTRUE);
  photpresel->SetShowerShapeType("2012ShowerShape");
  //photpresel->SetDoMCErrScaling(kTRUE);
  //photpresel->SetMCErrScale(1.07, 1.045); 
  //photpresel->SetMCErrScale(1, 1); 
  photpresel->SetJetsName(jetCorr->GetOutputName()); 
  //photpresel->SetRescaledBeamspotWidth(5.0);
  photpresel->SetIsData(isData);

  PhotonPairSelector         *photpreselinverteleveto = new PhotonPairSelector("PhotonPairSelectorPreselInvertEleVeto");
  photpreselinverteleveto->SetOutputName("GoodPhotonsPreselInvertEleVeto");
  photpreselinverteleveto->SetOutputVtxName("OutVtxPreselInvertEleVeto");    
  photpreselinverteleveto->SetPhotonSelType("MITPFSelection");
  photpreselinverteleveto->SetIdMVAType("2012IdMVA_globe");
  photpreselinverteleveto->SetVertexSelType("CiCMVA2012Selection");
  //photpreselinverteleveto->SetUseSingleLegConversions(kTRUE);
  photpreselinverteleveto->DoMCSmear(kTRUE);
  photpreselinverteleveto->DoDataEneCorr(kTRUE);
  photpreselinverteleveto->SetPhotonsFromBranch(kFALSE);
  photpreselinverteleveto->SetInputPhotonsName(photreg->GetOutputName());
  photpreselinverteleveto->SetMCSmearFactors(0.00825,0.01035,0.0116,0.0222,0.0138,0.0292,0.0275,0.0322,0.0310);
  photpreselinverteleveto->AddEnCorrPerRun(190450, 190781, 0.9944, 0.9944, 1.0006, 0.9948, 1.0083, 0.9881, 0.9956, 0.9878, 0.9991);
  photpreselinverteleveto->AddEnCorrPerRun(190782, 190949, 1.0085, 1.0085, 1.0146, 0.9858, 0.9993, 1.0047, 1.0120, 0.9962, 1.0074);
  photpreselinverteleveto->AddEnCorrPerRun(190950, 191833, 0.9939, 0.9939, 1.0001, 0.9879, 1.0014, 0.9906, 0.9981, 0.9874, 0.9986);
  photpreselinverteleveto->AddEnCorrPerRun(191834, 193686, 0.9916, 0.9916, 0.9978, 0.9855, 0.9990, 0.9907, 0.9982, 0.9893, 1.0005);
  photpreselinverteleveto->AddEnCorrPerRun(193746, 194210, 0.9921, 0.9921, 0.9983, 0.9834, 0.9970, 0.9921, 0.9996, 0.9883, 0.9995);
  photpreselinverteleveto->AddEnCorrPerRun(194211, 194479, 0.9926, 0.9926, 0.9988, 0.9845, 0.9980, 0.9928, 1.0003, 0.9895, 1.0007);
  photpreselinverteleveto->AddEnCorrPerRun(194480, 195147, 0.9927, 0.9927, 0.9989, 0.9862, 0.9998, 0.9954, 1.0028, 0.9891, 1.0003);
  photpreselinverteleveto->AddEnCorrPerRun(195148, 195350, 0.9928, 0.9928, 0.9990, 0.9852, 0.9988, 1.0007, 1.0081, 0.9943, 1.0054);
  photpreselinverteleveto->AddEnCorrPerRun(195351, 195395, 0.9931, 0.9931, 0.9992, 0.9861, 0.9996, 0.9933, 1.0008, 0.9900, 1.0012);
  photpreselinverteleveto->AddEnCorrPerRun(195396, 195530, 0.9932, 0.9932, 0.9994, 0.9861, 0.9996, 0.9951, 1.0025, 0.9869, 0.9982);
  photpreselinverteleveto->AddEnCorrPerRun(195531, 999999, 0.9924, 0.9924, 0.9986, 0.9857, 0.9992, 0.9962, 1.0036, 0.9885, 0.9997);
   photpreselinverteleveto->SetShowerShapeType("2012ShowerShape");
   photpreselinverteleveto->SetDoShowerShapeScaling(kTRUE);  
   //photpreselinverteleveto->SetDoMCErrScaling(kTRUE);
   //photpreselinverteleveto->SetMCErrScale(1.07, 1.045);  
   //photpreselinverteleveto->SetMCErrScale(1, 1);  
   photpreselinverteleveto->SetApplyEleVeto(kFALSE);
   photpreselinverteleveto->SetInvertElectronVeto(kFALSE);
   photpreselinverteleveto->SetJetsName(jetCorr->GetOutputName());  
   //photpreselinverteleveto->SetRescaledBeamspotWidth(5.0);   
   photpreselinverteleveto->SetIsData(isData);  
   
   PhotonPairSelector         *photpreselnosmear = new PhotonPairSelector("PhotonPairSelectorPreselNoSmear");
   photpreselnosmear->SetOutputName("GoodPhotonsPreselNoSmear");
   photpreselnosmear->SetPhotonSelType("MITPFSelection");
   photpreselnosmear->SetVertexSelType("CiCMVA2012Selection");
   //photpreselnosmear->SetUseSingleLegConversions(kTRUE);  
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

  
  jetCorr          ->Add(photcic);
  //jetCorr          ->Add(photcicnoeleveto);  
  jetCorr          ->Add(photpresel);  
  jetCorr          ->Add(photpreselinverteleveto);  
  jetCorr          ->Add(photpreselnosmear);  

  photcic         ->Add(phottreecic);
  photcicnoeleveto         ->Add(phottreecicnoeleveto);
  photpresel    ->Add(phottreepresel);
  photpreselinverteleveto    ->Add(phottreepreselinverteleveto);
  photpreselnosmear    ->Add(phottreepreselnosmear);


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
  //ana->AddFile("/scratch/bendavid/root/DA2D6D6F-BB95-E111-9564-001A92971B7E.root");

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
