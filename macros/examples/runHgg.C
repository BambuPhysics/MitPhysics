// $Id: runHgg.C,v 1.6 2012/04/03 08:23:23 mingyang Exp $
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
void runHgg(const char *fileset    = "0000",
            const char *skim       = "noskim",
	    const char *dataset = "f11--h123gg-gf-v14b-pu",
            //const char *dataset = "r11a-pho-j16-v1",   
            const char *book       = "local/filefi/025",
            const char *catalogDir = "/home/cmsprod/catalog",
            const char *outputName = "hgg",
            int         nEvents    = 2000)
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
//   hltModE->AddTrigger("HLT_Photon10_L1R",132440,137028);
//   hltModE->AddTrigger("HLT_Photon15_Cleaned_L1R",138564,140401);
//   hltModE->AddTrigger("HLT_Ele15_SW_CaloEleId_L1R",141956,144114);
//   hltModE->AddTrigger("HLT_Ele17_SW_CaloEleId_L1R",144115,147145);
//   hltModE->AddTrigger("HLT_Ele17_SW_TightEleId_L1R",147146,148102);
//   hltModE->AddTrigger("HLT_Ele27_SW_TightCaloEleIdTrack_L1R_v1",147146,148102);  
//   hltModE->AddTrigger("HLT_Ele22_SW_TighterCaloIdIsol_L1R_v1",148103,159999);  
//   hltModE->AddTrigger("HLT_Ele22_SW_TighterCaloIdIsol_L1R_v2",148103,159999);  
//   hltModE->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v1",160000,999999);
//   hltModE->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2",160000,999999);
//   hltModE->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v3",160000,999999);
//   hltModE->AddTrigger("HLT_Ele15_LW_L1R",1,1);
//   hltModE->AddTrigger("HLT_Ele17_SW_CaloEleId_L1R",1,1);
//   hltModE->AddTrigger("HLT_Ele22_SW_TighterCaloIdIsol_L1R_v2",1,1);


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
  hltModE->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v10",179890,999999);
  hltModE->SetTrigObjsName("MyHltElecObjs");
  hltModE->SetAbortIfNotAccepted(isData);

  HLTMod *hltModES = new HLTMod("HLTModES");
  hltModES->AddTrigger("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v*",160000,999999);
  hltModES->AddTrigger("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2",160000,999999);
  hltModES->AddTrigger("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3",160000,999999);
  hltModES->AddTrigger("HLT_Ele25_WP80_PFMT40_v*",160000,999999);
  hltModES->AddTrigger("HLT_Ele27_WP80_PFMT50_v*",160000,999999);
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
  
  hltModP->AddTrigger("HLT_Photon26_CaloIdXL_IsoXL_Photon18_CaloIdXL_IsoXL_Mass60_v*",178420,999999);
  hltModP->AddTrigger("HLT_Photon26_CaloIdXL_IsoXL_Photon18_R9IdT_Mass60_v*",178420,999999);
  hltModP->AddTrigger("HLT_Photon26_R9IdT_Photon18_CaloIdXL_IsoXL_Mass60_v*",178420,999999);
  hltModP->AddTrigger("HLT_Photon26_R9IdT_Photon18_R9IdT_Mass60_v*",178420,999999);
  hltModP->AddTrigger("HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL_v*",178420,999999);
  hltModP->AddTrigger("HLT_Photon36_CaloIdL_IsoVL_Photon22_R9Id_v*",178420,999999);
  hltModP->AddTrigger("HLT_Photon36_R9Id_Photon22_CaloIdL_IsoVL_v*",178420,999999);
  hltModP->AddTrigger("HLT_Photon36_R9Id_Photon22_R9Id_v*",178420,999999);
    
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

  ElectronIDMod *elecId = new ElectronIDMod;
  elecId->SetVertexName(goodPVFilterModE->GetOutputName());
  elecId->SetIDType                    ("VBTFWorkingPointLowPtId");
  elecId->SetIsoType                   ("PFIso");
  elecId->SetNExpectedHitsInnerCut     (0);
  elecId->SetEtaMax(999.);
  elecId->SetChargeFilter(kFALSE);
  elecId->SetApplySpikeRemoval(kFALSE);
  elecId->SetApplyEcalFiducial(kTRUE);
  elecId->SetApplyEcalSeeded(kTRUE);
  elecId->SetPtMin(20.0);
  
  ElectronIDMod *elecIdS = new ElectronIDMod("ElectronIDModS");
  elecIdS->SetVertexName(goodPVFilterModES->GetOutputName());
  elecIdS->SetIDType                    ("VBTFWorkingPoint70Id");
  elecIdS->SetIsoType                   ("VBTFWorkingPoint70Iso");
  elecIdS->SetNExpectedHitsInnerCut     (0);
  elecIdS->SetEtaMax(999.);
  elecIdS->SetChargeFilter(kFALSE);
  elecIdS->SetApplySpikeRemoval(kFALSE);
  elecIdS->SetApplyEcalFiducial(kTRUE);
  elecIdS->SetApplyEcalSeeded(kTRUE);
  elecIdS->SetPtMin(20.0);
  elecIdS->SetOutputName("GoodElectronsS");
  
  PublisherMod<PFJet,Jet> *pubJet = new PublisherMod<PFJet,Jet>("JetPub");
  pubJet->SetInputName("AKt5PFJets");
  pubJet->SetOutputName("PubAKt5PFJets");
  
  PublisherMod<PFJet,Jet> *pubJetOpen = new PublisherMod<PFJet,Jet>("JetPubOpen");
  pubJetOpen->SetInputName("AKt5PFJets");
  pubJetOpen->SetOutputName("PubAKt5PFJetsOpen");  

  JetCorrectionMod *jetCorr = new JetCorrectionMod;
  jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/START42_V12_AK5PF_L1FastJet.txt")).Data())); 
  jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/START42_V12_AK5PF_L2Relative.txt")).Data())); 
  jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/START42_V12_AK5PF_L3Absolute.txt")).Data())); 
  if(isData){ 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/START42_V12_AK5PF_L2L3Residual.txt")).Data())); 
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
  
  PhotonMvaMod *photreg = new PhotonMvaMod;
  photreg->SetOutputName("GoodPhotonsRegr");
  photreg->SetIsData(isData);

    
  PhotonPairSelector         *photcic = new PhotonPairSelector("PhotonPairSelectorCiC");
  photcic->SetOutputName("GoodPhotonsCIC");
  photcic->SetPhotonSelType("CiCSelection");
  photcic->SetVertexSelType("CiCMVASelection");
  photcic->DoMCSmear(kTRUE);
  photcic->DoDataEneCorr(kTRUE);
  photcic->SetPhotonsFromBranch(kFALSE);
  photcic->SetInputPhotonsName(photreg->GetOutputName());
  photcic->SetMCSmearFactors(0.0089, 0.0089, 0.0109, 0.0156, 0.0203,0.0303,0.0326,0.0318,0.0331);
  photcic->AddEnCorrPerRun(160000,167913,0.9905,0.9905,0.9971,0.9976,1.0094,0.9994,1.0044,0.9968,1.0079);
  photcic->AddEnCorrPerRun(170249,172619,0.9909,0.9909,0.9975,0.9994,1.0112,0.9962,1.0012,0.9962,1.0072);
  photcic->AddEnCorrPerRun(172620,173692,0.9909,0.9909,0.9975,0.9977,1.0096,0.9963,1.0013,0.9947,1.0057);
  photcic->AddEnCorrPerRun(175860,177139,0.9911,0.9911,0.9977,0.9990,1.0109,0.9922,0.9973,0.9967,1.0077);
  photcic->AddEnCorrPerRun(177140,178421,0.9910,0.9910,0.9975,0.9987,1.0105,0.9921,0.9972,0.9975,1.0085);
  photcic->AddEnCorrPerRun(178424,999999,0.9903,0.9903,0.9969,0.9976,1.0095,0.9889,0.9940,0.9976,1.0086);
  photcic->SetDoMCR9Scaling(kTRUE);
  photcic->SetMCR9Scale(1.0048, 1.00492);
  photcic->SetDoMCErrScaling(kTRUE);
  photcic->SetMCErrScale(1.07, 1.045);    
  photcic->SetIsData(isData);

  PhotonPairSelector         *photcicnoeleveto = new PhotonPairSelector("PhotonPairSelectorCiCNoEleVeto");
  photcicnoeleveto->SetOutputName("GoodPhotonsCICNoEleVeto");
  photcicnoeleveto->SetPhotonSelType("CiCSelection");
  photcicnoeleveto->SetVertexSelType("CiCMVASelection");
  photcicnoeleveto->DoMCSmear(kTRUE);
  photcicnoeleveto->DoDataEneCorr(kTRUE);
  photcicnoeleveto->SetPhotonsFromBranch(kFALSE);
  photcicnoeleveto->SetInputPhotonsName(photreg->GetOutputName());
  photcicnoeleveto->SetMCSmearFactors(0.0089, 0.0089, 0.0109, 0.0156, 0.0203,0.0303,0.0326,0.0318,0.0331);
  photcicnoeleveto->AddEnCorrPerRun(160000,167913,0.9905,0.9905,0.9971,0.9976,1.0094,0.9994,1.0044,0.9968,1.0079);
  photcicnoeleveto->AddEnCorrPerRun(170249,172619,0.9909,0.9909,0.9975,0.9994,1.0112,0.9962,1.0012,0.9962,1.0072);
  photcicnoeleveto->AddEnCorrPerRun(172620,173692,0.9909,0.9909,0.9975,0.9977,1.0096,0.9963,1.0013,0.9947,1.0057);
  photcicnoeleveto->AddEnCorrPerRun(175860,177139,0.9911,0.9911,0.9977,0.9990,1.0109,0.9922,0.9973,0.9967,1.0077);
  photcicnoeleveto->AddEnCorrPerRun(177140,178421,0.9910,0.9910,0.9975,0.9987,1.0105,0.9921,0.9972,0.9975,1.0085);
  photcicnoeleveto->AddEnCorrPerRun(178424,999999,0.9903,0.9903,0.9969,0.9976,1.0095,0.9889,0.9940,0.9976,1.0086);
  photcicnoeleveto->SetDoMCR9Scaling(kTRUE);
  photcicnoeleveto->SetMCR9Scale(1.0048, 1.00492);
  photcicnoeleveto->SetDoMCErrScaling(kTRUE);
  photcicnoeleveto->SetMCErrScale(1.07, 1.045);    
  photcicnoeleveto->SetApplyEleVeto(kFALSE);
  photcicnoeleveto->SetIsData(isData);
  
  
  PhotonPairSelector         *photpresel = new PhotonPairSelector("PhotonPairSelectorPresel");
  photpresel->SetOutputName("GoodPhotonsPresel");
  photpresel->SetPhotonSelType("MITSelection");
  photpresel->SetVertexSelType("CiCMVASelection");
  photpresel->DoMCSmear(kTRUE);
  photpresel->DoDataEneCorr(kTRUE);
  photpresel->SetPhotonsFromBranch(kFALSE);
  photpresel->SetInputPhotonsName(photreg->GetOutputName());
  //photpresel->SetMCSmearFactors(0.0045, 0.0084, 0.0109, 0.0156, 0.0203,0.0303,0.0326,0.0318,0.0331);
  photpresel->SetMCSmearFactors(0.0067,0.0077,0.0096,0.0141,0.0196,0.0268,0.0279,0.0293,0.0301);//ming:smear(sigE/E)
  //photpresel->AddEnCorrPerRun(160000,167913,0.9905,0.9905,0.9971,0.9976,1.0094,0.9994,1.0044,0.9968,1.0079);
  //photpresel->AddEnCorrPerRun(170249,172619,0.9909,0.9909,0.9975,0.9994,1.0112,0.9962,1.0012,0.9962,1.0072);
  //photpresel->AddEnCorrPerRun(172620,173692,0.9909,0.9909,0.9975,0.9977,1.0096,0.9963,1.0013,0.9947,1.0057);
  //photpresel->AddEnCorrPerRun(175860,177139,0.9911,0.9911,0.9977,0.9990,1.0109,0.9922,0.9973,0.9967,1.0077);
  //photpresel->AddEnCorrPerRun(177140,178421,0.9910,0.9910,0.9975,0.9987,1.0105,0.9921,0.9972,0.9975,1.0085);
  //photpresel->AddEnCorrPerRun(178424,999999,0.9903,0.9903,0.9969,0.9976,1.0095,0.9889,0.9940,0.9976,1.0086); 
  photpresel->AddEnCorrPerRun(160431,167913,0.9941,0.9941,1.0004,0.9916,1.0045,1.0033,1.0082,0.9958,1.0064);//ming:Emc/Edata
  photpresel->AddEnCorrPerRun(170000,172619,0.9954,0.9954,1.0016,0.9937,1.0066,0.9976,1.0025,0.9940,1.0046);
  photpresel->AddEnCorrPerRun(172620,173692,0.9955,0.9955,1.0017,0.9929,1.0058,0.9986,1.0035,0.9923,1.0029);
  photpresel->AddEnCorrPerRun(175830,177139,0.9958,0.9958,1.0021,0.9944,1.0073,0.9968,1.0017,0.9933,1.004);
  photpresel->AddEnCorrPerRun(177140,178421,0.9962,0.9962,1.0025,0.9946,1.0075,0.9960,1.0010,0.9944,1.005);
  photpresel->AddEnCorrPerRun(178424,180252,0.9961,0.9961,1.0024,0.9942,1.0071,0.9921,0.9970,0.9953,1.0059); 
  photpresel->SetDoMCR9Scaling(kTRUE);
  photpresel->SetMCR9Scale(1.0035, 1.0035);  
  photpresel->SetDoMCSigIEtaIEtaScaling(kTRUE);
  photpresel->SetDoMCWidthScaling(kTRUE);  
  photpresel->SetDoMCErrScaling(kTRUE);
  photpresel->SetMCErrScale(1.07, 1.045);    
  photpresel->SetIsData(isData);

  PhotonPairSelector         *photpreselinverteleveto = new PhotonPairSelector("PhotonPairSelectorPreselInvertEleVeto");
  photpreselinverteleveto->SetOutputName("GoodPhotonsPreselInvertEleVeto");
  photpreselinverteleveto->SetPhotonSelType("MITSelection");
  photpreselinverteleveto->SetVertexSelType("CiCMVASelection");
  photpreselinverteleveto->DoMCSmear(kTRUE);
  photpreselinverteleveto->DoDataEneCorr(kTRUE);
  photpreselinverteleveto->SetPhotonsFromBranch(kFALSE);
  photpreselinverteleveto->SetInputPhotonsName(photreg->GetOutputName());
  photpreselinverteleveto->SetMCSmearFactors(0.0045, 0.0084, 0.0109, 0.0156, 0.0203,0.0303,0.0326,0.0318,0.0331);
  photpreselinverteleveto->AddEnCorrPerRun(160000,167913,0.9905,0.9905,0.9971,0.9976,1.0094,0.9994,1.0044,0.9968,1.0079);
  photpreselinverteleveto->AddEnCorrPerRun(170249,172619,0.9909,0.9909,0.9975,0.9994,1.0112,0.9962,1.0012,0.9962,1.0072);
  photpreselinverteleveto->AddEnCorrPerRun(172620,173692,0.9909,0.9909,0.9975,0.9977,1.0096,0.9963,1.0013,0.9947,1.0057);
  photpreselinverteleveto->AddEnCorrPerRun(175860,177139,0.9911,0.9911,0.9977,0.9990,1.0109,0.9922,0.9973,0.9967,1.0077);
  photpreselinverteleveto->AddEnCorrPerRun(177140,178421,0.9910,0.9910,0.9975,0.9987,1.0105,0.9921,0.9972,0.9975,1.0085);
  photpreselinverteleveto->AddEnCorrPerRun(178424,999999,0.9903,0.9903,0.9969,0.9976,1.0095,0.9889,0.9940,0.9976,1.0086); 
  photpreselinverteleveto->SetDoMCR9Scaling(kTRUE);
  photpreselinverteleveto->SetMCR9Scale(1.0035, 1.0035);
  photpreselinverteleveto->SetDoMCSigIEtaIEtaScaling(kTRUE);
  photpreselinverteleveto->SetDoMCWidthScaling(kTRUE);  
  photpreselinverteleveto->SetDoMCErrScaling(kTRUE);
  photpreselinverteleveto->SetMCErrScale(1.07, 1.045);    
  photpreselinverteleveto->SetApplyEleVeto(kFALSE);
  photpreselinverteleveto->SetInvertElectronVeto(kTRUE);
  photpreselinverteleveto->SetIsData(isData);  
  
  PhotonPairSelector         *photpreselnosmear = new PhotonPairSelector("PhotonPairSelectorPreselNoSmear");
  photpreselnosmear->SetOutputName("GoodPhotonsPreselNoSmear");
  photpreselnosmear->SetPhotonSelType("MITSelection");
  photpreselnosmear->SetVertexSelType("CiCMVASelection");
  photpreselnosmear->SetPhotonsFromBranch(kFALSE);
  photpreselnosmear->SetInputPhotonsName(photreg->GetOutputName());
  photpreselnosmear->SetIsData(isData);  
  
  
  PhotonTreeWriter *phottreecic = new PhotonTreeWriter("PhotonTreeWriterCiC");
  phottreecic->SetPhotonsFromBranch(kFALSE);
  phottreecic->SetInputPhotonsName(photcic->GetOutputName());
  phottreecic->SetEnableJets(kTRUE);
  phottreecic->SetPFJetsFromBranch(kFALSE);
  phottreecic->SetPFJetName(jetCorr->GetOutputName());
  phottreecic->SetExcludeDoublePrompt(excludedoubleprompt);
  phottreecic->SetIsData(isData);

  PhotonTreeWriter *phottreecicnoeleveto = new PhotonTreeWriter("PhotonTreeWriterCiCNoEleVeto");
  phottreecicnoeleveto->SetPhotonsFromBranch(kFALSE);
  phottreecicnoeleveto->SetInputPhotonsName(photcicnoeleveto->GetOutputName());
  phottreecicnoeleveto->SetEnableJets(kTRUE);
  phottreecicnoeleveto->SetPFJetsFromBranch(kFALSE);
  phottreecicnoeleveto->SetPFJetName(jetCorr->GetOutputName());
  phottreecicnoeleveto->SetApplyElectronVeto(kFALSE);
  phottreecicnoeleveto->SetExcludeDoublePrompt(excludedoubleprompt);
  phottreecicnoeleveto->SetIsData(isData);  
  
  PhotonTreeWriter *phottreepresel = new PhotonTreeWriter("PhotonTreeWriterPresel");
  phottreepresel->SetPhotonsFromBranch(kFALSE);
  phottreepresel->SetInputPhotonsName(photpresel->GetOutputName());
  phottreepresel->SetEnableJets(kTRUE);
  phottreepresel->SetPFJetsFromBranch(kFALSE);
  phottreepresel->SetPFJetName(jetCorr->GetOutputName());  
  phottreepresel->SetExcludeDoublePrompt(excludedoubleprompt);  
  phottreepresel->SetIsData(isData);  
  
  PhotonTreeWriter *phottreepreselinverteleveto = new PhotonTreeWriter("PhotonTreeWriterPreselInvertEleVeto");
  phottreepreselinverteleveto->SetPhotonsFromBranch(kFALSE);
  phottreepreselinverteleveto->SetInputPhotonsName(photpreselinverteleveto->GetOutputName());
  phottreepreselinverteleveto->SetEnableJets(kTRUE);
  phottreepreselinverteleveto->SetPFJetsFromBranch(kFALSE);
  phottreepreselinverteleveto->SetPFJetName(jetCorr->GetOutputName());  
  phottreepreselinverteleveto->SetApplyElectronVeto(kFALSE);  
  phottreepreselinverteleveto->SetExcludeDoublePrompt(excludedoubleprompt);    
  phottreepreselinverteleveto->SetIsData(isData); 
  
  PhotonTreeWriter *phottreepreselnosmear = new PhotonTreeWriter("PhotonTreeWriterPreselNoSmear");
  phottreepreselnosmear->SetPhotonsFromBranch(kFALSE);
  phottreepreselnosmear->SetInputPhotonsName(photpreselnosmear->GetOutputName());
  phottreepreselnosmear->SetEnableJets(kTRUE);
  phottreepreselnosmear->SetPFJetsFromBranch(kFALSE);
  phottreepreselnosmear->SetPFJetName(jetCorr->GetOutputName());  
  phottreepreselnosmear->SetExcludeDoublePrompt(excludedoubleprompt);  
  phottreepreselnosmear->SetIsData(isData);    
  
  
  PhotonIDMod         *photidcic = new PhotonIDMod("PhotonIDModPresel");
  photidcic->SetPtMin(25.0);
  photidcic->SetOutputName("GoodPhotonsPreselid");
  photidcic->SetOutputName("MITSelection");
  photidcic->SetApplyElectronVeto(kTRUE);
  photidcic->SetIsData(isData);

  PhotonTreeWriter *phottreesingle = new PhotonTreeWriter("PhotonTreeWriterSingle");
  phottreesingle->SetWriteDiphotonTree(kFALSE);
  phottreesingle->SetPhotonsFromBranch(kFALSE);
  phottreesingle->SetInputPhotonsName(photidcic->GetOutputName());  
  phottreesingle->SetIsData(isData);
  
  PhotonTreeWriter *phottreeE = new PhotonTreeWriter("PhotonTreeWriterE");
  phottreeE->SetGoodElectronsFromBranch(kFALSE);
  phottreeE->SetGoodElectronName(elecId->GetOutputName());  
  phottreeE->SetLoopOnGoodElectrons(kTRUE);
  phottreeE->SetApplyElectronVeto(kFALSE);
  phottreeE->SetIsData(isData);

  PhotonTreeWriter *phottreeES = new PhotonTreeWriter("PhotonTreeWriterES");
  phottreeES->SetWriteDiphotonTree(kFALSE);
  phottreeES->SetGoodElectronsFromBranch(kFALSE);
  phottreeES->SetGoodElectronName(elecIdS->GetOutputName());  
  phottreeES->SetLoopOnGoodElectrons(kTRUE);
  phottreeES->SetApplyElectronVeto(kFALSE);
  phottreeES->SetIsData(isData);  
  

  //------------------------------------------------------------------------------------------------
  // making analysis chain
  //------------------------------------------------------------------------------------------------
  // this is how it always starts
  runLumiSel      ->Add(mcselmod);
    
  if (TString(dataset).Contains("-h")) {    
    mcselmod        ->Add(sysMod);
  }
  
  // high level trigger is always first
  //mcselmod         ->Add(hltModE);
  //mcselmod         ->Add(hltModES);
  mcselmod         ->Add(hltModP);

  hltModP         ->Add(goodPVFilterMod);
  //hltModE         ->Add(goodPVFilterModE);
  //hltModES        ->Add(goodPVFilterModES);
  
  //goodPVFilterMod ->Add(muonId);
  goodPVFilterMod->Add(photreg);
  photreg->Add(pubJet);
  pubJet->Add(jetCorr);
  
  // simple object id modules
  //goodPVFilterModE ->Add(elecId);
  //goodPVFilterModES ->Add(elecIdS);

  
  //jetCorr          ->Add(photcic);
  //jetCorr          ->Add(photcicnoeleveto);  
  jetCorr          ->Add(photpresel);  
  //jetCorr          ->Add(photpreselinverteleveto);  
  //jetCorr          ->Add(photpreselnosmear);  

  //photcic         ->Add(phottreecic);
  //photcicnoeleveto         ->Add(phottreecicnoeleveto);
  photpresel    ->Add(phottreepresel);
  //photpreselinverteleveto    ->Add(phottreepreselinverteleveto);
  //photpreselnosmear    ->Add(phottreepreselnosmear);


  //jetCorr          ->Add(photidcic);
  //photidcic       ->Add(phottreesingle);
  
  //elecId->Add(phottreeE);
  //elecIdS->Add(phottreeES);
  

  TFile::SetOpenTimeout(0);
  TFile::SetCacheFileDir("./rootfilecache",kTRUE,kTRUE);
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
  //ana->AddFile("/mnt/hadoop/cmsprod/filefi/025/r11b-pho-n30-v1/4863E9D1-BC1C-E111-99BA-001A92810AF4.root");

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
