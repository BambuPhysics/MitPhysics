 // $Id: FakeLeptonExampleAnaMod.cc,v 1.5 2009/08/11 09:16:01 loizides Exp $

#include "MitPhysics/FakeMods/interface/FakeLeptonExampleAnaMod.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/DataTree/interface/Names.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/GenJetCol.h"
#include "MitAna/DataTree/interface/MCParticleCol.h"
#include "MitAna/DataTree/interface/MetCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/TrackCol.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/FakeMods/interface/FakeObject.h"
#include "MitPhysics/FakeMods/interface/FakeEventHeader.h"
#include <TH1D.h>
#include <TH2D.h>

using namespace mithep;

ClassImp(mithep::FakeLeptonExampleAnaMod)

//--------------------------------------------------------------------------------------------------
FakeLeptonExampleAnaMod::FakeLeptonExampleAnaMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fUseMCFake(kFALSE),
  fPerformFakeMuonMetCorrection(kTRUE),
  fSampleName("NotSet"),
  fFakeEventHeaderName(ModNames::gkFakeEventHeadersName),
  fElectronFakeableObjectsName(ModNames::gkElFakeableObjsName),
  fMuonFakeableObjectsName(ModNames::gkMuFakeableObjsName),
  fMCPartBranchName(Names::gkMCPartBrn),
  fGenJetBranchName(Names::gkSC5GenJetBrn),
  fTrackBranchName(Names::gkTrackBrn),
  fMuonBranchName(Names::gkMuonBrn),
  fMetName("NotSet"),
  fCleanJetsName(ModNames::gkCleanJetsName),
  fTriggerObjectsName("NotSet"),
  fParticles(0),
  fGenJets(0),
  fTracks(0),
  fMuons(0),
  fMet(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void FakeLeptonExampleAnaMod::Process()
{
  // Process entries of the tree.
  LoadBranch(fTrackBranchName);
  LoadBranch(fMuonBranchName);

  //***********************************************************************************************
  //Import Collections
  //***********************************************************************************************

  //Obtain all cleaned objects
  MuonOArr *CleanMuons = dynamic_cast<MuonOArr* >
    (FindObjThisEvt(ModNames::gkCleanMuonsName));
  ParticleOArr *CleanLeptons = dynamic_cast<mithep::ParticleOArr*>
    (FindObjThisEvt(ModNames::gkMergedLeptonsName));

  //Get Met
  if (!fMetName.IsNull()) {
    fMet = GetObjThisEvt<MetCol>(fMetName);
  } 
  const Met *originalCaloMet = 0;
  if (fMet) {
    originalCaloMet = fMet->At(0);
  } else {
    std::cout << "Error: Met Collection " << fMetName << " could not be loaded." << std::endl;
  }
  ObjArray<Jet> *OriginalCleanJets = dynamic_cast<ObjArray<Jet>* > 
    (FindObjThisEvt(fCleanJetsName.Data()));

  //Obtain the collection of fake objects
  ElectronCol *ElectronFakeableObjects = 0;
  if(!fElectronFakeableObjectsName.IsNull())
    ElectronFakeableObjects = dynamic_cast<ElectronCol* >
      (FindObjThisEvt(fElectronFakeableObjectsName.Data()));
  MuonCol *MuonFakeableObjects = 0;
  if (!fMuonFakeableObjectsName.IsNull())
    MuonFakeableObjects = dynamic_cast<MuonCol* >
      (FindObjThisEvt(fMuonFakeableObjectsName.Data()));
  ChargedParticleOArr *FakeableObjects = new ChargedParticleOArr;
  if (ElectronFakeableObjects) {
    for (UInt_t i=0; i<ElectronFakeableObjects->GetEntries(); i++)
      FakeableObjects->Add(ElectronFakeableObjects->At(i));
  }
  if (MuonFakeableObjects) {
    for (UInt_t i=0; i<MuonFakeableObjects->GetEntries(); i++)
      FakeableObjects->Add(MuonFakeableObjects->At(i));  
  }
  Collection<FakeEventHeader> *FakeEventHeaders = 0;
  if (!fUseMCFake) {
    if (!fFakeEventHeaderName.IsNull()) {
      FakeEventHeaders = 
        dynamic_cast<Collection<FakeEventHeader>* >(FindObjThisEvt(fFakeEventHeaderName.Data()));
      if (!FakeEventHeaders) {
        std::string errorMsg = "Error: FakeEventHeader with name  " + 
          std::string(fFakeEventHeaderName.Data()) + "could not be loaded.\n";
        Fatal("Process()",errorMsg.c_str());        
      }
    } else 
      std::cout << "Error: FakeEventHeaders  " << fFakeEventHeaderName.Data() 
                << " could not be loaded." << std::endl;
  }  

  //***********************************************************************************************
  //If we use MC Fakes, then create a new FakeEventHeader containing no fakes with weight = 1
  //This ensures that in the loop over FakeEventHeaders we do the correct thing.
  //***********************************************************************************************
  if (fUseMCFake) {
    ObjArray <FakeEventHeader> *tmpFakeEventHeaders = new  ObjArray <FakeEventHeader> ;
    tmpFakeEventHeaders->SetOwner(kTRUE);

    FakeEventHeader *initialFakeEvent = new FakeEventHeader();
    for (UInt_t j=0;j<OriginalCleanJets->GetEntries();j++)
      initialFakeEvent->AddJets(OriginalCleanJets->At(j));

    tmpFakeEventHeaders->AddOwned(initialFakeEvent);
    FakeEventHeaders = dynamic_cast<Collection<FakeEventHeader>* > (tmpFakeEventHeaders);
  }

  //***********************************************************************************************
  //Loop over Fake Event Headers
  //***********************************************************************************************
  for (UInt_t i=0;i<FakeEventHeaders->GetEntries() ; i++) {

    //Create leptons collection containing real leptons and fake leptons
    ObjArray<Particle> *leptons = NULL;
    if (fUseMCFake) {
      leptons = CleanLeptons;
    } else {
      leptons = new ObjArray<Particle>;

      for (UInt_t j=0;j<CleanLeptons->GetEntries() ; j++) {
        leptons->Add(CleanLeptons->At(j));
      }
      for (UInt_t j=0;j<FakeEventHeaders->At(i)->FakeObjsSize() ; j++) {
        leptons->Add(FakeEventHeaders->At(i)->FakeObj(j)->FakeParticle());
      }
    }
    //we have to sort leptons
    leptons->Sort();


    //Construct the Clean Jet collection.
    ObjArray<Jet> *CleanJets = NULL;
    if (fUseMCFake) {
      CleanJets = OriginalCleanJets;
    } else {
      CleanJets = new ObjArray<Jet>;
      for (UInt_t j=0;j<FakeEventHeaders->At(i)->NJets() ; j++) {
        CleanJets->Add(FakeEventHeaders->At(i)->UnfakedJet(j));
      }
    }

    //Perform correction for potential fake muons
    //have to add fake muon momentum to originalCaloMet;
    const Met *caloMet = originalCaloMet;
    Double_t FakeMuonMetCorrection_X = 0.0;
    Double_t FakeMuonMetCorrection_Y = 0.0;
    for (UInt_t j=0;j<FakeEventHeaders->At(i)->FakeObjsSize() ; j++) {
      if (FakeEventHeaders->At(i)->FakeObj(j)->ObjType() == kMuon) {
        FakeMuonMetCorrection_X += FakeEventHeaders->At(i)->FakeObj(j)->Px();
        FakeMuonMetCorrection_Y += FakeEventHeaders->At(i)->FakeObj(j)->Py();
      }
    }
    
    if (!fUseMCFake && fPerformFakeMuonMetCorrection) {
      caloMet = new Met(originalCaloMet->Px()+FakeMuonMetCorrection_X,
                        originalCaloMet->Py()+FakeMuonMetCorrection_Y);
    }

    //*********************************************************************************************
    //Construct the event weight using fake rate and corrections
    //*********************************************************************************************
    //fake rate has to be corrected by the amount lost when those denominators 
    //became fakes in data. If a denominator fakes a lepton in data, it goes in the 2lepton
    //final state, and we don't count it in this prediction. So we have to add it back.
    Double_t eventweight = FakeEventHeaders->At(i)->Weight();
    if (FakeEventHeaders->At(i)->FakeObjsSize() > 0 && FakeEventHeaders->At(i)->Weight() < 1) {
      eventweight = eventweight / (1.0 - FakeEventHeaders->At(i)->Weight());
    }

    //*********************************************************************************************
    //another correction to account for events lost due to only the fake lepton firing the trigger
    //The numbers need to be changed for your analysis. 
    //Given numbers are for the 2 lepton final state.
    //*********************************************************************************************
    if (CleanLeptons->GetEntries() >= 1 && FakeEventHeaders->At(i)->FakeObjsSize() >= 1) {
      if (CleanLeptons->At(0)->ObjType() == kElectron && 
          FakeEventHeaders->At(i)->FakeObj(0)->FakeParticle()->ObjType() == kElectron) {
        eventweight = eventweight * 1.06;
      } else if (CleanLeptons->At(0)->ObjType() == kMuon && 
                 FakeEventHeaders->At(i)->FakeObj(0)->FakeParticle()->ObjType() == kMuon) {
        eventweight = eventweight * 1.12;
      } else if (CleanLeptons->At(0)->ObjType() == kElectron && 
                 FakeEventHeaders->At(i)->FakeObj(0)->FakeParticle()->ObjType() == kMuon) {
        eventweight = eventweight * 1.17;
      } else if (CleanLeptons->At(0)->ObjType() == kMuon && 
                FakeEventHeaders->At(i)->FakeObj(0)->FakeParticle()->ObjType() == kElectron) {
        eventweight = eventweight * 1.17;        
      }
    }

    //*********************************************************************************************
    //For FR method (fUseMCFake == kFALSE)
    //Make analysis specific cuts. 
    //For example for 2 lepton final state we require that the event contains
    //one and only one clean lepton with pt > 10 GeV. 
    //*********************************************************************************************
    if (!fUseMCFake) {
      if (CleanLeptons->GetEntries() != 1 || CleanLeptons->At(0)->Pt() <= 10.0)
        continue;
    }

    //*********************************************************************************************
    //Fill some distributions before preselection
    //*********************************************************************************************
    
    CompositeParticle *dilepton = NULL;

    if (leptons->GetEntries()>=2) {

      dilepton = new CompositeParticle();
      dilepton->AddDaughter(leptons->At(0));
      dilepton->AddDaughter(leptons->At(1));

      //Dilepton Charge will be filled like this
      // -2: -- , -1: -+, 1: +-, 2:++
      if (dilepton->Charge() == 0) {
        if (leptons->At(0)->Charge() == 1) {
          fDileptonCharge->Fill(1.0,eventweight);
        } else {
          fDileptonCharge->Fill(-1.0,eventweight);
        }
      } else {
        fDileptonCharge->Fill(dilepton->Charge(),eventweight);
      }
    }

    //*********************************************************************************************
    //Kinematic PreSelection
    //Example given is for the two lepton final state
    //*********************************************************************************************
    //make sure 2nd highest pt lepton has Pt > 10
    if (leptons->GetEntries() < 2 || leptons->At(1)->Pt() <= 10) continue;
    
    //make sure the 3rd highest pt lepton has pt <= 10.
    if (leptons->GetEntries() >= 3 && leptons->At(2)->Pt() > 10) continue;
   
    //charge of the leptons should be opposite
    if (dilepton->Charge() != 0) continue;

    //*********************************************************************************************
    //Get nonisolated soft muons
    //*********************************************************************************************
    ObjArray<Muon> *DirtyMuons = new ObjArray<Muon>;
    for (UInt_t m=0; m<fMuons->GetEntries(); ++m) {
      const Muon *mu = fMuons->At(m);
      if(!mu->GlobalTrk()) continue;
      if(mu->Pt() < 5.0)   continue;

      //remove the fake
       Bool_t isFakedMuon = kFALSE;
       for (UInt_t f=0;f<FakeEventHeaders->At(i)->FakeObjsSize() ; f++) {
         if (mu->HasTrackerTrk() && 
             (dynamic_cast<const mithep::ChargedParticle*>
              (FakeEventHeaders->At(i)->FakeObj(f)->FakeParticle()))->TrackerTrk() &&
             (dynamic_cast<const mithep::ChargedParticle*>
              (FakeEventHeaders->At(i)->FakeObj(f)->FakeParticle()))->TrackerTrk() ==
             mu->TrackerTrk()
           )
           isFakedMuon = kTRUE;
       }
      
      //remove clean muons
      Bool_t isCleanMuon = kFALSE;
      for (UInt_t j=0; j<CleanMuons->GetEntries(); j++) {
        if(fMuons->At(m) == CleanMuons->At(j)) isCleanMuon = kTRUE;
      }

      if(!isCleanMuon 
         && !(isFakedMuon && !fUseMCFake)
        ) DirtyMuons->Add(mu);
    }

    //*********************************************************************************************
    //Get Clean Tracks excluding the good leptons
    //*********************************************************************************************
    ObjArray<Track> *CleanExtraTracks = new ObjArray<Track>;
    Int_t nTracks = 0;

    Double_t z0Average = ( (dynamic_cast<const mithep::ChargedParticle*>(leptons->At(0)))->Trk()->Z0()
      + (dynamic_cast<const mithep::ChargedParticle*>(leptons->At(1)))->Trk()->Z0()) /2 ;
    
    for (UInt_t t=0; t<fTracks->GetEntries(); ++t) {
      Bool_t isLepton = kFALSE;
      
      if (MathUtils::DeltaR(fTracks->At(t)->Phi(),fTracks->At(t)->Eta(),leptons->At(0)->Phi(),
                            leptons->At(0)->Eta()) > 0.01 &&
          MathUtils::DeltaR(fTracks->At(t)->Phi(),fTracks->At(t)->Eta(),leptons->At(1)->Phi(),
                            leptons->At(1)->Eta()) > 0.01
        ) {
      } else {
        isLepton = kTRUE;
      }
      
      MDB(kAnalysis, 8) {
        std::cout << "Track " << t << " : "  << fTracks->At(t)->Pt() << " " << fTracks->At(t)->Eta() 
                  << " " << fTracks->At(t)->Phi() << " islepton=" << isLepton << std::endl;
      }             
      
      if ( !isLepton && fTracks->At(t)->Pt() > 3.0 
           && fTracks->At(t)->NHits() >= 8 
           && fabs(z0Average - fTracks->At(t)->Z0()) < 0.5 ) {
        CleanExtraTracks->Add(fTracks->At(t));
        nTracks++;
      }
    }

    //*********************************************************************************************
    //The code below is an example analysis for the HWW analysis. 
    //*********************************************************************************************


    //*********************************************************************************************
    //Define Event Variables
    //*********************************************************************************************
    //delta phi between the 2 leptons in degrees
    Double_t deltaPhiLeptons = MathUtils::DeltaPhi(leptons->At(0)->Phi(), 
                                                 leptons->At(1)->Phi())* 360.0 / 2 / TMath::Pi();
    
    Double_t deltaEtaLeptons = leptons->At(0)->Eta() - leptons->At(1)->Eta();
    
    Double_t deltaPhiDileptonMet = MathUtils::DeltaPhi(caloMet->Phi(), 
                                                     dilepton->Phi())*360.0 / 2 / TMath::Pi();
    
    Double_t mtHiggs = TMath::Sqrt(2.0*dilepton->Pt() * caloMet->Pt()*
                                 (1.0 - cos(deltaPhiDileptonMet * 2 * TMath::Pi() / 360.0)));
    
    //angle between MET and closest lepton
    Double_t deltaPhiMetLepton[2] = {MathUtils::DeltaPhi(caloMet->Phi(), leptons->At(0)->Phi()),
                                   MathUtils::DeltaPhi(caloMet->Phi(), leptons->At(1)->Phi())};
    
    Double_t mTW[2] = {TMath::Sqrt(2.0*leptons->At(0)->Pt()*caloMet->Pt()*
                                 (1.0 - cos(deltaPhiMetLepton[0]))),
                     TMath::Sqrt(2.0*leptons->At(1)->Pt()*caloMet->Pt()*
                                 (1.0 - cos(deltaPhiMetLepton[1])))};
    
    Double_t minDeltaPhiMetLepton = (deltaPhiMetLepton[0] < deltaPhiMetLepton[1])?
      deltaPhiMetLepton[0]:deltaPhiMetLepton[1];
    minDeltaPhiMetLepton = minDeltaPhiMetLepton * 360.0 / 2 / TMath::Pi();
    
    //count the number of central Jets for vetoing
    Int_t nCentralJets = 0;
    for (UInt_t j=0; j<CleanJets->GetEntries(); j++) {
      if (fabs(CleanJets->At(j)->Eta()) < 2.5)
        nCentralJets++;
    }
    
    //Lepton Type
    Int_t finalstateType = -1;
    if (leptons->At(0)->ObjType() == kMuon && 
        leptons->At(1)->ObjType() == kMuon ){ // mumu
      finalstateType = 10;
    } else if(leptons->At(0)->ObjType() == kElectron && 
              leptons->At(1)->ObjType() == kElectron ){ // ee
      finalstateType = 11;
    } else if((leptons->At(0)->ObjType() == kElectron && leptons->At(1)->ObjType() == kMuon) || 
              (leptons->At(0)->ObjType() == kMuon && leptons->At(1)->ObjType() == kElectron)) {
      finalstateType = 12;
    } else {
      std::cerr << "Error: finalstate lepton type not supported" << std::endl;
    }

    //*********************************************************************************************
    //Define Cuts
    //*********************************************************************************************
    const Int_t nCuts = 9;
    Bool_t passCut[nCuts] = {kFALSE, kFALSE, kFALSE, kFALSE,
                           kFALSE, kFALSE, kFALSE, kFALSE, kFALSE};
    
    if(leptons->At(0)->Pt() > 20.0 &&
       leptons->At(1)->Pt() > 10.0 &&
       caloMet->Pt()    > 30.0 &&
       dilepton->Mass() > 12.0
      )                              passCut[0] = kTRUE;
    //above cuts are for preselction to be fed into TMVA
    
    if(nCentralJets < 1)     passCut[1] = kTRUE;
    
    if (finalstateType == 10){ // mumu
      if(caloMet->Pt()	> 50.0 &&
         caloMet->Pt()	< 200.0)                 passCut[2] = kTRUE;
      if(deltaPhiLeptons	< 45.0)          passCut[3] = kTRUE;
      if(dilepton->Mass()	< 50.0)          passCut[4] = kTRUE;
      if(leptons->At(0)->Pt()	> 35.0   &&
         leptons->At(0)->Pt()	< 55.0)          passCut[5] = kTRUE;
      if(leptons->At(1)->Pt()	> 25.0)          passCut[6] = kTRUE;
    }
    else if(finalstateType == 11 ){ // ee
      if(caloMet->Pt()	  > 51.0   &&
         caloMet->Pt()	  < 200.0)               passCut[2] = kTRUE;
      if(deltaPhiLeptons	  < 45.0)        passCut[3] = kTRUE;
      if(dilepton->Mass()	  < 40.0)        passCut[4] = kTRUE;
      if(leptons->At(0)->Pt() > 25.0   &&
	 leptons->At(0)->Pt() < 49.0)            passCut[5] = kTRUE;
      if(leptons->At(1)->Pt() > 25.0)            passCut[6] = kTRUE;
    }      
    else if(finalstateType == 12) { //emu
      if(caloMet->Pt()			> 45.0 &&
         caloMet->Pt()			< 105.0) passCut[2] = kTRUE;
      if(deltaPhiLeptons		< 70.0)  passCut[3] = kTRUE;
      if(dilepton->Mass()		< 45.0)  passCut[4] = kTRUE;
      if(leptons->At(0)->Pt()	> 25.0   &&
         leptons->At(0)->Pt()	< 50.0)          passCut[5] = kTRUE;
      if(leptons->At(1)->Pt()	> 25.0)          passCut[6] = kTRUE;
    }
    
    if (DirtyMuons->GetEntries() < 1)            passCut[7] = kTRUE;
    if (CleanExtraTracks->GetEntries() < 4)      passCut[8] = kTRUE;

    //*********************************************************************************************
    //Final Decision
    //*********************************************************************************************
    Bool_t passAllCuts = kTRUE;
    for(Int_t c=0; c<nCuts; c++) passAllCuts = passAllCuts & passCut[c];
    
    //*****************************************************************************************
    //Histograms after no cuts
    //*****************************************************************************************
    fLeptonPtMax->Fill(leptons->At(0)->Pt(),eventweight);
    fLeptonPtMin->Fill(leptons->At(1)->Pt(),eventweight);
    fMetPtHist->Fill(caloMet->Pt(),eventweight);                             
    fDeltaPhiLeptons->Fill(deltaPhiLeptons,eventweight);
    fDeltaEtaLeptons->Fill(deltaEtaLeptons,eventweight);
    fDileptonMass->Fill(dilepton->Mass(),eventweight);    

    //*********************************************************************************************
    // N-1 Histograms
    //*********************************************************************************************

    //N Jet Veto  
    Bool_t pass = kTRUE;
    for (Int_t k=0;k<nCuts;k++) {
      if (k != 1) {
        pass = (pass && passCut[k]);      
      }
    }
    if (pass) {
      fNCentralJets_NMinusOne->Fill(nCentralJets,eventweight);
    }     

    //Met Cut
    pass = kTRUE;
    for (Int_t k=0;k<nCuts;k++) {
      if (k != 2) {
        pass = (pass && passCut[k]);      
      }
    }
    if (pass) {
      fMetPtHist_NMinusOne->Fill(caloMet->Pt(),eventweight);  
    }
    
    //DeltaPhiLeptons
    pass = kTRUE;
    for (Int_t k=0;k<nCuts;k++) {
      if (k != 3) {
        pass = (pass && passCut[k]);      
      }
    }
    if (pass) {
      fDeltaPhiLeptons_NMinusOne->Fill(deltaPhiLeptons,eventweight); 
    }
    
    //dilepton mass
    pass = kTRUE;
    for (Int_t k=0;k<nCuts;k++) {
      if (k != 4)
        pass = (pass && passCut[k]);    
    }
    if (pass) {
      fDileptonMass_NMinusOne->Fill(dilepton->Mass(),eventweight);
    }
    
    //Lepton Pt Max
    pass = kTRUE;
    for (Int_t k=0;k<nCuts;k++) {
      if (k != 5) {
        pass = (pass && passCut[k]);      
      }
    }
    if (pass) {
      fLeptonPtMax_NMinusOne->Fill(leptons->At(0)->Pt(),eventweight);
    }
    
    //Lepton Pt Min
    pass = kTRUE;
    for (Int_t k=0;k<nCuts;k++) {
      if (k != 6) {
        pass = (pass && passCut[k]);      
      }
    }
    if (pass) {
      fLeptonPtMin_NMinusOne->Fill(leptons->At(1)->Pt(),eventweight);
    }
    
    //NDirtyMuons
    pass = kTRUE;
    for (Int_t k=0;k<nCuts;k++) {
      if (k != 7)
        pass = (pass && passCut[k]);    
    }
    if (pass) {
      fNDirtyMuonsHist_NMinusOne->Fill(DirtyMuons->GetEntries(),eventweight);
    }
    
    //NCleanExtraTracks
    pass = kTRUE;
    for (Int_t k=0;k<nCuts;k++) {
      if (k != 8)
        pass = (pass && passCut[k]);    
    }
    if (pass) {
      fNCleanExtraTracksHist_NMinusOne->Fill(CleanExtraTracks->GetEntries(),
                                             eventweight);
    }     

    //*********************************************************************************************
    //Plots after all Cuts
    //*********************************************************************************************
    if (passAllCuts) {
      fMinDeltaPhiLeptonMet_afterCuts->Fill(minDeltaPhiMetLepton,eventweight);
      fMtLepton1_afterCuts->Fill(mTW[0],eventweight);
      fMtLepton2_afterCuts->Fill(mTW[1],eventweight);
      fMtHiggs_afterCuts->Fill(mtHiggs,eventweight);
    }

    if (!fUseMCFake) {
      delete leptons;
      delete CleanJets;
      delete caloMet;
    }
    delete dilepton;
    delete DirtyMuons;
    delete CleanExtraTracks;
  }

  delete FakeableObjects;
  if (fUseMCFake) {
    delete FakeEventHeaders;
  }
  return;
}

//--------------------------------------------------------------------------------------------------
void FakeLeptonExampleAnaMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we typically initialize histograms and other analysis objects and request
  // branches. For this module, we request a branch of the MitTree.

  ReqBranch(fTrackBranchName,  fTracks);
  ReqBranch(fMuonBranchName,   fMuons);

  //Create your histograms here


  //***********************************************************************************************
  // Before preselection
  //***********************************************************************************************
  AddTH1(fDileptonCharge, "hDileptonCharge", ";DileptonCharge;Number of Events", 5, -2.5, 2.5);
  TAxis *xa = fDileptonCharge->GetXaxis();
  xa->SetBinLabel(1,"--");
  xa->SetBinLabel(2,"-+");
  xa->SetBinLabel(4,"+-");
  xa->SetBinLabel(5,"++");

  //***********************************************************************************************
  // Histograms after preselection
  //***********************************************************************************************

  AddTH1(fLeptonPtMax        ,"hLeptonPtMax",";Lepton P_t Max;Number of Events",150,0.,150.);
  AddTH1(fLeptonPtMin        ,"hLeptonPtMin",";Lepton P_t Min;Number of Events",150,0.,150.);
  AddTH1(fMetPtHist          ,"hMetPtHist",";Met;Number of Events",150,0.,300.);  
  AddTH1(fDeltaPhiLeptons    ,"hDeltaPhiLeptons",";#Delta#phi_{ll};Number of Events",90,0,180);
  AddTH1(fDeltaEtaLeptons    ,"hDeltaEtaLeptons",";#Delta#eta_{ll};Number of Events",100,-50.,5.0);
  AddTH1(fDileptonMass       ,"hDileptonMass",";Mass_{ll};Number of Events",150,0.,300.);
 
  //***********************************************************************************************
  // N-1 Histograms
  //***********************************************************************************************
  //All events
  AddTH1(fLeptonPtMax_NMinusOne            ,"hLeptonPtMax_NMinusOne",
                                            ";Lepton P_t Max;Number of Events",150,0.,150.);
  AddTH1(fLeptonPtMin_NMinusOne            ,"hLeptonPtMin_NMinusOne",
                                            ";Lepton P_t Min;Number of Events",150,0.,150.);
  AddTH1(fMetPtHist_NMinusOne              ,"hMetPtHist_NMinusOne",
                                            ";Met;Number of Events",150,0.,300.);  
  AddTH1(fMetPhiHist_NMinusOne             ,"hMetPhiHist_NMinusOne",
                                            ";#phi;Number of Events",28,-3.5,3.5);
  AddTH1(fMETdeltaPhilEtHist_NMinusOne     ,"hMETdeltaPhilEtHist_NMinusOne",
                                            ";METdeltaPhilEtHist;Number of Events",150,0.,300.);
  AddTH1(fNCentralJets_NMinusOne           ,"hNCentralJets_NMinusOne",
                                            ";Number of Central Jets;Number of Events",6,-0.5,5.5);
  AddTH1(fNDirtyMuonsHist_NMinusOne        ,"hNDirtyMuonsHist_NMinusOne",
                                            ";Number of Dirty Muons; Number of Events",6,-0.5,5.5);
  AddTH1(fNCleanExtraTracksHist_NMinusOne  ,"hNCleanExtraTracksHist_NMinusOne",
                                            ";Number of Clean Extra Tracks; Number of Events",
                                            15,-0.5,14.5);
  AddTH1(fDeltaPhiLeptons_NMinusOne        ,"hDeltaPhiLeptons_NMinusOne",
                                            ";#Delta#phi_{ll};Number of Events",90,0,180);
  AddTH1(fDeltaEtaLeptons_NMinusOne        ,"hDeltaEtaLeptons_NMinusOne",
                                            ";#Delta#eta_{ll};Number of Events",100,-5.0,5.0);
  AddTH1(fDileptonMass_NMinusOne           ,"hDileptonMass_NMinusOne",
                                            ";Mass_{ll};Number of Events",150,0.,300.);
  AddTH1(fMinDeltaPhiLeptonMet_NMinusOne   ,"hMinDeltaPhiLeptonMet_NMinusOne", 
                                            ";Min #Delta#phi_{l,Met};Number of Events",90,0.,180);


  //***********************************************************************************************
  // After all cuts Histograms
  //***********************************************************************************************

  AddTH1(fMinDeltaPhiLeptonMet_afterCuts    ,"hMinDeltaPhiLeptonMet_afterCuts", 
                                             ";Min #Delta#phi_{l,Met};Number of Events",90,0.,180);
  AddTH1(fMtLepton1_afterCuts               ,"hMtLepton1_afterCuts",
                                             ";M_t (Lepton1,Met);Number of Events",100,0.,200.);
  AddTH1(fMtLepton2_afterCuts               ,"hMtLepton2_afterCuts",
                                             ";M_t (Lepton2,Met);Number of Events",100,0.,200.);
  AddTH1(fMtHiggs_afterCuts                 ,"hMtHiggs_afterCuts",
                                             ";M_t (l1+l2+Met);Number of Events",150,0.,300.);

}
