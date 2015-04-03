// $Id: GenFakesMod.cc,v 1.6 2009/08/11 11:19:40 phedex Exp $

#include "MitPhysics/FakeMods/interface/GenFakesMod.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/PhotonCol.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/ParticleCol.h"
#include "MitAna/DataTree/interface/MCParticleCol.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/FakeMods/interface/FakeEventHeader.h"
#include "MitPhysics/FakeMods/interface/FakeObject.h"
#include "MitPhysics/FakeMods/interface/FakeRate.h"

using namespace mithep;

ClassImp(mithep::GenFakesMod)

//--------------------------------------------------------------------------------------------------
GenFakesMod::GenFakesMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fElectronFRFilename("InputRequired"),
  fMuonFRFilename("InputRequired"),
  fUse2DFakeRate(false),
  fUseFitFunction(false),
  fElectronFRFunctionName("InputRequired"),
  fMuonFRFunctionName("InputRequired"),
  fElectronFRHistName("InputRequired"),
  fMuonFRHistName("InputRequired"),
  fCleanJetsName(ModNames::gkCleanJetsName),
  fMCLeptonsName(ModNames::gkMCLeptonsName),
  fMCTausName(ModNames::gkMCTausName),
  fElFakeableObjsName(ModNames::gkElFakeableObjsName),
  fMuFakeableObjsName(ModNames::gkMuFakeableObjsName),
  fFakeEventHeadersName(ModNames::gkFakeEventHeadersName)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void GenFakesMod::LoadFakeRate()
{ 
  //Load FakeRate Probabilities.
  fFakeRate = new FakeRate(fElectronFRFilename,fMuonFRFilename,fElectronFRFunctionName,
                           fMuonFRFunctionName,fElectronFRHistName,
                           fMuonFRHistName,fUse2DFakeRate, fUseFitFunction );
}

//--------------------------------------------------------------------------------------------------
void GenFakesMod::Process()
{
  // Process entries of the tree.

  // get input Fakeable object collections
   const ElectronCol *ElFakeableObjs = 0;
   if (!fElFakeableObjsName.IsNull())
     ElFakeableObjs = GetObjThisEvt<ElectronCol>(fElFakeableObjsName);
   const MuonCol *MuFakeableObjs = 0;
   if (!fMuFakeableObjsName.IsNull())
     MuFakeableObjs = GetObjThisEvt<MuonCol>(fMuFakeableObjsName);

  const JetCol      *CleanJets       = 0;
  if (!fCleanJetsName.IsNull())
    CleanJets = GetObjThisEvt<JetCol>(fCleanJetsName);
  mithep::ParticleOArr *CleanLeptons = dynamic_cast<mithep::ParticleOArr*>
    (FindObjThisEvt(ModNames::gkMergedLeptonsName));

  //get monte carlo collections
  const MCParticleCol *GenLeptons = 0;
  if (!fMCLeptonsName.IsNull())
    GenLeptons = GetObjThisEvt<MCParticleCol>(fMCLeptonsName);
  const MCParticleCol *GenTaus = 0;
  if (!fMCTausName.IsNull())
    GenTaus = GetObjThisEvt<MCParticleCol>(fMCTausName);
  ObjArray<MCParticle> *GenLeptonsAndTaus = new ObjArray<MCParticle>;
  if (GenLeptons) {
    for (UInt_t i=0; i<GenLeptons->GetEntries(); i++)
      GenLeptonsAndTaus->Add(GenLeptons->At(i));
  }
  if (GenTaus) {
    for (UInt_t i=0; i<GenTaus->GetEntries(); i++)
      GenLeptonsAndTaus->Add(GenTaus->At(i));
  }

  // create final output collection
  ObjArray <FakeEventHeader> *FakeEventHeaders = new  ObjArray <FakeEventHeader> ;
  FakeEventHeaders->SetOwner(kTRUE);

  //initialize with one fake event containing no fake objects and all jets.
  FakeEventHeader *initialFakeEvent = new FakeEventHeader();
  for (UInt_t j=0;j<CleanJets->GetEntries();j++)
    initialFakeEvent->AddJets(CleanJets->At(j));

  FakeEventHeaders->AddOwned(initialFakeEvent);
  
  // *****************************************************************************************
  // Fake into Muons
  // Loop through all Muon Fakeable objects and consider the fake possibility.
  // *****************************************************************************************
  for (UInt_t n = 0; n < MuFakeableObjs->GetEntries(); ++n) {

    //make temporary fake event headers array
    ObjArray <FakeEventHeader> *tempFakeEventHeaders = new ObjArray <FakeEventHeader> ;
    tempFakeEventHeaders->SetOwner(kTRUE);

    //loop over all fake events generated so far - and perform an additional fake if necessary
    for (UInt_t i=0; i<FakeEventHeaders->GetEntries();i++) {           

      // *****************************************************************************************
      // Determine if the fakeable object was a clean lepton
      // *****************************************************************************************
      Bool_t isCleanLepton = false;
      for (UInt_t j = 0; j < CleanLeptons->GetEntries() ; j++) {
        Double_t deltaR = MathUtils::DeltaR(MuFakeableObjs->At(n)->Phi(),
                                            MuFakeableObjs->At(n)->Eta(),
                                            CleanLeptons->At(j)->Phi(), CleanLeptons->At(j)->Eta());

        if (deltaR < 0.3) {
          isCleanLepton = true;
          break;
        }
      }

      // *****************************************************************************************
      // Determine if the fakeable object was a real muon from Monte Carlo        
      // *****************************************************************************************
      Bool_t isGenMuon = false;
      for (UInt_t l=0; l<GenLeptonsAndTaus->GetEntries(); l++) {
        if (MathUtils::DeltaR(MuFakeableObjs->At(n)->Phi(), MuFakeableObjs->At(n)->Eta(),
                              GenLeptonsAndTaus->At(l)->Phi(), 
                              GenLeptonsAndTaus->At(l)->Eta()) < 0.3) {
          isGenMuon = true;
          break;
        }
      }

      //this is used to determine the weight of the unfaked event.
      double totalCumulativeFakeProbability = 0.0;

      // *****************************************************************************************
      // Perform Muon Fake
      // *****************************************************************************************
      
      //match fake to one of the jets
      int fakeToJetMatch = -1; //index of the jet that matches to the fake
      double minDR = 5000;
      for (UInt_t jj=0;jj<FakeEventHeaders->At(i)->NJets();jj++) {
        Double_t deltaR = MathUtils::DeltaR(FakeEventHeaders->At(i)->UnfakedJet(jj)->Mom(),
                                            MuFakeableObjs->At(n)->Mom());
        if (deltaR < minDR) {
          minDR = deltaR;
          fakeToJetMatch = jj;
        }
      }
      if (!(minDR < 0.5)) {
        fakeToJetMatch = -1;
      }
      
      //Obtain the muon FakeRate 
      Double_t muonFakeProb = 0.0;
      Double_t muonFakeProbStatErrorLow = 0.0;
      Double_t muonFakeProbStatErrorHigh = 0.0;
      Double_t muonFakeProbSysErrorLow = 0.0;
      Double_t muonFakeProbSysErrorHigh = 0.0;
      Double_t muonFakeProbErrorLow = 0.0;
      Double_t muonFakeProbErrorHigh = 0.0;
      if(fFakeRate) {
        muonFakeProb = fFakeRate->MuonFakeRate(MuFakeableObjs->At(n)->Et(),
                                               MuFakeableObjs->At(n)->Eta(),
                                               MuFakeableObjs->At(n)->Phi());
        muonFakeProbStatErrorLow = fFakeRate->MuonFakeRateStatErrorLow(MuFakeableObjs->At(n)->Et(),
                                                            MuFakeableObjs->At(n)->Eta(),
                                                            MuFakeableObjs->At(n)->Phi());
        muonFakeProbStatErrorHigh = 
          fFakeRate->MuonFakeRateStatErrorHigh(MuFakeableObjs->At(n)->Et(),
                                               MuFakeableObjs->At(n)->Eta(),
                                               MuFakeableObjs->At(n)->Phi());        
        muonFakeProbSysErrorLow = fFakeRate->MuonFakeRateSysErrorLow(MuFakeableObjs->At(n)->Et(),
                                                            MuFakeableObjs->At(n)->Eta(),
                                                            MuFakeableObjs->At(n)->Phi());
        muonFakeProbSysErrorHigh = fFakeRate->MuonFakeRateSysErrorHigh(MuFakeableObjs->At(n)->Et(),
                                                                      MuFakeableObjs->At(n)->Eta(),
                                                                      MuFakeableObjs->At(n)->Phi());
        muonFakeProbErrorLow = TMath::Sqrt(muonFakeProbStatErrorLow*muonFakeProbStatErrorLow +
                                           muonFakeProbSysErrorLow*muonFakeProbSysErrorLow);
        muonFakeProbErrorHigh = TMath::Sqrt(muonFakeProbStatErrorHigh*muonFakeProbStatErrorHigh +
                                            muonFakeProbSysErrorHigh*muonFakeProbSysErrorHigh);

      } else {
        Fatal("Process()","Error: fFakeRate is a NULL pointer.");   
      }

      //only fake into a muon if the fakeable object did not match to a clean lepton
      if (!isCleanLepton) {

        //create new fake event header
        FakeEventHeader *fakeMuonEvent = new FakeEventHeader();
        fakeMuonEvent->SetWeight(FakeEventHeaders->At(i)->Weight() * muonFakeProb);
        Double_t weightLowError = 0;
        Double_t weightHighError = 0;
        if (muonFakeProb > 0) { 
          weightLowError = FakeEventHeaders->At(i)->Weight()*muonFakeProb*
            TMath::Sqrt((FakeEventHeaders->At(i)->WeightLowError()/
                         FakeEventHeaders->At(i)->Weight())*
                        (FakeEventHeaders->At(i)->WeightLowError()/
                         FakeEventHeaders->At(i)->Weight()) +
                        (muonFakeProbErrorLow/muonFakeProb)*
                        (muonFakeProbErrorLow/muonFakeProb));
          weightHighError = FakeEventHeaders->At(i)->Weight()*muonFakeProb*
            TMath::Sqrt((FakeEventHeaders->At(i)->WeightHighError()/
                         FakeEventHeaders->At(i)->Weight())*
                        (FakeEventHeaders->At(i)->WeightHighError()/
                         FakeEventHeaders->At(i)->Weight()) +
                        (muonFakeProbErrorHigh/muonFakeProb)*
                        (muonFakeProbErrorHigh/muonFakeProb));
        }
        fakeMuonEvent->SetWeightLowError(weightLowError);
        fakeMuonEvent->SetWeightHighError(weightHighError);

        //add all previous fakes
        for (UInt_t f=0;f<FakeEventHeaders->At(i)->FakeObjsSize();f++) {
          fakeMuonEvent->AddFakeObject(FakeEventHeaders->At(i)->FakeObj(f));        
        }
        //add new fake
        fakeMuonEvent->AddFakeObject(MuFakeableObjs->At(n),kMuon,isCleanLepton,isGenMuon);

        //add all previous jets except the one matching to the new fake
        for (UInt_t jj=0;jj<FakeEventHeaders->At(i)->NJets();jj++) {
          if ((int)jj != fakeToJetMatch)
            fakeMuonEvent->AddJets(FakeEventHeaders->At(i)->UnfakedJet(jj));
        }
        
        //add fake event to the temporary fake event header array
        tempFakeEventHeaders->AddOwned(fakeMuonEvent);
        
        //increase cumulative fake probability
        totalCumulativeFakeProbability += muonFakeProb;
      }

      // *****************************************************************************************
      // Do not fake into Muon
      // *****************************************************************************************
      //create new fake event header
      FakeEventHeader *notFakeEvent = new FakeEventHeader();
      notFakeEvent->SetWeight(FakeEventHeaders->At(i)->Weight() * 
                              (1-totalCumulativeFakeProbability));
      //add previous fakes
      for (UInt_t f=0;f<FakeEventHeaders->At(i)->FakeObjsSize();f++) {
        notFakeEvent->AddFakeObject(FakeEventHeaders->At(i)->FakeObj(f));        
      }
      //add previous jets
      for (UInt_t jj=0;jj<FakeEventHeaders->At(i)->NJets();jj++) {
          notFakeEvent->AddJets(FakeEventHeaders->At(i)->UnfakedJet(jj));
      }
      tempFakeEventHeaders->AddOwned(notFakeEvent);   

    } //loop over all current fake event headers

    //replace current fake event headers with the new temporary ones.
    delete FakeEventHeaders;
    FakeEventHeaders = tempFakeEventHeaders;
  } //loop over all muon fakeable objects


  // *****************************************************************************************
  // Fake into Electrons
  // Loop through all Electron Fakeable objects and consider the fake possibility.
  // *****************************************************************************************
  for (UInt_t n = 0; n < ElFakeableObjs->GetEntries();n++) {

    //make temporary fake event headers array
    ObjArray <FakeEventHeader> *tempFakeEventHeaders = new ObjArray <FakeEventHeader> ;
    tempFakeEventHeaders->SetOwner(kTRUE);

    //loop over all fake events generated so far - and perform an additional fake if necessary
    for (UInt_t i=0; i<FakeEventHeaders->GetEntries();i++) {      
      
      // *****************************************************************************************
      // Determine if the fakeable object was a clean lepton
      // *****************************************************************************************
      Bool_t isCleanLepton = false;
      for (UInt_t j = 0; j < CleanLeptons->GetEntries() ; j++) {
        Double_t deltaR = MathUtils::DeltaR(ElFakeableObjs->At(n)->Phi(),
                                            ElFakeableObjs->At(n)->Eta(),
                                            CleanLeptons->At(j)->Phi(), CleanLeptons->At(j)->Eta());

        if (deltaR < 0.3) {
          isCleanLepton = true;
          break;
        }
      }

      // *****************************************************************************************
      // Determine if the fakeable object was a real electron from Monte Carlo        
      // *****************************************************************************************
      Bool_t isGenElectron = false;
      for (UInt_t l=0; l<GenLeptonsAndTaus->GetEntries(); l++) {
        if (MathUtils::DeltaR(ElFakeableObjs->At(n)->Phi(), 
                              ElFakeableObjs->At(n)->Eta(),
                              GenLeptonsAndTaus->At(l)->Phi(), 
                              GenLeptonsAndTaus->At(l)->Eta()) < 0.3) {
          isGenElectron = true;
        }
      }


      //this is used to determine the weight of the unfaked event.
      double totalCumulativeFakeProbability = 0.0;
      
      // *****************************************************************************************
      // Determine if the current electron fakeable object already corresponds to one of the
      // fake muons already in the FakeEventHeader, determined based on deltaR proximity.
      // If the current electron fakeable object corresponds to one of the fake muon, then
      // we do not allow it to fake an electron, since one denominator cannot fake two leptons.
      // *****************************************************************************************
      Bool_t alreadyFaked = false;
      for (UInt_t f = 0; f < FakeEventHeaders->At(i)->FakeObjsSize() ; f++) {
        double deltaR = MathUtils::DeltaR(FakeEventHeaders->At(i)->FakeObj(f)->Mom(),
                                          ElFakeableObjs->At(n)->Mom());
        if (deltaR < 0.3)
          alreadyFaked = true;
      }
      if (!alreadyFaked) {

        // *****************************************************************************************
        // Perform Electron Fake
        // *****************************************************************************************
  
        //match fake to one of the jets
        int fakeToJetMatch = -1; //index of the jet that matches to the fake
        double minDR = 5000;
        for (UInt_t jj=0;jj<FakeEventHeaders->At(i)->NJets();jj++) {          
          Double_t deltaR = MathUtils::DeltaR(FakeEventHeaders->At(i)->UnfakedJet(jj)->Mom(),
                                              ElFakeableObjs->At(n)->Mom());
          if (deltaR < minDR) {
            minDR = deltaR;
            fakeToJetMatch = jj;
          }
        }
        if (!(minDR < 0.5)) { 
          fakeToJetMatch = -1;
        }

        //Obtain the electron FakeRate 
        Double_t electronFakeProb = 0.0;
        Double_t electronFakeProbStatErrorLow = 0.0;
        Double_t electronFakeProbStatErrorHigh = 0.0;
        Double_t electronFakeProbSysErrorLow = 0.0;
        Double_t electronFakeProbSysErrorHigh = 0.0;
        Double_t electronFakeProbErrorLow = 0.0;
        Double_t electronFakeProbErrorHigh = 0.0;
        if(fFakeRate) {
          electronFakeProb = fFakeRate->ElectronFakeRate(ElFakeableObjs->At(n)->Et(),
                                                         ElFakeableObjs->At(n)->Eta(),
                                                         ElFakeableObjs->At(n)->Phi());
          electronFakeProbStatErrorLow = 
            fFakeRate->ElectronFakeRateStatErrorLow(ElFakeableObjs->At(n)->Et(),
                                                    ElFakeableObjs->At(n)->Eta(),
                                                    ElFakeableObjs->At(n)->Phi());
          electronFakeProbStatErrorHigh = 
            fFakeRate->ElectronFakeRateStatErrorHigh(ElFakeableObjs->At(n)->Et(),
                                                     ElFakeableObjs->At(n)->Eta(),
                                                     ElFakeableObjs->At(n)->Phi());                
          electronFakeProbSysErrorLow = 
            fFakeRate->ElectronFakeRateSysErrorLow(ElFakeableObjs->At(n)->Et(),
                                                   ElFakeableObjs->At(n)->Eta(),
                                                   ElFakeableObjs->At(n)->Phi());
          electronFakeProbSysErrorHigh = 
            fFakeRate->ElectronFakeRateSysErrorHigh(ElFakeableObjs->At(n)->Et(),
                                                    ElFakeableObjs->At(n)->Eta(),
                                                    ElFakeableObjs->At(n)->Phi());                
          electronFakeProbErrorLow = 
            TMath::Sqrt(electronFakeProbStatErrorLow*electronFakeProbStatErrorLow +
                        electronFakeProbSysErrorLow*electronFakeProbSysErrorLow);
          electronFakeProbErrorHigh = 
            TMath::Sqrt(electronFakeProbStatErrorHigh*electronFakeProbStatErrorHigh +
                        electronFakeProbSysErrorHigh*electronFakeProbSysErrorHigh);
        } else {
          Fatal("Process()","Error: fFakeRate is a NULL pointer.");   
        }
        
        //only fake into an electron if the fakeable object did not match to a clean lepton
        if (!isCleanLepton) {
          //create new fake event header
          FakeEventHeader *fakeElectronEvent = new FakeEventHeader();
          fakeElectronEvent->SetWeight(FakeEventHeaders->At(i)->Weight() * electronFakeProb);
          Double_t weightLowError = 0;
          Double_t weightHighError = 0;
          if (electronFakeProb) {
            weightLowError = FakeEventHeaders->At(i)->Weight()*electronFakeProb*
              TMath::Sqrt((FakeEventHeaders->At(i)->WeightLowError()/
                           FakeEventHeaders->At(i)->Weight())*
                          (FakeEventHeaders->At(i)->WeightLowError()/
                           FakeEventHeaders->At(i)->Weight()) +
                          (electronFakeProbErrorLow/electronFakeProb)*
                          (electronFakeProbErrorLow/electronFakeProb));
            weightHighError = FakeEventHeaders->At(i)->Weight()*electronFakeProb*
              TMath::Sqrt((FakeEventHeaders->At(i)->WeightHighError()/
                           FakeEventHeaders->At(i)->Weight())*
                          (FakeEventHeaders->At(i)->WeightHighError()/
                           FakeEventHeaders->At(i)->Weight()) +
                          (electronFakeProbErrorHigh/electronFakeProb)*
                          (electronFakeProbErrorHigh/electronFakeProb));
          }
          fakeElectronEvent->SetWeightLowError(weightLowError);
          fakeElectronEvent->SetWeightHighError(weightHighError);
          
          //add previous fakes
          for (UInt_t f=0;f<FakeEventHeaders->At(i)->FakeObjsSize();f++) {
            fakeElectronEvent->AddFakeObject(FakeEventHeaders->At(i)->FakeObj(f));
          }
          //add the new fake
          fakeElectronEvent->AddFakeObject(ElFakeableObjs->At(n),
                                           kElectron,isCleanLepton,isGenElectron);
          //add previous jets that do not match to the new fake
          for (UInt_t jj=0;jj<FakeEventHeaders->At(i)->NJets();jj++) {
            if ((int)jj != fakeToJetMatch)
              fakeElectronEvent->AddJets(FakeEventHeaders->At(i)->UnfakedJet(jj));
          }
          
          //add fake event to the temporary fake event header array
          tempFakeEventHeaders->AddOwned(fakeElectronEvent);
          //increase cumulative fake probability
          totalCumulativeFakeProbability += electronFakeProb;
        }
      }
      
      // *****************************************************************************************
      // Do not fake into anything
      // *****************************************************************************************
      //create new fake event header
      FakeEventHeader *notFakeEvent = new FakeEventHeader();
      notFakeEvent->SetWeight(FakeEventHeaders->At(i)->Weight() * 
                              (1-totalCumulativeFakeProbability));
      //add previous fakes
      for (UInt_t f=0;f<FakeEventHeaders->At(i)->FakeObjsSize();f++) {
        notFakeEvent->AddFakeObject(FakeEventHeaders->At(i)->FakeObj(f));        
      }
      //add previous jets
      for (UInt_t jj=0;jj<FakeEventHeaders->At(i)->NJets();jj++) {
          notFakeEvent->AddJets(FakeEventHeaders->At(i)->UnfakedJet(jj));
      }
      tempFakeEventHeaders->AddOwned(notFakeEvent);      

    } //for all current fake event headers

    //replace current fake event headers with the new temporary ones.
    delete FakeEventHeaders;
    FakeEventHeaders = tempFakeEventHeaders;
  } //loop over all fakeable objects

  // export FakeEventHeaders for other modules to use
  FakeEventHeaders->SetName(fFakeEventHeadersName);
  AddObjThisEvt(FakeEventHeaders);

  //delete temporary collections
  delete GenLeptonsAndTaus;
}
