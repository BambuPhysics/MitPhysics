// $Id: GenFakeableObjsMod.cc,v 1.15 2011/02/17 13:44:54 bendavid Exp $

#include "MitPhysics/FakeMods/interface/GenFakeableObjsMod.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/SuperClusterCol.h"
#include "MitAna/DataTree/interface/TrackCol.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/DataTree/interface/DecayParticleCol.h"
#include "MitAna/DataTree/interface/StableData.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Utils/interface/IsolationTools.h"
#include "MitPhysics/Mods/interface/ElectronIDMod.h"

using namespace mithep;

ClassImp(mithep::GenFakeableObjsMod)

//--------------------------------------------------------------------------------------------------
GenFakeableObjsMod::GenFakeableObjsMod(const char *name, const char *title) : 
  BaseMod(name,title),

  fApplyConvFilter(kTRUE),
  fNWrongHitsMax(1),
  fApplyD0Cut(kTRUE),
  fChargeFilter(kTRUE),
  fD0Cut(0.02),
  fCombIsolationCut(0.5),
  fTrackIsolationCut(-1.0),
  fEcalIsolationCut(-1.0),
  fHcalIsolationCut(-1.0),
  fVetoTriggerJet(kFALSE),
  fVetoGenLeptons(kTRUE),
  fVetoCleanLeptons(kFALSE),
  fElectronFOType("Iso"),
  fMuonFOType("IsoTrack"),
  fTriggerName("NotSpecified"),
  fTriggerObjectsName("NotSpecified"),
  fElectronBranchName(Names::gkElectronBrn),
  fMuonBranchName(Names::gkMuonBrn),
  fTrackBranchName(Names::gkTrackBrn),
  fGsfTrackBranchName(Names::gkGsfTrackBrn),
  fBarrelSuperClusterBranchName(Names::gkBarrelSuperClusterBrn),
  fEndcapSuperClusterBranchName(Names::gkEndcapSuperClusterBrn),
  fVertexName(ModNames::gkGoodVertexesName),
  fConversionBranchName(Names::gkMvfConversionBrn),
  fMCLeptonsName(ModNames::gkMCLeptonsName),
  fMCTausName(ModNames::gkMCTausName),
  fElFakeableObjsName(ModNames::gkElFakeableObjsName),
  fMuFakeableObjsName(ModNames::gkMuFakeableObjsName),
  fElFOType(kElFOUndef),
  fMuFOType(kMuFOUndef),
  fBarrelSuperClusters(0),
  fEndcapSuperClusters(0),
  fTracks(0),
  fGsfTracks(0),
  fVertices(0),
  fConversions(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void GenFakeableObjsMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we typically initialize histograms and other analysis objects and request
  // branches. For this module, we request a branch of the MitTree.

  //------------------------------------------------------------------------------------------------
  // Request the branches (no significant time safed by not doing this)
  //------------------------------------------------------------------------------------------------
  ReqBranch(fElectronBranchName,              fElectrons);
  ReqBranch(fMuonBranchName,                  fMuons);
  ReqBranch(fTrackBranchName,                 fTracks);
  ReqBranch(fGsfTrackBranchName,              fGsfTracks);
  ReqBranch(fBarrelSuperClusterBranchName,    fBarrelSuperClusters);
  ReqBranch(fEndcapSuperClusterBranchName,    fEndcapSuperClusters);
  ReqBranch(fConversionBranchName,            fConversions);

  if (fElectronFOType.CompareTo("Iso") == 0) 
    fElFOType = kElFOIso;
  else if (fElectronFOType.CompareTo("LooseIdLooseIso") == 0) 
    fElFOType = kElFOLooseIdLooseIso;
  else {
    SendError(kAbortAnalysis, "SlaveBegin",
              "The specified electron fakeable object %s is not defined.",
              fElectronFOType.Data());
    return;
  }

  if (fMuonFOType.CompareTo("IsoTrack") == 0) 
    fMuFOType = kMuFOIsoTrack;
  else if (fMuonFOType.CompareTo("Global") == 0) 
    fMuFOType = kMuFOGlobal;
  else if (fMuonFOType.CompareTo("TrackerMuon") == 0) 
    fMuFOType = kMuFOTrackerMuon;
  else {
    SendError(kAbortAnalysis, "SlaveBegin",
              "The specified muon fakeable object %s is not defined.",
              fMuonFOType.Data());
    return;
  }

  electronID = new ElectronIDMod();
  electronID->SetApplyConversionFilterType1(kFALSE);
  electronID->SetApplyConversionFilterType1(fApplyConvFilter);    
  electronID->SetNWrongHitsMax(fNWrongHitsMax);    
  electronID->SetApplyD0Cut(fApplyD0Cut);    
  electronID->SetChargeFilter(fChargeFilter);    
  electronID->SetD0Cut(fD0Cut);

}

//--------------------------------------------------------------------------------------------------
void GenFakeableObjsMod::Process()
{
  // Process entries of the tree.
  LoadBranch(fElectronBranchName);
  LoadBranch(fMuonBranchName);
  LoadBranch(fTrackBranchName);
  LoadBranch(fGsfTrackBranchName);
  LoadBranch(fBarrelSuperClusterBranchName);
  LoadBranch(fEndcapSuperClusterBranchName);
  LoadBranch(fConversionBranchName);

  fVertices = GetObjThisEvt<VertexOArr>(fVertexName);

  //Load Trigger Objects
  const TriggerObjectCol *triggerObjects = GetHLTObjects(fTriggerObjectsName);
  if (!triggerObjects && fVetoTriggerJet) {
    cout << "Error: Could not load Trigger Object Collection with name " 
         << fTriggerObjectsName << endl;
  }

  // get input clean object collections
  mithep::ParticleOArr *CleanLeptons = dynamic_cast<mithep::ParticleOArr*>
    (FindObjThisEvt(ModNames::gkMergedLeptonsName));

  //get input MC collections
  const MCParticleCol *GenLeptons = 0;
  if (!fMCLeptonsName.IsNull())
    GenLeptons = GetObjThisEvt<MCParticleCol>(fMCLeptonsName);
  const MCParticleCol *GenTaus = 0;
  if (!fMCTausName.IsNull())
    GenTaus = GetObjThisEvt<MCParticleCol>(fMCTausName);
  ObjArray<MCParticle> *GenLeptonsAndTaus = new ObjArray<MCParticle>;
  for (UInt_t i=0; i<GenLeptons->GetEntries(); i++)
    GenLeptonsAndTaus->Add(GenLeptons->At(i));
  for (UInt_t i=0; i<GenTaus->GetEntries(); i++)
    GenLeptonsAndTaus->Add(GenTaus->At(i));

  //Combine Barrel and Endcap superclusters into the same ObjArray
  ObjArray <SuperCluster> *SuperClusters = new ObjArray <SuperCluster>;
  for (UInt_t i=0; i<fBarrelSuperClusters->GetEntries(); i++) {
    SuperClusters->Add(fBarrelSuperClusters->At(i));
  }
  for (UInt_t i=0; i<fEndcapSuperClusters->GetEntries(); i++) {
    SuperClusters->Add(fEndcapSuperClusters->At(i));
  }

  //collections for duplicate removed electrons
  ObjArray<Electron> *DuplicateRemovedElectrons = new ObjArray<Electron>;
  std::vector<const Electron*> tmpDuplicateRemovedElectrons;

  // create final output collection
  ElectronArr *ElFakeableObjs = new ElectronArr;
  MuonArr *MuFakeableObjs = new MuonArr;


  //***********************************************************************************************
  //First do duplicate electron removal
  //***********************************************************************************************
  for (UInt_t i=0; i<fElectrons->GetEntries(); ++i) {  


    const Electron *e = fElectrons->At(i);   
    Bool_t isElectronOverlap = kFALSE;

    for (UInt_t j=0; j<tmpDuplicateRemovedElectrons.size(); ++j) {
      // why does someone calculate it if it is not used?
      //Double_t deltaR = MathUtils::DeltaR(tmpDuplicateRemovedElectrons[j]->Mom(), e->Mom());
      if (e->SCluster() == tmpDuplicateRemovedElectrons[j]->SCluster() ||
          e->GsfTrk() == tmpDuplicateRemovedElectrons[j]->GsfTrk()) {
        isElectronOverlap = kTRUE;
      }
    
    
      if (isElectronOverlap) {
        if (TMath::Abs(tmpDuplicateRemovedElectrons[j]->ESuperClusterOverP() - 1) > 
            TMath::Abs(e->ESuperClusterOverP() - 1)) {	   
          tmpDuplicateRemovedElectrons[j] = e;
        }      
        break;
      }
    }
    
    if (!isElectronOverlap) {
      tmpDuplicateRemovedElectrons.push_back(e);
    }
  }
  for (UInt_t i=0; i<tmpDuplicateRemovedElectrons.size(); ++i) {  
    DuplicateRemovedElectrons->Add(tmpDuplicateRemovedElectrons[i]);
  }

  //***********************************************************************************************
  //Fakeable Objects for Electron Fakes
  //Reco electron with full isolation
  //***********************************************************************************************
  if (fElFOType == kElFOIso) {

    for (UInt_t i=0; i<DuplicateRemovedElectrons->GetEntries(); i++) {  
      const Electron *denominator = DuplicateRemovedElectrons->At(i);

      //Veto denominators matching to real electrons      
      Bool_t IsGenLepton = false;
      for (UInt_t l=0; l<GenLeptonsAndTaus->GetEntries(); l++) {
        if (MathUtils::DeltaR(denominator->Phi(), denominator->Eta(),
                              GenLeptonsAndTaus->At(l)->Phi(), 
                              GenLeptonsAndTaus->At(l)->Eta()) < 0.1) {
          IsGenLepton = true;
        }
      }

      //Veto denominators matching to clean leptons
      Bool_t IsCleanLepton = false;
      for (UInt_t l=0; l<CleanLeptons->GetEntries(); l++) {
        if (MathUtils::DeltaR(denominator->Phi(), denominator->Eta(),
                              CleanLeptons->At(l)->Phi(), 
                              CleanLeptons->At(l)->Eta()) < 0.1) {
          IsCleanLepton = true;
          }
      }
      
      //Veto on Trigger jet
      Bool_t IsTriggerJet = false;
      if (fVetoTriggerJet) {
        for (UInt_t l=0; l<triggerObjects->GetEntries(); l++) {      
          Double_t deltaR = MathUtils::DeltaR(denominator->Phi(), 
                                              denominator->Eta(),
                                              triggerObjects->At(l)->Phi(), 
                                              triggerObjects->At(l)->Eta());
          if (triggerObjects->At(l)->TrigName() == fTriggerName.Data() 
              && triggerObjects->At(l)->Type() == TriggerObject::TriggerJet
              && deltaR < 0.3
            ) {
            IsTriggerJet = true;
            break;
          }
        }
      }

      const Electron *tmpEle = denominator;
      //****************************************************************************************
      // Isolation Cut
      //****************************************************************************************
      Double_t combIso = 
        denominator->TrackIsolationDr03() + TMath::Max(denominator->EcalRecHitIsoDr03() - 1.0, 0.0) + denominator->HcalTowerSumEtDr03();
      if (fabs(denominator->Eta()) > 1.5) {
        combIso = denominator->TrackIsolationDr03() + denominator->EcalRecHitIsoDr03() + denominator->HcalTowerSumEtDr03();
      }

      Bool_t passIsolationCut = (combIso / denominator->Pt() < 0.1) ;
      
      //****************************************************************************************
      // conversion filter
      //****************************************************************************************
      Bool_t  passConversionFilter = (TMath::Abs(denominator->ConvPartnerDCotTheta()) >= 0.02 || 
                                      TMath::Abs(denominator->ConvPartnerDist()) >= 0.02);
      
      //****************************************************************************************
      // D0 Cut        
      //****************************************************************************************
      Bool_t passD0Cut = ElectronTools::PassD0Cut(tmpEle,fVertices, kTRUE);
      
      //****************************************************************************************
      // Make denominator object cuts
      //****************************************************************************************
      if( denominator->Pt() > 10.0  
          && passIsolationCut
          && (passConversionFilter || !fApplyConvFilter)
          && (passD0Cut || !fApplyD0Cut)
          && !(fVetoCleanLeptons && IsCleanLepton)
          && !(fVetoGenLeptons && IsGenLepton)
          && !(fVetoTriggerJet && IsTriggerJet)
        ) {        
        Electron *tmpElectron = ElFakeableObjs->AddNew();
        tmpElectron->SetPtEtaPhi(denominator->Pt(), denominator->Eta(),denominator->Phi());
        tmpElectron->SetGsfTrk(denominator->GsfTrk());
        tmpElectron->SetSuperCluster(denominator->SCluster());        
      } 
    }
  } else if (fElFOType == kElFOLooseIdLooseIso) {
    for (UInt_t i=0; i<DuplicateRemovedElectrons->GetEntries(); i++) {  
      const Electron *denominator = DuplicateRemovedElectrons->At(i);
      
      //Veto denominators matching to real electrons      
      Bool_t IsGenLepton = false;
      for (UInt_t l=0; l<GenLeptonsAndTaus->GetEntries(); l++) {
        if (MathUtils::DeltaR(denominator->Phi(), denominator->Eta(),
                              GenLeptonsAndTaus->At(l)->Phi(), 
                              GenLeptonsAndTaus->At(l)->Eta()) < 0.1) {
          IsGenLepton = true;
        }
      }

      //Veto denominators matching to clean leptons
      Bool_t IsCleanLepton = false;
      for (UInt_t l=0; l<CleanLeptons->GetEntries(); l++) {
        if (MathUtils::DeltaR(denominator->Phi(), denominator->Eta(),
                              CleanLeptons->At(l)->Phi(), 
                              CleanLeptons->At(l)->Eta()) < 0.1) {
          IsCleanLepton = true;
          }
      }

      //Veto on Trigger jet
      Bool_t IsTriggerJet = false;
      if (fVetoTriggerJet) {
        for (UInt_t l=0; l<triggerObjects->GetEntries(); l++) {      
          Double_t deltaR = MathUtils::DeltaR(denominator->Phi(), 
                                              denominator->Eta(),
                                              triggerObjects->At(l)->Phi(), 
                                              triggerObjects->At(l)->Eta());
          if (triggerObjects->At(l)->TrigName() == fTriggerName.Data() 
              && triggerObjects->At(l)->Type() == TriggerObject::TriggerJet
              && deltaR < 0.3
            ) {
            IsTriggerJet = true;
            break;
          }
        }
      }

      const Electron *tmpEle = denominator;
      //****************************************************************************************
      // Id Cuts
      //****************************************************************************************
      Bool_t passIdCut = ElectronTools::PassCustomID(denominator, ElectronTools::kVBTFWorkingPoint90Id);
 

     //****************************************************************************************
      // Isolation Cut
      //****************************************************************************************
      Double_t combIso = 
        denominator->TrackIsolationDr03() + TMath::Max(denominator->EcalRecHitIsoDr03() - 1.0, 0.0) + denominator->HcalTowerSumEtDr03();
      if (fabs(denominator->Eta()) > 1.5) {
        combIso = denominator->TrackIsolationDr03() + denominator->EcalRecHitIsoDr03() + denominator->HcalTowerSumEtDr03();
      }
      Bool_t passIsolationCut = (combIso / denominator->Pt() < 0.3) ;
  
      //****************************************************************************************
      // conversion filter
      //****************************************************************************************     
      Bool_t  passConversionFilter = (TMath::Abs(denominator->ConvPartnerDCotTheta()) >= 0.02 || 
                                      TMath::Abs(denominator->ConvPartnerDist()) >= 0.02);

      //****************************************************************************************
      // D0 Cut        
      //****************************************************************************************
      Bool_t passD0Cut = ElectronTools::PassD0Cut(tmpEle,fVertices, kTRUE);
      
      //****************************************************************************************
      // Make denominator object cuts
      //****************************************************************************************
      if( denominator->Pt() > 10.0  
          && passIdCut
          && passIsolationCut
          && (passConversionFilter || !fApplyConvFilter)
          && (passD0Cut || !fApplyD0Cut)
          && denominator->PassLooseID()
          && !(fVetoCleanLeptons && IsCleanLepton)
          && !(fVetoGenLeptons && IsGenLepton)
          && !(fVetoTriggerJet && IsTriggerJet)
        ) {        
        Electron *tmpElectron = ElFakeableObjs->AddNew();
        tmpElectron->SetPtEtaPhi(denominator->Pt(), denominator->Eta(),denominator->Phi());
        tmpElectron->SetGsfTrk(denominator->GsfTrk());
        tmpElectron->SetSuperCluster(denominator->SCluster());
      }
    }
  }

  //***********************************************************************************************
  //Fakeable Objects for Muon Fakes
  //***********************************************************************************************
  if (fMuFOType == kMuFOIsoTrack) {
    for (UInt_t i=0; i<fTracks->GetEntries(); i++) {
      const Track *track = fTracks->At(i);
      Double_t trackIsolation = IsolationTools::TrackIsolation(track, 0.4, 0.015, 1.0, 
                                                               0.2, fTracks);
      //Determine if muon fakeable object matches a gen lepton
      Bool_t IsGenLepton = false;
      for (UInt_t l=0; l<GenLeptonsAndTaus->GetEntries(); l++) {
        if (MathUtils::DeltaR(track->Phi(), track->Eta(),
                              GenLeptonsAndTaus->At(l)->Phi(), 
                              GenLeptonsAndTaus->At(l)->Eta()) < 0.3) {
          IsGenLepton = true;
        }
      }
      
      //Veto denominators matching to clean leptons
      Bool_t IsCleanLepton = false;
      for (UInt_t l=0; l<CleanLeptons->GetEntries(); l++) {
        if (MathUtils::DeltaR(track->Phi(), track->Eta(),
                              CleanLeptons->At(l)->Phi(), 
                              CleanLeptons->At(l)->Eta()) < 0.1) {
          IsCleanLepton = true;
          }
      }

      //Veto on Trigger jet
      Bool_t IsTriggerJet = false;
      if (fVetoTriggerJet) {
        for (UInt_t l=0; l<triggerObjects->GetEntries(); l++) {      
          Double_t deltaR = MathUtils::DeltaR(track->Phi(), track->Eta(),
                                              triggerObjects->At(l)->Phi(), 
                                              triggerObjects->At(l)->Eta());
          if (triggerObjects->At(l)->TrigName() == fTriggerName.Data() 
              && triggerObjects->At(l)->Type() == TriggerObject::TriggerJet
              && deltaR < 0.3
            ) {
            IsTriggerJet = true;
            break;
          }
        }
      }
      
      //****************************************************************************************
      // D0 Cut        
      //****************************************************************************************
      double d0Min = 99999;
      for(UInt_t i0 = 0; i0 < fVertices->GetEntries(); i0++) {
        double pD0 = track->D0Corrected(*fVertices->At(i0));
        if(TMath::Abs(pD0) < TMath::Abs(d0Min)) d0Min = TMath::Abs(pD0);
      }


      //define denominator cuts
      if (track->Pt() > 10.0 && trackIsolation < 10.0
          && d0Min < 0.025
          && !(fVetoGenLeptons && IsGenLepton)  
          && !(fVetoCleanLeptons && IsCleanLepton)  
          && !(fVetoTriggerJet && IsTriggerJet)
        ) {
        //add to fakeable objects
        Muon* tmpMuon = MuFakeableObjs->AddNew();
        tmpMuon->SetTrackerTrk(track);
      }
    }
  } else if (fMuFOType == kMuFOGlobal) {
    for (UInt_t i=0; i<fMuons->GetEntries(); i++) {
      const Muon *denominator = fMuons->At(i);
      Double_t totalIsolation = denominator->IsoR03SumPt() + denominator->IsoR03EmEt() + 
        denominator->IsoR03HadEt();

      //Determine if muon fakeable object matches a gen lepton
      Bool_t IsGenLepton = false;
      for (UInt_t l=0; l<GenLeptonsAndTaus->GetEntries(); l++) {
        if (MathUtils::DeltaR(denominator->Phi(), denominator->Eta(),
                              GenLeptonsAndTaus->At(l)->Phi(), 
                              GenLeptonsAndTaus->At(l)->Eta()) < 0.3) {
          IsGenLepton = true;
        }
      }
      
      //Veto denominators matching to clean leptons
      Bool_t IsCleanLepton = false;
      for (UInt_t l=0; l<CleanLeptons->GetEntries(); l++) {
        if (MathUtils::DeltaR(denominator->Phi(), denominator->Eta(),
                              CleanLeptons->At(l)->Phi(), 
                              CleanLeptons->At(l)->Eta()) < 0.1) {
          IsCleanLepton = true;
          }
      }

      //Veto on Trigger jet
      Bool_t IsTriggerJet = false;
      if (fVetoTriggerJet) {
        for (UInt_t l=0; l<triggerObjects->GetEntries(); l++) {      
          Double_t deltaR = MathUtils::DeltaR(denominator->Phi(), 
                                              denominator->Eta(),
                                              triggerObjects->At(l)->Phi(), 
                                              triggerObjects->At(l)->Eta());
          if (triggerObjects->At(l)->TrigName() == fTriggerName.Data() 
              && triggerObjects->At(l)->Type() == TriggerObject::TriggerJet
              && deltaR < 0.3
            ) {
            IsTriggerJet = true;
            break;
          }
        }
      }

      //****************************************************************************************
      // D0 Cut        
      //****************************************************************************************
      double d0Min = 99999;
      for(UInt_t i0 = 0; i0 < fVertices->GetEntries(); i0++) {
        double pD0 = denominator->Trk()->D0Corrected(*fVertices->At(i0));
        if(TMath::Abs(pD0) < TMath::Abs(d0Min)) d0Min = TMath::Abs(pD0);
      }
      
      if (denominator->Pt() > 10.0 && totalIsolation < 10.0 && denominator->HasGlobalTrk()
          && d0Min < 0.025
          && !(fVetoGenLeptons && IsGenLepton)
          && !(fVetoCleanLeptons && IsCleanLepton)
          && !(fVetoTriggerJet && IsTriggerJet)
        ) {
        Muon* tmpMuon = MuFakeableObjs->AddNew();
        tmpMuon->SetGlobalTrk(denominator->GlobalTrk());
        tmpMuon->SetTrackerTrk(denominator->TrackerTrk());
      }
    }
  } else if (fMuFOType == kMuFOTrackerMuon) {
    for (UInt_t i=0; i<fMuons->GetEntries(); i++) {
      const Muon *denominator = fMuons->At(i);
      Double_t totalIsolation = 
        denominator->IsoR03SumPt() + 
        denominator->IsoR03EmEt() + 
        denominator->IsoR03HadEt();

      //Determine if muon fakeable object matches a gen lepton
      Bool_t IsGenLepton = false;
      for (UInt_t l=0; l<GenLeptonsAndTaus->GetEntries(); l++) {
        if (MathUtils::DeltaR(denominator->Phi(), denominator->Eta(),
                              GenLeptonsAndTaus->At(l)->Phi(), 
                              GenLeptonsAndTaus->At(l)->Eta()) < 0.3) {
          IsGenLepton = true;
        }
      }

      //Veto denominators matching to clean leptons
      Bool_t IsCleanLepton = false;
      for (UInt_t l=0; l<CleanLeptons->GetEntries(); l++) {
        if (MathUtils::DeltaR(denominator->Phi(), denominator->Eta(),
                              CleanLeptons->At(l)->Phi(), 
                              CleanLeptons->At(l)->Eta()) < 0.1) {
          IsCleanLepton = true;
          }
      }

      //Veto on Trigger jet
      Bool_t IsTriggerJet = false;
      if (fVetoTriggerJet) {
        for (UInt_t l=0; l<triggerObjects->GetEntries(); l++) {      
          Double_t deltaR = MathUtils::DeltaR(denominator->Phi(), 
                                              denominator->Eta(),
                                              triggerObjects->At(l)->Phi(), 
                                              triggerObjects->At(l)->Eta());
          if (triggerObjects->At(l)->TrigName() == fTriggerName.Data() 
              && triggerObjects->At(l)->Type() == TriggerObject::TriggerJet
              && deltaR < 0.3
            ) {
            IsTriggerJet = true;
            break;
          }
        }
      }

      //****************************************************************************************
      // D0 Cut        
      //****************************************************************************************
      double d0Min = 99999;
      for(UInt_t i0 = 0; i0 < fVertices->GetEntries(); i0++) {
        double pD0 = denominator->Trk()->D0Corrected(*fVertices->At(i0));
        if(TMath::Abs(pD0) < TMath::Abs(d0Min)) d0Min = TMath::Abs(pD0);
      }

      if (denominator->Pt() > 10.0 && totalIsolation < 10.0 && denominator->HasTrackerTrk()
          && d0Min < 0.025
          && !(fVetoGenLeptons && IsGenLepton)
          && !(fVetoCleanLeptons && IsCleanLepton)
          && !(fVetoTriggerJet && IsTriggerJet)
        ) {
        Muon* tmpMuon = MuFakeableObjs->AddNew();
        tmpMuon->SetTrackerTrk(denominator->TrackerTrk());
      }
    }
  }

  //***********************************************************************************************
  //Export the fakeable object collections for other modules to use
  //***********************************************************************************************
  ElFakeableObjs->SetName(fElFakeableObjsName);
  AddObjThisEvt(ElFakeableObjs);
  MuFakeableObjs->SetName(fMuFakeableObjsName);
  AddObjThisEvt(MuFakeableObjs);

  delete GenLeptonsAndTaus;
  delete SuperClusters;
  delete DuplicateRemovedElectrons;
}
