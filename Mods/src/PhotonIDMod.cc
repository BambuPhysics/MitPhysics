#include "TDataMember.h"
#include "TTree.h"
#include "TRandom3.h"
#include "MitPhysics/Mods/interface/PhotonIDMod.h"
#include "MitAna/DataTree/interface/PhotonCol.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Utils/interface/IsolationTools.h"
#include "MitPhysics/Utils/interface/PhotonTools.h"

#include <TSystem.h>

using namespace mithep;

ClassImp(mithep::PhotonIDMod)
  
//--------------------------------------------------------------------------------------------------
PhotonIDMod::PhotonIDMod(const char *name, const char *title) : 
BaseMod(name,title),
  fPhotonBranchName  (Names::gkPhotonBrn),
  fGoodPhotonsName   (ModNames::gkGoodPhotonsName),
  fTrackBranchName   (Names::gkTrackBrn),
  fBeamspotBranchName(Names::gkBeamSpotBrn),
  fPileUpDenName     (Names::gkPileupEnergyDensityBrn),
  fConversionName    ("MergedConversions"),
  fElectronName      ("Electrons"),
  fGoodElectronName  ("Electrons"),
  fPVName            (Names::gkPVBeamSpotBrn),
// MC specific stuff...
  fMCParticleName    (Names::gkMCPartBrn),
  fPileUpName        (Names::gkPileupInfoBrn),
  fPFCandsName       ("PFCandidates"),
  
  fPhotonIDType      ("Custom"),
  fPhotonIsoType     ("Custom"),
  fPhotonPtMin       (15.0),
  fHadOverEmMax      (0.02),
  fApplySpikeRemoval (false),
  fApplyPixelSeed    (true),
  fApplyElectronVeto (false),
  fInvertElectronVeto(false),
  fApplyElectronVetoConvRecovery(false),
  fApplyConversionId (false),
  fApplyTriggerMatching(false),
  fPhotonR9Min       (0.5),
  fPhIdType          (kIdUndef),
  fPhIsoType         (kIsoUndef),
  fFiduciality       (true),
  fEtaWidthEB        (0.01),
  fEtaWidthEE        (0.028),
  fAbsEtaMax         (999.99),
  fApplyR9Min        (false),
  
  fEffAreaEcalEE     (0.089),
  fEffAreaHcalEE     (0.156),
  fEffAreaTrackEE    (0.261),
  
  fEffAreaEcalEB     (0.183),
  fEffAreaHcalEB     (0.062),
  fEffAreaTrackEB    (0.306),
  
  fPhotons(0),
  fTracks(0),
  fBeamspots(0),
  fPileUpDen(0),
  fConversions(0),
  fElectrons(0),
  fGoodElectrons(0),
  fPV(0),
  fMCParticles       (0),
  fPileUp            (0),
  fPFCands           (0),
  
// MVA ID Stuff
  fbdtCutBarrel      (0.0744), //cuts give the same effiiciency (relative to preselection) with cic
  fbdtCutEndcap      (0.0959), //cuts give the same effiiciency (relative to preselection) with cic  

// ------------------------------------------------------------------------------
// this stuff should go away ..... (fab)
  fVariableType       (10), //please use 4 which is the correct type
  fEndcapWeights      (gSystem->Getenv("CMSSW_BASE")+TString("/src/MitPhysics/data/TMVAClassificationPhotonID_Endcap_PassPreSel_Variable_10_BDTnCuts2000_BDT.weights.xml")),
  fBarrelWeights      (gSystem->Getenv("CMSSW_BASE")+TString("/src/MitPhysics/data/TMVAClassificationPhotonID_Barrel_PassPreSel_Variable_10_BDTnCuts2000_BDT.weights.xml")),
// ------------------------------------------------------------------------------  

  fIdMVATypeName     ("2011IdMVA"),
  fIdMVAType         (MVATools::k2011IdMVA),

//   fDoMCR9Scaling         (false),
//   fMCR9ScaleEB           (1.0),
//   fMCR9ScaleEE           (1.0),
//   fDoMCSigIEtaIEtaScaling(false),
//   fDoMCWidthScaling      (false),
  
  fDoMCErrScaling        (false),
  fMCErrScaleEB          (1.0),
  fMCErrScaleEE          (1.0),
  
  fPhotonsFromBranch(true),
  fPVFromBranch      (true),
  fGoodElectronsFromBranch (true),
  fIsData(false),

// ------------------------------------------------
// added by fab: smearing/scaling options...
  fDoDataEneCorr                 (false),
  fDoMCSmear                     (false),
  
  fDoShowerShapeScaling          (false),

// ---------------------------------------
  fDataEnCorr_EBlowEta_hR9       (0.),
  fDataEnCorr_EBlowEta_lR9       (0.),
  fDataEnCorr_EBhighEta_hR9      (0.),
  fDataEnCorr_EBhighEta_lR9      (0.),
  fDataEnCorr_EElowEta_hR9       (0.),
  fDataEnCorr_EElowEta_lR9       (0.),
  fDataEnCorr_EEhighEta_hR9      (0.),
  fDataEnCorr_EEhighEta_lR9      (0.),
  fRunStart                      (0),
  fRunEnd                        (0),

  fMCSmear_EBlowEta_hR9          (0.),
  fMCSmear_EBlowEta_lR9          (0.),
  fMCSmear_EBhighEta_hR9         (0.),
  fMCSmear_EBhighEta_lR9         (0.),
  fMCSmear_EElowEta_hR9          (0.),
  fMCSmear_EElowEta_lR9          (0.),
  fMCSmear_EEhighEta_hR9         (0.),
  fMCSmear_EEhighEta_lR9         (0.),
 
  fRng                           (new TRandom3()),
  fRhoAlgo                       (mithep::PileupEnergyDensity::kKt6PFJets)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void PhotonIDMod::Process()
{
  // Process entries of the tree. 

  LoadEventObject(fPhotonBranchName,   fPhotons);
  
  const BaseVertex *bsp = NULL;
  Double_t _tRho = -1.;
  Float_t rho2012 = -99;
  const TriggerObjectCol *trigObjs = 0;  
  if (fPhotons->GetEntries()>0) {
    LoadEventObject(fTrackBranchName,    fTracks);
    LoadEventObject(fBeamspotBranchName, fBeamspots);
    LoadEventObject(fPileUpDenName,      fPileUpDen);
    LoadEventObject(fConversionName,     fConversions);
    LoadEventObject(fElectronName,       fElectrons);
    LoadEventObject(fGoodElectronName,   fGoodElectrons);
    LoadEventObject(fPVName,             fPV);
    LoadEventObject(fPFCandsName,             fPFCands);

    
    if(fBeamspots->GetEntries() > 0)
      bsp = fBeamspots->At(0);

    if(fPileUpDen->GetEntries() > 0)
      _tRho = (Double_t) fPileUpDen->At(0)->RhoRandomLowEta();

    if (fPileUpDen->At(0)->RhoKt6PFJets()>0.) rho2012 = fPileUpDen->At(0)->RhoKt6PFJets();
    else rho2012 = fPileUpDen->At(0)->Rho();
    
    //get trigger object collection if trigger matching is enabled
    if (fApplyTriggerMatching) {
      trigObjs = GetHLTObjects(fTrigObjectsName);
    }    
  
  }

  Float_t theRho = _tRho;
  if (fRhoAlgo != mithep::PileupEnergyDensity::nAllAlgos) {
    if (fRhoAlgo == mithep::PileupEnergyDensity::kKt6PFJets)
      theRho = rho2012;
    else if(fPileUpDen->GetEntries() != 0)
      theRho = fPileUpDen->At(0)->Rho(fRhoAlgo);
  }

  // ------------------------------------------------------------
  // Get Event header for Run info etc.
  const EventHeader* evtHead   = this->GetEventHeader();
  unsigned int       evtNum    = evtHead->EvtNum();
  UInt_t             runNumber = evtHead->RunNum();
  Float_t            _runNum   = (Float_t) runNumber;
  Float_t            _lumiSec  = (Float_t) evtHead->LumiSec();

  PhotonOArr *GoodPhotons = new PhotonOArr;
  GoodPhotons->SetName(fGoodPhotonsName);
  GoodPhotons->SetOwner(true);

  for (UInt_t i=0; i<fPhotons->GetEntries(); ++i) {    
    // need to cpoy the photon in order to be able to scale R9 etc.
    Photon *ph = new Photon(*fPhotons->At(i));        
    
    if (fFiduciality == true &&
        (ph->SCluster()->AbsEta()>=2.5 || (ph->SCluster()->AbsEta()>=1.4442 && ph->SCluster()->AbsEta()<=1.566) ) ) 
      continue;
    
    if (fInvertElectronVeto && PhotonTools::PassElectronVeto(ph,fGoodElectrons) && false) {
      continue;
    }
    
    // -----------------------------------------------------------------------------------
    // Do all the scaling ....
    if (fDoMCErrScaling && !fIsData) {
      if (ph->SCluster()->AbsEta()<1.5) PhotonTools::ScalePhotonError(ph,fMCErrScaleEB);
      else PhotonTools::ScalePhotonError(ph,fMCErrScaleEE);
    }
    
    if (fDoShowerShapeScaling && !fIsData) {
      PhotonTools::ScalePhotonShowerShapes(ph,fSSType);
    }

    // fab: replaced by shower-shape scaling above...
    //     if (fDoMCR9Scaling && !fIsData) {
    //       if (ph->SCluster()->AbsEta()<1.5) PhotonTools::ScalePhotonR9(ph,fMCR9ScaleEB);
    //       else PhotonTools::ScalePhotonR9(ph,fMCR9ScaleEE);
    //     }

    //     if (fDoMCSigIEtaIEtaScaling && !fIsData) {
    //       if (ph->SCluster()->AbsEta()<1.5) ph->SetCoviEtaiEta(0.87*ph->CoviEtaiEta() + 0.0011);
    //       else ph->SetCoviEtaiEta(0.99*ph->CoviEtaiEta());
    //     }

    //     if (fDoMCWidthScaling && !fIsData) {
    //       ph->SetEtaWidth(0.99*ph->EtaWidth());
    //       ph->SetPhiWidth(0.99*ph->PhiWidth());
    //     }
    
    PhotonTools::eScaleCats escalecat;

    escalecat = PhotonTools::EScaleCat(ph);

    // now we dicide if we either scale (Data) or Smear (MC) the Photons
    if (fIsData) {
      if (fDoDataEneCorr) {
        // starting with scale = 1.
        double scaleFac = 1.;
        // checking the run Rangees ...
        Int_t runRange = FindRunRangeIdx(runNumber);
        if(runRange > -1) {
          scaleFac *= GetDataEnCorr(runRange, escalecat);
        }
        PhotonTools::ScalePhoton(ph, scaleFac);
      }
    }
    
    if (fDoMCSmear) {

      double width = 0;
      
      width = GetMCSmearFac(escalecat);
      
      if (!fIsData) {
        // get the seed to do deterministic smearing...
        UInt_t seedBase = (UInt_t) evtNum + (UInt_t) _runNum + (UInt_t) _lumiSec;
        UInt_t seed     = seedBase + (UInt_t) ph->E() +
          (UInt_t) (TMath::Abs(10.*ph->SCluster()->Eta()));
        // get the smearing for MC photons..
        PhotonTools::SmearPhoton(ph, fRng, width, seed);
      }
      PhotonTools::SmearPhotonError(ph, width);
    }
    
    // ---------------------------------------------------------------------
    // set the photonIdMVA value of requested...
    double idMvaVal = fTool.GetMVAbdtValue(ph,fPV->At(0),fTracks, fPV, theRho, fPFCands, fElectrons, fApplyElectronVeto);

    ph->SetIdMva(idMvaVal);

    // ---------------------------------------------------------------------
    // check if we use the Vgamma2011 Selection. If yes, bypass all the below...
    if( fPhIdType == kTrivialSelection ) {
      if( ph->Pt() > fPhotonPtMin && ph->SCluster()->AbsEta() <= fAbsEtaMax )
        GoodPhotons->AddOwned(ph);
      continue; // go to next Photons
    }
    
    // ---------------------------------------------------------------------
    // check if we use the Vgamma2011 Selection. If yes, bypass all the below...
    if( fPhIdType == kVgamma2011Selection ) {
      if( ph->Pt() > fPhotonPtMin && ph->SCluster()->AbsEta() <= fAbsEtaMax && PhotonTools::PassVgamma2011Selection(ph, _tRho) )
        GoodPhotons->AddOwned(ph);
      continue; // go to next Photons
    }

    // ---------------------------------------------------------------------
    // check if we use the CiC Selection. If yes, bypass all the below...
    if(fPhIdType == kBaseLineCiC) {
      if( PhotonTools::PassCiCSelection(ph, fPV->At(0), fTracks, fElectrons, fPV, _tRho, fPhotonPtMin, fApplyElectronVeto) )
        GoodPhotons->AddOwned(ph);
      continue; // go to next Photons
    }
    
    if(fPhIdType == kBaseLineCiCPF) {
      if( PhotonTools::PassSinglePhotonPreselPFISO(ph,fElectrons,fConversions,bsp,fTracks,fPV->At(0),rho2012,fPFCands,fApplyElectronVeto,fInvertElectronVeto) && PhotonTools::PassCiCPFIsoSelection(ph, fPV->At(0), fPFCands, fPV, rho2012, fPhotonPtMin, false,0.,1.,0.,1.) )
        GoodPhotons->AddOwned(ph);
      continue; // go to next Photons
    }    

    //add for mono photon: cic id without presel
    if(fPhIdType == kBaseLineCiCPFNoPresel) {
      if(PhotonTools::PassCiCPFIsoSelectionWithEleVeto(ph, fElectrons,fConversions,bsp,fPV->At(0), fPFCands, fPV, rho2012, fPhotonPtMin,fApplyElectronVeto,fInvertElectronVeto) )
        GoodPhotons->AddOwned(ph);
      continue; // go to next Photons
    }    
    // ---------------------------------------------------------------------

    //add for mono photon: egamma medium working point
    if(fPhIdType == kEgammaMedium) {
      if(PhotonTools::PassEgammaMediumSelectionWithEleVeto(ph, fElectrons,fConversions,bsp,fPV->At(0), fPFCands, fPV, rho2012, fPhotonPtMin,fApplyElectronVeto,fInvertElectronVeto) )
        GoodPhotons->AddOwned(ph);
      continue; // go to next Photons
    } 

    //loose photon preselection for subsequent mva
    if(fPhIdType == kMITPhSelection ) {
      if( ph->Pt()>fPhotonPtMin && PhotonTools::PassSinglePhotonPresel(ph,fElectrons,fConversions,bsp,fTracks,fPV->At(0),_tRho,fApplyElectronVeto,fInvertElectronVeto) ) {
        GoodPhotons->AddOwned(ph);
      }
      continue;
    }

    //loose photon preselection for subsequent mva
    if(fPhIdType == kMITPFPhSelection ) {
      if( ph->Pt()>fPhotonPtMin && PhotonTools::PassSinglePhotonPreselPFISO(ph,fElectrons,fConversions,bsp,fTracks,fPV->At(0),rho2012,fPFCands,fApplyElectronVeto,fInvertElectronVeto) ) {
        GoodPhotons->AddOwned(ph);
      }
      continue;
    }    
    
    //loose photon preselection for subsequent mva
    if(fPhIdType == kMITPFPhSelectionNoEcal) {
      if( ph->Pt()>fPhotonPtMin && PhotonTools::PassSinglePhotonPreselPFISONoEcal(ph,fElectrons,fConversions,bsp,fTracks,fPV->At(0),rho2012,fPFCands,fApplyElectronVeto,fInvertElectronVeto) ) {
        GoodPhotons->AddOwned(ph);
      }
      continue;
    }   
    
    //loose photon preselection for subsequent mva
    if(fPhIdType == kMITPFPhSelection_NoTrigger ) {
      if( ph->Pt()>fPhotonPtMin && PhotonTools::PassSinglePhotonPreselPFISO_NoTrigger(ph,fElectrons,fConversions,bsp,fTracks,fPV->At(0),rho2012,fPFCands,fApplyElectronVeto,fInvertElectronVeto) ) {
        GoodPhotons->AddOwned(ph);
      }
      continue;
    }    

    
    Bool_t isbarrel = ph->SCluster()->AbsEta()<1.5;
    
    // add MingMings MVA ID on single Photon level
    if(fPhIdType == kMITMVAId ) {
      // we compute the bdt val already before, so use it ...
      //if( ph->Pt()>fPhotonPtMin && PhotonTools::PassSinglePhotonPresel(ph,fElectrons,fConversions,bsp,fTracks,fPV->At(0),_tRho,fApplyElectronVeto) && fTool.PassMVASelection(ph, fPV->At(0) ,fTracks, fPV, _t_TRho,fbdtCutBarrel,fbdtCutEndcap, fElectrons, fApplyElectronVeto) ) {
        
      if( ph->Pt()>fPhotonPtMin &&
          PhotonTools::PassSinglePhotonPresel(
                                              ph, fElectrons, fConversions, bsp, fTracks, fPV->At(0), _tRho,
                                              fApplyElectronVeto
                                              ) &&
          ( ( isbarrel && ph->IdMva() > fbdtCutBarrel ) ||
            (             ph->IdMva() > fbdtCutEndcap ) )
          )
        GoodPhotons->AddOwned(ph);
      
      continue;
    }  // go to next Photon
    
    if (ph->Pt() <= fPhotonPtMin) 
      continue;    // add good electron
    
    Bool_t passSpikeRemovalFilter = true;
    
    if (ph->SCluster() && ph->SCluster()->Seed()) {
      if(ph->SCluster()->Seed()->Energy() > 5.0 && 
         ph->SCluster()->Seed()->EMax() / ph->SCluster()->Seed()->E3x3() > 0.95
         ) {
        passSpikeRemovalFilter = false;
      }
    }
    
    // For Now Only use the EMax/E3x3 prescription.
    //if(ph->SCluster()->Seed()->Energy() > 5.0 && 
    //   (1 - (ph->SCluster()->Seed()->E1x3() + ph->SCluster()->Seed()->E3x1() - 2*ph->SCluster()->Seed()->EMax())) > 0.95
    //  ) {
    //  passSpikeRemovalFilter = false;
    //}

    if (fApplySpikeRemoval && !passSpikeRemovalFilter) continue;
    
    if (ph->HadOverEm() >= fHadOverEmMax) 
      continue;

    if (fApplyPixelSeed == true &&
        ph->HasPixelSeed() == true) 
      continue;

    if (fApplyElectronVeto && !PhotonTools::PassElectronVeto(ph,fElectrons) ) continue;

    if (fApplyElectronVetoConvRecovery && !PhotonTools::PassElectronVetoConvRecovery(ph,fElectrons,fConversions,bsp) ) continue;

    if (fApplyConversionId && !PhotonTools::PassConversionId(ph,PhotonTools::MatchedConversion(ph,fConversions,bsp))) continue;

    if (fApplyTriggerMatching && !PhotonTools::PassTriggerMatching(ph,trigObjs)) continue;

    Bool_t idcut = false;
    switch (fPhIdType) {
    case kTight:
      idcut = ph->IsTightPhoton();
      break;
    case kLoose:
      idcut = ph->IsLoosePhoton();
      break;
    case kCustomId:
      idcut = true;
    default:
      break;
    }

    if (!idcut) 
      continue;

    Bool_t isocut = false;
    switch (fPhIsoType) {      
    case kNoIso:
      isocut = true;
      break;
    case kCombinedIso:
      {
        Double_t totalIso = ph->HollowConeTrkIsoDr04()+
          ph->EcalRecHitIsoDr04() +
          ph->HcalTowerSumEtDr04();
        if (totalIso/ph->Pt() < 0.25)
          isocut = true;
      }
      break;
    case kCustomIso:
      {
        if ( ph->HollowConeTrkIsoDr04() < (1.5 + 0.001*ph->Pt()) && ph->EcalRecHitIsoDr04()<(2.0+0.006*ph->Pt()) && ph->HcalTowerSumEtDr04()<(2.0+0.0025*ph->Pt()) )
          isocut = true;
      }
      break;
        
    case kMITPUCorrected:
      {
        // compute the PU corrections only if Rho is available
        // ... otherwise (_tRho = 0.0) it's the std isolation
        isocut = true;
        Double_t fEffAreaEcal = fEffAreaEcalEB;
        Double_t fEffAreaHcal = fEffAreaHcalEB;
        Double_t fEffAreaTrack = fEffAreaTrackEB;

        if( !isbarrel ) {
          fEffAreaEcal = fEffAreaEcalEE;
          fEffAreaHcal = fEffAreaHcalEE;
          fEffAreaTrack = fEffAreaTrackEE;
        }         

        Double_t EcalCorrISO =   ph->EcalRecHitIsoDr04();
        if(_tRho > -0.5 ) EcalCorrISO -= _tRho * fEffAreaEcal;
        if ( EcalCorrISO > (2.0+0.006*ph->Pt()) ) isocut = false; 
        if ( isocut || true ) {
          Double_t HcalCorrISO = ph->HcalTowerSumEtDr04(); 
          if(_tRho > -0.5 ) HcalCorrISO -= _tRho * fEffAreaHcal;
          if ( HcalCorrISO > (2.0+0.0025*ph->Pt()) ) isocut = false;
        }
        if ( isocut || true ) {
          Double_t TrackCorrISO = IsolationTools::TrackIsolationNoPV(ph, bsp, 0.4, 0.04, 0.0, 0.015, 0.1, TrackQuality::highPurity, fTracks);
          if(_tRho > -0.5 )
            TrackCorrISO -= _tRho * fEffAreaTrack;
          if ( TrackCorrISO > (1.5 + 0.001*ph->Pt()) ) isocut = false;
        }
        break;
      }
    default:
      break;
    }
    
    if (!isocut) 
      continue;
    
    if ( fApplyR9Min && ph->R9() <= fPhotonR9Min) 
      continue;

    if ((isbarrel && ph->CoviEtaiEta() >= fEtaWidthEB) ||
        (!isbarrel && ph->CoviEtaiEta() >= fEtaWidthEE))
      continue;

    if (ph->AbsEta() >= fAbsEtaMax) 
      continue;
    
    // add good electron
    GoodPhotons->AddOwned(ph);
  }

  // sort according to pt
  GoodPhotons->Sort();

  // add to event for other modules to use
  AddObjThisEvt(GoodPhotons);  
  
  return; 
  
}

//--------------------------------------------------------------------------------------------------
void PhotonIDMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we just request the photon collection branch.

  ReqEventObject(fPhotonBranchName,   fPhotons,     fPhotonsFromBranch);
  ReqEventObject(fTrackBranchName,    fTracks,      true);
  ReqEventObject(fBeamspotBranchName, fBeamspots,   true);
  ReqEventObject(fConversionName,     fConversions,  true);
  ReqEventObject(fElectronName,       fElectrons,   true);
  ReqEventObject(fGoodElectronName,       fGoodElectrons,   fGoodElectronsFromBranch);  
  ReqEventObject(fPVName, fPV, fPVFromBranch);
  ReqEventObject(fPileUpDenName,      fPileUpDen, true);
  ReqEventObject(fPFCandsName,      fPFCands, true);

  if (!fIsData) {
    ReqBranch(fPileUpName,            fPileUp);
    ReqBranch(fMCParticleName,        fMCParticles);
  }   

  if(fIdMVATypeName.CompareTo("2011IdMVA") == 0)
    fIdMVAType =       MVATools::k2011IdMVA;
  else if (fIdMVATypeName.CompareTo("2012IdMVA_globe") == 0)
    fIdMVAType =       MVATools::k2012IdMVA_globe;
  else if (fIdMVATypeName.CompareTo("2012IdMVA") == 0)
    fIdMVAType =       MVATools::k2012IdMVA;
  else if (fIdMVATypeName.CompareTo("2013FinalIdMVA_8TeV") == 0)
    fIdMVAType =       MVATools::k2013FinalIdMVA_8TeV;
  else if (fIdMVATypeName.CompareTo("2013FinalIdMVA_7TeV") == 0)
    fIdMVAType =       MVATools::k2013FinalIdMVA_7TeV;
  else if (fIdMVATypeName.CompareTo("2011IdMVA_HZg") == 0)
    fIdMVAType =       MVATools::k2011IdMVA_HZg;
  else if (fIdMVATypeName.CompareTo("None") == 0)
    fIdMVAType =       MVATools::kNone;
  else {
    std::cerr<<" Id MVA "<<fIdMVATypeName<<" not implemented."<<std::endl;
    return;
  }

  fTool.InitializeMVA(fIdMVAType);
  
  if (fPhotonIDType.CompareTo("Tight") == 0) 
    fPhIdType = kTight;
  else if (fPhotonIDType.CompareTo("Loose") == 0) 
    fPhIdType = kLoose;
  else if (fPhotonIDType.CompareTo("LooseEM") == 0) {
    throw std::runtime_error("Invalid congiruation LooseEM for PhotonIDMod::PhotonIDType");
  }
  else if (fPhotonIDType.CompareTo("Custom") == 0)
    fPhIdType = kCustomId;
  else if (fPhotonIDType.CompareTo("BaseLineCiC") == 0) {
    fPhIdType = kBaseLineCiC;
    fPhotonIsoType = "NoIso";
  }
  else if (fPhotonIDType.CompareTo("BaseLineCiCPF") == 0) {
    fPhIdType = kBaseLineCiCPF;
    fPhotonIsoType = "NoIso";
  } 
  else if (fPhotonIDType.CompareTo("BaseLineCiCPFNoPresel") == 0) {
    fPhIdType = kBaseLineCiCPFNoPresel;
    fPhotonIsoType = "NoIso";
  }  
  else if (fPhotonIDType.CompareTo("EgammaMedium") == 0) {
    fPhIdType = kEgammaMedium;
    fPhotonIsoType = "NoIso";
  } 
  else if (fPhotonIDType.CompareTo("MITMVAId") == 0) {
    fPhIdType = kMITMVAId;
    fPhotonIsoType = "NoIso";
    // this is now a 'generic' MVAId: set the MVAType using 'SetIdMVAType(const char*)'  (fab)
    //fTool.InitializeMVA(fVariableType,fEndcapWeights,fBarrelWeights);   // moved down for after we know which MVA type to load...
  }
  else if (fPhotonIDType.CompareTo("MITSelection") == 0)  {
    fPhIdType = kMITPhSelection;
    fPhotonIsoType = "NoIso";    
  }
  else if (fPhotonIDType.CompareTo("MITPFSelection") == 0)  {
    fPhIdType = kMITPFPhSelection;
    fPhotonIsoType = "NoIso";    
  }
  else if (fPhotonIDType.CompareTo("MITPFSelectionNoEcal") == 0){
    fPhIdType =       kMITPFPhSelectionNoEcal; 
    fPhotonIsoType = "NoIso";
  }
  else if (fPhotonIDType.CompareTo("MITPFSelection_NoTrigger") == 0)  {
    fPhIdType = kMITPFPhSelection_NoTrigger;
    fPhotonIsoType = "NoIso";
  }
  else if (fPhotonIDType.CompareTo("Vgamma2011Selection") == 0)  {
    fPhIdType = kVgamma2011Selection;
    fPhotonIsoType = "NoIso";    
  }
  else if (fPhotonIDType.CompareTo("TrivialSelection") == 0)  {
    fPhIdType = kTrivialSelection;
    fPhotonIsoType = "NoIso";    
  }
  else {
    SendError(kAbortAnalysis, "SlaveBegin",
              "The specified photon identification %s is not defined.",
              fPhotonIDType.Data());
    return;
  }

  if (fPhotonIsoType.CompareTo("NoIso") == 0 )
    fPhIsoType = kNoIso;
  else if (fPhotonIsoType.CompareTo("CombinedIso") == 0 )
    fPhIsoType = kCombinedIso;
  else if (fPhotonIsoType.CompareTo("Custom") == 0 )
    fPhIsoType = kCustomIso;
  else if (fPhotonIsoType.CompareTo("MITPUCorrected") == 0 )
    fPhIsoType = kMITPUCorrected;
  else {
    SendError(kAbortAnalysis, "SlaveBegin",
              "The specified photon isolation %s is not defined.",
              fPhotonIsoType.Data());
    return;
  }
  
  if      (fShowerShapeType.CompareTo("None")            == 0)
    fSSType =       PhotonTools::kNoShowerShapeScaling;
  else if (fShowerShapeType.CompareTo("2011ShowerShape") == 0)
    fSSType =       PhotonTools::k2011ShowerShape;
  else if (fShowerShapeType.CompareTo("2012ShowerShape") == 0)
    fSSType =       PhotonTools::k2012ShowerShape;
  else {
    std::cerr<<"shower shape scale "<<fShowerShapeType<<" not implemented."<<std::endl;
    return;
  }
  
}


//---------------------------------------------------------------------------------------------------
Int_t PhotonIDMod::FindRunRangeIdx(UInt_t run)
{
  // this routine looks for the idx of the run-range
  Int_t runRange=-1;
  for (UInt_t iRun = 0; iRun<fRunStart.size(); ++iRun) {
    if (run >= fRunStart[iRun] && run <= fRunEnd[iRun]) {
      runRange = (Int_t) iRun;
      return runRange;
    }
  }
  return runRange;
}

//---------------------------------------------------------------------------------------------------
Double_t PhotonIDMod::GetDataEnCorr(Int_t runRange, PhotonTools::eScaleCats cat)
{
  switch (cat) {
  case PhotonTools::kEBhighEtaGold:
    return fDataEnCorr_EBhighEta_hR9[runRange];
  case PhotonTools::kEBhighEtaBad:
    return fDataEnCorr_EBhighEta_lR9[runRange];
  case PhotonTools::kEBlowEtaGold:
    return fDataEnCorr_EBlowEta_hR9[runRange];
  case PhotonTools::kEBlowEtaBad:
    return fDataEnCorr_EBlowEta_lR9[runRange];
  case PhotonTools::kEEhighEtaGold:
    return fDataEnCorr_EEhighEta_hR9[runRange];
  case PhotonTools::kEEhighEtaBad:
    return fDataEnCorr_EEhighEta_lR9[runRange];
  case PhotonTools::kEElowEtaGold:
    return fDataEnCorr_EElowEta_hR9[runRange];
  case PhotonTools::kEElowEtaBad:
    return fDataEnCorr_EElowEta_lR9[runRange];
  default:
    return 1.;
  }
}

//---------------------------------------------------------------------------------------------------
Double_t PhotonIDMod::GetMCSmearFac(PhotonTools::eScaleCats cat)
{
  switch (cat) {
  case PhotonTools::kEBhighEtaGold:
    return fMCSmear_EBhighEta_hR9;
  case PhotonTools::kEBhighEtaBad:
    return fMCSmear_EBhighEta_lR9;
  case PhotonTools::kEBlowEtaGold:
    return fMCSmear_EBlowEta_hR9;
  case PhotonTools::kEBlowEtaBad:
    return fMCSmear_EBlowEta_lR9;
  case PhotonTools::kEEhighEtaGold:
    return fMCSmear_EEhighEta_hR9;
  case PhotonTools::kEEhighEtaBad:
    return fMCSmear_EEhighEta_lR9;
  case PhotonTools::kEElowEtaGold:
    return fMCSmear_EElowEta_hR9;
  case PhotonTools::kEElowEtaBad:
    return fMCSmear_EElowEta_lR9;
  default:
    return 1.;
  }
}

