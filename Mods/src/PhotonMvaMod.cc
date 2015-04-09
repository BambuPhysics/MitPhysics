#include <TSystem.h>
#include <TRandom3.h>
#include <TDataMember.h>
#include <TH1D.h>
#include <TNtuple.h>
#include "MitAna/DataTree/interface/PhotonCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/StableData.h"
#include "MitAna/DataTree/interface/StableParticle.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Utils/interface/IsolationTools.h"
#include "MitPhysics/Utils/interface/PhotonTools.h"
#include "MitPhysics/Utils/interface/VertexTools.h"
#include "MitPhysics/Utils/interface/PhotonFix.h"
#include "MitPhysics/Utils/interface/MVATools.h"
#include "MitPhysics/Mods/interface/PhotonMvaMod.h"

using namespace mithep;

ClassImp(mithep::PhotonMvaMod)

//--------------------------------------------------------------------------------------------------
PhotonMvaMod::PhotonMvaMod(const char *name, const char *title) :
  BaseMod                 (name,title),
  fPhotonBranchName       (Names::gkPhotonBrn),
  fElectronName           (Names::gkElectronBrn),
  fGoodElectronName       (Names::gkElectronBrn),
  fConversionName         (Names::gkMvfConversionBrn),
  fTrackBranchName        (Names::gkTrackBrn),
  fPileUpDenName          (Names::gkPileupEnergyDensityBrn),
  fPVName                 (Names::gkPVBeamSpotBrn),
  fBeamspotName           (Names::gkBeamSpotBrn),
  fPFCandName             (Names::gkPFCandidatesBrn),
  fMCParticleName         (Names::gkMCPartBrn),
  fPileUpName             (Names::gkPileupInfoBrn),
  fGoodPhotonsName        (ModNames::gkGoodPhotonsName),
  fPhotonPtMin            (20.0),
  fPhotonEtaMax           (2.5),
  fIsData                 (false),
  fApplyShowerRescaling   (false),
  fPhotonsFromBranch      (true),
  fPVFromBranch           (true),
  fGoodElectronsFromBranch(kTRUE),
  fPhotons                (0),
  fElectrons              (0),
  fConversions            (0),
  fTracks                 (0),
  fPileUpDen              (0),
  fPV                     (0),
  fBeamspot               (0),
  fPFCands                (0),
  fMCParticles            (0),
  fPileUp                 (0),
  fDoRegression           (kTRUE),
  fPhFixString            ("4_2"),
  fRegWeights             (gSystem->Getenv("MIT_DATA") + TString("gbrv2ph_52x.root")),
  fApplyEleVeto           (true),
  fRegressionVersion      (2),
  fMinNumPhotons          (2),
  fDoPreselection         (kTRUE)
{
  // Constructor.
}

PhotonMvaMod::~PhotonMvaMod()
{
  // Destructor.
}

//--------------------------------------------------------------------------------------------------
void PhotonMvaMod::Process()
{
  // ------------------------------------------------------------
  // Process entries of the tree.
  LoadEventObject(fPhotonBranchName,   fPhotons);

  // -----------------------------------------------------------
  // Output Photon Collection. Will ALWAYS contain 0 or 2 Photons
  PhotonOArr *GoodPhotons = new PhotonOArr;
  GoodPhotons->SetName(fGoodPhotonsName);
  GoodPhotons->SetOwner(kTRUE);
  // add to event for other modules to us
  AddObjThisEvt(GoodPhotons);

  if (fPhotons->GetEntries() < fMinNumPhotons)
    return;

  LoadEventObject(fPVName,       fPV);
  LoadEventObject(fPileUpDenName,fPileUpDen);

  // ------------------------------------------------------------
  // store preselected Photons (and which CiCCategory they are)
  PhotonOArr* preselPh  = new PhotonOArr;
  // 1. pre-selection; but keep non-passing photons in second container...
  for (UInt_t i=0; i<fPhotons->GetEntries(); ++i) {
    const Photon *ph = fPhotons->At(i);
    if (fDoPreselection) {
      if (ph->SCluster()->AbsEta() >= fPhotonEtaMax ||
          (ph->SCluster()->AbsEta()>=1.4442 && ph->SCluster()->AbsEta()<=1.566))
        continue;
      if (ph->Et()                <  fPhotonPtMin)
        continue;
      if (ph->HadOverEm()         >  0.15)
        continue;
      if(ph->IsEB()) {
        if (ph->CoviEtaiEta() > 0.015)
          continue;
      }
      else {
        if (ph->CoviEtaiEta() > 0.035)
          continue;
      }
    }
    preselPh->Add(ph);
  }
  assert(preselPh);

  if (preselPh->GetEntries() < fMinNumPhotons) {
    delete preselPh;
    return;
  }

  // second loop to sort & assign the right Categories..
  for (unsigned int iPh = 0; iPh <preselPh->GetEntries(); ++iPh) {
    const Photon *ph = preselPh->At(iPh);
    Photon       *outph = new Photon(*ph);
    if (fDoRegression) {

    //if (!egcor.IsInitialized())
    //  egcor.Initialize(fPhFixString,fPhFixFile,fRegWeights,fRegressionVersion);
    //if (fRegressionVersion>0)
    //  egcor.CorrectEnergyWithError(outph,fPV,fPileUpDen->At(0)->RhoKt6PFJets(),
    //				     fRegressionVersion, fApplyShowerRescaling && !fIsData);
          ThreeVectorC scpos = outph->SCluster()->Point();
      outph->SetCaloPosXYZ(scpos.X(),scpos.Y(),scpos.Z());
    }
    GoodPhotons->AddOwned(outph);
  }
  // sort according to pt
  GoodPhotons->Sort();

  // delete auxiliary photon collection...
  delete preselPh;

  return;
}

//--------------------------------------------------------------------------------------------------
void PhotonMvaMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here, we just request the
  // photon collection branch.

  ReqEventObject(fPhotonBranchName,   fPhotons,      fPhotonsFromBranch);
  ReqEventObject(fTrackBranchName,    fTracks,       true);
  ReqEventObject(fElectronName,       fElectrons,    true);
  ReqEventObject(fGoodElectronName,   fGoodElectrons,fGoodElectronsFromBranch);
  ReqEventObject(fPileUpDenName,      fPileUpDen,    true);
  ReqEventObject(fPVName,             fPV,           fPVFromBranch);
  ReqEventObject(fConversionName,     fConversions,  true);
  ReqEventObject(fBeamspotName,       fBeamspot,     true);
  ReqEventObject(fPFCandName,         fPFCands,      true);
  if (!fIsData) {
    ReqBranch(fPileUpName,            fPileUp);
    ReqBranch(fMCParticleName,        fMCParticles);
  }

  if (fIsData)
    fPhFixFile = gSystem->Getenv("MIT_DATA") + TString("/PhotonFixGRPV22.dat");
  else
    fPhFixFile = gSystem->Getenv("MIT_DATA") + TString("/PhotonFixSTART42V13.dat");
}
