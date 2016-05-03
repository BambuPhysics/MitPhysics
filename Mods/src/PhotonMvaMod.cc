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
#include "MitPhysics/Init/interface/Constants.h"
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
  fGoodPhotonsName        (ModNames::gkGoodPhotonsName),
  fPhotonPtMin            (20.0),
  fPhotonEtaMax           (2.5),
  fIsData                 (false),
  fApplyShowerRescaling   (false),
  fPhotonsFromBranch      (true),
  fDoRegression           (kTRUE),
  fPhFixString            ("4_2"),
  fRegWeights             (gSystem->Getenv("MIT_DATA") + TString("/gbrv2ph_52x.root")),
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
  auto* photons = GetObject<PhotonCol>(fPhotonBranchName);

  // -----------------------------------------------------------
  // Output Photon Collection. Will ALWAYS contain 0 or 2 Photons
  PhotonOArr *GoodPhotons = new PhotonOArr;
  GoodPhotons->SetName(fGoodPhotonsName);
  GoodPhotons->SetOwner(kTRUE);

  // add to event for other modules to us
  AddObjThisEvt(GoodPhotons);

  if (photons->GetEntries() < fMinNumPhotons)
    return;

  // ------------------------------------------------------------
  // store preselected Photons (and which CiCCategory they are)
  PhotonOArr preselPh;

  // 1. pre-selection; but keep non-passing photons in second container...
  for (UInt_t i=0; i<photons->GetEntries(); ++i) {
    const Photon *ph = photons->At(i);
    if (fDoPreselection) {
      if (ph->SCluster()->AbsEta() >= fPhotonEtaMax ||
          (ph->SCluster()->AbsEta()>=1.4442 && ph->SCluster()->AbsEta()<=1.566))
        continue;
      if (ph->Et()                <  fPhotonPtMin)
        continue;
      if (ph->HadOverEm()         >  0.15)
        continue;
      if (ph->SCluster()->AbsEta() < mithep::gkPhoEBEtaMax) {
        if (ph->CoviEtaiEta() > 0.015)
          continue;
      }
      else {
        if (ph->CoviEtaiEta() > 0.035)
          continue;
      }
    }
    preselPh.Add(ph);
  }

  if (preselPh.GetEntries() < fMinNumPhotons)
    return;

  // second loop to sort & assign the right Categories..
  for (unsigned int iPh = 0; iPh <preselPh.GetEntries(); ++iPh) {
    const Photon *ph = preselPh.At(iPh);
    Photon       *outph = new Photon(*ph);
    if (fDoRegression) {
      
      //if (!egcor.IsInitialized())
      //  egcor.Initialize(fPhFixString,fPhFixFile,fRegWeights,fRegressionVersion);
      //if (fRegressionVersion>0)
      //  egcor.CorrectEnergyWithError(outph,pv,fPileUpDen->At(0)->RhoKt6PFJets(),
      //				     fRegressionVersion, fApplyShowerRescaling && !fIsData);
    }
    GoodPhotons->AddOwned(outph);
  }
  // sort according to pt
  GoodPhotons->Sort();
}

//--------------------------------------------------------------------------------------------------
void PhotonMvaMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here, we just request the
  // photon collection branch.

  if (fIsData)
    fPhFixFile = gSystem->Getenv("MIT_DATA") + TString("/PhotonFixGRPV22.dat");
  else
    fPhFixFile = gSystem->Getenv("MIT_DATA") + TString("/PhotonFixSTART42V13.dat");
}
