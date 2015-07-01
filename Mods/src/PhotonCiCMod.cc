#include <TNtuple.h>
#include <TRandom3.h>
#include "MitAna/DataTree/interface/PhotonCol.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Utils/interface/IsolationTools.h"
#include "MitPhysics/Utils/interface/PhotonTools.h"
#include "MitPhysics/Mods/interface/PhotonCiCMod.h"

using namespace mithep;

ClassImp(mithep::PhotonCiCMod)

//--------------------------------------------------------------------------------------------------
PhotonCiCMod::PhotonCiCMod(const char *name, const char *title) :
  BaseMod           (name,title),
  fPhotonBranchName (Names::gkPhotonBrn),
  fGoodPhotonsName  (ModNames::gkGoodPhotonsName),
  fTrackBranchName  (Names::gkTrackBrn),
  fPileUpDenName    (Names::gkPileupEnergyDensityBrn),
  fElectronName     ("Electrons"),
  fPhotonPtMin      (20.0),
  fApplySpikeRemoval(kFALSE),
  fAbsEtaMax        (999.99),
  fPVName           (Names::gkPVBeamSpotBrn),
  fPVFromBranch     (kTRUE),
  fConversionName   (Names::gkMvfConversionBrn),
  fDataEnCorr_EB_hR9(0),
  fDataEnCorr_EB_lR9(0),
  fDataEnCorr_EE_hR9(0),
  fDataEnCorr_EE_lR9(0),
  fRunStart         (0),
  fRunEnd           (0),
  fMCSmear_EB_hR9   (0.0089),
  fMCSmear_EB_lR9   (0.0199),
  fMCSmear_EE_hR9   (0.0409),
  fMCSmear_EE_lR9   (0.0246),
  fIsData           (false),
  fRnd3             (new TRandom3()),
  fMCParticleName   (Names::gkMCPartBrn),
  fPileUpName       ("PileupInfo")
{
  // Constructor.
}

PhotonCiCMod::~PhotonCiCMod(){
  if (fRnd3)
    delete fRnd3;
}

//--------------------------------------------------------------------------------------------------
void PhotonCiCMod::Process()
{
  // Process entries of the tree.

  auto* photons = GetObject<PhotonCol>(fPhotonBranchName);

  PhotonOArr *GoodPhotons = new PhotonOArr;
  GoodPhotons->SetName(fGoodPhotonsName);
  GoodPhotons->SetOwner(kTRUE);

  Double_t _tRho = 0.;
  auto* tracks = GetObject<TrackCol>(fTrackBranchName);
  auto* pileupDen = GetObject<PileupEnergyDensityCol>(fPileUpDenName);
  auto* electrons = GetObject<ElectronCol>(fElectronName);
  auto* pv = GetObject<VertexCol>(fPVName);
  auto* conversions = GetObject<DecayParticleCol>(fConversionName);
  auto* beamspot = GetObject<BeamSpotCol>("BeamSpot");

  MCParticleCol const* mcParticles = 0;
  PileupInfoCol const* pileup = 0;
  if ( !fIsData ) {
    mcParticles = GetObject<MCParticleCol>(fMCParticleName);
    pileup = GetObject<PileupInfoCol>(fPileUpName);
  }

  Float_t numPU = -1.;
  if ( !fIsData )
    numPU = (Float_t) pileup->At(0)->GetPU_NumInteractions();

  if (pileupDen->GetEntries() > 0)
    _tRho = (Double_t) pileupDen->At(0)->RhoRandomLowEta();

  bool doVtxSelection = true;
  bool doMCSmear      = true;

  const EventHeader* evtHead = this->GetEventHeader();

  unsigned int evtNum = evtHead->EvtNum();
  Float_t _evtNum1 = (Float_t) ( (int) (evtNum/10000.) );
  Float_t _evtNum2 = (Float_t) ( (int) (evtNum % 10000)  );

  //double evtNumTest = (int) ( ( (double) _evtNum1 )*10000. + (double) _evtNum2 );

  UInt_t runNumber = evtHead->RunNum();
  Float_t _runNum  = (Float_t) runNumber;
  Float_t _lumiSec = (Float_t) evtHead->LumiSec();

  //unsigned int numVertices = pv->GetEntries();

  const BaseVertex *bsp = dynamic_cast<const BaseVertex*>(beamspot->At(0));

  PhotonOArr* preselPh  = new PhotonOArr;

  // 1. we do the pre-selection; but keep the non-passing photons in a second container...
  for (UInt_t i=0; i<photons->GetEntries(); ++i) {
    const Photon *ph = photons->At(i);

    if (ph->SCluster()->AbsEta()>=2.5 || (ph->SCluster()->AbsEta()>=1.4442 && ph->SCluster()->AbsEta()<=1.566))
      continue;
    if (ph->Et() < 20.)
      continue;
    if (ph->HadOverEm() > 0.15)
      continue;
    if (ph->AbsEta() < 1.5) {
      if (ph->CoviEtaiEta() > 0.013)
	continue;
    }
    else {
      if (ph->CoviEtaiEta() > 0.03)
	continue;
    }

    Bool_t passSpikeRemovalFilter = kTRUE;

    if (ph->SCluster() && ph->SCluster()->Seed()) {
      if (ph->SCluster()->Seed()->Energy() > 5.0 &&
	  ph->SCluster()->Seed()->EMax() / ph->SCluster()->Seed()->E3x3() > 0.95)
        passSpikeRemovalFilter = kFALSE;
   }
    if (fApplySpikeRemoval && !passSpikeRemovalFilter)
      continue;

    preselPh->Add(ph);
  }

  // sort both by pt... again ;)
  preselPh->Sort();

  unsigned int numPairs = 0;
  if (preselPh->GetEntries() > 0)
    numPairs = (preselPh->GetEntries()-1)*preselPh->GetEntries()/2;

  // create all possible pairs of pre-selected photons

  std::vector<unsigned int> idxFst;
  std::vector<unsigned int> idxSec;
  std::vector<bool> pairPasses;

  if (numPairs > 0) {
    for (unsigned int iFst = 0; iFst <preselPh->GetEntries() - 1; ++iFst) {
      for (unsigned int iSec = iFst + 1; iSec <preselPh->GetEntries(); ++iSec) {
        idxFst.push_back(iFst);
        idxSec.push_back(iSec);
        pairPasses.push_back(true);
      }
    }
  }

  // array to store the index of 'chosen Vtx' for each pair
  const Vertex** theVtx = new const Vertex*[numPairs];
  UInt_t* theVtxIdx     = new UInt_t[numPairs];

  // arays to store the Vtx 'fixed' photons
  Photon** fixPhFst = new Photon*[numPairs];
  Photon** fixPhSec = new Photon*[numPairs];

  // sotr pair-indices for pairs passing the selection
  std::vector<unsigned int> passPairs;
  passPairs.resize(0);

  unsigned int theChosenVtx = 0;

  float kinPh1[20];
  float kinPh2[20];

  for (int i=0; i<10;  ++i) {
    kinPh1[i] =-99.;
    kinPh2[i] =-99.;
  }

  float ptBefore1 = -99.;
  float ptBefore2 = -99.;

  bool print = false;
  if (evtNum == 17031 && false) {
    std::cout<<" ------------------------------------------- "<<std::endl;
    std::cout<<"   printing info for event #"<<evtNum<<std::endl;
    print = true;
  }

  for (unsigned int iPair = 0; iPair < numPairs; ++iPair) {

    // copy the photons for manipulation
    fixPhFst[iPair] = new Photon(*preselPh->At(idxFst[iPair]));
    fixPhSec[iPair] = new Photon(*preselPh->At(idxSec[iPair]));

    // if this is Data, scale the energy, if MC smear...
    FourVectorM scMomFst = fixPhFst[iPair]->Mom();
    FourVectorM scMomSec = fixPhSec[iPair]->Mom();
    double scaleFac1 = 1.;
    double scaleFac2 = 1.;
    if (fIsData) {
      if ( fRunStart.size() > 0) {
        // find run in rnage
        Int_t runRange=-1;
        for (UInt_t iRun = 0; iRun<fRunStart.size(); ++iRun) {
          if ( runNumber >= fRunStart[iRun] && runNumber <= fRunEnd[iRun]) {
            runRange = (Int_t) iRun;
            break;
          }
        }
        if (runRange > -1) {
          if (fixPhFst[iPair]->IsEB())
            if (fixPhFst[iPair]->R9() > 0.94)
              scaleFac1 += fDataEnCorr_EB_hR9[runRange];
            else
              scaleFac1 += fDataEnCorr_EB_lR9[runRange];
          else
            if (fixPhFst[iPair]->R9() > 0.94)
              scaleFac1 += fDataEnCorr_EE_hR9[runRange];
            else
              scaleFac1 += fDataEnCorr_EE_lR9[runRange];
          if (fixPhSec[iPair]->IsEB())
            if (fixPhSec[iPair]->R9() > 0.94)
              scaleFac2 += fDataEnCorr_EB_hR9[runRange];
            else
              scaleFac2 += fDataEnCorr_EB_lR9[runRange];
          else
            if (fixPhSec[iPair]->R9() > 0.94)
              scaleFac2 += fDataEnCorr_EE_hR9[runRange];
            else
              scaleFac2 += fDataEnCorr_EE_lR9[runRange];
        }
      }
    } else {
      // get the smearing for MC photons..
      UInt_t seedBase = (UInt_t) evtNum + (UInt_t) _runNum + (UInt_t) _lumiSec;
      UInt_t seed1    = seedBase + (UInt_t) fixPhFst[iPair]->E() + (UInt_t) (TMath::Abs(10.*fixPhFst[iPair]->SCluster()->Eta()));
      UInt_t seed2    = seedBase + (UInt_t) fixPhSec[iPair]->E() + (UInt_t) (TMath::Abs(10.*fixPhSec[iPair]->SCluster()->Eta()));

      double width1 = 0.;
      double width2 = 0.;
      if (fixPhFst[iPair]->IsEB())
        if (fixPhFst[iPair]->R9() > 0.94)
          width1 = fMCSmear_EB_hR9;
        else
          width1 = fMCSmear_EB_lR9;
      else
        if (fixPhFst[iPair]->R9() > 0.94)
          width1 = fMCSmear_EE_hR9;
        else
          width1 = fMCSmear_EE_lR9;


      if (fixPhSec[iPair]->IsEB())
        if (fixPhSec[iPair]->R9() > 0.94)
          width2 = fMCSmear_EB_hR9;
        else
          width2 = fMCSmear_EB_lR9;
      else
        if (fixPhSec[iPair]->R9() > 0.94)
          width2 = fMCSmear_EE_hR9;
        else
          width2 = fMCSmear_EE_lR9;

      if (doMCSmear) {
        fRnd3->SetSeed(seed1);
        scaleFac1 = fRnd3->Gaus(1.,width1);
        fRnd3->SetSeed(seed2);
        scaleFac2 = fRnd3->Gaus(1.,width2);
      }
    }

    if (iPair==0) {
      ptBefore1 = fixPhFst[iPair]->Pt();
      ptBefore2 = fixPhSec[iPair]->Pt();
    }

    if (print && false) {
      std::cout<<" Photon Pair #"<<iPair+1<<std::endl;
      std::cout<<"      Ph1 px = "<<fixPhFst[iPair]->Mom().X()<<std::endl;
      std::cout<<"          py = "<<fixPhFst[iPair]->Mom().Y()<<std::endl;
      std::cout<<"          pz = "<<fixPhFst[iPair]->Mom().Z()<<std::endl;
      std::cout<<"           E = "<<fixPhFst[iPair]->Mom().E()<<std::endl;
      std::cout<<"           M = "<<fixPhFst[iPair]->Mom().M()<<std::endl;
      std::cout<<"      Ph2 px = "<<fixPhSec[iPair]->Mom().X()<<std::endl;
      std::cout<<"          py = "<<fixPhSec[iPair]->Mom().Y()<<std::endl;
      std::cout<<"          pz = "<<fixPhSec[iPair]->Mom().Z()<<std::endl;
      std::cout<<"           E = "<<fixPhSec[iPair]->Mom().E()<<std::endl;
      std::cout<<"           M = "<<fixPhSec[iPair]->Mom().M()<<std::endl;
    }

    fixPhFst[iPair]->SetMom(scaleFac1*scMomFst.X(), scaleFac1*scMomFst.Y(), scaleFac1*scMomFst.Z(), scaleFac1*scMomFst.E());
    fixPhSec[iPair]->SetMom(scaleFac2*scMomSec.X(), scaleFac2*scMomSec.Y(), scaleFac2*scMomSec.Z(), scaleFac2*scMomSec.E());

    if (print && false) {
      std::cout<<"      SF 1 = "<<scaleFac1<<std::endl;
      std::cout<<"      SF 2 = "<<scaleFac2<<std::endl;

      std::cout<<" Photon Pair #"<<iPair+1<<std::endl;
      std::cout<<"      Ph1 px = "<<fixPhFst[iPair]->Mom().X()<<std::endl;
      std::cout<<"          py = "<<fixPhFst[iPair]->Mom().Y()<<std::endl;
      std::cout<<"          pz = "<<fixPhFst[iPair]->Mom().Z()<<std::endl;
      std::cout<<"           E = "<<fixPhFst[iPair]->Mom().E()<<std::endl;
      std::cout<<"           M = "<<fixPhFst[iPair]->Mom().M()<<std::endl;
      std::cout<<"      Ph2 px = "<<fixPhSec[iPair]->Mom().X()<<std::endl;
      std::cout<<"          py = "<<fixPhSec[iPair]->Mom().Y()<<std::endl;
      std::cout<<"          pz = "<<fixPhSec[iPair]->Mom().Z()<<std::endl;
      std::cout<<"           E = "<<fixPhSec[iPair]->Mom().E()<<std::endl;
      std::cout<<"           M = "<<fixPhSec[iPair]->Mom().M()<<std::endl;
    }

    // store the vertex for this pair
    if (doVtxSelection) {
      unsigned int iVtx = FindBestVertex(fixPhFst[iPair],fixPhSec[iPair], pv, bsp, conversions, print);
      theVtx[iPair]    =  pv->At(iVtx);
      theVtxIdx[iPair] =  iVtx;
      if (iPair == 0) theChosenVtx = iVtx;
    } else
      theVtx[iPair] =  pv->At(0);


    // fix the kinematics for both events
    FourVectorM newMomFst = fixPhFst[iPair]->MomVtx(theVtx[iPair]->Position());
    FourVectorM newMomSec = fixPhSec[iPair]->MomVtx(theVtx[iPair]->Position());
    fixPhFst[iPair]->SetMom(newMomFst.X(), newMomFst.Y(), newMomFst.Z(), newMomFst.E());
    fixPhSec[iPair]->SetMom(newMomSec.X(), newMomSec.Y(), newMomSec.Z(), newMomSec.E());


    if (print && false) {
      std::cout<<"         Vtx = "<<theChosenVtx<<std::endl;
      std::cout<<" Photon Pair #"<<iPair+1<<std::endl;
      std::cout<<"      Ph1 px = "<<fixPhFst[iPair]->Mom().X()<<std::endl;
      std::cout<<"          py = "<<fixPhFst[iPair]->Mom().Y()<<std::endl;
      std::cout<<"          pz = "<<fixPhFst[iPair]->Mom().Z()<<std::endl;
      std::cout<<"           E = "<<fixPhFst[iPair]->Mom().E()<<std::endl;
      std::cout<<"           M = "<<fixPhFst[iPair]->Mom().M()<<std::endl;
      std::cout<<"      Ph2 px = "<<fixPhSec[iPair]->Mom().X()<<std::endl;
      std::cout<<"          py = "<<fixPhSec[iPair]->Mom().Y()<<std::endl;
      std::cout<<"          pz = "<<fixPhSec[iPair]->Mom().Z()<<std::endl;
      std::cout<<"           E = "<<fixPhSec[iPair]->Mom().E()<<std::endl;
      std::cout<<"           M = "<<fixPhSec[iPair]->Mom().M()<<std::endl;
    }



    if (iPair != 0) {
      // check if both photons pass the CiC selection
      bool pass1 = PhotonTools::PassCiCSelection(fixPhFst[iPair], theVtx[iPair], tracks, electrons, pv, _tRho, 40., true, false);
      bool pass2 = PhotonTools::PassCiCSelection(fixPhSec[iPair], theVtx[iPair], tracks, electrons, pv, _tRho, 30., true, false);
      if ( pass1 && pass2 )
        passPairs.push_back(iPair);
    } else {
      bool pass1 = PhotonTools::PassCiCSelection(fixPhFst[iPair], theVtx[iPair], tracks, electrons, pv, _tRho, 40., true, false, kinPh1);
      bool pass2 = PhotonTools::PassCiCSelection(fixPhSec[iPair], theVtx[iPair], tracks, electrons, pv, _tRho, 30., true, false, kinPh2);

      if ( pass1 && pass2 )
        passPairs.push_back(iPair);
    }
  }

  // loop over all passing pairs and find the one with the highest sum Et
  Photon* phHard = NULL;
  Photon* phSoft = NULL;

  const Vertex* _theVtx = NULL;

  double maxSumEt = 0.;
  for (unsigned int iPair=0; iPair<passPairs.size(); ++iPair){
    double sumEt = fixPhFst[passPairs[iPair]]->Et();
    sumEt += fixPhSec[passPairs[iPair]]->Et();
    if ( sumEt > maxSumEt ) {
      maxSumEt = sumEt;
      phHard = fixPhFst[passPairs[iPair]];
      phSoft = fixPhSec[passPairs[iPair]];
      _theVtx = theVtx[iPair];
      theChosenVtx = theVtxIdx[iPair];
    }
  }

  Float_t _theVtxZ = -999.;
  Float_t catPh1 = 0.;
  Float_t catPh2 = 0.;
  Float_t catEvt = -1.;

  bool ph1IsLowR9 = false;
  bool ph2IsLowR9 = false;

  if (phHard && phSoft) {
    GoodPhotons->AddOwned(phHard);
    GoodPhotons->AddOwned(phSoft);
    if (_theVtx)
      _theVtxZ=_theVtx->Position().Z();

    catPh1=1.;
    catPh2=1.;
    catEvt=0.;

    if (phHard->SCluster()->AbsEta()>1.5)
      catPh1=3;
    if (phHard->R9() < 0.94) {
      catPh1=catPh1+1;
      ph1IsLowR9 = true;
    }

    if (phSoft->SCluster()->AbsEta()>1.5)
      catPh2=3;

    if (phSoft->R9() < 0.94) {
      catPh2=catPh2+1;
      ph2IsLowR9 = true;
    }

    if (catPh1 > 2.5 || catPh2 > 2.5)
      catEvt=2.;
    if (ph1IsLowR9 || ph2IsLowR9)
      catEvt=catEvt+1.;
  }

  // sort according to pt
  GoodPhotons->Sort();

  // add to event for other modules to use
  AddObjThisEvt(GoodPhotons);

  delete preselPh;

  bool doFill = (phHard && phSoft);
  Float_t _mass = ( doFill ? (phHard->Mom()+phSoft->Mom()).M() : -100.);
  Float_t _ptgg = ( doFill ? (phHard->Mom()+phSoft->Mom()).Pt() : -100.);
  if (_ptgg < 40. && doFill)
    catEvt = catEvt+4.;

  Float_t _pth    = -100.;
  Float_t _decayZ = -100.;
  if (!fIsData)
    FindHiggsPtAndZ(mcParticles, _pth, _decayZ);

  Float_t fillEvent[] = { (float)_tRho,
                          _pth,
                          _decayZ,
                          _theVtxZ,
                          numPU,
                          _mass,
                          _ptgg,
                          _evtNum1,
                          _evtNum2,
                          _runNum,
                          _lumiSec,
                          (float) theChosenVtx,
                          (float) numPairs,
                          catPh1,
                          catPh2,
                          catEvt,
                          kinPh1[0],
                          kinPh1[1],
                          kinPh1[2],
                          kinPh1[3],
                          kinPh1[4],
                          kinPh1[5],
                          kinPh1[6],
                          kinPh1[7],
                          kinPh1[8],
                          kinPh1[9],
                          kinPh1[10],
                          kinPh1[11],
                          kinPh1[12],
                          kinPh1[13],
                          kinPh1[14],
                          kinPh1[15],
                          kinPh1[16],
                          kinPh1[17],
                          kinPh1[18],
                          kinPh1[19],
                          kinPh2[0],
                          kinPh2[1],
                          kinPh2[2],
                          kinPh2[3],
                          kinPh2[4],
                          kinPh2[5],
                          kinPh2[6],
                          kinPh2[7],
                          kinPh2[8],
                          kinPh2[9],
                          kinPh2[10],
                          kinPh2[11],
                          kinPh2[12],
                          kinPh2[13],
                          kinPh2[14],
                          kinPh2[15],
                          kinPh2[16],
                          kinPh2[17],
                          kinPh2[18],
                          kinPh2[19],
                          ptBefore1,
                          ptBefore2
  };



  if (_mass > 0.) {
    //std::cout<<catPh1<<"  "<<catPh2<<"  "<<catEvt<<"  "<<_mass<<std::endl;
    hCiCTuple->Fill(fillEvent);
  }

  delete[] theVtx;
  delete[] theVtxIdx;

  return;

}

//--------------------------------------------------------------------------------------------------
void PhotonCiCMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we just request the photon collection branch.

  hCiCTuple = new TNtuple("hCiCTuple","hCiCTuple","rho:higgspt:higgsZ:vtxZ:numPU:mass:ptgg:evtnum1:evtnum2:runnum:lumisec:ivtx:npairs:ph1Cat:ph2Cat:evtCat:ph1Iso1:ph1Iso2:ph1Iso3:ph1Cov:ph1HoE:ph1R9:ph1DR:ph1Pt:ph1Eta:ph1Phi:ph1Eiso3:ph1Eiso4:ph1Hiso4:ph1TisoA:ph1TisoW:ph1Tiso:ph1Et:ph1E:ph1Pass:ph1CatDebug:ph2Iso1:ph2Iso2:ph2Iso3:ph2Cov:ph2HoE:ph2R9:ph2DR:ph2Pt:ph2Eta:ph2Phi:ph2Eiso3:ph2Eiso4:ph2Hiso4:ph2TisoA:ph2TisoW:ph2Tiso:ph2Et:ph2E:ph2Pass:ph2CatDebug:ph1UPt:ph2UPt");

  AddOutput(hCiCTuple);

}

// return the index of the bext vertex
unsigned int PhotonCiCMod::FindBestVertex(Photon* ph1, Photon* ph2, VertexCol const* pv, const BaseVertex *bsp, DecayParticleCol const* conversions, bool print) {

  // loop over all vertices and assigne the ranks
  int* ptbal_rank  = new int[pv->GetEntries()];
  int* ptasym_rank = new int[pv->GetEntries()];
  int* total_rank  = new int[pv->GetEntries()];
  double* ptbal = new double[pv->GetEntries()];
  double* ptasym = new double[pv->GetEntries()];

  unsigned int numVertices = pv->GetEntries();

  double ptgg = 0.;    // stored for later in the conversion

  if (print && false) {
    std::cout<<" --------------------------------- "<<std::endl;
    std::cout<<"   looking for Vtx for photon pair "<<std::endl;
    std::cout<<"                            pt 1 = "<<ph1->Et()<<std::endl;
    std::cout<<"                            pt 2 = "<<ph2->Et()<<std::endl;
    std::cout<<"                    among #"<<numVertices<<" Vertices."<<std::endl;
  }


  for (unsigned int iVtx = 0; iVtx < numVertices; ++iVtx) {
    if (print)
      std::cout<<std::endl<<"       Vertex #"<<iVtx<<std::endl;

    const Vertex* tVtx = pv->At(iVtx);
    ptbal [iVtx] = 0.0;
    ptasym[iVtx] = 0.0;
    ptbal_rank [iVtx] = 1;
    ptasym_rank[iVtx] = 1;

    // compute the photon momenta with respect to this Vtx
    FourVectorM newMomFst = ph1->MomVtx(tVtx->Position());
    FourVectorM newMomSec = ph2->MomVtx(tVtx->Position());

    FourVectorM higgsMom = newMomFst+newMomSec;

    double ph1Eta = newMomFst.Eta();
    double ph2Eta = newMomSec.Eta();

    double ph1Phi = newMomFst.Phi();
    double ph2Phi = newMomSec.Phi();

    if (print && iVtx == 0 && false ) {
      std::cout<<"     new momenta Et1 = "<<newMomFst.Et()<<std::endl;
      std::cout<<"                 Eta = "<<newMomFst.Eta()<<std::endl;
      std::cout<<"                 Phi = "<<newMomFst.Phi()<<std::endl;
      std::cout<<"                  Px = "<<newMomFst.Px()<<std::endl;
      std::cout<<"                  Py = "<<newMomFst.Py()<<std::endl;
      std::cout<<"     new momenta Et2 = "<<newMomSec.Et()<<std::endl;
      std::cout<<"                 Eta = "<<newMomSec.Eta()<<std::endl;
      std::cout<<"                 Phi = "<<newMomSec.Phi()<<std::endl;
      std::cout<<"                  Px = "<<newMomSec.Px()<<std::endl;
      std::cout<<"                  Py = "<<newMomSec.Py()<<std::endl;
    }

    FourVectorM totTrkMom;
    for (unsigned int iTrk = 0; iTrk < tVtx->NTracks(); ++iTrk) {
      const Track* tTrk = tVtx->Trk(iTrk);
      //if (tTrk->Pt()<1.) continue;
      double tEta = tTrk->Eta();
      double tPhi = tTrk->Phi();
      double dEta1 = TMath::Abs(tEta-ph1Eta);
      double dEta2 = TMath::Abs(tEta-ph2Eta);
      double dPhi1 = TMath::Abs(tPhi-ph1Phi);
      double dPhi2 = TMath::Abs(tPhi-ph2Phi);
      if (dPhi1 > M_PI) dPhi1 = 2*M_PI - dPhi1;
      if (dPhi2 > M_PI) dPhi2 = 2*M_PI - dPhi2;

      double dR1 = TMath::Sqrt(dEta1*dEta1+dPhi1*dPhi1);
      double dR2 = TMath::Sqrt(dEta2*dEta2+dPhi2*dPhi2);

      if ( ( iVtx == 0 || iVtx == 1 ) && print && false) {
        std::cout<<"  Track #"<<iTrk<<std::endl;
        std::cout<<"      pt = "<<tTrk->Pt()<<std::endl;
        std::cout<<"     eta = "<<tTrk->Eta()<<std::endl;
        std::cout<<"     phi = "<<tTrk->Phi()<<std::endl;
        std::cout<<"      px = "<<tTrk->Px()<<std::endl;
        std::cout<<"      py = "<<tTrk->Py()<<std::endl;
        std::cout<<"     dR1 = "<<dR1<<std::endl;
        std::cout<<"     dR2 = "<<dR2<<std::endl;
      }

      if (dR1 < 0.05 || dR2 < 0.05) continue;

      if (iTrk == 0) totTrkMom = tTrk->Mom4(0);
      else totTrkMom = totTrkMom + tTrk->Mom4(0);

    }

    if (iVtx ==0 && print && false) {
      std::cout<<"   Trk passes cuts: "<<std::endl;
      std::cout<<"            px tot = "<<totTrkMom.Px()<<std::endl;
      std::cout<<"            py tot = "<<totTrkMom.Py()<<std::endl;
    }

    double ptvtx   = totTrkMom.Pt();
    if (iVtx ==0 && print && false)
      std::cout<<"    Total TkPt = "<<ptvtx<<std::endl;
    double pthiggs = higgsMom.Pt();
    if (iVtx ==0 && print && false)
      std::cout<<"    Total H Pt = "<<pthiggs<<std::endl;
    if (iVtx == 0) ptgg = pthiggs;
    ptbal [iVtx]  = (totTrkMom.Px()*(newMomFst.Px()+newMomSec.Px()));
    ptbal [iVtx] += (totTrkMom.Py()*(newMomFst.Py()+newMomSec.Py()));
    ptbal [iVtx]  = -ptbal[iVtx]/pthiggs;
    ptasym[iVtx] = (ptvtx - pthiggs)/(ptvtx + pthiggs);

    if (iVtx ==0 && print && false) {
      std::cout<<"    Results: ptbal = "<<ptbal [iVtx]<<std::endl;
      std::cout<<"            ptasym = "<<ptasym[iVtx]<<std::endl;
    }

    for (unsigned int cVtx =0; cVtx < iVtx; ++cVtx) {
      if (ptbal [iVtx] > ptbal [cVtx])
        ptbal_rank[cVtx]++;
      else
        ptbal_rank[iVtx]++;
      if (ptasym [iVtx] > ptasym [cVtx])
        ptasym_rank[cVtx]++;
      else
        ptasym_rank[iVtx]++;
    }
  }


  // compute the total rank
  for (unsigned int iVtx = 0; iVtx < numVertices; ++iVtx) {
    if (print && false) {
      std::cout<<"     Vertex #"<<iVtx<<"  has rank PTB "<<ptbal_rank[iVtx]<<"   ("<<ptbal[iVtx]<<")"<<std::endl;
      std::cout<<"     Vertex #"<<iVtx<<"  has rank PTSYM "<<ptasym_rank[iVtx]<<"   ("<<ptasym[iVtx]<<")"<<std::endl;
    }
    ptasym_rank [iVtx] = ptbal_rank [iVtx]*ptasym_rank [iVtx]*(iVtx+1);
    total_rank [iVtx] = 0;
    for (unsigned int cVtx =0; cVtx < iVtx; ++cVtx) {
      if (ptasym_rank [iVtx] > ptasym_rank [cVtx])
        total_rank[iVtx]++;
      else if (ptasym_rank [iVtx] == ptasym_rank [cVtx]) {
        if (ptbal_rank [iVtx] > ptbal_rank [cVtx])
          total_rank[iVtx]++;
        else
          total_rank[cVtx]++;
      }
      else
        total_rank[cVtx]++;
    }
  }

  if (print && false) std::cout<<std::endl;

  unsigned int bestIdx = 0;
  for (unsigned int iVtx = 0; iVtx < numVertices; ++iVtx) {
    if (total_rank[iVtx] == 0) bestIdx = iVtx;
    if (print && false) {
      std::cout<<"     Vertex #"<<iVtx<<"  has rank "<<total_rank[iVtx]<<std::endl;
    }
  }

  delete[] ptbal_rank  ;
  delete[] ptasym_rank ;
  delete[] ptbal       ;
  delete[] ptasym      ;

  //return bestIdx;

  // check if there's a conversion among the pre-selected photons
  const DecayParticle* conv1 = PhotonTools::MatchedCiCConversion(ph1, conversions, 0.1, 0.1, 0.1, print);
  const DecayParticle* conv2 = PhotonTools::MatchedCiCConversion(ph2, conversions, 0.1, 0.1, 0.1, print);

  if (print && false) {
    if (conv1) {
      std::cout<<" Photon 1 has has conversion with P = "<<conv1->Prob()<<std::endl;
      std::cout<<"                                Rho = "<<conv1->Position().Rho()<<std::endl;
      std::cout<<"                                  Z = "<<conv1->Position().Z()<<std::endl;
    }
    if (conv2) {
      std::cout<<" Photon 2 has has conversion with P = "<<conv2->Prob()<<std::endl;
      std::cout<<"                                Rho = "<<conv2->Position().Rho()<<std::endl;
      std::cout<<"                                  Z = "<<conv2->Position().Z()<<std::endl;
    }
  }

  if (conv1 && ( conv1->Prob() < 0.0005))
    conv1 = NULL;
  if (conv2 && ( conv2->Prob() < 0.0005))
    conv2 = NULL;

  double zconv  = 0.;
  double dzconv = 0.;

  if (conv1 || conv2) {
    if      (!conv2) {
      //const mithep::ThreeVector caloPos1(ph1->SCluster()->Point());
      const mithep::ThreeVector caloPos1(ph1->CaloPos());
      zconv  = conv1->Z0EcalVtx(bsp->Position(), caloPos1);
      if (ph1->IsEB()) {
        double rho = conv1->Position().Rho();
        if      (rho < 15.)
          dzconv = 0.06;
        else if (rho < 60.)
          dzconv = 0.67;
        else
          dzconv = 2.04;
      }
      else {
        double z = conv1->Position().Z();
        if      (TMath::Abs(z) < 50. )
          dzconv = 0.18;
        else if (TMath::Abs(z) < 100.)
          dzconv = 0.61;
        else
          dzconv = 0.99;
      }
    }
    else if (!conv1) {
      //const mithep::ThreeVector caloPos2(ph2->SCluster()->Point());
      const mithep::ThreeVector caloPos2(ph2->CaloPos());
      zconv  = conv2->Z0EcalVtx(bsp->Position(), caloPos2);
      if (ph2->IsEB()) {
        double rho = conv2->Position().Rho();
        if      (rho < 15.)
          dzconv = 0.06;
        else if (rho < 60.)
          dzconv = 0.67;
        else
          dzconv = 2.04;
      }
      else {
        double z = conv2->Position().Z();
        if      (TMath::Abs(z) < 50. )
          dzconv = 0.18;
        else if (TMath::Abs(z) < 100.)
          dzconv = 0.61;
        else
          dzconv = 0.99;
      }
    }
    else {
      //const mithep::ThreeVector caloPos1(ph1->SCluster()->Point());
      const mithep::ThreeVector caloPos1(ph1->CaloPos());
      double z1  = conv1->Z0EcalVtx(bsp->Position(), caloPos1);
      double dz1 = 0.;
      if (ph1->IsEB()) {
        double rho = conv1->Position().Rho();
        if      (rho < 15.)
          dz1 = 0.06;
        else if (rho < 60.)
          dz1 = 0.67;
        else
          dz1 = 2.04;
      }
      else {
        double z = conv1->Position().Z();
        if      (TMath::Abs(z) < 50. )
          dz1 = 0.18;
        else if (TMath::Abs(z) < 100.)
          dz1 = 0.61;
        else
          dz1 = 0.99;
      }
      //const mithep::ThreeVector caloPos2(ph2->SCluster()->Point());
      const mithep::ThreeVector caloPos2(ph2->CaloPos());
      double z2  = conv2->Z0EcalVtx(bsp->Position(), caloPos2);
      double dz2 = 0.;
      if ( ph2->IsEB() ) {
        double rho = conv2->Position().Rho();
        if      (rho < 15.)
          dz2 = 0.06;
        else if (rho < 60.)
          dz2 = 0.67;
        else
          dz2 = 2.04;
      }
      else {
        double z = conv2->Position().Z();
        if      (TMath::Abs(z) < 50. )
          dz2 = 0.18;
        else if (TMath::Abs(z) < 100.)
          dz2 = 0.61;
        else
          dz2 = 0.99;
      }

      if (print) {
        std::cout<<"  z1 = "<<z1<<std::endl;
        std::cout<<" dz1 = "<<dz1<<std::endl;
        std::cout<<"  z2 = "<<z2<<std::endl;
        std::cout<<" dz2 = "<<dz2<<std::endl;
      }

      zconv  = ( 1./(1./dz1/dz1 + 1./dz2/dz2 )*(z1/dz1/dz1 + z2/dz2/dz2) ) ;  // weighted average
      dzconv = TMath::Sqrt( 1./(1./dz1/dz1 + 1./dz2/dz2)) ;
    }

    if (print) {
      std::cout<<"  Conversion Z = "<<zconv<<std::endl;
      std::cout<<"            dZ = "<<dzconv<<std::endl;
    }

    if (true) {
      // loop over all ranked Vertices and choose the closest to the Conversion one
      int maxVertices = ( ptgg > 30 ? 3 : 5);
      double minDz = -1.;

      if (print)
        std::cout<<std::endl<<"  looping over vertices... "<<ptgg<<"  "<<maxVertices<<std::endl;

      for (unsigned int iVtx =0; iVtx < numVertices; ++iVtx) {

        if (print)
          std::cout<<"     "<<iVtx<<"   has rank  "<<total_rank[iVtx]<<std::endl;

        if (total_rank[iVtx] < maxVertices) {
          const Vertex* tVtx = pv->At(iVtx);
          double tDz = TMath::Abs(zconv - tVtx->Z());
          if (print)
            std::cout<<"     is considered with tDz = "<<tDz<<std::endl;
          if ((minDz < 0. || tDz < minDz) && ( tDz < dzconv)) {
            minDz = tDz;
            bestIdx = iVtx;
            if (print) std::cout<<"      and is the best now."<<std::endl;
          }
        }
      }
    }
    else {
      unsigned int bestIdxTmp = bestIdx;

      // loop over all ranked Vertices and choose the closest to the Conversion one
      double minDz = -1.;
      int maxVertices = ( ptgg > 30 ? 3 : 5);

      if (print)
        std::cout<<std::endl<<"  looping over vertices... "<<ptgg<<"  "<<maxVertices<<std::endl;

      for (unsigned int iVtx =0; iVtx < numVertices; ++iVtx) {

        if (print)
          std::cout<<"     "<<iVtx<<"   has rank  "<<total_rank[iVtx]<<std::endl;

        const Vertex* tVtx = pv->At(iVtx);
        double tDz = TMath::Abs(zconv - tVtx->Z());
        if (print)
          std::cout<<"     is considered with tDz = "<<tDz<<std::endl;
        if ((minDz < 0. || tDz < minDz) && (tDz < dzconv)) {
          minDz = tDz;
          bestIdxTmp = iVtx;
          if (print)
	    std::cout<<"      and is the best now."<<std::endl;
        }
      }

      // check if best Vtx is among higest ranked ones
      if (total_rank[bestIdxTmp] < maxVertices)
        bestIdx = bestIdxTmp;
    }
  }

  delete[] total_rank  ;
  return bestIdx;
}

void PhotonCiCMod::FindHiggsPtAndZ(MCParticleCol const* mcParticles, Float_t& pt, Float_t& decayZ) {

  pt = -999.;
  decayZ = -999.;

  // loop over all GEN particles and look for status 1 photons
  for (UInt_t i=0; i<mcParticles->GetEntries(); ++i) {
    const MCParticle* p = mcParticles->At(i);
    if (!(p->Is(MCParticle::kH)))
      continue;
    pt = p->Pt();
    decayZ = p->DecayVertex().Z();
    break;
  }

  return;
 }
