#include "MitPhysics/Utils/interface/PhotonTools.h"
#include "MitPhysics/Utils/interface/ElectronTools.h"
#include "MitPhysics/Utils/interface/IsolationTools.h"
#include "MitPhysics/Init/interface/Constants.h"
#include "MitAna/DataTree/interface/StableData.h"
#include <TRandom3.h>

#include <limits>

ClassImp(mithep::PhotonTools)

using namespace mithep;

//--------------------------------------------------------------------------------------------------
PhotonTools::PhotonTools()  
{
  // Constructor.
}

Double_t
mithep::PhotonTools::PhotonEffectiveArea(EPhotonEffectiveAreaType type, Double_t SCEta, EPhotonEffectiveAreaTarget target)
{
  if (target == kPhoEANoCorr)
    return 0.;

  double etaBinning1[] = {0., 1., 1.479, 2., 2.2, 2.3, 2.4, std::numeric_limits<double>::max()};
  double etaBinning2[] = {0., 0.9, 1.479, 2., 2.2, std::numeric_limits<double>::max()};

  double* etaBinning = 0;
  std::vector<double> areas;

  switch (target) {
  case kPhoEAPhys14:
    etaBinning = etaBinning1;
    switch (type) {
    case kPhoChargedHadron03:
      areas = {0.0234, 0.0189, 0.0171, 0.0129, 0.011, 0.0074, 0.0034};
      break;
    case kPhoNeutralHadron03:
      areas = {0.0053, 0.0103, 0.0057, 0.007, 0.0152, 0.0232, 0.1709};
      break;
    case kPhoPhoton03:
      areas = {0.078, 0.0629, 0.0264, 0.0462, 0.074, 0.0924, 0.1484};
      break;
    default:
      return 0.;
    }

  case kPhoEAHighPtV2:
    etaBinning = etaBinning2;
    switch (type) {
    case kPhoPhoton03:
      areas = {0.17, 0.14, 0.11, 0.14, 0.22};
      break;
    default:
      return 0.;
    }

  case kPhoEASpring1550ns:
    etaBinning = etaBinning1;
    switch (type) {
    case kPhoChargedHadron03:
      areas = {0.0157, 0.0143, 0.0115, 0.0094, 0.0095, 0.0068, 0.0053};
      break;
    case kPhoNeutralHadron03:
      areas = {0.0143, 0.0210, 0.0147, 0.0082, 0.0124, 0.0186, 0.0320};
      break;
    case kPhoPhoton03:
      areas = {0.0725, 0.0604, 0.0320, 0.0512, 0.0766, 0.0949, 0.1160};
      break;
    default:
      return 0.;
    }

 case kPhoEASpring15:
    etaBinning = etaBinning1;
    switch (type) {
    case kPhoNeutralHadron03:
      areas = {0.0599, 0.0819, 0.0696, 0.0360, 0.0360, 0.0462, 0.0656};
      break;
    case kPhoPhoton03:
      areas = {0.1271, 0.1101, 0.0756, 0.1175, 0.1498, 0.1857, 0.2183};
      break;
    default:
      return 0.;
    }

  default:
    return 0.;
  }

  double absEta = std::abs(SCEta);
  unsigned etaBin = 0;
  while (absEta >= etaBinning[etaBin + 1])
    ++etaBin;
  return areas.at(etaBin);
}

Bool_t
mithep::PhotonTools::PassID(Photon const* pho, EPhIdType type)
{
  double hOverECut, sigmaIEtaIEtaCut;
  bool isEB = pho->SCluster()->AbsEta() < gkPhoEBEtaMax;

  switch (type) {
  case kPhys14Loose:
    hOverECut        = isEB ? 0.028  : 0.093;
    sigmaIEtaIEtaCut = isEB ? 0.0107 : 0.0272;
    break;
  case kPhys14Medium:
    hOverECut        = isEB ? 0.012  : 0.023;
    sigmaIEtaIEtaCut = isEB ? 0.0100 : 0.0267;
    break;
  case kPhys14Tight:
    hOverECut        = isEB ? 0.010  : 0.015;
    sigmaIEtaIEtaCut = isEB ? 0.0100 : 0.0265;
    break;
  case kHighPtV2:
    hOverECut        = isEB ? 0.05 : 0.05;
    sigmaIEtaIEtaCut = isEB ? 0.0105 : 0.028;
    break;
  case kSpring15Loose50ns:
    hOverECut        = isEB ? 0.05 : 0.05;
    sigmaIEtaIEtaCut = isEB ? 0.0103 : 0.0277;
    break;
  case kSpring15Medium50ns:
    hOverECut        = isEB ? 0.05 : 0.05;
    sigmaIEtaIEtaCut = isEB ? 0.0100 : 0.0267;
    break;
  case kSpring15Tight50ns:
    hOverECut        = isEB ? 0.05 : 0.05;
    sigmaIEtaIEtaCut = isEB ? 0.0100 : 0.0267;
    break;
  case kSpring15Loose:
    hOverECut        = isEB ? 0.05 : 0.05;
    sigmaIEtaIEtaCut = isEB ? 0.0102 : 0.0274;
    break;
  case kSpring15Medium:
    hOverECut        = isEB ? 0.05 : 0.05;
    sigmaIEtaIEtaCut = isEB ? 0.0102 : 0.0268;
    break;
  case kSpring15Tight:
    hOverECut        = isEB ? 0.05 : 0.05;
    sigmaIEtaIEtaCut = isEB ? 0.0100 : 0.0268;
    break;
  default:
    return false;
  }

  if (pho->HadOverEmTow() > hOverECut)
    return false;

  //This is for backwards compatibility
  double ietaieta = -1;
  if (pho->CoviEtaiEta5x5() < 0) 
    ietaieta = pho->CoviEtaiEta();
  else 
    ietaieta = pho->CoviEtaiEta5x5();

  //check the cut
  if (ietaieta > sigmaIEtaIEtaCut)
    return false;

  return true;
}

Bool_t
mithep::PhotonTools::PassIsoFootprintRhoCorr(Photon const* pho, EPhIsoType isoType, Vertex const* pv, PFCandidateCol const* pfCands, Double_t rho)
{
  double chIso = 0.;
  double nhIso = 0.;
  double phIso = 0.;

  IsolationTools::PFEGIsoFootprintRemoved(pho, pv, pfCands, 0.3, chIso, nhIso, phIso);

  return PassIsoRhoCorr(pho, isoType, rho, chIso, nhIso, phIso);
}

Bool_t
mithep::PhotonTools::PassIsoRhoCorr(Photon const* pho, EPhIsoType isoType, Double_t rho)
{
  return PassIsoRhoCorr(pho, isoType, rho, pho->PFChargedHadronIso(), pho->PFNeutralHadronIso(), pho->PFPhotonIso());
}

Bool_t
mithep::PhotonTools::PassIsoRhoCorr(Photon const* pho, EPhIsoType isoType, Double_t rho, Double_t chIso, Double_t nhIso, Double_t phIso)
{
  double scEta = pho->SCluster()->AbsEta();
  bool isEB = scEta < gkPhoEBEtaMax;
  double pEt = pho->Et();

  double chIsoCut = 0.;
  double nhIsoCut = 0.;
  double phIsoCut = 0.;

  EPhotonEffectiveAreaTarget eaTarget;

  switch (isoType){
  case kPhys14LooseIso:
    chIsoCut = isEB ? 2.67 : 1.79;
    nhIsoCut = isEB ? (7.23 + TMath::Exp(0.0028 * pEt + 0.5408)) : (8.89 + 0.01725 * pEt);
    phIsoCut = isEB ? (2.11 + 0.0014 * pEt) : (3.09 + 0.0091 * pEt);
    eaTarget = kPhoEAPhys14;
    break;
  case kPhys14MediumIso:
    chIsoCut = isEB ? 1.79 : 1.09;
    nhIsoCut = isEB ? (0.16 + TMath::Exp(0.0028 * pEt + 0.5408)) : (4.31 + 0.0172 * pEt);
    phIsoCut = isEB ? (1.90 + 0.0014 * pEt) : (1.90 + 0.0091 * pEt);
    eaTarget = kPhoEAPhys14;
    break;
  case kPhys14TightIso:
    chIsoCut = isEB ? 1.66 : 1.04;
    nhIsoCut = isEB ? (0.14 + TMath::Exp(0.0028 * pEt + 0.5408)) : (3.89 + 0.0172 * pEt);
    phIsoCut = isEB ? (1.40 + 0.0014 * pEt) : (1.40 + 0.0091 * pEt);
    eaTarget = kPhoEAPhys14;
    break;
  case kHighPtV2Iso:
    chIsoCut = 5.;
    nhIsoCut = std::numeric_limits<double>::max();
    phIsoCut = 2.5 + (isEB ? 2.75 : 2.) - (scEta < 2. ? 0.0045 : 0.003) * pEt;
    eaTarget = kPhoEAPhys14;
    break;
  case kSpring15Loose50nsIso:
    chIsoCut = isEB ? 2.44 : 1.84;
    nhIsoCut = isEB ? (2.57 + TMath::Exp(0.0044 * pEt + 0.5809)) : (4. + TMath::Exp(0.004 * pEt + 0.9402));
    phIsoCut = isEB ? (1.92 + 0.0043 * pEt) : (2.15 + 0.0041 * pEt);
    eaTarget = kPhoEASpring1550ns;
    break;
  case kSpring15Medium50nsIso:
    chIsoCut = isEB ? 1.31 : 1.25;
    nhIsoCut = isEB ? (0.60 + TMath::Exp(0.0044 * pEt + 0.5809)) : (1.65 + TMath::Exp(0.004 * pEt + 0.9402));
    phIsoCut = isEB ? (1.33 + 0.0043 * pEt) : (1.02 + 0.0041 * pEt);
    eaTarget = kPhoEASpring1550ns;
    break;
  case kSpring15Tight50nsIso:
    chIsoCut = isEB ? 0.91 : 0.65;
    nhIsoCut = isEB ? (0.33 + TMath::Exp(0.0044 * pEt + 0.5809)) : (0.93 + TMath::Exp(0.004 * pEt + 0.9402));
    phIsoCut = isEB ? (0.61 + 0.0043 * pEt) : (0.54 + 0.0041 * pEt);
    eaTarget = kPhoEASpring1550ns;
    break;
  case kSpring15LooseIso:
    chIsoCut = isEB ? 3.32 : 1.97;
    nhIsoCut = isEB ? (1.92 + 0.014 * pEt + 0.000019 * pEt * pEt) : (11.86 + 0.0139 * pEt + 0.000025 * pEt * pEt);
    phIsoCut = isEB ? (0.81 + 0.0053 * pEt) : (0.83 + 0.0034 * pEt);
    eaTarget = kPhoEASpring15;
    break;
  case kSpring15MediumIso:
    chIsoCut = isEB ? 1.37 : 1.10;
    nhIsoCut = isEB ? (1.06 + 0.014 * pEt + 0.000019 * pEt * pEt) : (2.69 + 0.0139 * pEt + 0.000025 * pEt * pEt);
    phIsoCut = isEB ? (0.28 + 0.0053 * pEt) : (0.39 + 0.0034 * pEt);
    eaTarget = kPhoEASpring15;
    break;
  case kSpring15TightIso:
    chIsoCut = isEB ? 1.66 : 1.04;
    nhIsoCut = isEB ? (0.97 + 0.014 * pEt + 0.000019 * pEt * pEt) : (2.09 + 0.0139 * pEt + 0.000025 * pEt * pEt);
    phIsoCut = isEB ? (0.08 + 0.0053 * pEt) : (0.16 + 0.0034 * pEt);
    eaTarget = kPhoEASpring15;
    break;
  default:
    return false;
  }

  if (IsolationTools::PFPhotonIsolationRhoCorr(scEta, chIso, rho, eaTarget, kPhoChargedHadron03) > chIsoCut)
    return false;
  if (IsolationTools::PFPhotonIsolationRhoCorr(scEta, nhIso, rho, eaTarget, kPhoNeutralHadron03) > nhIsoCut)
    return false;
  if (IsolationTools::PFPhotonIsolationRhoCorr(scEta, phIso, rho, eaTarget, kPhoPhoton03) > phIsoCut)
    return false;

  return true;
}

//--------------------------------------------------------------------------------------------------
PhotonTools::eScaleCats PhotonTools::EScaleCat(const Photon *p)
{
  if (p->SCluster()->AbsEta()<1.0) {
    if (p->R9()>0.94) {
      return kEBlowEtaGold;
    }
    else return kEBlowEtaBad;
  }
  else if (p->SCluster()->AbsEta()<1.5) {
    if (p->R9()>0.94) return kEBhighEtaGold;
    else return kEBhighEtaBad;    
  }
  else if (p->SCluster()->AbsEta()<2.0) {
    if (p->R9()>0.94) return kEElowEtaGold;
    else return kEElowEtaBad;    
  }
  else {
    if (p->R9()>0.94) return kEEhighEtaGold;
    else return kEEhighEtaBad;    
  }  
  
}


void PhotonTools::ScalePhoton(Photon* p, Double_t scale) {
  if( !p ) return;
  FourVectorM mom = p->Mom();
  p->SetMom(scale*mom.X(), scale*mom.Y(), scale*mom.Z(), scale*mom.E());
  
  double oldscale = std::max(0.,p->EnergyScale());
  
  p->SetEnergyScale(oldscale*scale);
  
}

void PhotonTools::SmearPhoton(Photon* p, TRandom3* rng, Double_t width, UInt_t iSeed) {
  
  if( !p  )       return;
  if( !rng)       return;
  if( width < 0.) return;

  rng->SetSeed(iSeed);
  FourVectorM mom = p->Mom();
  Double_t scale = -1.;

  while (scale<0.) {
    scale = rng->Gaus(1.,width);
  }
  
  p->SetMom(scale*mom.X(), scale*mom.Y(), scale*mom.Z(), scale*mom.E());

  // moved here from SmearPhotonError, in order to being able to
  // use different smearing for Error and actual Momentum (any reason not to? Josh?)
  Double_t smear = p->EnergySmearing();
  p->SetEnergySmearing(TMath::Sqrt(smear*smear+width*width*p->E()*p->E()));
    

  return;
}

void PhotonTools::SmearPhotonError(Photon* p, Double_t width) {
  
  if( !p  )       return;
  if( width < 0.) return;  
  
  // WARNING: This is not set correctly if different smearing is used for error and actual momentum...
  // therefore move this into 'SmearPhoton'.... any reason not to ? Josh ?
  //Double_t smear = p->EnergySmearing();
  //p->SetEnergySmearing(TMath::Sqrt(smear*smear+width*width*p->E()*p->E()));
  
  Double_t err = p->EnergyErrSmeared();
  if (err>=0.0) {
    p->SetEnergyErrSmeared(TMath::Sqrt(err*err + width*width));
    
    Double_t smear = p->EnergyErrSmearing();    
    p->SetEnergyErrSmearing(TMath::Sqrt(smear*smear+width*width));    
  }
  
}

void PhotonTools::ScalePhotonR9(Photon* p, Double_t scale) {
 p->SetR9(scale*p->R9()); 
}


void PhotonTools::ScalePhotonError(Photon* p, Double_t scale) {
 p->SetEnergyErrSmeared(scale*p->EnergyErrSmeared()); 
}

//--------------------------------------------------------------------------------------------------
Bool_t PhotonTools::PassConversionId(const Photon *p, const DecayParticle *c) {

  if (!c) return kTRUE;
  
  ThreeVector dirconvsc = ThreeVector(p->SCluster()->Point()) - c->Position();
  Double_t deta = c->Eta()-dirconvsc.Eta();
  Double_t dphi = MathUtils::DeltaPhi(c->Phi(),dirconvsc.Phi());
  Double_t eoverp = p->SCluster()->Energy()/c->P();
  
  if (p->IsEB() && eoverp>2.0) return kFALSE;
  if (p->IsEE() && eoverp>3.0) return kFALSE;
  
  if (p->IsEE() && TMath::Abs(deta)>0.01) return kFALSE;
  if (p->IsEE() && TMath::Abs(dphi)>0.01) return kFALSE;

  return kTRUE;
    
}

//--------------------------------------------------------------------------------------------------
Bool_t PhotonTools::PassElectronVeto(const Photon *p, const ElectronCol *els) {

  Bool_t pass = kTRUE;
  for (UInt_t i=0; i<els->GetEntries(); ++i) {
    const Electron *e = els->At(i);
    if (e->SCluster()==p->SCluster() && e->GsfTrk()->NExpectedHitsInner()==0) {
      pass = kFALSE;
    }
  }
  
  return pass;
}

//--------------------------------------------------------------------------------------------------
const Electron *PhotonTools::MatchedElectron(const Photon *p, const ElectronCol *els) {

  for (UInt_t i=0; i<els->GetEntries(); ++i) {
    const Electron *e = els->At(i);
    if ( e->SCluster()==p->SCluster() ) {
      return e;
    }
  }
  
  return 0;
}

//--------------------------------------------------------------------------------------------------
const Photon *PhotonTools::MatchedPhoton(const Electron *e, const PhotonCol *phs) {

  for (UInt_t i=0; i<phs->GetEntries(); ++i) {
    const Photon *p = phs->At(i);
    if ( p->SCluster()==e->SCluster() ) {
      return p;
    }
  }
  
  return 0;
}

//--------------------------------------------------------------------------------------------------
const SuperCluster *PhotonTools::MatchedSC(const SuperCluster *psc, const SuperClusterCol *scs, Double_t drMin) {

  Double_t drsmallest = 999.;
  const SuperCluster *match = 0;
  for (UInt_t i=0; i<scs->GetEntries(); ++i) {
    const SuperCluster *sc = scs->At(i);
    Double_t dr = MathUtils::DeltaR(*sc,*psc);
    if ( dr<drsmallest && dr<drMin ) {
      drsmallest = dr;
      match = sc;
    }
  }
  
  return match;
}

//--------------------------------------------------------------------------------------------------
const SuperCluster *PhotonTools::MatchedPFSC(const SuperCluster *psc, const PhotonCol *pfphos, const ElectronCol *eles, Double_t drMin) {

  
  if (pfphos) {
    for (UInt_t i=0; i<pfphos->GetEntries(); ++i) {
      const Photon *pfpho = pfphos->At(i);
      if (psc == pfpho->SCluster() && pfpho->PFSCluster()) return pfpho->PFSCluster();
    }
  }
  
  for (UInt_t i=0; i<eles->GetEntries(); ++i) {
    const Electron *ele = eles->At(i);
    if (psc == ele->SCluster() && ele->PFSCluster()) return ele->PFSCluster();
  }  

  Double_t drsmallest = 999.;
  const SuperCluster *match = 0;
  for (UInt_t i=0; i<eles->GetEntries(); ++i) {
    const Electron *ele = eles->At(i);
    //if (psc == ele->SCluster() && ele->HasPFSCluster()) return ele->PFSCluster();
    if (!ele->PFSCluster()) continue;
    
    Double_t dr = MathUtils::DeltaR(*ele->PFSCluster(),*psc);
    if ( dr<drsmallest && dr<drMin ) {
      drsmallest = dr;
      match = ele->PFSCluster();
    }    
  }  
    
  return match;
}

//--------------------------------------------------------------------------------------------------
Double_t PhotonTools::ElectronVetoCiC(const Photon *p, const ElectronCol *els) {
  
  for (UInt_t i=0; i<els->GetEntries(); ++i) {
    const Electron *e = els->At(i);
    //if (MathUtils::DeltaR(*e->SCluster(),*p->SCluster())<1) {
    if (e->SCluster()==p->SCluster() ) {
      //if( e->GsfTrk()->NExpectedHitsInner()==0 && e->GsfTrk()->Pt() > 2.5 ) {
      if( e->GsfTrk()->NExpectedHitsInner()==0 ) {
        double dEta = e->DeltaEtaSuperClusterTrackAtVtx();
        double dPhi = e->DeltaPhiSuperClusterTrackAtVtx();
        double dR = TMath::Sqrt(dEta*dEta+dPhi*dPhi);
        return dR;
      }
    }
  }  
  return 99.;
}

//--------------------------------------------------------------------------------------------------
Bool_t PhotonTools::PassElectronVetoConvRecovery(const Photon *p, const ElectronCol *els, const DecayParticleCol *conversions, const BaseVertex *v) {

  Bool_t pass = kTRUE;
  for (UInt_t i=0; i<els->GetEntries(); ++i) {
    const Electron *e = els->At(i);

    // HACVK to match CMSSW bug...
    if (e->SCluster()==p->SCluster() && e->GsfTrk()->NExpectedHitsInner()==0 && ElectronTools::PassConversionFilter(e, conversions, 
                                                                                                                    v, 0, 0., 1e-6, kTRUE, kFALSE) ) {
      //                                                         v, 0, 1e-6, 2.0, kFALSE, kFALSE) ) {
      pass = kFALSE;
    }
  }
  
  return pass;
}

//--------------------------------------------------------------------------------------------------
Bool_t PhotonTools::PassTriggerMatching(const Photon *p, const TriggerObjectCol *trigobjs)
{
  
  for (UInt_t i=0; i<trigobjs->GetEntries(); ++i) {
    const TriggerObject *trigobj = trigobjs->At(i);
    if (trigobj->TriggerType()==TriggerObject::TriggerCluster || trigobj->TriggerType()==TriggerObject::TriggerElectron || trigobj->TriggerType()==TriggerObject::TriggerPhoton) {
      if (MathUtils::DeltaR(p->SCluster(),trigobj)<0.3) {
        return kTRUE;
      }
    }
  }
  
  return kFALSE;
  
  
}

//--------------------------------------------------------------------------------------------------
const DecayParticle *PhotonTools::MatchedConversion(const Photon *p, const DecayParticleCol *conversions, 
                                               const BaseVertex *vtx, Int_t nWrongHitsMax, Double_t probMin,
                                               Double_t lxyMin, Double_t dRMin) {
  
  return MatchedConversion(p->SCluster(), conversions, vtx, nWrongHitsMax, probMin, lxyMin, dRMin);
  
}

//--------------------------------------------------------------------------------------------------
const DecayParticle *PhotonTools::MatchedConversion(const SuperCluster *sc, const DecayParticleCol *conversions, 
                                               const BaseVertex *vtx, Int_t nWrongHitsMax, Double_t probMin,
                                               Double_t lxyMin, Double_t dRMin) {
  
  const DecayParticle *match = 0;
  Double_t rhosmallest = 999.;
  for (UInt_t i=0; i<conversions->GetEntries(); ++i) {
    const DecayParticle *c = conversions->At(i);
    ThreeVector dirconvsc = ThreeVector(sc->Point()) - c->Position();
    Double_t dr = MathUtils::DeltaR(*c,dirconvsc);
    Double_t rho = c->Position().Rho();
    if (dr<dRMin && rho<rhosmallest && c->Prob()>probMin && c->LxyCorrected(vtx)>lxyMin) {
      Int_t nhb1 = dynamic_cast<const StableData*>(c->DaughterDat(0))->NHitsBeforeVtx();
      Int_t nhb2 = dynamic_cast<const StableData*>(c->DaughterDat(1))->NHitsBeforeVtx();
      if (TMath::Max(nhb1,nhb2)<=nWrongHitsMax) {
        rhosmallest = rho;
        match = c;
      }
    }
    
  }
  
  return match;
  
}

//--------------------------------------------------------------------------------------------------
const DecayParticle *PhotonTools::MatchedConversion(const Track *t, const DecayParticleCol *conversions, 
                                               const BaseVertex *vtx, Int_t nWrongHitsMax, Double_t probMin,
                                               Double_t lxyMin) {
  
  for (UInt_t i=0; i<conversions->GetEntries(); ++i) {
    const DecayParticle *c = conversions->At(i);
    if (c->Prob()>probMin && c->LxyCorrected(vtx)>lxyMin) {
      Int_t nhb1 = dynamic_cast<const StableData*>(c->DaughterDat(0))->NHitsBeforeVtx();
      Int_t nhb2 = dynamic_cast<const StableData*>(c->DaughterDat(1))->NHitsBeforeVtx();
      if (TMath::Max(nhb1,nhb2)<=nWrongHitsMax) {
        const Track *ct1 = dynamic_cast<const ChargedParticle*>(c->Daughter(0))->Trk();
        const Track *ct2 = dynamic_cast<const ChargedParticle*>(c->Daughter(1))->Trk();
        if (t==ct1 || t==ct2) return c;
      }
    }
    
  }
  
  return 0;
  
}

PhotonTools::DiphotonR9EtaCats PhotonTools::DiphotonR9EtaCat(const Photon *p1, const Photon *p2) {
  
  if (p1->IsEB() && p2->IsEB()) {
    if (p1->R9()>0.93 && p2->R9()>0.93) return kCat1;
    else return kCat2;
    
  }
  else {
    if (p1->R9()>0.93 && p2->R9()>0.93) return kCat3;
    else return kCat4; 
  }
  
}

PhotonTools::DiphotonR9EtaConversionCats PhotonTools::DiphotonR9EtaConversionCat(const Photon *p1, const Photon *p2, const DecayParticleCol *conversions, const BaseVertex *v) {
  
  const DecayParticle *conv1 = MatchedConversion(p1, conversions, v);
  const DecayParticle *conv2 = MatchedConversion(p2, conversions, v);
    
  if (p1->IsEB() && p2->IsEB()) {
    if (p1->R9()>0.93 && p2->R9()>0.93) return kNewCat1;
    else if (conv1||conv2) return kNewCat2;
    else return kNewCat3;
    
  }
  else {
    if (p1->R9()>0.93 && p2->R9()>0.93) return kNewCat4;
    else if (conv1||conv2) return kNewCat5;
    else return kNewCat6;
  }
  
}

PhotonTools::CiCBaseLineCats PhotonTools::CiCBaseLineCat(const Photon *p) {
  if( p->SCluster()->AbsEta()<1.5 ) {
    if ( p->R9() > 0.94 ) return kCiCCat1;
    else return kCiCCat2;
  } else {
    if ( p->R9() > 0.94 ) return kCiCCat3;
    else return kCiCCat4;
  }

  // shoud NEVER happen, but you never know...
  return kCiCNoCat;
}

//--------------------------------------------------------------------------------------------------
const DecayParticle *PhotonTools::MatchedCiCConversion(const Photon *p, const DecayParticleCol *conversions, 
                                                       Double_t dPhiMin,
                                                       Double_t dEtaMin,
                                                       Double_t dRMin,
                                                       bool     print,
                                                       int*                 numLegs,
                                                       int*                 convIdx) {
  
  // if there are no conversons, return
  if ( !p || !conversions)  return NULL;

  const DecayParticle *match = NULL;
  int                  matchIdx = -1;

  double minDR   = 999.;

  double phPhi = p->SCluster()->Phi();
  double phEta = p->SCluster()->Eta();
  
  if(print)
    std::cout<<"  --- conv match photon eta = "<<phEta<<"  phi = "<<phPhi<<std::endl;  
  

  for (UInt_t i=0; i<conversions->GetEntries(); ++i) {
    const DecayParticle *c = conversions->At(i);
    
    if(print)
      std::cout<< "   c "<<i+1<<"  pt = "<<c->Pt()<<"  with Nlegs = "<<c->NDaughters()<<"   "<<c->Position().X()<<"   "<<c->Position().Y()<<"   "<<c->Position().Z()<<"   "<<std::endl;
    
    if (c->Pt()   < 10. )                       continue; // is this refittedPairMomentum?
    if (c->NDaughters()==2 && c->Prob() < 1e-6) continue;

    ThreeVector dirconvsc = ThreeVector(p->SCluster()->Point()) - c->Position();      
    Double_t dR = MathUtils::DeltaR(dirconvsc,c->Mom());  
    
    if(dR < minDR) {
      minDR = dR;
      match = c;
      matchIdx = (int) i;

      if(print) {
        std::cout<<" conv "<<i+1<<" matches with dR   = "<<minDR<<std::endl;
      }
    } 
  }
  
  if( minDR < dRMin && match ) {
    if ( numLegs ) (*numLegs) = match->NDaughters();
    if ( convIdx ) (*convIdx) = matchIdx;
    if(print)
      std::cout<<"    best conversion is chosen"<<std::endl;
    return match;
  } else {
    if(print)
      std::cout<<"    NO conversion is chosen"<<std::endl;
    return NULL;
  }
}

bool PhotonTools::PassCiCSelection(const Photon* ph, const Vertex* vtx, 
                                   const TrackCol* trackCol,
                                   const ElectronCol* eleCol,
                                   const VertexCol* vtxCol,
                                   double rho, double ptmin, 
                                   bool applyEleVeto,
                                   bool print, float* kin) {

  
  // these values are taken from the H2GGlobe code... (actually from Marco/s mail)
  float cic4_allcuts_temp_sublead[] = { 
    3.8,         2.2,         1.77,        1.29,
    11.7,        3.4,         3.9,         1.84,
    3.5,         2.2,         2.3,         1.45,
    0.0106,      0.0097,      0.028,       0.027,
    0.082,       0.062,       0.065,       0.048,
    0.94,        0.36,        0.94,        0.32,
    1.,          0.062,       0.97,        0.97,
    1.5,         1.5,         1.5,         1.5 };  // the last line is PixelmatchVeto and un-used
  
  // cut on Et instead of Pt???    
  Bool_t isbarrel = ph->SCluster()->AbsEta()<1.5;
    
  // compute all relevant observables first
  double ecalIso3 = ph->EcalRecHitIsoDr03();
  double ecalIso4 = ph->EcalRecHitIsoDr04();
  double hcalIso4 = ph->HcalTowerSumEtDr04();

  unsigned int wVtxInd = 0;

  double trackIso1 = IsolationTools::CiCTrackIsolation(ph, vtx, 0.3, 0.02, 0.0, 0.0, 0.1, 1.0, trackCol);

  // track iso only
  double trackIso3 = IsolationTools::CiCTrackIsolation(ph, vtx, 0.3, 0.02, 0.0, 0.0, 0.1, 1.0, trackCol);

  // track iso worst vtx
  double trackIso2 = IsolationTools::CiCTrackIsolation(ph, vtx, 0.4, 0.02, 0.0, 0.0, 0.1, 1.0, trackCol, &wVtxInd, vtxCol);
  
  double combIso1 = ecalIso3+hcalIso4+trackIso1 - 0.17*rho;
  double combIso2 = ecalIso4+hcalIso4+trackIso2 - 0.52*rho;
  
  double tIso1 = (combIso1) *50./ph->Et();
  double tIso2 = (combIso2) *50./(ph->MomVtx(vtxCol->At(wVtxInd)->Position()).Pt());
  //double tIso2 = (combIso2) *50./ph->Et();
  double tIso3 = (trackIso3)*50./ph->Et();
  
  double covIEtaIEta  =ph->CoviEtaiEta();
  double HoE = ph->HadOverEm();
  
  double R9 = ph->R9();
  
  double dRTrack = PhotonTools::ElectronVetoCiC(ph, eleCol);

  // check which category it is ...
  int _tCat = 1;
  if ( !isbarrel ) _tCat = 3;
  if ( R9 < 0.94 ) _tCat++;
  
  if(print) {
    std::cout<<" -------------------------- "<<std::endl;

    std::cout<<" trackIso1  = "<<trackIso1<<std::endl;
    std::cout<<" trackIso2  = "<<trackIso2<<std::endl;
    std::cout<<" trackIso3  = "<<trackIso3<<std::endl;
    std::cout<<" ecalIso3   = "<<ecalIso3<<std::endl;
    std::cout<<" ecalIso4   = "<<ecalIso4<<std::endl;
    std::cout<<" hcalIso4   = "<<hcalIso4<<std::endl;

    std::cout<<" photon Et  = "<<ph->Et()<<std::endl;
    std::cout<<"        Eta = "<<ph->SCluster()->Eta()<<std::endl;
    std::cout<<"        HoE = "<<HoE<<"   ("<<cic4_allcuts_temp_sublead[_tCat-1+4*4]<<")"<<std::endl;
    std::cout<<"         R9 = "<<R9<<"   ("<<cic4_allcuts_temp_sublead[_tCat-1+5*4]<<")"<<std::endl;
    std::cout<<"         dR = "<<dRTrack<<"   ("<<cic4_allcuts_temp_sublead[_tCat-1+6*4]<<")"<<std::endl;
    std::cout<<"       iso1 = "<<tIso1<<"   ("<<cic4_allcuts_temp_sublead[_tCat-1+0*4]<<")"<<std::endl; 
    std::cout<<"       iso2 = "<<tIso2<<"   ("<<cic4_allcuts_temp_sublead[_tCat-1+1*4]<<")"<<std::endl; 
    std::cout<<"       iso3 = "<<tIso3<<"   ("<<cic4_allcuts_temp_sublead[_tCat-1+2*4]<<")"<<std::endl; 
  }

  if(kin) {
    kin[0] = tIso1;
    kin[1] = tIso2;
    kin[2] = tIso3;
    kin[3] = covIEtaIEta;
    kin[4] = HoE;
    kin[5] = R9;
    kin[6] = dRTrack;
    kin[7] = (float) ph->Pt();
    kin[8] = (float) ph->Eta();
    kin[9] = (float) ph->Phi();

    // iso quantities separate
    kin[10] = ecalIso3;
    kin[11] = ecalIso4;
    kin[12] = hcalIso4;
    kin[13] = trackIso1;
    kin[14] = trackIso2;
    kin[15] = trackIso3;

    kin[16] = (float) ph->Et();
    kin[17] = (float) ph->E();

    kin[18] = (float) ph->SCluster()->Eta();
    kin[19] = (float) ph->SCluster()->Phi();
  }

  float passCuts = 1.;

  if ( ph->Pt()     <= ptmin      ) passCuts = -1.;

  // not needed anymore, do in pre-selection...
  //if (  ph->SCluster()->AbsEta()>=2.5 || (ph->SCluster()->AbsEta()>=gkPhoEBEtaMax && ph->SCluster()->AbsEta()<=1.566)) passCuts = -1.;
  
  if(   ! (    tIso1                          < cic4_allcuts_temp_sublead[_tCat-1+0*4]
               && tIso2                       < cic4_allcuts_temp_sublead[_tCat-1+1*4]
               && tIso3                       < cic4_allcuts_temp_sublead[_tCat-1+2*4]
               && covIEtaIEta                 < cic4_allcuts_temp_sublead[_tCat-1+3*4]
               && HoE                         < cic4_allcuts_temp_sublead[_tCat-1+4*4]
               && R9                          > cic4_allcuts_temp_sublead[_tCat-1+5*4]
               && ( dRTrack > cic4_allcuts_temp_sublead[_tCat-1+6*4] || !applyEleVeto ) ) )   passCuts = -1.;
  
  if(print) std::cout<<"   ---> "<<passCuts<<std::endl;

  if(kin) {    
    kin[20] = passCuts;
    kin[21] = (float) _tCat;
  }    

  if(passCuts > 0.) return true;
  return false;
}

bool PhotonTools::PassCiCPFIsoSelection(const Photon* ph, 
                                        const Vertex* vtx, 
                                        const PFCandidateCol*    pfCol,
                                        const VertexCol*   vtxCol,
                                        double rho, double ptmin,bool dor9rescale, double p0b, double p1b,double p0e, double p1e, 
                                        std::vector<double>* kin  // store variables for debugging...
                                        ) {
  
  // these values are taken from the H2GGlobe code... (actually from Marco/s mail)
  float cic4_allcuts_temp_sublead[] = { 
    6.0,         4.7,         5.6,         3.6,
    10.0,        6.5,         5.6,         4.4,
    3.8,         2.5,         3.1,         2.2,
    0.0108,      0.0102,      0.028,       0.028,
    0.124,       0.092,       0.142,       0.063,
    0.94,        0.28,        0.94,        0.24 };  
  
  // cut on Et instead of Pt???    
  Bool_t isbarrel = ph->SCluster()->AbsEta()<1.5;
      
  
  // compute all relevant observables first
  double ecalIso3 = IsolationTools::PFGammaIsolation(ph, 0.3, 0.0, pfCol);
  double ecalIso4 = IsolationTools::PFGammaIsolation(ph, 0.4, 0.0, pfCol);

  unsigned int wVtxInd = 0;

  double trackIsoSel03 = IsolationTools::PFChargedIsolation(ph, vtx, 0.3, 0.0, pfCol);

  // track iso worst vtx
  double trackIsoWorst04 = IsolationTools::PFChargedIsolation(ph, vtx, 0.4, 0.00, pfCol, &wVtxInd, vtxCol);
  
  double combIso1 = ecalIso3+trackIsoSel03   + 2.5 - 0.09*rho;
  double combIso2 = ecalIso4+trackIsoWorst04 + 2.5 - 0.23*rho;
  
  double tIso1 = (combIso1) *50./ph->Et();
  double tIso2 = (combIso2) *50./(ph->MomVtx(vtxCol->At(wVtxInd)->Position()).Pt());
  //double tIso2 = (combIso2) *50./ph->Et();
  double tIso3 = (trackIsoSel03)*50./ph->Et();
  
  double covIEtaIEta  =ph->CoviEtaiEta();
  double HoE = ph->HadOverEm();
  
  double R9 = ph->R9();

  if(dor9rescale){
    if(isbarrel){
      R9 = p0b + p1b * R9;
    }else{
      R9 = p0e + p1e * R9;
    }
  }
  
  // check which category it is ...
  int _tCat = 1;
  if ( !isbarrel ) _tCat = 3;
  if ( R9 < 0.94 ) _tCat++;
  
  float passCuts = 1.;

  if( kin ) {
    kin->resize(0);

    kin->push_back(tIso1);
    kin->push_back(tIso2);
    kin->push_back(tIso3);
    kin->push_back(covIEtaIEta);
    kin->push_back(HoE);
    kin->push_back(R9);

    kin->push_back( (double) wVtxInd );
    kin->push_back( ecalIso3 );
    kin->push_back( ecalIso4 );

    kin->push_back( trackIsoSel03 );
    kin->push_back( trackIsoWorst04 );

    kin->push_back( combIso1 );
    kin->push_back( combIso2 );
  }


  if ( ph->Pt()     <= ptmin      ) passCuts = -1.;

  // not needed anymore, do in pre-selection...
  //if (  ph->SCluster()->AbsEta()>=2.5 || (ph->SCluster()->AbsEta()>=gkPhoEBEtaMax && ph->SCluster()->AbsEta()<=1.566)) passCuts = -1.;
  
  if(   ! (    tIso1                          < cic4_allcuts_temp_sublead[_tCat-1+0*4]
               && tIso2                       < cic4_allcuts_temp_sublead[_tCat-1+1*4]
               && tIso3                       < cic4_allcuts_temp_sublead[_tCat-1+2*4]
               && covIEtaIEta                 < cic4_allcuts_temp_sublead[_tCat-1+3*4]
               && HoE                         < cic4_allcuts_temp_sublead[_tCat-1+4*4]
               && R9                          > cic4_allcuts_temp_sublead[_tCat-1+5*4] ) )   passCuts = -1.;
  

  if(passCuts > 0.) return true;
  return false;
}

//for mono photon: cic photon id with conversion safe eleveto
bool PhotonTools::PassCiCPFIsoSelectionWithEleVeto(const Photon* ph, 
                                                   const ElectronCol *els,
                                                   const DecayParticleCol *conversions, const BaseVertex *bs,
                                                   const Vertex* vtx, 
                                                   const PFCandidateCol*    pfCol,
                                                   const VertexCol*   vtxCol,
                                                   double rho, double ptmin,
                                                   Bool_t applyElectronVeto, Bool_t invertElectronVeto,
                                                   std::vector<double>* kin  // store variables for debugging...
                                                   ) {

  Bool_t PassEleVetoRaw = PhotonTools::PassElectronVetoConvRecovery(ph, els, conversions, bs);  
  Bool_t PassEleVeto = (!applyElectronVeto && !invertElectronVeto) || (applyElectronVeto && !invertElectronVeto && PassEleVetoRaw) || (!applyElectronVeto && invertElectronVeto && !PassEleVetoRaw);
  if(!PassEleVeto){
    return false;
  }
 
  // these values are taken from the H2GGlobe code... (actually from Marco/s mail)
  float cic4_allcuts_temp_sublead[] = { 
    6.0,         4.7,         5.6,         3.6,
    10.0,        6.5,         5.6,         4.4,
    3.8,         2.5,         3.1,         2.2,
    0.0108,      0.0102,      0.028,       0.028,
    0.124,       0.092,       0.142,       0.063,
    0.94,        0.28,        0.94,        0.24 };  
  
  // cut on Et instead of Pt???    
  Bool_t isbarrel = ph->SCluster()->AbsEta()<1.5;
      
  
  // compute all relevant observables first
  double ecalIso3 = IsolationTools::PFGammaIsolation(ph, 0.3, 0.0, pfCol);
  double ecalIso4 = IsolationTools::PFGammaIsolation(ph, 0.4, 0.0, pfCol);

  unsigned int wVtxInd = 0;

  double trackIsoSel03 = IsolationTools::PFChargedIsolation(ph, vtx, 0.3, 0.0, pfCol);

  // track iso worst vtx
  double trackIsoWorst04 = IsolationTools::PFChargedIsolation(ph, vtx, 0.4, 0.00, pfCol, &wVtxInd, vtxCol);
  
  double combIso1 = ecalIso3+trackIsoSel03   + 2.5 - 0.09*rho;
  double combIso2 = ecalIso4+trackIsoWorst04 + 2.5 - 0.23*rho;
  
  double tIso1 = (combIso1) *50./ph->Et();
  double tIso2 = (combIso2) *50./(ph->MomVtx(vtxCol->At(wVtxInd)->Position()).Pt());
  //double tIso2 = (combIso2) *50./ph->Et();
  double tIso3 = (trackIsoSel03)*50./ph->Et();
  
  double covIEtaIEta  =ph->CoviEtaiEta();
  double HoE = ph->HadOverEm();
  
  double R9 = ph->R9();
  
  // check which category it is ...
  int _tCat = 1;
  if ( !isbarrel ) _tCat = 3;
  if ( R9 < 0.94 ) _tCat++;
  
  float passCuts = 1.;

  if( kin ) {
    kin->resize(0);

    kin->push_back(tIso1);
    kin->push_back(tIso2);
    kin->push_back(tIso3);
    kin->push_back(covIEtaIEta);
    kin->push_back(HoE);
    kin->push_back(R9);

    kin->push_back( (double) wVtxInd );
    kin->push_back( ecalIso3 );
    kin->push_back( ecalIso4 );

    kin->push_back( trackIsoSel03 );
    kin->push_back( trackIsoWorst04 );

    kin->push_back( combIso1 );
    kin->push_back( combIso2 );
  }


  if ( ph->Pt()     <= ptmin      ) passCuts = -1.;

  // not needed anymore, do in pre-selection...
  //if (  ph->SCluster()->AbsEta()>=2.5 || (ph->SCluster()->AbsEta()>=gkPhoEBEtaMax && ph->SCluster()->AbsEta()<=1.566)) passCuts = -1.;
  
  if(   ! (    tIso1                          < cic4_allcuts_temp_sublead[_tCat-1+0*4]
               && tIso2                       < cic4_allcuts_temp_sublead[_tCat-1+1*4]
               && tIso3                       < cic4_allcuts_temp_sublead[_tCat-1+2*4]
               && covIEtaIEta                 < cic4_allcuts_temp_sublead[_tCat-1+3*4]
               && HoE                         < cic4_allcuts_temp_sublead[_tCat-1+4*4]
               && R9                          > cic4_allcuts_temp_sublead[_tCat-1+5*4] ) )   passCuts = -1.;
  

  if(passCuts > 0.) return true;
  return false;
}


bool PhotonTools::PassEgammaMediumSelectionWithEleVeto(const Photon* ph,
                                                       const ElectronCol *els,
                                                       const DecayParticleCol *conversions, const BaseVertex *bs,
                                                       const Vertex* vtx,
                                                       const PFCandidateCol* pfCol,
                                                       const VertexCol* vtxCol,
                                                       double rho, double ptmin,
                                                       Bool_t applyElectronVeto, Bool_t invertElectronVeto,
                                                       std::vector<double>* kin // store variables for debugging...
                                                       ) {
  Bool_t PassEleVetoRaw = PhotonTools::PassElectronVetoConvRecovery(ph, els, conversions, bs);
  Bool_t PassEleVeto = (!applyElectronVeto && !invertElectronVeto) || (applyElectronVeto && !invertElectronVeto && PassEleVetoRaw) || (!applyElectronVeto && invertElectronVeto && !PassEleVetoRaw);
  
  ////debug
  //std::cout
  //<< " ph->Pt() " << ph->Pt()
  //<< " ph->Eta() " << ph->Eta()
  //<< " ph->IsEB() " << ph->IsEB()
  //<< " ph->IsEE() " << ph->IsEE()
  //<< " PassEleVeto " << PassEleVeto
  //<< std::endl;
  
  if(!PassEleVeto){
    return false;
  }
  // prepare the values of the cuts (cat 0/1 is barrel/endcap)
  float egmedium_all_cuts[] = {
    0.05, 0.05,
    0.011, 0.033,
    1.5, 1.2,
    1.0, 1.5,
    0.7, 1.0};
  // prepare the values of the EA (cat 0/1/2 is CH/NH/Photons)
  float effective_area[] = {
    0.012, 0.030, 0.148,
    0.010, 0.057, 0.130,
    0.014, 0.039, 0.112,
    0.012, 0.015, 0.216,
    0.016, 0.024, 0.262,
    0.020, 0.039, 0.260,
    0.012, 0.072, 0.266};
    
  // determine the eta category (needed for the EA)
  double absEta = ph->SCluster()->AbsEta();
  int _etaCat = 0;
  if (absEta >= 1 && absEta < 1.479) _etaCat = 1;
  if (absEta >= 1.479 && absEta < 2.0)_etaCat = 2;
  if (absEta >= 2.0 && absEta < 2.2) _etaCat = 3;
  if (absEta >= 2.2 && absEta < 2.3) _etaCat = 4;
  if (absEta >= 2.3 && absEta < 2.4) _etaCat = 5;
  if (absEta >= 2.4) _etaCat = 6;
  // determine EB/EE category
  Bool_t isbarrel = true;
  if (_etaCat >= 2) 
    isbarrel = false;
    
  // compute all relevant observables first
  double CHIso03 = IsolationTools::PFChargedIsolation(ph, vtx, 0.3, 0.0, pfCol);
  double NHIso03 = IsolationTools::PFNeutralHadronIsolation(ph, 0.3, 0.0, pfCol);
  double PHIso03 = IsolationTools::PFGammaIsolation(ph, 0.3, 0.0, pfCol);
  double combIso1 = TMath::Max(CHIso03 - rho*effective_area[0+_etaCat*3], 0.);
  double combIso2 = TMath::Max(NHIso03 - rho*effective_area[1+_etaCat*3], 0.);
  double combIso3 = TMath::Max(PHIso03 - rho*effective_area[2+_etaCat*3], 0.);
  double covIEtaIEta =ph->CoviEtaiEta();
  double HoE = ph->HadOverEmTow(); //single tower H/E
  // check which category it is ...
  int _tCat = 1;
  if ( !isbarrel ) 
    _tCat = 2;
  float passCuts = 1.;
  if (ph->Pt() <= ptmin) 
    passCuts = -1.;
  if (!ph->IsEB() && !ph->IsEE()) 
    passCuts = -1.;
  if( !(
      HoE < egmedium_all_cuts[_tCat-1+0*2]
      &&covIEtaIEta < egmedium_all_cuts[_tCat-1+1*2]
      &&combIso1 < egmedium_all_cuts[_tCat-1+2*2]
      &&combIso2 < (egmedium_all_cuts[_tCat-1+3*2] + 0.04 * ph->Pt())
      &&combIso3 < (egmedium_all_cuts[_tCat-1+4*2] + 0.005 * ph->Pt()) 
      )) 
      passCuts = -1.;

  //// debug
  //std::cout
  //<< " rho " << rho
  //<< " CHIso03 " << CHIso03
  //<< " NHIso03 " << NHIso03
  //<< " PHIso03 " << PHIso03
  //<< " combIso1 " << combIso1
  //<< " combIso2 " << combIso2
  //<< " combIso3 " << combIso3
  //<< " covIEtaIEta " << covIEtaIEta
  //<< " HoE " << HoE
  //<< " Decision " << (bool) (
      //HoE < egmedium_all_cuts[_tCat-1+0*2]
      //&&covIEtaIEta < egmedium_all_cuts[_tCat-1+1*2]
      //&&combIso1 < egmedium_all_cuts[_tCat-1+2*2]
      //&&combIso2 < (egmedium_all_cuts[_tCat-1+3*2] + 0.04 * ph->Pt())
      //&&combIso3 < (egmedium_all_cuts[_tCat-1+4*2] + 0.005 * ph->Pt()) 
      //)
  //<< std::endl;

  if(passCuts > 0.) 
    return true;
  return false;
}

const MCParticle *PhotonTools::MatchMC(const Particle *ph, const MCParticleCol *c, Bool_t matchElectrons) {

//  printf("Start loop\n");
  for (UInt_t i=0; i<c->GetEntries(); ++i) {
    const MCParticle *p = c->At(i);
//     if (p->IsGenerated() && p->AbsPdgId()==11 && (p->DistinctMother()->AbsPdgId()==23|| p->DistinctMother()->AbsPdgId()==24)) {
//       printf("pdgid = %i, status = %i, pt = %5f, mother pdgid = %i\n",p->PdgId(),p->Status(),p->Pt(),p->Mother()->PdgId());
//     }
    if (matchElectrons && p->AbsPdgId()==11 && p->IsGenerated() && p->Status()==3 && MathUtils::DeltaR(*ph,*p) < 0.3 && p->Mother() && (p->Mother()->AbsPdgId()==23 || p->Mother()->AbsPdgId()==24 || p->Mother()->AbsPdgId()==22)) {
      return p;
    }
    if ( p->AbsPdgId()==22 && p->IsGenerated() && MathUtils::DeltaR(*ph,*p) < 0.3 && p->Mother() && (p->Mother()->AbsPdgId()==25 || p->Mother()->AbsPdgId()<=21) ) {
      return p;
    }
  }
  return 0;
}

PhotonTools::DiphotonR9EtaPtCats PhotonTools::DiphotonR9EtaPtCat(const Photon *p1, const Photon *p2) {

  PhotonTools::CiCBaseLineCats cat1 = CiCBaseLineCat(p1);
  PhotonTools::CiCBaseLineCats cat2 = CiCBaseLineCat(p2);
  
  PhotonTools::DiphotonR9EtaPtCats evtcat = PhotonTools::kOctCat0;
  
  bool ph1IsEB = (cat1 ==  PhotonTools::kCiCCat1 || cat1 == PhotonTools::kCiCCat2);
  bool ph2IsEB = (cat2 ==  PhotonTools::kCiCCat1 || cat2 == PhotonTools::kCiCCat2);

  bool ph1IsHR9 = (cat1 ==  PhotonTools::kCiCCat1 || cat1 == PhotonTools::kCiCCat3);
  bool ph2IsHR9 = (cat2 ==  PhotonTools::kCiCCat1 || cat2 == PhotonTools::kCiCCat3);
  
  if( ph1IsEB && ph2IsEB )
    evtcat = ( ph1IsHR9 && ph2IsHR9 ? PhotonTools::kOctCat0 : PhotonTools::kOctCat1);
  else
    evtcat = ( ph1IsHR9 && ph2IsHR9 ? PhotonTools::kOctCat2 : PhotonTools::kOctCat3);
  
  float ptgg = (p1->Mom()+p2->Mom()).Pt();
  if (ptgg<40.0) evtcat = PhotonTools::DiphotonR9EtaPtCats(evtcat + 4);
  
  
  return evtcat;
}

Bool_t  PhotonTools::PassSinglePhotonPresel(const Photon *p,const ElectronCol *els, const DecayParticleCol *conversions, const BaseVertex *bs, const TrackCol* trackCol,const Vertex *vtx, double rho, Bool_t applyElectronVeto, Bool_t invertElectronVeto) {
  float ScEta=p->SCluster()->Eta();
  float Et=p->Et();
  float R9=p->R9();
  float HoE=p->HadOverEm();
  float CovIEtaIEta=p->CoviEtaiEta();
  float EcalIsoDr03=p->EcalRecHitIsoDr03();
  float HcalIsoDr03=p->HcalTowerSumEtDr03();
  float TrkIsoHollowDr03=p->HollowConeTrkIsoDr03();
  float NewEcalIso=EcalIsoDr03-0.012*Et;
  float NewHcalIso=HcalIsoDr03-0.005*Et;
  float NewTrkIsoHollowDr03=TrkIsoHollowDr03-0.002*Et;
  Bool_t IsBarrel=kFALSE;
  Bool_t IsEndcap=kFALSE;
  Bool_t PassEleVetoRaw = PhotonTools::PassElectronVetoConvRecovery(p, els, conversions, bs);  
  Bool_t PassEleVeto = (!applyElectronVeto && !invertElectronVeto) || (applyElectronVeto && !invertElectronVeto && PassEleVetoRaw) || (!applyElectronVeto && invertElectronVeto && !PassEleVetoRaw);
  float AbsTrackIsoCIC=IsolationTools::CiCTrackIsolation(p,vtx, 0.3, 0.02, 0.0, 0.0, 0.1, 1.0,trackCol, NULL, NULL, (!applyElectronVeto ? els : NULL) );
  float HcalEcalPUCorr=EcalIsoDr03+HcalIsoDr03-0.17*rho;

  if(fabs(ScEta)<gkPhoEBEtaMax){IsBarrel=kTRUE;}
  if(fabs(ScEta)>1.566 && fabs(ScEta)<2.5){IsEndcap=kTRUE;}
  if((!IsBarrel) && (!IsEndcap)){
    return kFALSE;
  }
  if(!PassEleVeto){
    return kFALSE;
  }
  if(R9<=0.9){
    if(HoE<0.075 && ((IsBarrel && CovIEtaIEta<0.014) || (IsEndcap && CovIEtaIEta<0.034)) && NewEcalIso<4 && NewHcalIso<4 && NewTrkIsoHollowDr03<4 && HcalEcalPUCorr<3 && AbsTrackIsoCIC<2.8 && TrkIsoHollowDr03<4.0) {
      return kTRUE;
    }
  }
  if(R9>0.9){
    if(((IsBarrel && HoE<0.082 && CovIEtaIEta<0.014) || (IsEndcap && HoE <0.075 && CovIEtaIEta<0.034)) && NewEcalIso<50 && NewHcalIso<50 && NewTrkIsoHollowDr03<50 && HcalEcalPUCorr<3 && AbsTrackIsoCIC<2.8 && TrkIsoHollowDr03<4.0) {
      return kTRUE;  
    }
  }
  return kFALSE;
}


Bool_t  PhotonTools::PassSinglePhotonPreselPFISO_NoTrigger(const Photon *p,const ElectronCol *els, const DecayParticleCol *conversions, const BaseVertex *bs, const TrackCol* trackCol,const Vertex *vtx, double rho, const PFCandidateCol *fPFCands, Bool_t applyElectronVeto, Bool_t invertElectronVeto) {
  
  float ScEta=p->SCluster()->Eta();
  float CovIEtaIEta=p->CoviEtaiEta();

  Bool_t IsBarrel=kFALSE;
  Bool_t IsEndcap=kFALSE;
  Bool_t PassEleVetoRaw = PhotonTools::PassElectronVetoConvRecovery(p, els, conversions, bs);  
  Bool_t PassEleVeto = (!applyElectronVeto && !invertElectronVeto) || (applyElectronVeto && !invertElectronVeto && PassEleVetoRaw) || (!applyElectronVeto && invertElectronVeto && !PassEleVetoRaw);
  float ChargedIso_selvtx_DR002To0p02=IsolationTools::PFChargedIsolation(p,vtx, 0.2, 0.,fPFCands);
  if(fabs(ScEta)<gkPhoEBEtaMax){IsBarrel=kTRUE;}
  if(fabs(ScEta)>1.566 && fabs(ScEta)<2.5){IsEndcap=kTRUE;}
  if((!IsBarrel) && (!IsEndcap)){
    return kFALSE;
  }
  if(!PassEleVeto){
    return kFALSE;
  }
  if( ((IsBarrel && CovIEtaIEta<0.02) || (IsEndcap && CovIEtaIEta<0.06)) && ChargedIso_selvtx_DR002To0p02<4) {
    return kTRUE;
  }

  return kFALSE;
}



Bool_t  PhotonTools::PassSinglePhotonPreselPFISO(const Photon *p,const ElectronCol *els, const DecayParticleCol *conversions, const BaseVertex *bs, const TrackCol* trackCol,const Vertex *vtx, double rho, const PFCandidateCol *fPFCands, Bool_t applyElectronVeto, Bool_t invertElectronVeto) {


  float ScEta=p->SCluster()->Eta();
  float Et=p->Et();
  float R9=p->R9();
  float HoE=p->HadOverEm();
  float CovIEtaIEta=p->CoviEtaiEta();
  float EcalIsoDr03=p->EcalRecHitIsoDr03();
  float HcalIsoDr03=p->HcalTowerSumEtDr03();
  float TrkIsoHollowDr03=p->HollowConeTrkIsoDr03();

  float NewEcalIso=EcalIsoDr03-0.012*Et;
  float NewHcalIso=HcalIsoDr03-0.005*Et;
  float NewTrkIsoHollowDr03=TrkIsoHollowDr03-0.002*Et;


  Bool_t IsBarrel=kFALSE;
  Bool_t IsEndcap=kFALSE;
  Bool_t PassEleVetoRaw = PhotonTools::PassElectronVetoConvRecovery(p, els, conversions, bs);  
  Bool_t PassEleVeto = (!applyElectronVeto && !invertElectronVeto) || (applyElectronVeto && !invertElectronVeto && PassEleVetoRaw) || (!applyElectronVeto && invertElectronVeto && !PassEleVetoRaw);
  float ChargedIso_selvtx_DR002To0p02=IsolationTools::PFChargedIsolation(p,vtx, 0.2, 0.,fPFCands);
  if(fabs(ScEta)<gkPhoEBEtaMax){IsBarrel=kTRUE;}
  if(fabs(ScEta)>1.566 && fabs(ScEta)<2.5){IsEndcap=kTRUE;}
  if((!IsBarrel) && (!IsEndcap)){
    return kFALSE;
  }
  if(!PassEleVeto){
    return kFALSE;
  }
  if(R9<=0.9){
    if(HoE<0.075 && ((IsBarrel && CovIEtaIEta<0.014) || (IsEndcap && CovIEtaIEta<0.034)) && NewEcalIso<4 && NewHcalIso<4 && NewTrkIsoHollowDr03<4 && ChargedIso_selvtx_DR002To0p02<4) {
      return kTRUE;
    }
  }
  if(R9>0.9){
    if(((IsBarrel && HoE<0.082 && CovIEtaIEta<0.014) || (IsEndcap && HoE <0.075 && CovIEtaIEta<0.034)) && NewEcalIso<50 && NewHcalIso<50 && NewTrkIsoHollowDr03<50 && ChargedIso_selvtx_DR002To0p02<4) {
      return kTRUE;  
    }
  }
  return kFALSE;
}

bool PhotonTools::PassVgamma2011Selection(const Photon* ph, double rho) {

  bool isEB = (ph->SCluster()->AbsEta()<1.5);
  
  if (ph->HasPixelSeed())                                           return false;  // ? is this what we want ?
  if (ph->HadOverEm()            > 0.05)                            return false;
  if (ph->CoviEtaiEta()          > (isEB ? 0.011 : 0.03) )          return false;
  if (ph->HollowConeTrkIsoDr04() > ( 2.0 + 0.001 *ph->Pt() + 
                                     (isEB ? 0.0167 : 0.032)*rho )) return false; 
  if (ph->EcalRecHitIsoDr04()    > ( 4.2 + 0.006 *ph->Pt() + 
                                     (isEB ? 0.183  : 0.090)*rho )) return false; 
  if (ph->HcalTowerSumEtDr04()   > ( 2.2 + 0.0025*ph->Pt() + 
                                     (isEB ? 0.062  : 0.180)*rho )) return false; 

  // spike cleaning...
  if ( ph->CoviEtaiEta() < 0.001 && TMath::Sqrt(TMath::Abs(ph->SCluster()->Seed()->CoviPhiiPhi())) < 0.001)
    return false;

  return true;
}


Bool_t  PhotonTools::PassSinglePhotonPreselPFISONoEcal(const Photon *p,const ElectronCol *els, const DecayParticleCol *conversions, const BaseVertex *bs, const TrackCol* trackCol,const Vertex *vtx, double rho, const PFCandidateCol *fPFCands, Bool_t applyElectronVeto, Bool_t invertElectronVeto) {


  float ScEta=p->SCluster()->Eta();
  float Et=p->Et();
  float R9=p->R9();
  float HoE=p->HadOverEm();
  float CovIEtaIEta=p->CoviEtaiEta();
  float HcalIsoDr03=p->HcalTowerSumEtDr03();
  float TrkIsoHollowDr03=p->HollowConeTrkIsoDr03();

  float NewHcalIso=HcalIsoDr03-0.005*Et;
  float NewTrkIsoHollowDr03=TrkIsoHollowDr03-0.002*Et;


  Bool_t IsBarrel=kFALSE;
  Bool_t IsEndcap=kFALSE;
  Bool_t PassEleVetoRaw = PhotonTools::PassElectronVetoConvRecovery(p, els, conversions, bs);  
  Bool_t PassEleVeto = (!applyElectronVeto && !invertElectronVeto) || (applyElectronVeto && !invertElectronVeto && PassEleVetoRaw) || (!applyElectronVeto && invertElectronVeto && !PassEleVetoRaw);
  float ChargedIso_selvtx_DR002To0p02=IsolationTools::PFChargedIsolation(p,vtx, 0.2, 0.,fPFCands);
  if(fabs(ScEta)<gkPhoEBEtaMax){IsBarrel=kTRUE;}
  if(fabs(ScEta)>1.566 && fabs(ScEta)<2.5){IsEndcap=kTRUE;}
  if((!IsBarrel) && (!IsEndcap)){
    return kFALSE;
  }
  if(!PassEleVeto){
    return kFALSE;
  }
  if(R9<=0.9){
    if(HoE<0.075 && ((IsBarrel && CovIEtaIEta<0.014) || (IsEndcap && CovIEtaIEta<0.034)) && NewHcalIso<4 && NewTrkIsoHollowDr03<4 && ChargedIso_selvtx_DR002To0p02<4) {
      return kTRUE;
    }
  }
  if(R9>0.9){
    if(((IsBarrel && HoE<0.082 && CovIEtaIEta<0.014) || (IsEndcap && HoE <0.075 && CovIEtaIEta<0.034)) && NewHcalIso<50 && NewTrkIsoHollowDr03<50 && ChargedIso_selvtx_DR002To0p02<4) {
      return kTRUE;  
    }
  }
  return kFALSE;
}

Bool_t  PhotonTools::PassSinglePhotonPreselPFISONoEcalNoPFChargedIso(const Photon *p,const ElectronCol *els, const DecayParticleCol *conversions, const BaseVertex *bs, const TrackCol* trackCol,const Vertex *vtx, double rho, const PFCandidateCol *fPFCands, Bool_t applyElectronVeto, Bool_t invertElectronVeto) {


  float ScEta=p->SCluster()->Eta();
  float Et=p->Et();
  float R9=p->R9();
  float HoE=p->HadOverEm();
  float CovIEtaIEta=p->CoviEtaiEta();
  float HcalIsoDr03=p->HcalTowerSumEtDr03();
  float TrkIsoHollowDr03=p->HollowConeTrkIsoDr03();

  float NewHcalIso=HcalIsoDr03-0.005*Et;
  float NewTrkIsoHollowDr03=TrkIsoHollowDr03-0.002*Et;


  Bool_t IsBarrel=kFALSE;
  Bool_t IsEndcap=kFALSE;
  Bool_t PassEleVetoRaw = PhotonTools::PassElectronVetoConvRecovery(p, els, conversions, bs);  
  Bool_t PassEleVeto = (!applyElectronVeto && !invertElectronVeto) || (applyElectronVeto && !invertElectronVeto && PassEleVetoRaw) || (!applyElectronVeto && invertElectronVeto && !PassEleVetoRaw);
 
  if(fabs(ScEta)<gkPhoEBEtaMax){IsBarrel=kTRUE;}
  if(fabs(ScEta)>1.566 && fabs(ScEta)<2.5){IsEndcap=kTRUE;}
  if((!IsBarrel) && (!IsEndcap)){
    return kFALSE;
  }
  if(!PassEleVeto){
    return kFALSE;
  }
  if(R9<=0.9){
    if(HoE<0.075 && ((IsBarrel && CovIEtaIEta<0.014) || (IsEndcap && CovIEtaIEta<0.034)) && NewHcalIso<4 && NewTrkIsoHollowDr03<4) {
      return kTRUE;
    }
  }
  if(R9>0.9){
    if(((IsBarrel && HoE<0.082 && CovIEtaIEta<0.014) || (IsEndcap && HoE <0.075 && CovIEtaIEta<0.034)) && NewHcalIso<50 && NewTrkIsoHollowDr03<50) {
      return kTRUE;  
    }
  }
  return kFALSE;
}

void PhotonTools::ScalePhotonShowerShapes(Photon* p, PhotonTools::ShowerShapeScales scale) {

  switch(scale) {
    
  case kNoShowerShapeScaling:
    break;

  case k2011ShowerShape:
    //regression sigmaE
    if (p->SCluster()->AbsEta()<1.5)
      PhotonTools::ScalePhotonError(p,1.07);
    else
      PhotonTools::ScalePhotonError(p,1.045);
    
    //R9
    if (p->SCluster()->AbsEta()<1.5) p->SetR9(1.0035*p->R9());
    else  p->SetR9(1.0035*p->R9());
    //CoviEtaiEta(SigiEtaiEta)
    if (p->SCluster()->AbsEta()<1.5) p->SetCoviEtaiEta(0.87*p->CoviEtaiEta()+0.0011);
    else  p->SetCoviEtaiEta(0.99*p->CoviEtaiEta());
    //EtaWidth
    if (p->SCluster()->AbsEta()<1.5) p->SetEtaWidth(0.99*p->EtaWidth());
    else  p->SetEtaWidth(0.99*p->EtaWidth());
    //PhiWidth
    if (p->SCluster()->AbsEta()<1.5) p->SetPhiWidth(0.99*p->PhiWidth());
    else  p->SetPhiWidth(0.99*p->PhiWidth());
    break;
    
  case k2012ShowerShape:        
    //R9
    if (p->SCluster()->AbsEta()<1.5) p->SetR9(1.0045*p->R9()+0.001);
    else  p->SetR9(1.0086*p->R9()-0.0007);
    //CoviEtaiEta(SigiEtaiEta)
    if (p->SCluster()->AbsEta()<1.5) p->SetCoviEtaiEta(0.891832*p->CoviEtaiEta()+0.0009133);
    else  p->SetCoviEtaiEta(0.9947*p->CoviEtaiEta()+0.00003);
    //EtaWidth
    if (p->SCluster()->AbsEta()<1.5) p->SetEtaWidth(1.04302*p->EtaWidth()- 0.000618);
    else  p->SetEtaWidth(0.903254*p->EtaWidth()+ 0.001346);
    //PhiWidth
    if (p->SCluster()->AbsEta()<1.5) p->SetPhiWidth(1.00002*p->PhiWidth()- 0.000371);
    else  p->SetPhiWidth(0.99992*p->PhiWidth()+ 4.8e-07);
    //s4Ratio
    if (p->SCluster()->AbsEta()<1.5) p->SetS4Ratio(1.01894*p->S4Ratio()-0.01034);
    else  p->SetS4Ratio(1.04969*p->S4Ratio()-0.03642);
    break;
  default:
    break;
    
  }
  return;
}
    
