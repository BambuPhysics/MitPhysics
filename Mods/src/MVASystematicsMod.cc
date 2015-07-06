#include <iostream>
#include <sstream>
#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TNtuple.h>
#include <TRandom3.h>
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/DataTree/interface/Names.h"
#include "MitAna/DataTree/interface/BeamSpotCol.h"
#include "MitAna/DataTree/interface/TrackCol.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/DataTree/interface/PhotonCol.h"
#include "MitPhysics/Utils/interface/IsolationTools.h"
#include "MitPhysics/Utils/interface/PhotonTools.h"
#include "MitPhysics/Mods/interface/MVASystematicsMod.h"


using namespace mithep;

ClassImp(mithep::MVASystematicsMod)

//--------------------------------------------------------------------------------------------------
MVASystematicsMod::MVASystematicsMod(const char *name, const char *title) :
  BaseMod        (name,title),
  // define all the Branches to load
  fMCParticleName(Names::gkMCPartBrn),
  fPVName        (Names::gkPVBeamSpotBrn),
  fEBSCName      (Names::gkBarrelSuperClusterBrn),
  fEESCName      (Names::gkEndcapSuperClusterBrn),  
  // -----------------------------------
  fMCR9ScaleEB   (1.0),
  fMCR9ScaleEE   (1.0),
  fIsData        (false),
  fTupleName     ("hMVAtuple")
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void MVASystematicsMod::Begin()
{
  // Run startup code on the client machine. For this module, we dont do anything here.
}

//--------------------------------------------------------------------------------------------------
void MVASystematicsMod::Process()
{  
  IncNEventsProcessed();
  MCParticleCol* mcParticles = 0;
  if (!fIsData)
    mcParticles = GetObject<MCParticleCol>(fMCParticleName);

  auto* pv = GetObject<VertexCol>(fPVName);
  auto* ebSC = GetObject<SuperClusterCol>(fEBSCName);
  auto* eeSC = GetObject<SuperClusterCol>(fEESCName);
  
  const Vertex *vtx = 0;
  if (pv->GetEntries()>0) vtx = pv->At(0);
  
  const MCParticle *h = 0;
  const MCParticle *p1 = 0;
  const MCParticle *p2 = 0;
  
  const SuperCluster *sc1 = 0;
  const SuperCluster *sc2 = 0;
  
  Float_t _pth    = -100.;
  Float_t _y      = -100.;
  Float_t _genmass = -100.;
  if( !fIsData ) h = FindHiggsPtAndY(_pth, _y, _genmass, mcParticles);

  if (h && h->NDaughters()>=2) {
    p1 = h->Daughter(0);
    p2 = h->Daughter(1);
  }
  
  Float_t pt1 = -100;
  Float_t eta1 = -100;
  Float_t phi1 = -100;

  Float_t pt2 = -100;
  Float_t eta2 = -100;
  Float_t phi2 = -100;  
  
  Float_t scet1 = -100;
  Float_t sceta1 = -100;
  Float_t scphi1 = -100;
  Float_t scr91 = -100;
  bool iseb1 = kFALSE;

  Float_t scet2 = -100;
  Float_t sceta2 = -100;
  Float_t scphi2 = -100;  
  Float_t scr92 = -100;  
  bool iseb2 = kFALSE;
  
  if (p1) {
    pt1 = p1->Pt();
    eta1 = p1->Eta();
    phi1 = p1->Phi();
    sc1 = MatchSC(p1,iseb1, ebSC, eeSC);
  }

  if (p2) {
    pt2 = p2->Pt();
    eta2 = p2->Eta();
    phi2 = p2->Phi();
    sc2 = MatchSC(p2,iseb2, ebSC, eeSC);
  }
  
  if (sc1) {
    double r9scale;
    if   (iseb1)
      r9scale = fMCR9ScaleEB;
    else
      r9scale = fMCR9ScaleEE;
      
    // compute the probe momentum wr to the chosen Vtx
    FourVectorM scmom;  
    ThreeVectorC scp;
    if (vtx)
      scp = sc1->Point() - vtx->Position();
    else
      scp = sc1->Point();
    scp = scp/scp.R();
    scmom.SetXYZT(sc1->Energy()*scp.X(), sc1->Energy()*scp.Y(), sc1->Energy()*scp.Z(), sc1->Energy());
      
    scet1 = scmom.Pt();
    sceta1 = sc1->Eta();
    scphi1 = sc1->Phi();
    scr91 = r9scale*sc1->R9();
    
  }


  if (sc2) {
    double r9scale;
    if (iseb2)
      r9scale = fMCR9ScaleEB;
    else
      r9scale = fMCR9ScaleEE;
      
    // compute the probe momentum wr to the chosen Vtx
    FourVectorM scmom;  
    ThreeVectorC scp;
    if (vtx)
      scp = sc2->Point() - vtx->Position();
    else
      scp = sc2->Point();
    scp = scp/scp.R();
    scmom.SetXYZT(sc2->Energy()*scp.X(), sc2->Energy()*scp.Y(), sc2->Energy()*scp.Z(), sc2->Energy());
      
    scet2 = scmom.Pt();
    sceta2 = sc2->Eta();
    scphi2 = sc2->Phi();
    scr92 = r9scale*sc2->R9();
    
  }

  Float_t fill[] = {_pth, _y, _genmass,pt1,eta1,phi1,pt2,eta2,phi2,scet1,sceta1,scphi1,scr91,scet2,sceta2,scphi2,scr92};

  hMVAtuple->Fill(fill);


  return;
}

//--------------------------------------------------------------------------------------------------
void MVASystematicsMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we typically initialize histograms and other analysis objects and request
  // branches or objects created by earlier modules.

  hMVAtuple = new TNtuple(fTupleName.Data(),fTupleName.Data(),
			  "hpt:hy:hm:pt1:eta1:phi1:pt2:eta2:phi2:scet1:sceta1:scphi1:scr91:scet2:sceta2:scphi2:scr92");

  AddOutput(hMVAtuple);
}

//--------------------------------------------------------------------------------------------------
void MVASystematicsMod::SlaveTerminate()
{
  // Run finishing code on the computer (slave) that did the analysis. For this module, we dont do
  // anything here.
}

//--------------------------------------------------------------------------------------------------
void MVASystematicsMod::Terminate()
{
  // Run finishing code on the client computer. For this module, we dont do anything here.
}

//--------------------------------------------------------------------------------------------------
const MCParticle *MVASystematicsMod::FindHiggsPtAndY(Float_t& pt, Float_t& Y, Float_t& mass, MCParticleCol const* mcParticles)
{
  const MCParticle *h = 0;
  
  pt = -999.;
  Y  = -999;
  mass = -999.;

  // loop over all GEN particles and look for status 1 photons
  for (UInt_t i=0; i<mcParticles->GetEntries(); ++i) {
    const MCParticle* p = mcParticles->At(i);
    if (p->Is(MCParticle::kH)) {
      pt=p->Pt();
      Y = p->Rapidity();
      mass = p->Mass();
      h = p;
      break;
    }
  }
  
  return h;
}

//--------------------------------------------------------------------------------------------------
const SuperCluster *MVASystematicsMod::MatchSC(const MCParticle *p, bool &iseb, SuperClusterCol const* ebSC, SuperClusterCol const* eeSC)
{

  iseb = kFALSE;
  
  for (UInt_t i=0; i<ebSC->GetEntries(); ++i) {
    iseb = kTRUE;
    const SuperCluster *sc = ebSC->At(i);
    if (MathUtils::DeltaR(sc,p)<0.2) return sc;
  }

  for (UInt_t i=0; i<eeSC->GetEntries(); ++i) {
    iseb = kFALSE;
    const SuperCluster *sc = eeSC->At(i);
    if (MathUtils::DeltaR(sc,p)<0.2) return sc;
  }

  return 0;
}
