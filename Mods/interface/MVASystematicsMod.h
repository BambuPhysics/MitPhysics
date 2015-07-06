//--------------------------------------------------------------------------------------------------
// MVASystematicsMod
//
// This module compues photon eff from Z->mumugamma
//
// Authors: F,.Stoeckli
//--------------------------------------------------------------------------------------------------

#ifndef MITHGG_MODS_MVASYSTEMATICSMOD_H
#define MITHGG_MODS_MVASYSTEMATICSMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/CollectionsFwd.h"
#include "MitAna/DataTree/interface/PileupInfoCol.h"
#include "MitAna/DataCont/interface/Types.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityCol.h"
#include "MitAna/DataTree/interface/MCParticleCol.h"

class TH1D;
class TH2D;
class TNtuple;
class MCEventInfo;
class PileupEnergyDensityCol;
class PFCandidateCol;
class Vertex;
class BaseVertex;
class BeamSpotCol;

class TRandom3;

namespace mithep 
{
  class MVASystematicsMod : public BaseMod
  {
  public:
    MVASystematicsMod(const char *name  = "MVASystematicsMod", 
		const char *title = "Photon Efficiency Analysis");


    // setting all the input Names
    void                SetIsData(bool b) { fIsData = b; }
    void                SetTupleName(const char* c)    { fTupleName = c; }
    const MCParticle   *FindHiggsPtAndY(Float_t& pt, Float_t& Y, Float_t& mass, MCParticleCol const*);
    const SuperCluster *MatchSC(const MCParticle *p, bool &iseb, SuperClusterCol const*, SuperClusterCol const*);
    
    void                SetMCR9Scale(Double_t ebscale, Double_t eescale) { fMCR9ScaleEB = ebscale; fMCR9ScaleEE = eescale; }
    

  protected:

    void                     Begin();
    void                     Process();
    void                     SlaveBegin();
    void                     SlaveTerminate();
    void                     Terminate();

    // Names for the input Collections
    TString             fMCParticleName;
    TString             fPVName;
    TString             fEBSCName;
    TString             fEESCName;
    
   // is it Data or MC?
    Double_t fMCR9ScaleEB;
    Double_t fMCR9ScaleEE;
   
    Bool_t              fIsData;
    
    TString fTupleName;

    // The output Ntuple
    TNtuple* hMVAtuple;
    ClassDef(MVASystematicsMod, 1)
  };
}
#endif
