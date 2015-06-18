//--------------------------------------------------------------------------------------------------
// M.Yang 2011/10/12
//
// PhotonMvaMod
//
// Precompute regression energy corrections and id bdt output to save memory and cpu time in
// subsequent module chains.
//
// Authors: J.Bendavid
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_PHOTONMVAMOD_H
#define MITPHYSICS_MODS_PHOTONMVAMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/PhotonFwd.h"
#include "MitAna/DataTree/interface/TrackCol.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/DataTree/interface/BeamSpotCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityCol.h"
#include "MitAna/DataTree/interface/DecayParticleCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/DecayParticleCol.h"
#include "MitAna/DataTree/interface/PileupInfoCol.h"
#include "MitAna/DataTree/interface/MCParticleCol.h"

#include "MitPhysics/Utils/interface/PhotonTools.h"
#include "MitPhysics/Utils/interface/EGEnergyCorrector.h"

#include "MitPhysics/Utils/interface/MVATools.h"
#include "MitPhysics/Utils/interface/VertexTools.h"

class TNtuple;
class TRandom3;
class TH1D;

namespace mithep 
{  
  class PhotonMvaMod : public BaseMod
  {
  public:
    PhotonMvaMod(const char *name ="PhotonMvaMod", 
		 const char *title="Selecting PhotonPairs");
    
    ~PhotonMvaMod();

    // setting all the input Names
    void                SetInputPhotonsName(const char *n){ fPhotonBranchName= n;        }
    void                SetPhotonsFromBranch(bool b)      { fPhotonsFromBranch = b;      }
    void                SetTrackName(const char *n)       { fTrackBranchName = n;        }
    void                SetElectronName(const char *n)    { fElectronName = n;           }
    void                SetConversionName(const char *n)  { fConversionName = n;         }
    void                SetPUDensityName(const char *n)   { fPileUpDenName = n;          }
    void                SetPVName(const char *n)          { fPVName = n;                 }
    void                SetPVFromBranch(bool b)           { fPVFromBranch = b;           }
    void                SetMCParticle(const char *n)      { fMCParticleName = n;         }
    void                SetPUInfoName(const char *n)      { fPileUpName = n;             }
    void                SetBeamspotName(const char *n)    { fBeamspotName = n;           }
    void                SetPFCandName(const char *n)      { fPFCandName = n;             }


    // get/set the Names for the output Photon Collection
    const char         *GetOutputName()             const { return fGoodPhotonsName;     }   
    void                SetOutputName(const char *n)      { fGoodPhotonsName=n;          }    

    // set basic Cut variables (FOR PRE-SELECTION)
    void                SetPtMin(Double_t pt)             { fPhotonPtMin     = pt;       }
    void                SetAbsEtaMax(Double_t eta)        { fPhotonEtaMax    = eta;	 }

    // is Data Or Not?
    void                SetIsData (Bool_t b) { fIsData = b;};
    
    void                SetApplyShowerRescaling(Bool_t b) { fApplyShowerRescaling = b; }

    void                ApplyEleVeto(bool a)            { fApplyEleVeto  = a; }

    void                SetDoRegression(Bool_t b)         { fDoRegression = b; }

    void                SetRegressionVersion(UInt_t v)     { fRegressionVersion = v; }
    void                SetRegressionWeights(TString f)    { fRegWeights = f; }


    void                SetMinNumPhotons ( UInt_t i ) { fMinNumPhotons = i;}
    void                SetDoPreselection( Bool_t b ) { fDoPreselection = b;}
    
  protected:
    void                Process();
    void                SlaveBegin();

    // private auxiliary methods...
    // Names for the input Collections
    TString                       fPhotonBranchName;
    TString                       fElectronName;
    TString                       fGoodElectronName;
    TString                       fConversionName;
    TString                       fTrackBranchName;
    TString                       fPileUpDenName;    
    TString                       fPVName;
    TString                       fBeamspotName;
    TString                       fPFCandName;
    TString                       fMCParticleName;
    TString                       fPileUpName;

    TString                       fGoodPhotonsName;      // name of exported "good photon" collection
    

    // Basic Pre-Selection kinematics
    Double_t                      fPhotonPtMin;          // min pt cut fro PRE-SELECTION!
    Double_t                      fPhotonEtaMax;         // max eta cut for PRE-SELECTION!
    
    // is it Data or MC?
    Bool_t                        fIsData;
    Bool_t                        fApplyShowerRescaling;
    
    // in case there's some PV pre-selection
    Bool_t                        fPhotonsFromBranch;
    Bool_t                        fPVFromBranch;
    Bool_t                        fGoodElectronsFromBranch;

    const PhotonCol              *fPhotons;
    const ElectronCol            *fElectrons;
    const ElectronCol            *fGoodElectrons;    
    const DecayParticleCol       *fConversions;
    const TrackCol               *fTracks;
    const PileupEnergyDensityCol *fPileUpDen;
    const VertexCol              *fPV;
    const BeamSpotCol            *fBeamspot;
    const PFCandidateCol         *fPFCands;
    const MCParticleCol          *fMCParticles;
    const PileupInfoCol          *fPileUp;    

    // EGEnergyCorrector egcor;
    Bool_t                        fDoRegression;
    TString                       fPhFixString;
    TString                       fPhFixFile;
    TString                       fRegWeights; 
  
    const TH1D                   *fEtaCorrections;
    
    // --------------------------------
    bool                          fApplyEleVeto;
    UInt_t                        fRegressionVersion;
    UInt_t                        fMinNumPhotons;
    Bool_t                        fDoPreselection;

    ClassDef(PhotonMvaMod, 1) // Photon identification module
  };
}
#endif
