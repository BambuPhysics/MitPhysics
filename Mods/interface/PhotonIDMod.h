//--------------------------------------------------------------------------------------------------
// PhotonIDMod
//
// This module applies photon identification criteria and exports a pointer to a collection
// of "good photons" according to the specified identification scheme.
//
// Authors: S.Xie, C.Loizides
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_PHOTONIDMOD_H
#define MITPHYSICS_MODS_PHOTONIDMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/PhotonFwd.h"
#include "MitAna/DataTree/interface/TrackCol.h"
#include "MitAna/DataTree/interface/BeamSpotCol.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityCol.h"
#include "MitAna/DataTree/interface/DecayParticleCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/DataTree/interface/DecayParticleCol.h"
#include "MitAna/DataTree/interface/PileupInfoCol.h"
#include "MitAna/DataTree/interface/MCParticleCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"

#include "MitPhysics/Utils/interface/MVATools.h"
#include "MitPhysics/Utils/interface/PhotonTools.h"
#include "MitPhysics/Utils/interface/RhoUtilities.h"

#include "TFormula.h"

class TRandom3;

namespace mithep 
{
  class PhotonIDMod : public BaseMod
  {
  public:
    PhotonIDMod(const char *name="PhotonIDMod", 
                const char *title="Photon identification module");
    ~PhotonIDMod();

    Bool_t              GetApplySpikeRemoval()      const { return fApplySpikeRemoval;   }
    Bool_t              GetApplyPixelSeed()         const { return fApplyPixelSeed;      }
    const char         *GetGoodName()               const { return GetGoodPhotonsName(); }   
    const char         *GetGoodPhotonsName()        const { return fGoodPhotonsName;     }   
    Double_t            GetHadOverEmMax()           const { return fHadOverEmMax;        }
    const char         *GetIDType()                 const { return fPhotonIDType;        }
    const char         *GetInputName()              const { return fPhotonBranchName;    }   
    const char         *GetIsoType()                const { return fPhotonIsoType;       }
    const char         *GetOutputName()             const { return GetGoodPhotonsName(); }   
    Double_t            GetPtMin()                  const { return fPhotonPtMin;         }
    Bool_t              GetApplyFiduciality()       const { return fFiduciality;         }
    Double_t            GetEtaWidthEB()	      const { return fEtaWidthEB;	   }
    Double_t            GetEtaWidthEE()	      const { return fEtaWidthEE;	   }
    Double_t            GetAbsEtaMax()	      const { return fAbsEtaMax;	   }
    void                SetApplySpikeRemoval(Bool_t b)    { fApplySpikeRemoval  = b;     }
    void                SetApplyPixelSeed(Bool_t b)       { fApplyPixelSeed  = b;        }
    void                SetApplyElectronVeto(Bool_t b)    { fApplyElectronVeto = b;      }
    void                SetInvertElectronVeto(Bool_t b)   { fInvertElectronVeto = b;     }      
    void                SetApplyElectronVetoConvRecovery(Bool_t b) { fApplyElectronVetoConvRecovery = b; }
    void                SetApplyConversionId(Bool_t b)    { fApplyConversionId = b;      }
    void                SetApplyTriggerMatching(Bool_t b)      { fApplyTriggerMatching = b;  }      
    void                SetGoodName(const char *n)        { SetGoodPhotonsName(n);       }   
    void                SetGoodPhotonsName(const char *n) { fGoodPhotonsName = n;        }   
    void                SetHadOverEmMax(Double_t hoe)     { fHadOverEmMax    = hoe;      }
    void                SetIDType(const char *type)       { fPhotonIDType    = type;     }
    void                SetInputName(const char *n)       { fPhotonBranchName= n;        }   
    void                SetTrackName(const char *n)       { fTrackBranchName = n;        }   
    void                SetBeamspotName(const char *n)    { fBeamspotBranchName = n;     }   
    void                SetIsoType(const char *type)      { fPhotonIsoType   = type;     }
    void                SetOutputName(const char *n)      { SetGoodPhotonsName(n);       }    
    void                SetPtMin(Double_t pt)             { fPhotonPtMin     = pt;       }
    void                SetR9Min(Double_t x)              { fPhotonR9Min     = x;        }
    void                SetEtaWidthEB(Double_t x)	    { fEtaWidthEB      = x;	   }
    void                SetEtaWidthEE(Double_t x)         { fEtaWidthEE      = x;	   }
    void                SetAbsEtaMax(Double_t x)          { fAbsEtaMax       = x;	   }
    void                SetAbsEtaMin(Double_t x)          { fAbsEtaMin       = x;	   }
    void                SetApplyR9Min(Bool_t b)           { fApplyR9Min      = b;        }
    void                SetApplyFiduciality(Bool_t b)     { fFiduciality = b;            }      
    void                SetEffAreas(Double_t ecal, Double_t hcal, Double_t track) { 
      fEffAreaEcalEE = ecal; fEffAreaHcalEE = hcal; fEffAreaTrackEE = track;
      fEffAreaEcalEB = ecal; fEffAreaHcalEB = hcal; fEffAreaTrackEB = track;
    }
    void                SetEffAreasEEEB(Double_t ecalEE, Double_t hcalEE, Double_t trackEE,
                                        Double_t ecalEB, Double_t hcalEB, Double_t trackEB) { 
      fEffAreaEcalEE = ecalEE; fEffAreaHcalEE = hcalEE; fEffAreaTrackEE = trackEE;
      fEffAreaEcalEB = ecalEB; fEffAreaHcalEB = hcalEB; fEffAreaTrackEB = trackEB;
    }
    void                SetChargedHadronIsoMaxEB(const char* expr) { fIsoCalcCHEB.SetTitle(expr); }
    void                SetChargedHadronIsoMaxEE(const char* expr) { fIsoCalcCHEE.SetTitle(expr); }
    void                SetNeutralHadronIsoMaxEB(const char* expr) { fIsoCalcNHEB.SetTitle(expr); }
    void                SetNeutralHadronIsoMaxEE(const char* expr) { fIsoCalcNHEE.SetTitle(expr); }
    void                SetPhotonIsoMaxEB(const char* expr) { fIsoCalcPhEB.SetTitle(expr); }
    void                SetPhotonIsoMaxEE(const char* expr) { fIsoCalcPhEE.SetTitle(expr); }

    void                SetTriggerObjectsName(const char *n)   { fTrigObjectsName = n;       }
    

    void                SetPhotonsFromBranch(bool b)           { fPhotonsFromBranch = b;           }

    void                SetPVName(const char *n)          { fPVName = n;                 }
    void                SetPVFromBranch(bool b)           { fPVFromBranch = b;           }
    void                SetIsData (Bool_t b) { fIsData = b;}

    void                SetShowerShapeType(const char *type) { fShowerShapeType        = type;    } 

    // methods to set the MC smearing/energy correction values
    void                AddEnCorrPerRun( UInt_t sRun, UInt_t eRun,
                                         Double_t corr_EBlowEta_hR9,
                                         Double_t corr_EBlowEta_lR9,
                                         Double_t corr_EBhighEta_hR9,
                                         Double_t corr_EBhighEta_lR9,                                         
                                         Double_t corr_EElowEta_hR9,
                                         Double_t corr_EElowEta_lR9,
                                         Double_t corr_EEhighEta_hR9,
                                         Double_t corr_EEhighEta_lR9) {

      fDataEnCorr_EBlowEta_hR9.push_back(corr_EBlowEta_hR9);
      fDataEnCorr_EBlowEta_lR9.push_back(corr_EBlowEta_lR9);
      fDataEnCorr_EBhighEta_hR9.push_back(corr_EBhighEta_hR9);
      fDataEnCorr_EBhighEta_lR9.push_back(corr_EBhighEta_lR9);      
      fDataEnCorr_EElowEta_hR9.push_back(corr_EElowEta_hR9);
      fDataEnCorr_EElowEta_lR9.push_back(corr_EElowEta_lR9);      
      fDataEnCorr_EEhighEta_hR9.push_back(corr_EEhighEta_hR9);
      fDataEnCorr_EEhighEta_lR9.push_back(corr_EEhighEta_lR9);       
      fRunStart.push_back         (sRun);
      fRunEnd.push_back           (eRun);
    };

    void                SetMCSmearFactors(Double_t _EBlowEta_hR9,
                                          Double_t _EBlowEta_lR9,
                                          Double_t _EBhighEta_hR9, 
                                          Double_t _EBhighEta_lR9,
                                          Double_t _EElowEta_hR9, 
                                          Double_t _EElowEta_lR9,
                                          Double_t _EEhighEta_hR9,
                                          Double_t _EEhighEta_lR9) {
      fMCSmear_EBlowEta_hR9 = _EBlowEta_hR9;
      fMCSmear_EBlowEta_lR9 = _EBlowEta_lR9;
      fMCSmear_EBhighEta_hR9 = _EBhighEta_hR9;
      fMCSmear_EBhighEta_lR9 = _EBhighEta_lR9;      
      fMCSmear_EElowEta_hR9 = _EElowEta_hR9;
      fMCSmear_EElowEta_lR9 = _EElowEta_lR9;
      fMCSmear_EEhighEta_hR9 = _EEhighEta_hR9;
      fMCSmear_EEhighEta_lR9 = _EEhighEta_lR9;      
    };




    void                SetGoodElectronsFromBranch(Bool_t b) { fGoodElectronsFromBranch = b; }
    void                SetGoodElectronName(TString name) { fGoodElectronName = name; }

    void                SetBdtCutBarrel(double a) {fbdtCutBarrel = a; }
    void                SetBdtCutEndcap(double a) {fbdtCutEndcap = a; }


    void                DoDataEneCorr(bool a)              { fDoDataEneCorr = a; }
    void                DoMCSmear(bool a)                  { fDoMCSmear     = a; }
    void                SetDoShowerShapeScaling(Bool_t b)  { fDoShowerShapeScaling = b; }

    // replaced by ShowerShapeScaling (fab)
    /*     void                SetDoMCR9Scaling(Bool_t b)        { fDoMCR9Scaling = b; } */
    /*     void                SetMCR9Scale(Double_t ebscale, Double_t eescale) { fMCR9ScaleEB = ebscale; fMCR9ScaleEE = eescale; } */
    /*     void                SetDoMCSigIEtaIEtaScaling(Bool_t b)        { fDoMCSigIEtaIEtaScaling = b; } */
    /*     void                SetDoMCWidthScaling(Bool_t b)        { fDoMCWidthScaling = b; } */

    void                SetDoMCErrScaling(Bool_t b)        { fDoMCErrScaling = b; }
    void                SetMCErrScale(Double_t ebscale, Double_t eescale) { fMCErrScaleEB = ebscale; fMCErrScaleEE = eescale; }

    void                SetIdMVAType(const char *type)    { fIdMVATypeName        = type;    } 

    void                SetRhoType(RhoUtilities::RhoType); /*DEPRECATED*/
    void                SetRhoAlgo(UInt_t algo) { fRhoAlgo = algo; }
    
    enum EPhIdType {
      kIdUndef = 0,       //not defined
      kTight,             //"Tight"
      kLoose,             //"Loose"
      kLooseEM,           //"LooseEM"
      kBaseLineCiC,          //"2011" Hgg BaseLine CiC
      kBaseLineCiCPF,        //"2012" Hgg BaseLine CiC
      kBaseLineCiCPFNoPresel,//"2012" Hgg BaseLine CiC plus eleveto -- for mono photon
      kMITMVAId,             // MingMing MVA ID
      kMITPhSelection,       //MIT loose preselection (for mva)
      kMITPFPhSelection,     //MIT loose preselection (for mva)
      kMITPFPhSelectionNoEcal,
      kMITPFPhSelection_NoTrigger,    //MIT loose preselection (for mva, no Trigger)
      kVgamma2011Selection,  // Vgamma 2011 Photon ID
      kTrivialSelection,     // only pt & eta cuts
      kEgammaMedium,
      kCustomId             //"Custom"
    };

    enum EPhIsoType {
      kIsoUndef = 0,      //not defined        
      kNoIso,             //"NoIso"
      kCombinedIso,       //"CombinedIso"
      kCustomIso,         //"Custom"
      kMITPUCorrected,     //PileUp Corrected Hgg Isolation
      kPFPUCorrected
    };

  protected:
    void                Process();
    void                SlaveBegin();

    Int_t               FindRunRangeIdx(UInt_t run);
    Double_t            GetDataEnCorr(Int_t runRange, PhotonTools::eScaleCats cat);
    Double_t            GetMCSmearFac(PhotonTools::eScaleCats cat);   // last flag in case of special smearing for error computation
      
    TString             fPhotonBranchName;     //name of photon collection (input)
    TString             fGoodPhotonsName;      //name of exported "good photon" collection
    TString             fTrackBranchName;      // name of the track collection (only needed for PU corrected isolation)
    TString             fBeamspotBranchName;   //name of the Beamspot collection (only needed for PU corrected isolation)
    TString             fPileUpDenName;        //name of the PU density collection      
    TString             fConversionName;       //name of conversion branch
    TString             fElectronName;
    TString             fGoodElectronName;
    TString             fTrigObjectsName;        //name of trigger object collection
    TString             fPVName;
    TString             fMCParticleName;
    TString             fPileUpName;
    TString             fPFCandsName;
    TString             fPhotonIDType;         //type of photon identification we impose
    TString             fPhotonIsoType;        //type of photon isolation we impose
    Double_t            fPhotonPtMin;          //min pt cut
    Double_t            fHadOverEmMax;         //maximum of hadronic/em energy
    Bool_t              fApplySpikeRemoval;    //whether apply spike removal      
    Bool_t              fApplyPixelSeed;       //=true then apply pixel seed constraint
    Bool_t              fApplyElectronVeto;    //=true then apply electron veto (with no conversion recovery)
    Bool_t              fInvertElectronVeto;    //=true then invert electron veto (for cic selection only atm)      
    Bool_t              fApplyElectronVetoConvRecovery; //=true then apply electron veto with conversion recovery
    Bool_t              fApplyConversionId;    //=true then apply conversion id cuts
    Bool_t              fApplyTriggerMatching;   //match to hlt photon (default=0)      
    Double_t            fPhotonR9Min;          //min R9 value
    EPhIdType           fPhIdType;             //!identification scheme
    EPhIsoType          fPhIsoType;            //!isolation scheme
    Bool_t              fFiduciality;          //=true then apply fiducual requirement

    Double_t            fEtaWidthEB;  	 //max Eta Width in ECAL Barrel
    Double_t            fEtaWidthEE;  	 //max Eta Width in ECAL End Cap
    Double_t            fAbsEtaMax;  	         //max Abs Eta
    Double_t            fAbsEtaMin;  	         //min Abs Eta
    Bool_t              fApplyR9Min;           //apply R9 min
    Double_t            fEffAreaEcalEE;
    Double_t            fEffAreaHcalEE;
    Double_t            fEffAreaTrackEE;
    Double_t            fEffAreaEcalEB;
    Double_t            fEffAreaHcalEB;
    Double_t            fEffAreaTrackEB;

    TFormula            fIsoCalcCHEB;
    TFormula            fIsoCalcCHEE;
    TFormula            fIsoCalcNHEB;
    TFormula            fIsoCalcNHEE;
    TFormula            fIsoCalcPhEB;
    TFormula            fIsoCalcPhEE;

    const PhotonCol    *fPhotons;              //!photon branch
    const TrackCol     *fTracks;               //!track branch
    const BeamSpotCol  *fBeamspots;            //!beamspot branch    
    const PileupEnergyDensityCol *fPileUpDen;  //!rho branch
    const DecayParticleCol *fConversions;      //!conversion branch
    const ElectronCol  *fElectrons;            //!electron branch
    const ElectronCol  *fGoodElectrons;        //!electron branch
    const VertexCol*    fPV;                   //!
    const MCParticleCol          *fMCParticles;//!
    const PileupInfoCol          *fPileUp;     //!  
    const PFCandidateCol          *fPFCands;     //!  

    Double_t fbdtCutBarrel;
    Double_t fbdtCutEndcap;

    // ----------------------------------------------------------------
    // these guys should go away.... (let MVATools handle this)   (fab)
    int                         fVariableType;
    TString                     fEndcapWeights;
    TString                     fBarrelWeights;
    // ----------------------------------------------------------------
    MVATools*                   fTool;         //!
    TString                     fIdMVATypeName;

    // ----------------------------------------------------------------

    Bool_t fDoMCR9Scaling;
    Double_t fMCR9ScaleEB;
    Double_t fMCR9ScaleEE;
      
    Bool_t fDoMCSigIEtaIEtaScaling;
    Bool_t fDoMCWidthScaling;
      
    Bool_t fDoMCErrScaling;
    Double_t fMCErrScaleEB;
    Double_t fMCErrScaleEE;    

    Bool_t              fPhotonsFromBranch;
    Bool_t              fPVFromBranch;
    Bool_t              fGoodElectronsFromBranch;
    Bool_t              fIsData;

    // showershape
    TString                                fShowerShapeType;
    PhotonTools::ShowerShapeScales         fSSType;
      
    bool                  fDoDataEneCorr;
    bool                  fDoMCSmear;
    Bool_t                fDoShowerShapeScaling; 


    // Vectroes to hols smeraring/correction factors
    std::vector<Double_t> fDataEnCorr_EBlowEta_hR9;
    std::vector<Double_t> fDataEnCorr_EBlowEta_lR9;
    std::vector<Double_t> fDataEnCorr_EBhighEta_hR9;
    std::vector<Double_t> fDataEnCorr_EBhighEta_lR9;    
    std::vector<Double_t> fDataEnCorr_EElowEta_hR9;
    std::vector<Double_t> fDataEnCorr_EElowEta_lR9;
    std::vector<Double_t> fDataEnCorr_EEhighEta_hR9;
    std::vector<Double_t> fDataEnCorr_EEhighEta_lR9;
    
    std::vector<UInt_t>   fRunStart;
    std::vector<UInt_t>   fRunEnd;
      
    Double_t              fMCSmear_EBlowEta_hR9;
    Double_t              fMCSmear_EBlowEta_lR9;
    Double_t              fMCSmear_EBhighEta_hR9;
    Double_t              fMCSmear_EBhighEta_lR9;    
    Double_t              fMCSmear_EElowEta_hR9;
    Double_t              fMCSmear_EElowEta_lR9;
    Double_t              fMCSmear_EEhighEta_hR9;
    Double_t              fMCSmear_EEhighEta_lR9;    
      
    TRandom3* fRng;       //!

    UInt_t    fRhoAlgo;

    ClassDef(PhotonIDMod, 1) // Photon identification module
  };
}
#endif
