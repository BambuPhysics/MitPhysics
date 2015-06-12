//--------------------------------------------------------------------------------------------------
// MuonIDMod
//
// This module applies muon identification criteria and exports a pointer to a collection
// of "good muons" according to the specified ID scheme.
//
// See http://indico.cern.ch/contributionDisplay.py?contribId=1&confId=45945
// See http://indico.cern.ch/getFile.py/access?contribId=1&resId=0&materialId=slides&confId=42229
//
// Authors: S.Xie
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_MUONIDMOD_H
#define MITPHYSICS_MODS_MUONIDMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/MuonFwd.h"
#include "MitAna/DataTree/interface/VertexFwd.h"
#include "MitAna/DataTree/interface/TrackFwd.h"
#include "MitAna/DataTree/interface/PFCandidateFwd.h"
#include "MitPhysics/Utils/interface/MuonTools.h"
#include "MitPhysics/Utils/interface/MuonIDMVA.h"
#include "MitPhysics/Utils/interface/IsolationTools.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityCol.h"

#include "MitPhysics/Utils/interface/RhoUtilities.h"

namespace mithep 
{
  class MuonIDMod : public BaseMod
  {
  public:
    MuonIDMod(const char *name="MuonIDMod", 
              const char *title="Muon identification module");

    Double_t           GetCaloIsoCut()                const { return fCaloIsolationCut;   }
    const char        *GetClassType()                 const { return fMuonClassType;      }
    const char        *GetCleanName()                 const { return GetCleanMuonsName(); }   
    const char        *GetCleanMuonsName()            const { return fCleanMuonsName;     }   
    Double_t           GetCombIsoCut()                const { return fCombIsolationCut;   }
    const char        *GetIDType()                    const { return fMuonIDType;         }
    const char        *GetInputName()                 const { return fMuonBranchName;     }   
    const char        *GetIsoType()                   const { return fMuonIsoType;        }
    const char        *GetOutputName()                const { return GetCleanMuonsName(); }   
    Double_t           GetPtMin()                     const { return fMuonPtMin;          }
    Double_t           GetTrackIsoCut()               const { return fTrackIsolationCut;  }
    Bool_t             PassMuonMVA_BDTG_IdIso(const Muon *mu, const Vertex *vertex, 
                                              const PileupEnergyDensityCol *PileupEnergyDensity) const;
    Bool_t             PassMuonIsoRingsV0_BDTG_Iso(const Muon *mu, const Vertex *vertex, 
                                                   const PileupEnergyDensityCol *PileupEnergyDensity) const;
    Bool_t             PassMuonIsoDeltaR(const Muon *mu, const Vertex *vertex, 
                                         const PileupEnergyDensityCol *PileupEnergyDensity) const;
    void               SetPrintMVADebugInfo(Bool_t b)       { fPrintMVADebugInfo = b;     }
    void               SetApplyD0Cut(Bool_t b)              { fApplyD0Cut        = b;     }
    void               SetApplyDZCut(Bool_t b)              { fApplyDZCut        = b;     }
    void               SetCaloIsoCut(Double_t cut)          { fCaloIsolationCut  = cut;   }
    void               SetClassType(const char *type)       { fMuonClassType     = type;  }
    void               SetCleanMuonsName(const char *name)  { fCleanMuonsName    = name;  }   
    void               SetOldMuonsName(const char *n)       { fNonIsolatedMuonsName  = n; }  
    void               SetOldElectronsName(const char *n)   { fNonIsolatedElectronsName  = n; }  
    void               SetCleanName(const char *name)       { SetCleanMuonsName(name);    }   
    void               SetCombIsoCut(Double_t cut)          { fCombIsolationCut  = cut;   }
    void               SetCombRelativeIsoCut(Double_t cut)  { fCombRelativeIsolationCut  = cut; }
    void               SetPFIsoCut(Double_t cut)            { fPFIsolationCut  = cut;     }
    void               SetD0Cut(Double_t cut)               { fD0Cut             = cut;   }
    void               SetDZCut(Double_t cut)               { fDZCut             = cut;   }
    void               SetWhichVertex(Int_t d)              { fWhichVertex = d;           }
    void               SetEtaCut(Double_t cut)              { fEtaCut            = cut;   }
    void               SetIDType(const char *type)          { fMuonIDType        = type;  }
    void               SetInputName(const char *name)       { fMuonBranchName    = name;  }   
    void               SetIsoType(const char *type)         { fMuonIsoType       = type;  }
    void               SetOutputName(const char *name)      { SetCleanMuonsName(name);    }   
    void               SetPtMin(Double_t pt)                { fMuonPtMin         = pt;    }
    void               SetTrackIsoCut(Double_t cut)         { fTrackIsolationCut = cut;   }
    void               SetIntRadius(Double_t dr)            { fIntRadius = dr;            }
    void               SetPFNoPileUpName(const char *n)     { fPFNoPileUpName  = n;       } 
    void               SetPFPileUpName(const char *n)       { fPFPileUpName  = n;         }
    void               SetRhoType(RhoUtilities::RhoType);
    void               SetRhoAlgo(UInt_t algo)              { fRhoAlgo = algo;            }

    void               SetVertexName(const char* name) { fVertexName = name; }
    void               SetPVName(const char *n)             { fPVName = n;                }

    enum EMuIdType {
      kIdUndef = 0,       //not defined
      kWMuId,             //"WMuId"
      kZMuId,             //"ZMuId"
      kTight,             //"Tight"
      kmuonPOG2012CutBasedIDTight,             //"muonPOG2012CutBasedIDTight"
      kLoose,             //"Loose"
      kWWMuIdV1,          //"WWMuIdV1"
      kWWMuIdV2,          //"WWMuIdV2"
      kWWMuIdV3,          //"WWMuIdV3"
      kWWMuIdV4,          //"WWMuIdV4"
      kNoId,              //"NoId"        
      kCustomId,          //"Custom"
      kMVAID_BDTG_IDIso   //"BDTG ID + Iso03, Iso04 Combined"
    };
    enum EMuIsoType {
      kIsoUndef = 0,                      //"not defined"
      kTrackCalo,                         //"TrackCalo"
      kTrackCaloCombined,                 //"TrackCaloCombined"
      kTrackCaloSliding,                  //"TrackCaloSliding"
      kTrackCaloSlidingNoCorrection,      //"TrackCaloSlidingNoCorrection"
      kCombinedRelativeConeAreaCorrected, //"CombinedRelativeConeAreaCorrected"        
      kCombinedRelativeEffectiveAreaCorrected,
      kCustomIso,                         //"Custom"
      kPFIso,                             //"PFIso"
      kPFRadialIso,                       //"PFRadialIso"
      kPFIsoBetaPUCorrected,              //"PFISo with PUcorrection using delta Beta
      kPFIsoEffectiveAreaCorrected,       //"PFIso with EffectiveArea Pileup Correction"
      kPFIsoNoL,                          //"PFIsoNoL"
      kNoIso,                             //"NoIso"
      kMVAIso_BDTG_IDIso,                 //"BDTG ID + Iso03, Iso04 Combined"
      kIsoRingsV0_BDTG_Iso,               //"BDTG Iso Rings"
      kIsoDeltaR                          //"BGDT Iso dR"              
    };
    enum EMuClassType {
      kClassUndef = 0,    //not defined
      kAll,               //"All"
      kGlobal,            //"Global"
      kGlobalorTracker,   //"Global or Tracker Muon"
      kGlobalTracker,     //"GlobalTracker"
      kSta,               //"Standalone"
      kTrackerMuon,       //"TrackerMuon"
      kCaloMuon,          //"CaloMuon"
      kTrackerBased,      //"TrackerMuon or CaloMuon"
      kGlobalOnly         //"GlobalOnly"
    };

  protected:
    void               Process();
    void               SlaveBegin();
    void               Terminate();

    Bool_t             fPrintMVADebugInfo;   //print MVA debug information
    TString            fMuonBranchName;      //name of muon collection (input)
    TString            fCleanMuonsName;      //name of exported "good muon" collection
    TString            fNonIsolatedMuonsName;    //name of imported "old muon" collection
    TString            fNonIsolatedElectronsName;//name of imported "old electron" collection
    TString            fVertexName;          //name of vertex collection
    TString            fBeamSpotName;        //name of beamspot collection
    TString            fTrackName;           //name of track collection
    TString            fPFCandidatesName;    //name of pfcandidates collection
    TString            fPFNoPileUpName;  //name of pfnpu collection
    TString            fPFPileUpName;    //name of pfpu collection
    TString            fMuonIDType;          //type of muon id scheme we impose
    TString            fMuonIsoType;         //type of muon isolations scheme we impose
    TString            fMuonClassType;       //type of muon class we impose
    Double_t           fTrackIsolationCut;   //cut value for track isolation
    Double_t           fCaloIsolationCut;    //cut value for calo isolation
    Double_t           fCombIsolationCut;    //cut value for combined isolation
    Double_t           fCombRelativeIsolationCut; //cut value for combined relative isolation
    Double_t           fPFIsolationCut;      //cut value for combined isolation
    Double_t           fMuonPtMin;           //min muon pt
    Bool_t             fApplyD0Cut;          //=true then apply d0 cut (def=1)
    Bool_t             fApplyDZCut;          //=true then apply dz cut (def=1)
    Double_t           fD0Cut;               //max d0
    Double_t           fDZCut;               //max dz
    Int_t              fWhichVertex;         //vertex to use (-2: beamspot, -1: closest in Z)
    Double_t           fEtaCut;              //max eta, absolute value
    EMuIdType          fMuIDType;            //!muon id type (imposed)
    EMuIsoType         fMuIsoType;           //!muon iso type (imposed)
    EMuClassType       fMuClassType;         //!muon class type (imposed)
    const MuonCol      *fMuons;               //!muon collection
    const VertexCol   *fVertices;            //!vertices branch
    const BeamSpotCol *fBeamSpot;            //!beamspot branch
    const TrackCol    *fTracks;              //!track branch     
    const PFCandidateCol *fPFCandidates;     //!pfcandidate branch
    const PFCandidateCol *fPFNoPileUpCands;  //!pfnpu collection
    const PFCandidateCol *fPFPileUpCands;    //!pfpu collection
    Double_t           fIntRadius;           //!min IntRadius cut in pf isolation
    MuonCol            *fNonIsolatedMuons;    //!pointer to old muon collection 
    ElectronCol        *fNonIsolatedElectrons;//!pointer to old electron collection
    TString             fPileupEnergyDensityName;
    const PileupEnergyDensityCol *fPileupEnergyDensity;
    MuonTools          *fMuonTools;           // interface to tools for muon ID
    MuonIDMVA          *fMuonIDMVA;           // helper class for MuonMVA
    TString             fPVName;

    UInt_t              fRhoAlgo;
      
    ClassDef(MuonIDMod, 1) // Muon identification module
      };
}
#endif
