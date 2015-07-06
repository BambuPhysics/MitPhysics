//--------------------------------------------------------------------------------------------------
// GenFakeableObjsMod
//
// This Module generates a collection of Electron and Muon Fakeable Objects. The exact definition
// of the electron and muon fakeable objects can be specified according to the ElectronFOType 
// and MuonFOType enum types.
//
// Authors: S.Xie
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_FAKEMODS_GENFAKEABLEOBJSMOD_H
#define MITPHYSICS_FAKEMODS_GENFAKEABLEOBJSMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/CollectionsFwd.h"

namespace mithep 
{
  class GenFakeableObjsMod : public BaseMod
  {
    public:
      GenFakeableObjsMod(const char *name="GenFakeableObjsMod", 
                         const char *title="Fake Object Generation Module");

      Bool_t        GetApplyConversionFilter()       const { return fApplyConvFilter;              }
      Bool_t        GetApplyD0Cut()                  const { return fApplyD0Cut;                   }
      Bool_t        GetChargeFilter()                const { return fChargeFilter;                 }
      Bool_t        GetNWrongHitsMax()               const { return fNWrongHitsMax;                }
      Double_t      GetD0Cut()                       const { return fD0Cut;                        }
      Double_t      GetCombIsolationCut()            const { return fCombIsolationCut;             }
      Double_t      GetTrackIsolationCut()           const { return fTrackIsolationCut;            }
      Double_t      GetEcalIsolationCut()            const { return fEcalIsolationCut;             }
      Double_t      GetHcalIsolationCut()            const { return fHcalIsolationCut;             }
      Bool_t        GetVetoTriggerJet()              const { return fVetoTriggerJet;               }
      Bool_t        GetVetoGenLeptons()              const { return fVetoGenLeptons;               }
      Bool_t        GetVetoCleanLeptons()            const { return fVetoCleanLeptons;             }
      const char   *GetElectronsName()               const { return fElectronBranchName;           }
      const char   *GetMuonsName()                   const { return fMuonBranchName;               }
      const char   *GetElectronFOType()              const { return fElectronFOType;               }
      const char   *GetMuonFOType()                  const { return fMuonFOType;                   }
      const char   *GetTriggerName()                 const { return fTriggerName;                  }
      const char   *GetTriggerObjectsName()          const { return fTriggerObjectsName;           }
      const char   *GetTrackName()                   const { return fTrackBranchName;              }
      const char   *GetGsfTrackBranchName()          const { return fGsfTrackBranchName;           }
      const char   *GetBarrelSuperClustersName()     const { return fBarrelSuperClusterBranchName; }
      const char   *GetEndcapSuperClustersName()     const { return fEndcapSuperClusterBranchName; }
      const char   *GetJetBranchName()               const { return fJetBranchName;                }
      const char   *GetVertexName()                  const { return fVertexName;                   }
      const char   *GetConversionBranchName()        const { return fConversionBranchName;         }
      const char   *GetCleanPhotonsName()            const { return fCleanPhotonsName;             }
      const char   *GetCleanJetsName()               const { return fCleanJetsName;                }
      const char   *GetElFakeableObjsName()          const { return fElFakeableObjsName;           }
      const char   *GetMuFakeableObjsName()          const { return fMuFakeableObjsName;           }
      const char   *GetMCLeptonsName()               const { return fMCLeptonsName;                }
      const char   *GetMCTausName()                  const { return fMCTausName;                   }

      void SetApplyConversionFilter(Bool_t b)                { fApplyConvFilter            = b;    }
      void SetApplyD0Cut(Bool_t b)                           { fApplyD0Cut                 = b;    }
      void SetD0Cut(Double_t cut)                            { fD0Cut                      = cut;  }
      void SetCombIsolationCut(Double_t cut)                 { fCombIsolationCut           = cut;  }
      void SetTrackIsolationCut(Double_t cut)                { fTrackIsolationCut          = cut;  }
      void SetEcalIsolationCut(Double_t cut)                 { fEcalIsolationCut           = cut;  }
      void SetHcalIsolationCut(Double_t cut)                 { fHcalIsolationCut           = cut;  }
      void SetChargeFilter(Bool_t b)                         { fChargeFilter               = b;    }
      void SetNWrongHitsMax(UInt_t n)                        { fNWrongHitsMax              = n;    }
      void SetVetoTriggerJet(Bool_t b)                       { fVetoTriggerJet             = b;    }
      void SetVetoGenLeptons(Bool_t b)                       { fVetoGenLeptons             = b;    }
      void SetVetoCleanLeptons(Bool_t b)                     { fVetoCleanLeptons           = b;    }
      void SetElectronFOType(const char *name)               { fElectronFOType             = name; }
      void SetMuonFOType(const char *name)                   { fMuonFOType                 = name; }
      void SetTriggerName(const char *name)                  { fTriggerName                = name; }
      void SetTriggerObjectsName(const char *name)           { fTriggerObjectsName         = name; }
      void SetElectronBranchName(const char *name)           { fElectronBranchName         = name; }
      void SetMuonBranchName(const char *name)               { fMuonBranchName             = name; }
      void SetTrackBranchName(const char *name)              { fTrackBranchName            = name; }
      void SetGsfTrackBranchName(const char *name)           { fGsfTrackBranchName         = name; }
      void SetBarrelSuperClusterBranchName(const char *name) { fBarrelSuperClusterBranchName=name; }
      void SetEndcapSuperClusterBranchName(const char *name) { fEndcapSuperClusterBranchName=name; }
      void SetJetBranchName(const char *name)                { fJetBranchName              = name; }
      void SetVertexName(const char *name)                   { fVertexName                 = name; }
      void SetConversionBranchName(const char *name)         { fConversionBranchName       = name; }
      void SetCleanPhotonsName(const char *name)             { fCleanPhotonsName           = name; }
      void SetCleanJetsName(const char *name)                { fCleanJetsName              = name; }
      void SetElFakeableObjsName(const char *name)           { fElFakeableObjsName         = name; }
      void SetMuFakeableObjsName(const char *name)           { fMuFakeableObjsName         = name; }
      void SetMCLeptonsName(const char *name)                { fMCLeptonsName              = name; }
      void SetMCTausName(const char *name)                   { fMCTausName                 = name; }

     enum ElectronFOType {
        kElFOUndef = 0,    //not defined
        kElFOIso,          //"Reco Electron with full isolation"
        kElFOLooseIdLooseIso         //"Reco Electron with loose id and isolation
      };
      enum MuonFOType {
        kMuFOUndef = 0,    //not defined
        kMuFOIsoTrack,     //"Loosely Isolated Track"
        kMuFOGlobal,       //"GlobalMuon with loose isolation"
        kMuFOTrackerMuon   //"TrackerMuon with loose isolation"
      };

    protected:
      void             Process();
      void             SlaveBegin();

      Bool_t           fApplyConvFilter;               //whether remove conversions
      UInt_t           fNWrongHitsMax;                 //whether to use wrong hits req for conversion removal
      Bool_t           fApplyD0Cut;                    //whether apply d0 cut
      Bool_t           fChargeFilter;                  //whether apply GSF and CFT equal requirement
      Double_t         fD0Cut;                         //max d0
      Double_t         fCombIsolationCut;                  //max isolation
      Double_t         fTrackIsolationCut;             //max isolation
      Double_t         fEcalIsolationCut;              //max isolation
      Double_t         fHcalIsolationCut;              //max isolation
      Bool_t           fVetoTriggerJet;                //whether to veto on the leading jet
      Bool_t           fVetoGenLeptons;                //whether we exclude gen leptons
      Bool_t           fVetoCleanLeptons;              //whether we exclude clean leptons
      TString          fElectronFOType;                //type of electron Fakeable object
      TString          fMuonFOType;                    //type of muon Fakeable object
      TString          fTriggerName;                   //name of trigger 
      TString          fTriggerObjectsName;            //name of trigger objects
      TString          fElectronBranchName;            //name of electron brach              (input)
      TString          fMuonBranchName;                //name of muon brach                  (input)
      TString          fTrackBranchName;               //name of track brach                 (input)
      TString          fGsfTrackBranchName;            //name of track collection            (input)
      TString          fBarrelSuperClusterBranchName;  //name of barrel supercluster branch  (input)
      TString          fEndcapSuperClusterBranchName;  //name of endcap supercluster branch  (input)
      TString          fJetBranchName;                 //name of jet branch                  (input)
      TString          fVertexName;                    //name of vertex branch               (input)
      TString          fConversionBranchName;          //name of conversion collection       (input)
      TString          fCleanPhotonsName;              //name of clean photons               (input)
      TString          fCleanJetsName;                 //name of clean jets                  (input)
      TString          fMCLeptonsName;                 //name of MC leptons                  (input)
      TString          fMCTausName;                    //name of MC taus                     (input)
      TString          fElFakeableObjsName;            //name of fakeable objects           (output)
      TString          fMuFakeableObjsName;            //name of fakeable objects           (output)
      ElectronFOType   fElFOType;                      //Electron Fakeable Object type
      MuonFOType       fMuFOType;                      //Muon Fakeable Object type
      const ElectronCol      *fElectrons;             //!Electron branch
      const MuonCol          *fMuons;                 //!Muon branch
      const SuperClusterCol  *fBarrelSuperClusters;   //!Barrel Supercluster branch
      const SuperClusterCol  *fEndcapSuperClusters;   //!Endcap Supercluster branch
      const TrackCol         *fTracks;                //!Track branch
      const TrackCol	     *fGsfTracks;	      //!GsfTrack branch
      const JetCol           *fJets;                  //!Jet branch
      const VertexCol        *fVertices;              //!Vertex branch
      const DecayParticleCol *fConversions;           //!conversion collection       

      ClassDef(GenFakeableObjsMod, 1) // Fakeable objects generation module
  };
}
#endif
