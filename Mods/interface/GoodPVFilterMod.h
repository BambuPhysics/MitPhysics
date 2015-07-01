//--------------------------------------------------------------------------------------------------
// GoodPVFilterMod
//
// This module selects events with a good reconstructed Primary Vertex according to
// configurable cuts on the number of tracks and vertex position.
//
// Authors: J.Bendavid
//--------------------------------------------------------------------------------------------------

#ifndef MITMODS_MODS_GOODPVFILTERMOD_H
#define MITMODS_MODS_GOODPVFILTERMOD_H

#include <string>
#include <TString.h>
#include <TH1F.h>
#include "MitAna/DataTree/interface/VertexFwd.h" 
#include "MitAna/DataTree/interface/PileupInfoFwd.h"
#include "MitAna/TreeMod/interface/BaseMod.h" 

namespace mithep 
{
  class GoodPVFilterMod : public BaseMod {
    public:
      
      enum ECuts {
        eNTracks,
        eNDof,
        eZ,
        eRho
      };
      
      GoodPVFilterMod(const char *name="GoodPVFilterMod", const char *title="Good PV Filter Module");
      ~GoodPVFilterMod();

      Int_t                       GetNEvents()      const { return fNEvents;       }
      Int_t                       GetNAccepted()    const { return fNAcceped;      }
      Int_t                       GetNFailed()      const { return fNFailed;       }
      const char                 *GetOutputName()   const { return fGoodVertexesName; }
      const char                 *GetGoodVertexesName() const { return GetOutputName(); }
      void                        SetAbortIfNotAccepted(Bool_t b)   { fAbort = b;           }
      void                        SetIsMC(Bool_t b)                 { fIsMC = b;            }
      void                        SetMinVertexNTracks(UInt_t n)     { fMinVertexNTracks = n;}
      void                        SetMinNDof(UInt_t n)              { fMinNDof = n; 	    }
      void                        SetMaxAbsZ(Double_t x)  	    { fMaxAbsZ = x; 	    }
      void                        SetMaxRho(Double_t x)   	    { fMaxRho = x;          }
      void                        SetVertexesName(TString s)        { fVertexesName = s;    }
      void                        SetGoodVertexesName(TString s)    { SetOutputName(s);     }
      void                        SetOutputName(TString s)          { fGoodVertexesName = s; }
      

    protected:
      void                        BeginRun();
      const BitMask8              FailedCuts(const mithep::Vertex *v) const;
      virtual void                OnAccepted()  {/*could be implemented in derived classes*/}
      virtual void                OnFailed()    {/*could be implemented in derived classes*/}
      void                        Process();
      void                        SlaveBegin();
      void                        SlaveTerminate();

      Bool_t                      fAbort;         //=true then abort (sub-)modules if not accepted
      Bool_t                      fIsMC;
      UInt_t                      fMinVertexNTracks; //minimum number of tracks for the vertex
      UInt_t                      fMinNDof;       //minimum number of degrees of freedom
      Double_t                    fMaxAbsZ;       //maximum abs(z) of the vertex
      Double_t                    fMaxRho;        //maximum rho of the vertex
      TString                     fVertexesName;  //Name of PV collection
      TString                     fGoodVertexesName; //Name of newPV collection
      TString                     fPileupInfoName;
      Int_t                       fNEvents;       //!number of processed events
      Int_t                       fNAcceped;      //!number of accepted events
      Int_t                       fNFailed;       //!number of failed events
      TH1F                       *hVertexNTracks;
      TH1F                       *hVertexNDof;
      TH1F                       *hVertexRho;
      TH1F                       *hVertexZ;
      TH1F                       *hNVtx;
      TH1F                       *hNGoodVtx;
      TH1D                       *hNGenVtx;
      TH1D                       *hNPU;

    ClassDef(GoodPVFilterMod, 1) // L1 TAM module
  };
}
#endif
