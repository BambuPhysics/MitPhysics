//--------------------------------------------------------------------------------------------------
// TrackingPurityFilterMod
//
// This module implements the "fraction of high purity tracks" filter to remove
// scraping events.
//
// Authors: J.Bendavid
//--------------------------------------------------------------------------------------------------

#ifndef MITMODS_MODS_TRACKINGPURITYFILTERMOD_H
#define MITMODS_MODS_TRACKINGPURITYFILTERMOD_H

#include <string>
#include <TString.h>
#include <TH1F.h>
#include "MitAna/DataTree/interface/TrackFwd.h" 
#include "MitAna/TreeMod/interface/BaseMod.h" 

namespace mithep 
{
  class TrackingPurityFilterMod : public BaseMod {
    public:
      
      TrackingPurityFilterMod(const char *name="TrackingPurityFilterMod", const char *title="Good PV Filter Module");
      ~TrackingPurityFilterMod();

      Int_t                       GetNEvents()      const { return fNEvents;       }
      Int_t                       GetNAccepted()    const { return fNAcceped;      }
      Int_t                       GetNFailed()      const { return fNFailed;       }
      void                        SetAbortIfNotAccepted(Bool_t b)   { fAbort         = b; }
      void                        SetMinNTracksCut(UInt_t n) { fMinNTracksCut = n; }
      void                        SetMinHighPurFraction(Double_t x)  { fMinHighPurFraction = x; }

    protected:
      void                        BeginRun();
      virtual void                OnAccepted()  {/*could be implemented in derived classes*/}
      virtual void                OnFailed()    {/*could be implemented in derived classes*/}
      void                        Process();
      void                        SlaveBegin();
      void                        SlaveTerminate();

      Bool_t                      fAbort;         //=true then abort (sub-)modules if not accepted
      UInt_t                      fMinNTracksCut; //minimum number of tracks for the vertex
      Double_t                    fMinHighPurFraction;       //maximum abs(z) of the vertex
      TString                     fTracksName;    //Name of tracks collection
      Int_t                       fNEvents;       //!number of processed events
      Int_t                       fNAcceped;      //!number of accepted events
      Int_t                       fNFailed;       //!number of failed events
      const TrackCol             *fTracks;        //!track collection
      TH1F                       *hNTracksBefore;
      TH1F                       *hNTracks;
      TH1F                       *hNHighPurityTracks;
      TH1F                       *hHighPurityFraction;

    ClassDef(TrackingPurityFilterMod, 1) // L1 TAM module
  };
}
#endif
