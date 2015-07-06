#include <TFile.h>
#include <TTree.h>
#include "MitAna/DataTree/interface/Names.h"
#include "MitAna/DataTree/interface/Track.h"
#include "MitPhysics/Mods/interface/TrackingPurityFilterMod.h"

using namespace mithep;

ClassImp(mithep::TrackingPurityFilterMod)

//--------------------------------------------------------------------------------------------------
TrackingPurityFilterMod::TrackingPurityFilterMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fAbort(kTRUE),
  fMinNTracksCut(10),
  fMinHighPurFraction(0.2),
  fTracksName(Names::gkTrackBrn),
  fNEvents(0),
  fNAcceped(0),
  fNFailed(0),
  fTracks(0),
  hNTracksBefore(0),
  hNTracks(0),
  hNHighPurityTracks(0),
  hHighPurityFraction(0)
{
  // Constructor. 
}

//--------------------------------------------------------------------------------------------------
TrackingPurityFilterMod::~TrackingPurityFilterMod() 
{
  // Destructor.
}


//--------------------------------------------------------------------------------------------------
void TrackingPurityFilterMod::BeginRun()
{
  
}

//--------------------------------------------------------------------------------------------------
void TrackingPurityFilterMod::Process()
{
  // Increment counters and stop further processing of an event if current run is excluded

  auto* tracks = GetObject<TrackCol>(fTracksName);

  ++fNEvents; 
  
  UInt_t nTracks = 0;
  UInt_t nHighPurityTracks = 0;
  
  for (UInt_t i=0; i<tracks->GetEntries(); ++i) {
    const Track *t = tracks->At(i);
    ++nTracks;
    if (t->Quality().Quality(TrackQuality::highPurity))
      ++nHighPurityTracks;
  }
  
  Double_t highPurityFraction;
  
  if (nHighPurityTracks==0)
    highPurityFraction = 0.0;
  else
    highPurityFraction = static_cast<Double_t>(nHighPurityTracks)/static_cast<Double_t>(nTracks);
  
  hNTracksBefore->Fill(nTracks);
  hNHighPurityTracks->Fill(nHighPurityTracks);
  hHighPurityFraction->Fill(highPurityFraction);
  
  // take action if failed
  if (nTracks >= fMinNTracksCut && highPurityFraction<fMinHighPurFraction) {
    ++fNFailed;
    OnFailed();
    if (fAbort) {
      SkipEvent(); // abort processing of this event by sub-modules
    }
    return;
  } 

  // take action if accepted
  hNTracks->Fill(nTracks);
  
  ++fNAcceped;
  IncNEventsProcessed();
  OnAccepted();
}

//--------------------------------------------------------------------------------------------------
void TrackingPurityFilterMod::SlaveBegin()
{
  hNTracksBefore = new TH1F("hNTracksBefore", "hNTracksBefore", 401, -0.5,400.5);
  AddOutput(hNTracksBefore);
  
  hNTracks = new TH1F("hNTracks", "hNTracks", 401, -0.5,400.5);
  AddOutput(hNTracks);
  
  hNHighPurityTracks = new TH1F("hNHighPurityTracks", "hNHighPurityTracks", 401, -0.5,400.5);
  AddOutput(hNHighPurityTracks);
  
  hHighPurityFraction = new TH1F("hHighPurityFraction", "hHighPurityFraction", 100, 0.0, 1.0);
  AddOutput(hHighPurityFraction);
  
}

//--------------------------------------------------------------------------------------------------
void TrackingPurityFilterMod::SlaveTerminate()
{
  // Save number of accepted events.

  SaveNEventsProcessed();
}
