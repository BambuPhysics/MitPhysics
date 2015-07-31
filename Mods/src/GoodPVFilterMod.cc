#include "MitPhysics/Mods/interface/GoodPVFilterMod.h"
#include "MitAna/DataTree/interface/VertexCol.h"

using namespace mithep;

ClassImp(mithep::GoodPVFilterMod)

//--------------------------------------------------------------------------------------------------
void
GoodPVFilterMod::Process()
{
  auto* vertices = GetObject<VertexCol>(fVertexesName);

  VertexOArr *GoodVertexes = new VertexOArr;
  GoodVertexes->SetName(fGoodVertexesName);

  // Increment counters and stop further processing of an event if current run is excluded

  for (unsigned iV = 0; iV != vertices->GetEntries(); ++iV) {
    auto& vertex(*vertices->At(iV));

    unsigned nTracks = vertex.NTracksFit();
    unsigned nDof = vertex.Ndof();
    double z = vertex.Position().Z();
    double rho = vertex.Position().Rho();

    BitMask8 failed;

    if (nTracks < fMinVertexNTracks)
      failed.SetBit(eNTracks);

    if (nDof < fMinNDof)
      failed.SetBit(eNDof);

    if (std::abs(z) > fMaxAbsZ)
      failed.SetBit(eZ);

    if (rho > fMaxRho)
      failed.SetBit(eRho);

    if (failed.NBitsSet() > 1)
      continue;

    BitMask8 failedNTracks(failed);
    failedNTracks.ClearBit(eNTracks);
    if (failedNTracks.NBitsSet() == 0)
      hVertexNTracks->Fill(nTracks);

    BitMask8 failedNDof(failed);
    failedNDof.ClearBit(eNDof);
    if (failedNDof.NBitsSet() == 0)
      hVertexNDof->Fill(nDof);

    BitMask8 failedZ(failed);
    failedZ.ClearBit(eZ);
    if (failedZ.NBitsSet() == 0)
      hVertexZ->Fill(z);

    BitMask8 failedRho(failed);
    failedRho.ClearBit(eRho);
    if (failedRho.NBitsSet() == 0)
      hVertexRho->Fill(rho);

    if (failed.NBitsSet() == 0) {
      vertex.Mark();
      GoodVertexes->Add(&vertex);
    }
  }

  //fill histograms
  hNVtx->Fill(vertices->GetEntries());
  hNGoodVtx->Fill(GoodVertexes->GetEntries());

  // add objects for other modules to use
  AddObjThisEvt(GoodVertexes);

  // take action if failed
  if (GoodVertexes->GetEntries() == 0) {
    OnFailed();
    if (fAbort)
      SkipEvent(); // abort processing of this event by sub-modules

    return;
  }

  // take action if accepted
  OnAccepted();

  IncNEventsProcessed();
}

//--------------------------------------------------------------------------------------------------
void
GoodPVFilterMod::SlaveBegin()
{
  hVertexNTracks = new TH1F("hVertexNTracks", "hVertexNTracks", 401, -0.5,400.5);
  AddOutput(hVertexNTracks);

  hVertexNDof = new TH1F("hVertexNDof", "hVertexNDof", 401, -0.5,400.5);
  AddOutput(hVertexNDof);

  hVertexZ = new TH1F("hVertexZ", "hVertexZ", 100, -100.0, 100.0);
  AddOutput(hVertexZ);

  hVertexRho = new TH1F("hVertexRho", "hVertexRho", 100, 0.0, 20.0);
  AddOutput(hVertexRho);

  hNVtx = new TH1F("hNVtx", "hNVtx", 51, -0.5, 50.5);
  AddOutput(hNVtx);

  hNGoodVtx = new TH1F("hNGoodVtx", "hNGoodVtx", 51, -0.5, 50.5);
  AddOutput(hNGoodVtx);
}

//--------------------------------------------------------------------------------------------------
void
GoodPVFilterMod::SlaveTerminate()
{
  // Save number of accepted events.
  SaveNEventsProcessed();
}
