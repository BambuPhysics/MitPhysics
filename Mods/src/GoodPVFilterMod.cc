#include "MitPhysics/Mods/interface/GoodPVFilterMod.h"
#include <TFile.h>
#include <TTree.h>
#include "MitAna/DataTree/interface/Names.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitAna/DataTree/interface/Vertex.h"
#include "MitAna/DataTree/interface/PileupInfo.h"


using namespace mithep;

ClassImp(mithep::GoodPVFilterMod)

//--------------------------------------------------------------------------------------------------
GoodPVFilterMod::GoodPVFilterMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fAbort(kTRUE),
  fIsMC(kFALSE),
  fMinVertexNTracks(0),
  fMinNDof(5),
  fMaxAbsZ(15.0),
  fMaxRho(2.0),
  fGoodVertexesName(ModNames::gkGoodVertexesName),
  fPileupInfoName("PileupInfo"),
  fNEvents(0),
  fNAcceped(0),
  fNFailed(0),
  hVertexNTracks(0),
  hVertexRho(0),
  hVertexZ(0)
{
  // Constructor. 
}

//--------------------------------------------------------------------------------------------------
GoodPVFilterMod::~GoodPVFilterMod() 
{
  // Destructor.
}


//--------------------------------------------------------------------------------------------------
void GoodPVFilterMod::BeginRun()
{
  
}

//--------------------------------------------------------------------------------------------------
const BitMask8 GoodPVFilterMod::FailedCuts(const Vertex *v) const
{
  BitMask8 failedCuts;
  
  if (v->NTracksFit() < fMinVertexNTracks)
    failedCuts.SetBit(eNTracks);
  
  if (v->Ndof() < fMinNDof)
    failedCuts.SetBit(eNDof);
  
  if (TMath::Abs(v->Position().Z()) > fMaxAbsZ)
    failedCuts.SetBit(eZ);
  
  if (v->Position().Rho() > fMaxRho)
    failedCuts.SetBit(eRho);
  
  return failedCuts;
  
}

//--------------------------------------------------------------------------------------------------
void GoodPVFilterMod::Process()
{
  auto* vertices = GetObject<VertexCol>(fVertexesName);
  PileupInfoCol const* pileupInfo = 0;
  if (fIsMC)
    pileupInfo = GetObject<PileupInfoCol>(fPileupInfoName);
  
  VertexOArr *GoodVertexes = new VertexOArr;
  GoodVertexes->SetName(fGoodVertexesName);
  
  // Increment counters and stop further processing of an event if current run is excluded

  ++fNEvents; 
  Bool_t goodVertex = kFALSE;
  
  for (UInt_t i=0; i<vertices->GetEntries(); ++i) {
    const Vertex *v = vertices->At(i);
    BitMask8 failed = FailedCuts(v);
    
    if (failed.NBitsSet() > 1)
      continue;
       
    BitMask8 failedNTracks = failed;
    failedNTracks.ClearBit(eNTracks);
    if (!failedNTracks.NBitsSet())
      hVertexNTracks->Fill(v->NTracksFit());
    
    BitMask8 failedNDof = failed;
    failedNTracks.ClearBit(eNDof);
    if (!failedNDof.NBitsSet())
      hVertexNDof->Fill(v->Ndof());
    
    BitMask8 failedZ = failed;
    failedZ.ClearBit(eZ);
    if (!failedZ.NBitsSet())
      hVertexZ->Fill(v->Position().Z());
    
    BitMask8 failedRho = failed;
    failedRho.ClearBit(eRho);
    if (!failedRho.NBitsSet())
      hVertexRho->Fill(v->Position().Rho());
    
    if (!failed.NBitsSet()) {
      goodVertex = kTRUE;
      v->Mark();
      GoodVertexes->Add(v);
    }
  }
  
  //fill histograms
  hNVtx->Fill(vertices->GetEntries());
  hNGoodVtx->Fill(GoodVertexes->GetEntries());
  if (fIsMC) {
    Int_t npu = -99;
    for (UInt_t i=0; i<pileupInfo->GetEntries(); ++i) {
      const PileupInfo *puinfo = pileupInfo->At(i);
      if (puinfo->GetBunchCrossing()==0) {
        npu = puinfo->GetPU_NumInteractions();
        break;
      }
    }
    hNGenVtx->Fill(1 + npu);
    hNPU->Fill(npu);
  }

  // add objects for other modules to use
  AddObjThisEvt(GoodVertexes);  
  
  // take action if failed
  if (!goodVertex) {
    ++fNFailed;
    OnFailed();
    if (fAbort) {
      SkipEvent(); // abort processing of this event by sub-modules
    }
    return;
  } 
  //printf("goodvertex check 0\n");
  // take action if accepted
  ++fNAcceped;
  IncNEventsProcessed();
  OnAccepted();
}

//--------------------------------------------------------------------------------------------------
void GoodPVFilterMod::SlaveBegin()
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
  
  hNGenVtx = new TH1D("hNGenVtx", "hNGenVtx", 51, -0.5, 50.5);
  AddOutput(hNGenVtx);      

  hNPU = new TH1D("hNPU", "hNPU", 51, -0.5, 50.5);
  AddOutput(hNPU);      
}

//--------------------------------------------------------------------------------------------------
void GoodPVFilterMod::SlaveTerminate()
{
  // Save number of accepted events.

  SaveNEventsProcessed();
}
