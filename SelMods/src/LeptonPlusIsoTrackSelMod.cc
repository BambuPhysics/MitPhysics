// $Id: LeptonPlusIsoTrackSelMod.cc,v 1.3 2009/06/17 14:52:59 loizides Exp $

#include "MitPhysics/SelMods/interface/LeptonPlusIsoTrackSelMod.h"
#include "MitPhysics/Utils/interface/IsolationTools.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataTree/interface/PhotonCol.h" 
#include "MitAna/DataTree/interface/TrackCol.h" 
#include <TH1D.h>

using namespace mithep;

ClassImp(mithep::LeptonPlusIsoTrackSelMod)

//--------------------------------------------------------------------------------------------------
mithep::LeptonPlusIsoTrackSelMod::LeptonPlusIsoTrackSelMod(const char *name, const char *title) : 
  BaseSelMod(name,title),
  fLeptonColName("SetMe"),
  fTrackerTrackColName("SetMe"),
  fGsfTrackColName("SetMe"),
  fLeptonPtMin(0),
  fLeptonPtMax(5000),
  fLeptonEtaMin(-10),
  fLeptonEtaMax(10),
  fTrackPtMin(0),
  fTrackPtMax(5000),
  fTrackEtaMin(-10),
  fTrackEtaMax(10),
  fLeptonCol(0),
  fTrackerTrackCol(0),
  fGsfTrackCol(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void mithep::LeptonPlusIsoTrackSelMod::Process()
{
  // Process entries of the tree.

  //load the track branches
  LoadBranch(GetTrackerTrackColName());
  LoadBranch(GetGsfTrackColName());

  fNAccCounters->Fill(0);

  fLeptonCol = GetObjThisEvt<Collection<Particle> >(GetLeptonColName());
  if (!fLeptonCol) {
    this->SendError(kAbortModule, "Process", 
                    "Could not obtain collection with name %s!", GetLeptonColName());
    return;
  }

  if (!fLeptonCol && !fGsfTrackCol) {
    this->SendError(kAbortModule, "Process", 
                    "Could not obtain either collections with names %s , %s!", 
                    GetTrackerTrackColName(), GetGsfTrackColName());
    return;
  }

  fNAccCounters->Fill(1);

  UInt_t LeptonCounter = 0;
  for (UInt_t i=0;i<fLeptonCol->GetEntries();++i) {
    if (fLeptonCol->At(i)->Pt() >= fLeptonPtMin   && 
        fLeptonCol->At(i)->Pt() <= fLeptonPtMax   &&
        fLeptonCol->At(i)->Eta() >= fLeptonEtaMin && 
        fLeptonCol->At(i)->Eta() <= fLeptonEtaMax)
      LeptonCounter++;
  }
  if (LeptonCounter == 0) {
    this->SkipEvent();
    return;
  }

  fNAccCounters->Fill(2);

  UInt_t TrackCounter = 0;
  for(UInt_t i=0;i<fTrackerTrackCol->GetEntries();++i) {
    const Track *trk = fTrackerTrackCol->At(i);
    if (trk->Pt() >= fTrackPtMin   && 
        trk->Pt() <= fTrackPtMax   &&
        trk->Eta() >= fTrackEtaMin && 
        trk->Eta() <= fTrackEtaMax) {
      Double_t iso = IsolationTools::TrackIsolation(trk,0.3, 0.015,1.0,1000.0,fTrackerTrackCol);
      if (iso < 10.0) {
        //require that the track is not the same object as one of the leptons
        if (MathUtils::DeltaR(trk->Phi(), trk->Eta(),
                              fLeptonCol->At(0)->Phi(), 
                              fLeptonCol->At(0)->Eta()) >= 0.3)
          TrackCounter++;
      }        
    }
  }
  for(UInt_t i=0;i<fGsfTrackCol->GetEntries();++i) {
    const Track *trk = fGsfTrackCol->At(i);
    if (trk->Pt() >= fTrackPtMin   && 
        trk->Pt() <= fTrackPtMax   &&
        trk->Eta() >= fTrackEtaMin && 
        trk->Eta() <= fTrackEtaMax) {
      Double_t iso = IsolationTools::TrackIsolation(trk,0.3, 0.015,1.0,1000.0,fTrackerTrackCol);
      if (iso < 10.0) {
        if (MathUtils::DeltaR(trk->Phi(), trk->Eta(),
                              fLeptonCol->At(0)->Phi(), fLeptonCol->At(0)->Eta()) >= 0.3)
          TrackCounter++;
      }
    }
  }

  if (TrackCounter == 0) {
    this->SkipEvent();
    return;
  }

  fNAccCounters->Fill(3);
}

//--------------------------------------------------------------------------------------------------
void mithep::LeptonPlusIsoTrackSelMod::SlaveBegin()
{
  // Setup acceptence histogram.
  ReqBranch(GetTrackerTrackColName(),          fTrackerTrackCol);
  ReqBranch(GetGsfTrackColName(),              fGsfTrackCol);

  AddTH1(fNAccCounters,"hNAccCounters",";cut;#",5,-0.5,4.5);
  if (1) {
    TAxis *xa = fNAccCounters->GetXaxis();
    for(Int_t i=1;i<=fNAccCounters->GetNbinsX();++i)
      xa->SetBinLabel(i,"unused");
    xa->SetBinLabel(1,"Enter");
    xa->SetBinLabel(2,"Objs");
    xa->SetBinLabel(3,"AtLeastOneLepton");
    xa->SetBinLabel(4,"IsolatedTrack");
    xa->SetRangeUser(0,3);
  }
}
