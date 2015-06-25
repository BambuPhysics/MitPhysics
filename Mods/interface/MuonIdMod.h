//--------------------------------------------------------------------------------------------------
// MuonIdMod
//
// Authors: P.Dominguez, Y.Iiyama
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_MUONIDMOD_H
#define MITPHYSICS_MODS_MUONIDMOD_H

#include "MitPhysics/Mods/interface/IdMod.h"
#include "MitPhysics/Utils/interface/MuonTools.h"
#include "MitAna/DataTree/interface/Muon.h"
#include "MitAna/DataTree/interface/PileupEnergyDensity.h"

namespace mithep {

  class MuonIdMod : public IdMod<mithep::Muon> {
  public:
    MuonIdMod(char const* name = "MuonIdMod",
	      char const* title = "Muon Identification");

    ~MuonIdMod();

    UInt_t GetMuonClassType() const { return fMuonClassType; }
    Bool_t GetApplyD0Cut() const    { return fApplyD0Cut; }
    Bool_t GetApplyDZCut() const    { return fApplyDZCut; }
    Int_t  GetWhichVertex() const   { return fWhichVertex; }
    UInt_t GetRhoAlgo() const       { return fRhoAlgo; }

    void SetMuonClassType(UInt_t t) { fMuonClassType = t; }
    void SetApplyD0Cut(Bool_t b)    { fApplyD0Cut = b; }
    void SetApplyDZCut(Bool_t b)    { fApplyDZCut = b; }
    void SetWhichVertex(Int_t d)    { fWhichVertex = d; }
    void SetRhoAlgo(UInt_t algo)    { fRhoAlgo = algo; }

  protected:
    enum CutFlow {
      cAll,
      cMuonClass,
      cPt,
      cEta,
      cD0,
      cDZ,
      cId,
      cIsolation,
      nCuts
    };

    Bool_t IsGood(mithep::Muon const&) override;
    void IdBegin() override;

    Bool_t PassIsolation(mithep::Muon const&);

    UInt_t fMuonClassType = mithep::MuonTools::kGlobal;
    Bool_t fApplyD0Cut = kTRUE;
    Bool_t fApplyDZCut = kTRUE;
    Int_t  fWhichVertex = -1;
    UInt_t fRhoAlgo = mithep::PileupEnergyDensity::kFixedGridFastjetAll;

    ClassDef(MuonIdMod, 0)
  };

}

#endif
