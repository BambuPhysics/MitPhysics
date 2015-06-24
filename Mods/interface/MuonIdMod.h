//--------------------------------------------------------------------------------------------------
// MuonIdMod
//
// Authors: P.Dominguez, Y.Iiyama
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_MUONIDMOD_H
#define MITPHYSICS_MODS_MUONIDMOD_H

#include "MitPhysics/Mods/interface/IdMod.h"
#include "MitAna/DataTree/interface/Muon.h"

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

    UInt_t fMuonClassType;
    Bool_t fApplyD0Cut;          //=true then apply d0 cut (def=1)
    Bool_t fApplyDZCut;          //=true then apply dz cut (def=1)
    Int_t  fWhichVertex;         //vertex to use (-2: beamspot, -1: closest in Z)
    UInt_t fRhoAlgo;

    ClassDef(MuonIdMod, 0)
  };

}

#endif
