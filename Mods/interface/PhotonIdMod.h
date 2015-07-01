//--------------------------------------------------------------------------------------------------
// PhotonIdMod
//
// Authors: Z.Demiragli, Y.Iiyama
//--------------------------------------------------------------------------------------------------


#ifndef MITPHYSICS_MODS_PHOTONIDMOD_H
#define MITPHYSICS_MODS_PHOTONIDMOD_H

#include "MitPhysics/Mods/interface/IdMod.h"
#include "MitAna/DataTree/interface/Photon.h"
#include "MitAna/DataTree/interface/PileupEnergyDensity.h"

namespace mithep {

  class PhotonIdMod : public IdMod<mithep::Photon> {
  public:
    PhotonIdMod(char const* name = "PhotonIdMod", char const* title = "Photon Identification");
    ~PhotonIdMod();

    Bool_t GetApplyTriggerMatching() const { return fApplyTriggerMatching; }
    UInt_t GetRhoAlgo() const              { return fRhoAlgo; }

    void SetApplyTriggerMatching(Bool_t b) { fApplyTriggerMatching = b; }
    void SetRhoAlgo(UInt_t algo)           { fRhoAlgo = algo; }

  protected:
    enum CutFlow {
      cAll,
      cPt,
      cEta,
      cTriggerMatching,
      cId,
      cIsolation,
      nCuts
    };

    Bool_t IsGood(mithep::Photon const&) override;
    void IdBegin() override;

    Bool_t PassIdCut(Photon const&);
    Bool_t PassIsolationCut(Photon const&);

    Bool_t fApplyTriggerMatching = kFALSE;
    UInt_t fRhoAlgo = mithep::PileupEnergyDensity::kFixedGridFastjetAll;

    ClassDef(PhotonIdMod, 0)
  };

}

#endif
