//--------------------------------------------------------------------------------------------------
// PhotonIdMod
//
// Authors: Z.Demiragli, Y.Iiyama
//--------------------------------------------------------------------------------------------------


#ifndef MITPHYSICS_MODS_PHOTONIDMOD_H
#define MITPHYSICS_MODS_PHOTONIDMOD_H

#include "MitPhysics/Mods/interface/IdMod.h"
#include "MitAna/DataTree/interface/PhotonCol.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityCol.h"

namespace mithep {

  class PhotonIdMod : public IdMod {
  public:
    PhotonIdMod(char const* name = "PhotonIdMod", char const* title = "Photon Identification");
    ~PhotonIdMod();

    const char* GetTriggerObjectsName() const      { return fAuxInputNames[kTrigObjects]; }
    const char* GetPileupEnergyDensityName() const { return fAuxInputNames[kPileupEnergyDensity]; }

    Bool_t GetApplyTriggerMatching() const { return fApplyTriggerMatching; }
    UInt_t GetRhoAlgo() const              { return fRhoAlgo; }

    void SetTriggerObjectsName(const char* n)      { fAuxInputNames[kTrigObjects] = n; }
    void SetPileupEnergyDensityName(const char* n) { fAuxInputNames[kPileupEnergyDensity] = n; }

    void SetApplyTriggerMatching(Bool_t b) { fApplyTriggerMatching = b; }
    void SetRhoAlgo(UInt_t algo)           { fRhoAlgo = algo; }

  protected:
    void Process() override;
    void IdBegin() override;

    enum AuxInput {
      kTrigObjects,
      kPileupEnergyDensity,
      nAuxInputs
    };

    enum CutFlow {
      cAll,
      cPt,
      cEta,
      cTriggerMatching,
      cId,
      cIsolation,
      nCuts
    };

    Bool_t PassIdCut(Photon const&, TObject const**);
    Bool_t PassIsolationCut(Photon const&, TObject const**);

    template<class T> void GetAuxInput(AuxInput, TObject const**);

    TString fAuxInputNames[nAuxInputs] = {};

    Bool_t fApplyTriggerMatching = kFALSE;
    UInt_t fRhoAlgo = mithep::PileupEnergyDensity::kFixedGridFastjetAll;

    ClassDef(PhotonIdMod, 0)
  };

}

#endif
