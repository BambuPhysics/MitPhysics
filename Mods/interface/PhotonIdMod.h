//--------------------------------------------------------------------------------------------------
// PhotonIdMod
//
// Authors: Z.Demiragli, Y.Iiyama
//--------------------------------------------------------------------------------------------------


#ifndef MITPHYSICS_MODS_PHOTONIDMOD_H
#define MITPHYSICS_MODS_PHOTONIDMOD_H

#include "MitPhysics/Mods/interface/IdMod.h"
#include "MitPhysics/Utils/interface/PhotonTools.h"
#include "MitAna/DataTree/interface/Photon.h"
#include "MitAna/DataTree/interface/PileupEnergyDensity.h"

namespace mithep {

  class PhotonIdMod : public IdMod<mithep::Photon> {
  public:
    PhotonIdMod(char const* name = "PhotonIdMod", char const* title = "Photon Identification");
    ~PhotonIdMod();

    Bool_t GetApplyTriggerMatching() const   { return fApplyTriggerMatching; }
    Bool_t GetApplyPixelVeto() const         { return (fElectronVeto & (1 << PhotonTools::kPixelVeto)) == 1; }
    Bool_t GetApplyElectronVeto() const      { return (fElectronVeto & (1 << PhotonTools::kElectronVeto)) == 1; }
    Bool_t GetApplyCSafeElectronVeto() const { return (fElectronVeto & (1 << PhotonTools::kCSafeElectronVeto)) == 1; }
    UInt_t GetRhoAlgo() const                { return fRhoAlgo; }

    void SetApplyTriggerMatching(Bool_t b)   { fApplyTriggerMatching = b; }
    void SetApplyPixelVeto(Bool_t b)         { fElectronVeto |= (1 << PhotonTools::kPixelVeto); }
    void SetApplyElectronVeto(Bool_t b)      { fElectronVeto |= (1 << PhotonTools::kElectronVeto); }
    void SetApplyCSafeElectronVeto(Bool_t b) { fElectronVeto |= (1 << PhotonTools::kCSafeElectronVeto); }
    void SetRhoAlgo(UInt_t algo)             { fRhoAlgo = algo; }

  protected:
    enum CutFlow {
      cAll,
      cPt,
      cEta,
      cElectronVeto,
      cTriggerMatching,
      cId,
      cIsolation,
      nCuts
    };

    Bool_t IsGood(mithep::Photon const&) override;
    void IdBegin() override;

    Bool_t PassIdCut(Photon const&);
    Bool_t PassIsolationCut(Photon const&);

    Bool_t   fApplyTriggerMatching = kFALSE;
    UInt_t   fRhoAlgo = mithep::PileupEnergyDensity::kFixedGridFastjetAll;
    BitMask8 fElectronVeto = {};

    ClassDef(PhotonIdMod, 0)
  };

}

#endif
