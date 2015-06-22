//--------------------------------------------------------------------------------------------------
// IDMod
//
// Base class for object ID modules
//
// Authors: Y.Iiyama
//--------------------------------------------------------------------------------------------------


#ifndef MITPHYSICS_MODS_IDMOD_H
#define MITPHYSICS_MODS_IDMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h"
#include "MitAna/DataCont/interface/NamedFastArrayBasic.h"

#include "TString.h"
#include "TH1D.h"

#include <limits>

namespace mithep {

  class IDMod : public BaseMod {
  public:
    IDMod(char const* name, char const* title);
    ~IDMod();

    char const* GetInputName() const { return fInputName; }
    char const* GetOutputName() const;
    Bool_t GetIsFilterMode() const { return fIsFilterMode; }
    UInt_t GetIDType() const { return fIDType; }
    UInt_t GetIsoType() const { return fIsoType; }
    Double_t GetPtMin() const { return fPtMin; }
    Double_t GetEtaMax() const { return fEtaMax; }

    void SetInputName(char const* n) { fInputName = n; }
    void SetOutputName(char const* n);
    void SetIsFilterMode(Bool_t b) { fIsFilterMode = b; }
    void SetIDType(UInt_t t) { fIDType = t; }
    void SetIsoType(UInt_t t) { fIsoType = t; }
    void SetPtMin(Double_t m) { fPtMin = m; }
    void SetEtaMax(Double_t m) { fEtaMax = m; }

  protected:
    void SlaveBegin() override;
    void SlaveTerminate() override;
    virtual void IdBegin() {}
    virtual void IdTerminate() {}
    
    // fOutput:
    //  Collection of ID'ed objects (filter mode) or all objects (flag mode).
    //  A derived class must new the pointer to the desired output in the constructor
    BaseCollection* fOutput = 0;
    // fFlags: Good / bad flags published in flag mode
    NamedFastArrayBasic<Bool_t> fFlags;

    TH1D* fCutFlow = 0;

    TString fInputName = ""; // input collection of objects to be ID'ed

    Bool_t fIsFilterMode = kTRUE;
    UInt_t fIDType = 0xffffffff;
    UInt_t fIsoType = 0xffffffff;
    Double_t fPtMin = 0.;
    Double_t fEtaMax = std::numeric_limits<double>::max();

    ClassDef(IDMod, 0)
  };

  inline
  char const*
  mithep::IDMod::GetOutputName() const
  {
    if (fIsFilterMode)
      return fOutput->GetName();
    else
      return fFlags.GetName();
  }

  inline
  void
  mithep::IDMod::SetOutputName(char const* n)
  {
    if (fIsFilterMode)
      return fOutput->SetName(n);
    else
      return fFlags.SetName(n);
  }
}

#endif
