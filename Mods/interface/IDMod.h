//--------------------------------------------------------------------------------------------------
// IDMod
//
// Base class for object ID modules
//
// Authors: S.Xie, C.Loizides
//--------------------------------------------------------------------------------------------------


#ifndef MITPHYSICS_MODS_IDMOD_H
#define MITPHYSICS_MODS_IDMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h"
#include "MitAna/DataCont/interface/NamedFastArrayBasic.h"

#include <limits>

namespace mithep {

  class IDMod : public BaseMod {
  public:
    IDMod(char const* name, char const* title);
    ~IDMod();

    char const* GetInputName() const { return fInputName; }
    char const* GetOutputName() const { return fOutput->GetName(); }
    char const* GetFlagsName() const { return fFlags.GetName(); }
    Bool_t GetIsFilterMode() const { return fIsFilterMode; }
    UInt_t GetIDType() const { return fIDType; }
    UInt_t GetIsoType() const { return fIsoType; }
    Double_t GetPtMin() const { return fPtMin; }
    Double_t GetEtaMax() const { return fEtaMax; }

    void SetInputName(char const* n) { fInputName = n; }
    virtual void SetOutputName(char const*)  = 0; // BaseCollection does not have SetName() interface
    void SetFlagsName(char const* n) { fFlags.SetName(n); }
    void SetIsFilterMode(Bool_t b) { fIsFilterMode = b; }
    void SetIDType(UInt_t t) { fIDType = t; }
    void SetIsoType(UInt_t t) { fIsoType = t; }
    void SetPtMin(Double_t m) { fPtMin = m; }
    void SetEtaMax(Double_t m) { fEtaMax = m; }

  protected:
    Bool_t PublishOutput();
    void RetractOutput();
    
    // fOutput:
    //  Collection of ID'ed objects (filter mode) or all objects (flag mode).
    //  Need to:
    //   new the pointer to the desired output
    //   call PublishOutput() in SlaveBegin
    //   call RetractOutput() in SlaveTerminate
    BaseCollection* fOutput = 0;
    // fFlags: Good / bad flags published in flag mode
    NamedFastArrayBasic<Bool_t> fFlags;

    TString fInputName; // input collection of objects to be ID'ed

    Bool_t fIsFilterMode = kTRUE;
    UInt_t fIDType = 0xffffffff;
    UInt_t fIsoType = 0xffffffff;
    Double_t fPtMin = 0.;
    Double_t fEtaMax = std::numeric_limits<double>::max();

    ClassDef(IDMod, 0)
  };

}

#endif
