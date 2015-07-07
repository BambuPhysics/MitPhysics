//--------------------------------------------------------------------------------------------------
// MetCorrectionMod
//
// This module applies Type0/1 and XY shift MET corrections on the fly at analysis level.
// The methods are synchronized with JetMET POG 2012 studies, documented in this twiki
// https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMetAnalysis#7_7_6_MET_Corrections
//
// Authors: L.Di Matteo
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_METCORRECTIONMOD_H
#define MITPHYSICS_MODS_METCORRECTIONMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h"
#include "MitPhysics/Utils/interface/JetCorrector.h"
#include "MitAna/DataTree/interface/MetCol.h"
#include "MitAna/DataTree/interface/PileupEnergyDensity.h"

#include "MitAna/DataTree/interface/Names.h"
#include "MitPhysics/Init/interface/ModNames.h"

#include "TFormula.h"

namespace mithep {

  // (Get|Set)ExprShift(Data|MC)P(x|y) is removed. Use IsData flag + (Get|Set)ExprShiftP(x|y) instead.

  class MetCorrectionMod : public BaseMod {
  public:
    MetCorrectionMod(const char* name="MetCorrectionMod",
                     const char* title="Met correction module");
    ~MetCorrectionMod() {}

    const char*   GetInputName() const       { return fMetName; }
    const char*   GetOutputName() const      { return fOutput.GetName(); }
    const char*   GetCorrectedName() const   { return GetOutputName(); }
    const char*   GetJetsName() const        { return fJetsName; }
    Double_t      GetMinDz() const           { return fMinDz; }
    const char*   GetExprType0();
    const char*   GetExprShiftPx();
    const char*   GetExprShiftPy();
    JetCorrector* GetJetCorrector() const    { return fJetCorrector; }
    UInt_t        GetRhoAlgo() const         { return fRhoAlgo; }
    Double_t      GetMaxEMFraction() const   { return fMaxEMFraction; }
    Bool_t        GetSkipMuons() const       { return fSkipMuons; }

    void SetInputName(const char *name)        { fMetName = name; }
    void SetOutputName(const char *name)       { fOutput.SetName(name); }
    void SetCorrectedName(const char *name)    { SetOutputName(name); }
    void SetJetsName(const char *name)         { fJetsName = name; }
    void SetMinDz(Double_t d)                  { fMinDz = d; }
    void ApplyType0(bool b)                    { fApplyType0 = b; }
    void ApplyType1(bool b)                    { fApplyType1 = b; }
    void ApplyShift(bool b)                    { fApplyShift = b; }
    void SetExprType0(const char *expr)        { MakeFormula(0, expr); }
    void SetPFCandidatesName(const char *name) { fPFCandidatesName = name; }
    void SetExprShiftPx(const char *expr)  { MakeFormula(1, expr); }
    void SetExprShiftPy(const char *expr)  { MakeFormula(2, expr); }
    void AddJetCorrectionFromFile(char const* file);
    void SetJetCorrector(JetCorrector*);
    void SetRhoAlgo(UInt_t a)                  { fRhoAlgo = a; }
    void SetMaxEMFraction(Double_t m)          { fMaxEMFraction = m; }
    void SetSkipMuons(Bool_t s)                { fSkipMuons = s; }
    void IsData(bool b)                        { fIsData = b; }
    void SetPrint(bool b)                      { fPrint = b; }

  protected:
    void SlaveBegin() override;
    void SlaveTerminate() override;
    void Process() override;

    void MakeJetCorrector();
    void MakeFormula(UInt_t idx, char const* expr = "");

    TString       fMetName{"PFMet"};                           //name of met collection (input)
    TString       fJetsName{Names::gkPFJetBrn};                //name of uncorrected jet collection (input)
    TString       fPFCandidatesName{Names::gkPFCandidatesBrn}; //name of PF candidates collection (input)
    TString       fVertexName{ModNames::gkGoodVertexesName};   //name of vertex collection (input)
    Double_t      fMinDz = 0.2;                                //delta Z for separating PU verteces from PV
                  
    Bool_t        fApplyType0 = kTRUE; //switch on type 0 correction
    Bool_t        fApplyType1 = kTRUE; //switch on type 1 correction
    Bool_t        fApplyShift = kTRUE; //switch on XY shift correction
                  
    TFormula*     fFormulaType0 = 0;
    TFormula*     fFormulaShiftPx = 0;
    TFormula*     fFormulaShiftPy = 0;

    Bool_t        fOwnJetCorrector = kFALSE;
    JetCorrector* fJetCorrector = 0; //For type 1 correction. Can be created internally or set externally
    UInt_t        fRhoAlgo = PileupEnergyDensity::kFixedGridFastjetAll;

    Double_t      fMaxEMFraction = -1.; //maximum charged + neutral EM energy fraction
                                        //for jet to be in Type1 corr (< 0 -> does not skip jets)
    Bool_t        fSkipMuons = kFALSE; //remove muon P4 from jet P4 in Type 1 corr
                  
    Bool_t        fIsData = kTRUE; //flag for data/MC distinction
    Bool_t        fPrint = kFALSE; //flag for debug print out

    MetOArr       fOutput{1, "PFMetT0T1Shift"}; //using ObjArray to accomodate different MET types

    ClassDef(MetCorrectionMod, 1) // met correction module
  };

}

inline
char const*
mithep::MetCorrectionMod::GetExprType0()
{
  if (!fFormulaType0)
    MakeFormula(0);
  return fFormulaType0->GetExpFormula();
}

inline
char const*
mithep::MetCorrectionMod::GetExprShiftPx()
{
  if (!fFormulaShiftPx)
    MakeFormula(1);
  return fFormulaShiftPx->GetExpFormula();
}

inline
char const*
mithep::MetCorrectionMod::GetExprShiftPy()
{
  if (!fFormulaShiftPy)
    MakeFormula(2);
  return fFormulaShiftPy->GetExpFormula();
}

#endif
