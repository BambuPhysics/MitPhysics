//--------------------------------------------------------------------------------------------------
// JetIdMod
//
// This module applies jet identification criteria and exports a pointer to a collection
// of "good jet" according to the specified identification scheme.
//
// Authors: S.Xie, Y.Iiyama, S.Narayanan
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_JETIDMOD_H
#define MITPHYSICS_MODS_JETIDMOD_H
#include "MitPhysics/Mods/interface/IdMod.h"

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitPhysics/Utils/interface/JetIDMVA.h"


namespace mithep {
  class JetIdMVA; 

  class JetIdMod : public IdMod {
  public:
    JetIdMod(const char* name="JetIdMod",
                  const char* title="Jet identification module");

    const char* GetVertexName()                     { return fAuxInputNames[kVertices]; }

    Double_t GetPtCut()                             { return GetPtMin(); }
    Double_t GetEtaMaxCut()                         { return GetEtaMax(); }
    Bool_t   GetUseCorrection() const               { return fUseJetCorrection; }
    Double_t GetJetEEMFractionMinCut()  const       { return fJetEEMFractionMinCut;}
    Double_t   GetEtaMaxCut()                 const { return fJetEtaMaxCut;        }
    Bool_t     GetApplyBeta()              const    { return fApplyBetaCut;        }
    Bool_t     GetApplyMVA()               const    { return fApplyMVACut;         }
    Bool_t     GetApplyMVACHS()               const { return fApplyMVACHS;         }
    Bool_t     GetApplyPFLooseId()          const   { return fApplyPFLooseId; }


    void SetVertexName(const char *s)              { fAuxInputNames[kVertices] = s; }
    void              SetMinNJets(UInt_t n)                { fMinNJets = n;              }
    void              SetPtCut(Double_t cut)               { fJetPtCut = cut;            }
    void              SetUseCorrection(Bool_t b)           { fUseJetCorrection = b;      }
    void              SetEtaMaxCut(Double_t cut)           { fJetEtaMaxCut = cut;        }
    void              SetJetEEMFractionMinCut(Double_t cut){ fJetEEMFractionMinCut = cut;}
    void              SetChargedHadronFractionRange(Double_t min, Double_t max)
                                                           { fMinChargedHadronFraction = min; fMaxChargedHadronFraction = max; }
    void              SetNeutralHadronFraction(Double_t min, Double_t max)
                                                           { fMinNeutralHadronFraction = min; fMaxNeutralHadronFraction = max; }
    void              SetChargedEMFraction(Double_t min, Double_t max)
                                                           { fMinChargedEMFraction = min; fMaxChargedEMFraction = max; }
    void              SetNeutralEMFraction(Double_t min, Double_t max)
                                                           { fMinNeutralEMFraction = min; fMaxNeutralEMFraction = max; }
    void              SetApplyBetaCut(Bool_t b)            { fApplyBetaCut = b;          }
    void              SetApplyPFLooseId(Bool_t b)          { fApplyPFLooseId = b;        }
    void              SetApplyMVACut(Bool_t b)             { fApplyMVACut = b;           }
    void              SetApplyMVACHS(Bool_t b)             { fApplyMVACHS = b;           }

  protected:
    enum AuxInput {
      kVertices,
      nAuxInputs
    };

    enum CutFlow {
      cAll,
      cEta,
      cPt,
      cEEMFraction,
      cBeta,
      cPFLooseId,
      cMVA,
      cChargedHFrac,
      cNeutralHFrac,
      cChargedEMFrac,
      cNeutralEMFrac,
      nCuts
    };

    template<class T> void GetAuxInput(AuxInput, TObject const**);

    void Process() override;
    void IdBegin() override;
    
    TString  fAuxInputNames[nAuxInputs];

    Bool_t            fUseJetCorrection = kTRUE;        //=true then use corrected energy
    UInt_t            fMinNJets;              //minimum required number of jets passing the ID
    Double_t          fJetPtCut;              //jet pt cut
    Double_t          fJetEtaMaxCut;          //jet eta max cut
    Double_t          fJetEEMFractionMinCut = 0.01;  //jet Eem fraction min cut
    Double_t          fMinChargedHadronFraction = 0;
    Double_t          fMaxChargedHadronFraction = 1;
    Double_t          fMinNeutralHadronFraction = 0;
    Double_t          fMaxNeutralHadronFraction = 1;
    Double_t          fMinChargedEMFraction = 0;
    Double_t          fMaxChargedEMFraction = 1;
    Double_t          fMinNeutralEMFraction = 0;
    Double_t          fMaxNeutralEMFraction = 1;
    Bool_t            fApplyBetaCut = kFALSE;          //=true then apply beta cut
    Bool_t            fApplyPFLooseId = kFALSE;        //=true then apply PF loose ID
    Bool_t            fApplyMVACut = kFALSE;           //=true then apply MVA cut
    Bool_t            fApplyMVACHS = kFALSE;           //=true then apply MVA for CHS

    JetIDMVA         *fJetIDMVA;

    ClassDef(JetIdMod, 0) // Jet identification module
  };
}

#endif
