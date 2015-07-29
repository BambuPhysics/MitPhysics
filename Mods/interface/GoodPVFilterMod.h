//--------------------------------------------------------------------------------------------------
// GoodPVFilterMod
//
// This module selects events with a good reconstructed Primary Vertex according to
// configurable cuts on the number of tracks and vertex position.
//
// Authors: J.Bendavid
//--------------------------------------------------------------------------------------------------

#ifndef MITMODS_MODS_GOODPVFILTERMOD_H
#define MITMODS_MODS_GOODPVFILTERMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h"
#include "MitAna/DataTree/interface/Names.h"
#include "MitPhysics/Init/interface/ModNames.h"

#include "TString.h"
#include "TH1F.h"

namespace mithep {

  class GoodPVFilterMod : public BaseMod {
  public:
    enum ECuts {
      eNTracks,
      eNDof,
      eZ,
      eRho,
      nCuts
    };

    GoodPVFilterMod(const char* name = "GoodPVFilterMod", const char* title = "Good PV Filter Module") : BaseMod(name, title) {}
    ~GoodPVFilterMod() {}

    char const* GetOutputName() const { return fGoodVertexesName; }
    char const* GetGoodVertexesName() const { return GetOutputName(); }

    void SetAbortIfNotAccepted(Bool_t b)    { fAbort = b; }
    void SetMinVertexNTracks(UInt_t n)      { fMinVertexNTracks = n; }
    void SetMinNDof(UInt_t n)               { fMinNDof = n; }
    void SetMaxAbsZ(Double_t x)  	    { fMaxAbsZ = x; }
    void SetMaxRho(Double_t x)   	    { fMaxRho = x; }
    void SetInputName(char const* s)        { fVertexesName = s; }
    void SetVertexesName(char const* s)     { SetInputName(s); }
    void SetOutputName(char const* s)       { fGoodVertexesName = s; }
    void SetGoodVertexesName(char const* s) { SetOutputName(s); }

  protected:
    void Process() override;
    void SlaveBegin() override;
    void SlaveTerminate() override;

    virtual void OnAccepted()  {/*could be implemented in derived classes*/}
    virtual void OnFailed()    {/*could be implemented in derived classes*/}

    Bool_t   fAbort{kTRUE};         //=true then abort (sub-)modules if not accepted
    UInt_t   fMinVertexNTracks{0}; //minimum number of tracks for the vertex
    UInt_t   fMinNDof{5};       //minimum number of degrees of freedom
    Double_t fMaxAbsZ{15.};       //maximum abs(z) of the vertex
    Double_t fMaxRho{2.};        //maximum rho of the vertex
    TString  fVertexesName{Names::gkPVBrn};  //Name of PV collection
    TString  fGoodVertexesName{ModNames::gkGoodVertexesName}; //Name of newPV collection

    TH1F* hVertexNTracks = 0;
    TH1F* hVertexNDof = 0;
    TH1F* hVertexRho = 0;
    TH1F* hVertexZ = 0;
    TH1F* hNVtx = 0;
    TH1F* hNGoodVtx = 0;

    ClassDef(GoodPVFilterMod, 1)
  };

}
#endif
