//--------------------------------------------------------------------------------------------------
// SeparatePileUpMod
//
// This module applies PFNoPU selection
//
// Authors: G.Ceballos
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_SEPARATEPILEUPMOD_H
#define MITPHYSICS_MODS_SEPARATEPILEUPMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 

namespace mithep 
{
  class SeparatePileUpMod : public BaseMod
  {
    public:
      SeparatePileUpMod(const char *name="SeparatePileUpMod", 
                        const char *title="PFNoPU identification module");

      char const* GetPFCandidatesName() const    { return fPFCandidatesName; }  
      char const* GetPFPileUpName() const        { return fPFPileUpName; }  
      char const* GetPFNoPileUpName() const      { return fPFNoPileUpName; }  
      char const* GetAllVertexName() const       { return fAllVertexName; }
      char const* GetVertexName() const          { return fVertexName; }
      Bool_t      GetCheckClosestZVertex() const { return fCheckClosestZVertex; }
      Bool_t      GetUseAllVerteces() const      { return fUseAllVertices; }

      void SetPFCandidatesName(const char *n) { fPFCandidatesName = n; }  
      void SetPFPileUpName(const char *n)     { fPFPileUpName = n; }  
      void SetPFNoPileUpName(const char *n)   { fPFNoPileUpName = n; }
      void SetAllVertexName(const char *n)    { fAllVertexName = n; }
      void SetVertexName(const char *n)       { fVertexName = n; }
      void SetCheckClosestZVertex(Bool_t b)   { fCheckClosestZVertex = b; }
      void SetUseAllVerteces(Bool_t b)        { fUseAllVertices = b; }

    protected:
      void Process() override;

      TString fPFCandidatesName;    //name of PF collection (input)
      TString fPFPileUpName;        //name of exported PFPileUp collection (output)
      TString fPFNoPileUpName;      //name of exported PFNoPileUp collection (output)
      TString fAllVertexName;       //name of all vertex collection
      TString fVertexName;          //name of good vertex collection
      Bool_t  fCheckClosestZVertex; //boolean to use the closest vertex approach
      Bool_t  fUseAllVertices;      //boolean to use all vertices

    ClassDef(SeparatePileUpMod, 1) // PFNoPU identification module
  };
}
#endif
