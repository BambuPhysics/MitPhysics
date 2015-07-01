//--------------------------------------------------------------------------------------------------
// HKFactorProducer
//
// Produces the k factors from LO to NNLO.
//
// Authors: G. Gomez-Ceballos
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_HKFACTORPRODUCER_H
#define MITPHYSICS_MODS_HKFACTORPRODUCER_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/MCEventInfo.h"
#include "MitAna/DataTree/interface/EmbedWeightCol.h"
#include "MitPhysics/Mods/interface/HWWKFactorList.h"
#include "MitPhysics/Mods/interface/PwhgWrapper.h"

class TH1D;
class TH2D;

namespace mithep 
{
  class HKFactorProducer : public BaseMod
  {
    public:
      HKFactorProducer(const char *name="HKFactorProducer", 
                   const char *title="KFactor information module");
      ~HKFactorProducer();

      void               SetProcessID(Int_t d)	           { fProcessID            = d; }
      void               SetInputFilename(const char *s)   { fInputFileName        = s; }
      void               SetMCBosonsName(const char *s)    { fMCBosonsName         = s; }
      void               SetMCEventInfoName(const char *s) { fMCEvInfoName         = s; }
      void               SetEmbedWeightName(const char *s) { fEmbedWeightName      = s; }
      void               SetIsData(Bool_t b)               { fIsData               = b; }	 
      void               SetMakePDFNtuple(Bool_t b)        { fMakePDFNtuple        = b; }	 
      void               SetDoHiggsMhReweighting(Bool_t b) { fDoHiggsMhReweighting = b; }	 
      void               SetOutputName(const char *f)      { fOutputName           = f; }
      void               SetMh(Double_t x)                 { fMh                   = x; }	 
      void               SetWidth(Double_t x)              { fWidth                = x; }	 
      void               SetBWflag(Int_t d)                { fBWflag               = d; }	 

    protected:
      void               Process();
      void               SlaveBegin();
      void               SlaveTerminate();

      Int_t              fProcessID;            //process id (from pythia)
      TString            fInputFileName;        //input file name
      TString            fMCBosonsName;         //boson collection input name
      TString            fMCEvInfoName;         //event info branch name
      TString            fEmbedWeightName;      //tau embedding weight branch name
      Bool_t             fIsData;               //=true then it does nothing (def=0)
      Bool_t             fMakePDFNtuple;        //=true then it does nothing (def=0)
      Bool_t             fDoHiggsMhReweighting; //=true then it does reweighting
      Double_t           fMh;                   //fixed Higgs mass
      Double_t           fWidth;                //fixed Higgs width
      Int_t              fBWflag;               //running or fixed width
      TH1D              *hDHKFactor[10];        //!output histograms
      TFile             *fOutputFile; 	        //output file handle
      TString            fOutputName; 	        //output file name

      float              fTreeVariables[8];     //Ntuple variables
      TTree             *fTree;                 //ntuple tree
      pwhegwrapper       fWeightAlgo;           //Higgs reweighting algorithm

    ClassDef(HKFactorProducer, 1) // Module to produce k factors
  };
}
#endif
