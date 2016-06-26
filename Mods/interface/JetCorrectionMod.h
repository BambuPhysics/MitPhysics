//--------------------------------------------------------------------------------------------------
// JetCorrectionMod
//
// This module applies jet energy corrections on the fly at analysis level, using directly the
// FWLite oriented classes from CMSSW and the jet correction txt files from the release and/or
// checked out tags.
//
// Authors: J.Bendavid, Y.Iiyama
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_JETCORRECTIONMOD_H
#define MITPHYSICS_MODS_JETCORRECTIONMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h"
#include "MitAna/DataTree/interface/JetFwd.h"
#include "MitAna/DataTree/interface/Jet.h"
#include "MitAna/DataTree/interface/PileupEnergyDensity.h"
#include "MitAna/DataTree/interface/Names.h"

#include "MitPhysics/Utils/interface/JetCorrector.h"
#include "MitPhysics/Utils/interface/RhoUtilities.h"
#include "MitPhysics/Init/interface/ModNames.h"

class FactorizedJetCorrector;
namespace mithep
{
  class JetCorrectionMod : public BaseMod
  {
    public:
      JetCorrectionMod(const char *name="JetCorrectionMod",
		       const char *title="Jet correction module");
      ~JetCorrectionMod() {}

      char const*   GetInputName() const         { return fJetsName; }
      char const*   GetCorrectedName() const     { return GetCorrectedJetsName(); }
      char const*   GetCorrectedJetsName() const { return fCorrectedJets.GetName(); }
      char const*   GetOutputName() const        { return GetCorrectedJetsName(); }
      JetCorrector* GetCorrector() const         { return fCorrector; }

      void AddCorrectionFromFile(char const* file, JetCorrector::Corrector::FactorType = JetCorrector::Corrector::nFactorTypes);
      void SetCorrectedJetsName(char const* name)    { fCorrectedJets.SetName(name); }
      void SetCorrectedName(char const* name)        { SetCorrectedJetsName(name); }
      void SetInputName(char const* name)            { fJetsName = name; }
      void SetRhoType(RhoUtilities::RhoType); /*DEPRECATED*/
      void SetRhoAlgo(unsigned algo)                 { fRhoAlgo = algo; }
      void SetGenJetsName(char const* name)          { fGenJetsName = name; }
      void SetUncertaintySigma(Double_t s)           { fSigma = s; }
      void SetSmearingSeed(Int_t s)                  { fSmearingSeed = s; }
      void SetCorrector(JetCorrector*);

    protected:
      void SlaveBegin() override;
      void SlaveTerminate() override;
      void Process() override;

      void MakeCorrector();

      TString       fJetsName{Names::gkPFJetBrn}; //name of jet collection (input))
      TString       fRhoBranchName{Names::gkPileupEnergyDensityBrn}; //name of pileup energy density collection
      UInt_t        fRhoAlgo = PileupEnergyDensity::kFixedGridFastjetAll;
      TString       fGenJetsName{""};
      Double_t      fSigma = 0.; //shift JES by fSigma * Uncertainty (needs uncertainty source)
      Int_t         fSmearingSeed = 1234;
      Bool_t        fOwnCorrector = kFALSE;
      JetCorrector* fCorrector = 0;

      mithep::JetOArr fCorrectedJets{32, ModNames::gkCorrectedJetsName}; //output collection. Using ObjArray to allow output of derived classes

      ClassDef(JetCorrectionMod, 2) // Jet identification module
  };
}
#endif
