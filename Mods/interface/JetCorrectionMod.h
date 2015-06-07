//--------------------------------------------------------------------------------------------------
// JetCorrectionMod
//
// This module applies jet energy corrections on the fly at analysis level, using directly the
// FWLite oriented classes from CMSSW and the jet correction txt files from the release and/or
// checked out tags.  At the moment this is hardcoded to apply or re-apply L2+L3 corrections,
// but it is intended that this will become more configurable in the future.
//
// Authors: J.Bendavid
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_JETCORRECTIONMOD_H
#define MITPHYSICS_MODS_JETCORRECTIONMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/Jet.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataCont/interface/Types.h"

#include "MitPhysics/Utils/interface/RhoUtilities.h"

class FactorizedJetCorrector;
namespace mithep 
{
  class JetCorrectionMod : public BaseMod
  {
    public:
      JetCorrectionMod(const char *name="JetCorrectionMod", 
		       const char *title="Jet correction module");
     ~JetCorrectionMod();

      const char       *GetInputName()                 const      { return fJetsName;               }   
      const char       *GetCorrectedName()             const      { return GetCorrectedJetsName();  }     
      const char       *GetCorrectedJetsName()         const      { return fCorrectedJetsName;      }     
      const char       *GetOutputName()                const      { return GetCorrectedJetsName();  }
      void              AddCorrectionFromRelease(const std::string &path);
      void              AddCorrectionFromFile(const std::string &file);    
      void              ApplyL1FastJetCorrection(float maxEta=2.5, bool useFixedGrid=false); 
      void              ApplyL1FastJetCorrection(Jet * jet);
      void              SetCorrectedJetsName(const char *name)    { fCorrectedJetsName = name;      }     
      void              SetCorrectedName(const char *name)        { SetCorrectedJetsName(name);     }    
      void              SetInputName(const char *name)            { fJetsName = name;               }
      void              SetRhoType(RhoUtilities::RhoType); /*DEPRECATED*/
      void              SetRhoAlgo(mithep::PileupEnergyDensity::Algo algo) { fRhoAlgo = algo; }

    protected:
      void              SlaveBegin();
      void              Process();

      TString           fJetsName;              //name of jet collection (input)
      TString           fCorrectedJetsName;     //name of good jets collection (output)
      TString           fRhoBranchName;         //name of pileup energy density collection
      TString           fPFCandidatesName;      //name of PF candidates colleciont
      bool              fEnabledL1Correction; //switch on L1 correction
      float             rhoEtaMax; //parameter to choose which rho to use for L1 correction
      FactorizedJetCorrector *fJetCorrector;      //!CMSSW/FWLite jet corrections module
      TString           fEvtHdrName;	          // name of event header branch
      const EventHeader *fEventHeader;            // event header for current event

      std::vector<std::string> fCorrectionFiles;   //list of jet correction files
      const PileupEnergyDensityCol *fRho;          // collection of pileup energy density collection
      const PFCandidateCol         *fPFCandidates; // particle flow candidates collection handle
      
      BitMask8          fEnabledCorrectionMask;    //bitmask of enabled corrections
      std::vector<Jet::ECorr> fEnabledCorrections; //vector of enabled corrections

      bool              fUseFixedGrid;             // flag to use fixed grid method to compute energy density 

      mithep::PileupEnergyDensity::Algo fRhoAlgo;

      ClassDef(JetCorrectionMod, 2) // Jet identification module
  };
}
#endif
