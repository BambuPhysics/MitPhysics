//--------------------------------------------------------------------------------------------------
// $Id: PuppiJetMod.h,v 1.9 2011/03/01 17:27:22 mzanetti Exp $
//
// PuppiJetMod
//
// This mod processes an input collection of particles (i.e. from PuppiMod) and returns
// a collection of jets (of type JETTYPE)
//
// Authors:S.Narayanan
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_PUPPIJETMOD_H
#define MITPHYSICS_MODS_PUPPIJETMOD_H

#include <TVector2.h>

#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/ClusterSequenceArea.hh"

#include "MitAna/TreeMod/interface/BaseMod.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/FatJetCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/PFJetCol.h"

#include "MitCommon/DataFormats/interface/Vect4M.h"
#include "MitCommon/DataFormats/interface/Vect3.h"
#include "MitCommon/DataFormats/interface/Types.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataTree/interface/Names.h"



#include "TStopwatch.h"

namespace mithep
{
  template<typename JETTYPE>
   class PuppiJetMod : public BaseMod
  {
    public:
      enum JetAlgorithm {
        kAntiKT,
        kCambridgeAachen,
        kKT, // why?
        nJetAlgorithm
      };

      PuppiJetMod(const char *name = "PuppiJetMod",
                   const char *title = "Puppi jet module");
      ~PuppiJetMod();

      void IsData(Bool_t b)                { fIsData = b;           }
      void PublishOutput(Bool_t b)         { fPublishOutput = b;    }
      
      void SetProcessNJets(UInt_t n)       { fProcessNJets = n;     }

      void SetInputName(const char *n)      { fPFCandidatesName = n;         }

      void SetOutputName(const char *n)   { fJetsName = n;    }
      void SetR0(double d)                { fR0 = d;  }
      void SetBeVerbose(Bool_t b)          { fBeVerbose = b;  }
      void SetDebugFlag(int i)   { fDebugFlag = i; }

      void SetDoMatching(Bool_t b)              { fDoMatching = b;  }
      void SetMatchingJetsName(const char *n)   { fMatchingJetsName = n;  }

      void SetJetAlgorithm(JetAlgorithm j)      { fJetAlgorithm = j; }

      const char *GetOutputName()               { return fJetsName; }

    protected:
      void Process();
      void SlaveBegin();
      void SlaveTerminate();

      // Jet collection helpers
      std::vector <fastjet::PseudoJet>
            Sorted_by_pt_min_pt(std::vector <fastjet::PseudoJet> &jets,
                                float jetPtMin);

      Vect4M GetCorrectedMomentum(fastjet::PseudoJet fj_tmp, double thisJEC);

      void RunMatching(PFJet*);
      void RunMatching(FatJet*);


    private:
      Bool_t fIsData;                      //is this data or MC?
      Bool_t fPublishOutput;               //=true if output collection are published

      TString fPFCandidatesName;           //(i) name of PF candidates coll
      const PFCandidateCol *fPFCandidates; //particle flow candidates coll handle

      TString fJetsName;              //name of output jets collection
      Array<JETTYPE> *fJets;             //array of jets

      // Objects from fastjet we want to use
      double fR0;                    //fastjet clustering radius
      fastjet::JetDefinition *fJetDef;   //fastjet clustering definition
      fastjet::GhostedAreaSpec *fActiveArea;
      fastjet::AreaDefinition *fAreaDefinition;

      Bool_t fDoMatching;
      TString fMatchingJetsName;
      const Collection<Jet> *fMatchingJets;

      UInt_t fProcessNJets;

      Bool_t fBeVerbose;

      JetAlgorithm fJetAlgorithm;
      
      TStopwatch *fStopwatch;

      int fDebugFlag = -1;

      ClassDef(PuppiJetMod, 0)         //Puppi jets filler
  };

typedef PuppiJetMod<FatJet> PuppiFatJetMod;
typedef PuppiJetMod<PFJet>  PuppiPFJetMod;

#include "MitPhysics/Mods/src/PuppiJetMod.icc"

}


#endif
