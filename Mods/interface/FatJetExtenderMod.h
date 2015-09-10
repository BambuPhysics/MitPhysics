// FatJetExtender
//
// This module processes a collection of input FatJets, compute the substrucure
// and fill a output collections of fXlFatJets
//
// Authors: L.DiMatteo, S.Narayanan
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_FATJETEXTENDER_H
#define MITPHYSICS_MODS_FATJETEXTENDER_H

#include <TVector2.h>

#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/contrib/Njettiness.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/contrib/EnergyCorrelator.hh"

#include "MitAna/DataTree/interface/XlFatJetFwd.h"
#include "MitAna/DataTree/interface/XlFatJet.h"
#include "MitAna/DataTree/interface/XlSubJetFwd.h"

#include "MitAna/TreeMod/interface/BaseMod.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/FatJetCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/PFJet.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityCol.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/PhysicsUtils/interface/QGTagger.h"

#include "MitPhysics/SDAlgorithm/interface/AnalysisParameters.h"
#include "MitPhysics/SDAlgorithm/interface/Exception.h"
#include "MitPhysics/SDAlgorithm/interface/Message.h"
#include "MitPhysics/SDAlgorithm/interface/TopGluonModel.h"
#include "MitPhysics/SDAlgorithm/interface/BackgroundModel.h"
#include "MitPhysics/SDAlgorithm/interface/ISRModel.h"
#include "MitPhysics/SDAlgorithm/interface/Deconstruct.h"
#include "MitPhysics/SDAlgorithm/interface/ParseUtils.h"
#include "MitAna/PhysicsUtils/interface/CMSTopTagger.h"
#include "TStopwatch.h"

namespace mithep
{
  class FatJetExtenderMod : public BaseMod
  {
    public:
      typedef XlSubJet::ESubJetType ESubJetType;
      FatJetExtenderMod(const char *name = "FatJetExtenderMod",
                   const char *title = "XlFatJets Filler module");
      ~FatJetExtenderMod();

      void IsData(Bool_t b)                { fIsData = b;           }
      void SetQGTaggingOn(Bool_t b)        { fQGTaggingActive = b;  }
      void SetSubJetTypeOn(ESubJetType t) { fSubJetFlags |= (1<<t); }
      void SetSubJetTypeOff(ESubJetType t) { fSubJetFlags &= ~(1<<t); }
      void SetQGTaggerCHS(bool b)          { fQGTaggerCHS = b;      }
      void PublishOutput(Bool_t b)         { fPublishOutput = b;    }

      void SetProcessNJets(UInt_t n)       { fProcessNJets = n;     }

      void SetInputName(const char *n)      { fFatJetsName = n;         }
      void SetInputFromBranch(Bool_t b)     { fFatJetsFromBranch = b;   }

      void SetOutputName(const char *n)   { fXlFatJetsName = n;    }
      const char * GetOutputName()        { return fXlFatJetsName;  }

      void SetSoftDropZCut(double d)       { fSoftDropZCut = d;     }
      void SetSoftDropR0(double d)      { fSoftDropR0 = d;    }
      void SetPruneZCut(double d)          { fPruneZCut = d;        }
      void SetPruneDistCut(double d)       { fPruneDistCut = d;     }
      void SetFilterN(int n)               { fFilterN = n;          }
      void SetFilterRad(double d)          { fFilterRad = d;        }
      void SetTrimRad(double d)            { fTrimRad = d;          }
      void SetTrimPtFrac(double d)         { fTrimPtFrac = d;       }
      void SetConeSize(double d)           { fConeSize = d;         }
      void SetPFCandsName(const char *n)   { fPFCandidatesName = n; }
      void SetPileUpDenName(const char *n) { fPileUpDenName = n;    }
      void SetVertexesName(const char *n)  { fVertexesName = n;     }
      void SetDoShowerDeconstruction(Bool_t b) { fDoShowerDeconstruction = b; }
      void SetBeVerbose(Bool_t b)          { fBeVerbose = b;  }
      void SetDoECF(Bool_t b)              { fDoECF = b; }
      void SetDoQjets(Bool_t b)            { fDoQjets = b; }
      void SetNMaxMicrojets(unsigned int n)         { fNMaxMicrojets = n; }
      void SetDebugFlag(int i)   { fDebugFlag = i; }
      void SetSDInputCard(const char *s)   { fInputCard = s;  }
    protected:
      typedef std::vector<fastjet::PseudoJet> VPseudoJet;
      typedef std::vector<fastjet::PseudoJet const*> VPseudoJetPtr;

      void Process();
      void SlaveBegin();
      void SlaveTerminate();

      void FillXlFatJet (const FatJet *fatJet);
      void FillXlSubJets(VPseudoJetPtr const& fjSubJets, XlFatJet* pFatJet, XlSubJet::ESubJetType subJetType);

      // Jet collection helpers
      VPseudoJetPtr Sorted_by_pt_min_pt(VPseudoJet const& jets, float jetPtMin);
      // Subjet QGTagging helpers
      void   FillSubjetQGTagging(fastjet::PseudoJet const& jet, float constitsPtMin,
                                 XlSubJet *pSubJet, XlFatJet const* pFatJet);
      // Color pull helpers
      TVector2 GetPull(fastjet::PseudoJet const&, float constitsPtMin);
      double GetPullAngle(VPseudoJetPtr const& fjSubJets, float constitsPtMin);
      double fMicrojetConeSize = -1.0;

      double GetQjetVolatility (VPseudoJet const& constits, int QJetsN = 25, int seed = 12345);
      VPseudoJet FilterJetsByPt(VPseudoJet const&, double ptMin);
      double FindRMS(std::vector<float>);
      double FindMean(std::vector<float>);


      Vect4M GetCorrectedMomentum(fastjet::PseudoJet const&, double thisJEC);


    private:
      Bool_t fIsData;                      //is this data or MC?
      Bool_t fQGTaggingActive;             //=true if QGTagging info is filled
      Bool_t fQGTaggerCHS;                 //=true if QGTagging weights are taken from CHS
      Bool_t fPublishOutput;               //=true if output collection are published

      TString fFatJetsName;                   //(i) name of input jets
      Bool_t fFatJetsFromBranch;              //are input jets from Branch?
      const JetCol *fFatJets;                 //input jets

      TString fPFCandidatesName;           //(i) name of PF candidates coll
      Bool_t fPFCandidatesFromBranch;
      const PFCandidateCol *fPFCandidates; //particle flow candidates coll handle

      TString fPileUpDenName;              //(i) name of PU energy density coll
      Bool_t fPileUpDenFromBranch;         //is PU energy density from Branch?
      const PileupEnergyDensityCol *fPileUpDen; //PU energy density coll handle

      TString fVertexesName;               //(i) name of vertex coll
      Bool_t fVertexesFromBranch;          //are vertexex from Branch?
      const VertexCol *fVertexes;          //vertex coll handle

      TString fXlFatJetsName;              //name of output fXlFatJets collection
      XlFatJetArr* fXlFatJets{0};             //array of fXlFatJets

      // used for nsubjettiness calculation
      fastjet::contrib::Njettiness* fNJettiness{0};

      fastjet::contrib::EnergyCorrelatorDoubleRatio* fECR2b0{0};
      fastjet::contrib::EnergyCorrelatorDoubleRatio* fECR2b0p2{0};
      fastjet::contrib::EnergyCorrelatorDoubleRatio* fECR2b0p5{0};
      fastjet::contrib::EnergyCorrelatorDoubleRatio* fECR2b1{0};
      fastjet::contrib::EnergyCorrelatorDoubleRatio* fECR2b2{0};

      // Objects from fastjet we want to use
      fastjet::Pruner* fPruner{0};
      fastjet::Filter* fFilterer{0};
      fastjet::Filter* fTrimmer{0};           //no this is not a typo trimmers belong to fastjet Filter class

      double fSoftDropZCut;                //soft-drop Z cut
      double fSoftDropR0;               //soft-drop angular distance normalisation
      fastjet::contrib::SoftDrop* fSoftDropSDb0{0};
      fastjet::contrib::SoftDrop* fSoftDropSDb1{0};
      fastjet::contrib::SoftDrop* fSoftDropSDb2{0};
      fastjet::contrib::SoftDrop* fSoftDropSDbm1{0};

      double fPruneZCut;                   //pruning Z cut
      double fPruneDistCut;                //pruning distance cut
      int fFilterN;                        //number of subjets after filtering
      double fFilterRad;                   //filtered subjets radius
      double fTrimRad;                     //trimmed subjet radius
      double fTrimPtFrac;                  //trimmed subjet pt fraction
      double fConeSize;                    //fastjet clustering radius
      fastjet::JetDefinition* fCAJetDef{0};   //fastjet clustering definition
      fastjet::GhostedAreaSpec* fActiveArea{0};
      fastjet::AreaDefinition* fAreaDefinition{0};
      fastjet::CMSTopTagger* fCMSTopTagger{0};

      unsigned short fSubJetFlags = 1;    // flags turning on subjet types

      TString fXlSubJetsName[XlSubJet::nSubJetTypes];              //name of output fXlSubJets collection
      XlSubJetArr* fXlSubJets[XlSubJet::nSubJetTypes];             //array of fXlSubJets

      Deconstruction::Deconstruct* fDeconstruct{0};
      AnalysisParameters* fParam{0};
      Deconstruction::TopGluonModel* fSignal{0};
      Deconstruction::BackgroundModel* fBackground{0};
      Deconstruction::ISRModel* fISR{0};
      TString fInputCard;

      UInt_t fProcessNJets;

      Bool_t fDoShowerDeconstruction;

      QGTagger* fQGTagger{0};                 //QGTagger calculator

      Bool_t fBeVerbose;
      Bool_t fDoECF;
      Bool_t fDoQjets;
      unsigned int fNMaxMicrojets;

      TStopwatch* fStopwatch{0};

      // Counters : used to initialize seed for QJets volatility
      Long64_t fCounter;

      int fDebugFlag = -1;

      ClassDef(FatJetExtenderMod, 0)         //XlJets, Fat and Sub, filler
  };
}
#endif
