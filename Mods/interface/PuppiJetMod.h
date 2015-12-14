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
#include "fastjet/ClusterSequence.hh"

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
      void SetR0(float d)                { fR0 = d; fR0Squared = d*d;  }
      void SetBeVerbose(Bool_t b)          { fBeVerbose = b;  }
      void SetDebugFlag(int i)   { fDebugFlag = i; }

      void SetDoMatching(Bool_t b)              { fDoMatching = b;  }
      void SetMatchingJetsName(const char *n)   { fMatchingJetsName = n;  }

      void SetJetAlgorithm(JetAlgorithm j)      { fJetAlgorithm = j; }
      void SetDoGhosts(Bool_t b)                { fDoGhosts = b; }
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

      Bool_t PassJet(fastjet::PseudoJet&);  // decide if we want to keep a jet before allocating space for it

    private:
      Bool_t fIsData;                      //is this data or MC?
      Bool_t fPublishOutput;               //=true if output collection are published

      TString fPFCandidatesName;           //(i) name of PF candidates coll
      const PFCandidateCol *fPFCandidates; //particle flow candidates coll handle

      TString fJetsName;              //name of output jets collection
      Array<JETTYPE> *fJets;             //array of jets

      // Objects from fastjet we want to use
      float fR0;                    //fastjet clustering radius
      float fR0Squared;             //fastjet clustering radius, squared 
      fastjet::JetDefinition *fJetDef;   //fastjet clustering definition
      fastjet::GhostedAreaSpec *fActiveArea;
      fastjet::AreaDefinition *fAreaDefinition;

      Bool_t fDoMatching;
      Bool_t fDoGhosts;
      TString fMatchingJetsName;
      const Collection<Jet> *fMatchingJets;

      UInt_t fProcessNJets;

      Bool_t fBeVerbose;

      JetAlgorithm fJetAlgorithm;
      
      TStopwatch *fStopwatch;

      int fDebugFlag = -1;
      int globalJetCounter=0, globalMatchedJetCounter=0;

      ClassDef(PuppiJetMod, 0)         //Puppi jets filler
  };

typedef PuppiJetMod<FatJet> PuppiFatJetMod;
typedef PuppiJetMod<PFJet>  PuppiPFJetMod;


//--------------------------------------------------------------------------------------------------
template<typename JETTYPE> PuppiJetMod<JETTYPE>::PuppiJetMod(const char *name, const char *title) :
  BaseMod (name,title),
  fIsData (kFALSE),
  fPublishOutput (kTRUE),
  fPFCandidatesName ("PuppiParticles"),
  fPFCandidates (0),
  fJetsName ("PuppiJets"),
  fR0(0.4),
  fR0Squared(0.16),
  fDoMatching(kTRUE),
  fDoGhosts(kFALSE),
  fMatchingJetsName("AKt4PFJetsCHS"),
  fMatchingJets(0),
  fProcessNJets (4),
  fJetAlgorithm(kAntiKT)
{
  // Constructor.

}

template<typename JETTYPE> PuppiJetMod<JETTYPE>::~PuppiJetMod()
{
}

template<> Bool_t PuppiJetMod<FatJet>::PassJet(fastjet::PseudoJet&);
template<> Bool_t PuppiJetMod<PFJet>::PassJet(fastjet::PseudoJet&);

//--------------------------------------------------------------------------------------------------
template<typename JETTYPE> void PuppiJetMod<JETTYPE>::Process()
{
  if (fBeVerbose) fStopwatch->Start(kTRUE);

  // make sure the out collections are empty before starting
  fJets->Delete();
  
  fPFCandidates = GetObject<PFCandidateCol>(fPFCandidatesName);
  if (fDoMatching)
    fMatchingJets = GetObject<Collection<Jet>>(fMatchingJetsName);
  
  int nPFCandidates = fPFCandidates->GetEntries();
  std::vector<fastjet::PseudoJet> inPseudoJets;

  // Loop over PFCandidates and unmark them : necessary for skimming
  // Convert PFCandidates into fastjet::PseudoJets
  std::vector<double> pts;
  for (int iC=0; iC!=nPFCandidates; ++iC) {
    const PFCandidate *pfCand = fPFCandidates->At((unsigned int)iC);
    pfCand->UnmarkMe();
    fastjet::PseudoJet newPseudoJet(pfCand->Px(),pfCand->Py(),pfCand->Pz(),pfCand->E());
    newPseudoJet.set_user_index(iC);
    pts.push_back(newPseudoJet.perp());
    inPseudoJets.push_back(newPseudoJet);
  }
  std::sort(pts.begin(),pts.end());
  fprintf(stderr,"PTs:\n");
  for (auto pt : pts)
    fprintf(stderr,"%.3f \n",pt);
  // cluster puppi jets

  fastjet::ClusterSequence *fjClusteringPtr;
  if (fDoGhosts) 
    fjClusteringPtr = new fastjet::ClusterSequenceArea(inPseudoJets,*fJetDef,*fAreaDefinition);
  else
    fjClusteringPtr = new fastjet::ClusterSequence(inPseudoJets,*fJetDef);

  fastjet::ClusterSequence fjClustering = *fjClusteringPtr;
  std::vector<fastjet::PseudoJet>baseJets = sorted_by_pt(fjClustering.inclusive_jets(10.)); 

  unsigned int nClusteredJets = baseJets.size();
  for (unsigned int iJ=0; iJ!=nClusteredJets; ++iJ) {
    // convert from pseudojets to Bambu jets
    if (iJ>fProcessNJets)
      break;
    fastjet::PseudoJet baseJet = baseJets[iJ];
    if (!PassJet(baseJet))
      continue;
    ++globalJetCounter;
    JETTYPE *newJet = fJets->AddNew();
    newJet->SetRawPtEtaPhiM(baseJet.pt(),baseJet.eta(),baseJet.phi(),baseJet.m());
    std::vector<fastjet::PseudoJet> baseConstituents = baseJet.constituents();
    // add PFCandidates to the jet and keep track of energy fractions
    const PFCandidate *pfCand = 0;
    float chargedEMEnergy=0;
    float chargedHadEnergy=0;
    float neutralEMEnergy=0;
    float neutralHadEnergy=0;
    float muonEnergy=0;
    int chargedMult=0;
    int neutralMult=0;
    int muonMult=0;
    for (fastjet::PseudoJet &constituent : baseConstituents) {
      if (constituent.user_index()<0) 
        continue; // ghost particle
      pfCand = fPFCandidates->At((unsigned int)constituent.user_index());
      newJet->AddPFCand(pfCand);
      switch (pfCand->PFType()) {
        case PFCandidate::eHadron:
          chargedHadEnergy+=pfCand->E();
          ++chargedMult;
          break;
        case PFCandidate::eElectron:
          chargedEMEnergy+=pfCand->E();
          ++chargedMult;
          break;
        case PFCandidate::eMuon:
          muonEnergy+=pfCand->E();
          ++muonMult;
          break;
        case PFCandidate::eGamma:
        case PFCandidate::eEGammaHF:
          neutralEMEnergy+=pfCand->E();
          ++neutralMult;
          break;
        case PFCandidate::eNeutralHadron:
        case PFCandidate::eHadronHF:
          neutralHadEnergy+=pfCand->E();
          ++neutralMult;
          break;
        default:
          break;
      }
    }
    newJet->SetChargedEmEnergy(chargedEMEnergy);
    newJet->SetChargedHadronEnergy(chargedHadEnergy);
    newJet->SetNeutralEmEnergy(neutralEMEnergy);
    newJet->SetNeutralHadronEnergy(neutralHadEnergy);
    newJet->SetChargedMuEnergy(muonEnergy);
    newJet->SetChargedMultiplicity(chargedMult);
    newJet->SetNeutralMultiplicity(neutralMult);
    newJet->SetMuonMultiplicity(muonMult);

    if (fDoMatching) {
      assert(fMatchingJets!=NULL);
      RunMatching(newJet);
    }
  }

  fJets->Trim();

  if (fjClusteringPtr)
    delete fjClusteringPtr;

  if (fBeVerbose) fStopwatch->Stop();
  return;
}


template <typename JETTYPE> void PuppiJetMod<JETTYPE>::RunMatching(PFJet *newJet) {
  // dR matches a new jet with an existing collection
  const PFJet *matchedJet = 0;
  float bestDR2 = 999;
  for (unsigned int iJ=0; iJ!=fMatchingJets->GetEntries(); ++iJ) {
    // do dR Matching
    const PFJet *tmpJet = dynamic_cast<const PFJet*>(fMatchingJets->At(iJ));
    float dR2 = MathUtils::DeltaR2<const PFJet, PFJet>(tmpJet,newJet);
    if (dR2<bestDR2 && dR2<fR0Squared) {
      matchedJet = tmpJet;
      bestDR2 = dR2;
    }
  }
  if (matchedJet==0)
    return;
  ++globalMatchedJetCounter;
  // copy btags
  for (unsigned int iBtag=0; iBtag!=(unsigned int)Jet::nBTagAlgos; ++iBtag) {
    newJet->SetBJetTagsDisc(matchedJet->BJetTagsDisc(iBtag),iBtag);
  }
  return;
}


template<typename JETTYPE> void PuppiJetMod<JETTYPE>::RunMatching(FatJet *newJet) {
  // dR matches a new jet with an existing collection
  const FatJet *matchedJet = 0;
  float bestDR2 = 999;
  for (unsigned int iJ=0; iJ!=fMatchingJets->GetEntries(); ++iJ) {
    // do dR Matching
    const FatJet *tmpJet = dynamic_cast<const FatJet*>(fMatchingJets->At(iJ));
    float dR2 = MathUtils::DeltaR2<const FatJet, FatJet>(tmpJet,newJet);
    if (dR2<bestDR2 && dR2<fR0Squared) {
      matchedJet = tmpJet;
      bestDR2 = dR2;
    }
  }
  if (matchedJet==0)
    return;
  ++globalMatchedJetCounter;
  // copy btags
  for (unsigned int iBtag=0; iBtag!=(unsigned int)Jet::nBTagAlgos; ++iBtag) {
    newJet->SetBJetTagsDisc(matchedJet->BJetTagsDisc(iBtag),iBtag);
  }
  // copy fatjet b-tagging infos
  newJet->SetTau1IVF(matchedJet->Tau1IVF());
  newJet->SetTau2IVF(matchedJet->Tau2IVF());
  newJet->SetTau3IVF(matchedJet->Tau3IVF());
  newJet->SetTauDot(matchedJet->tauDot());
  newJet->SetZRatio(matchedJet->zRatio());

  for (int i=0; i!=3; ++i)  
    newJet->AddTauIVFAxis(matchedJet->GetTauIVFAxis(i));
  
  for (FatJet::TrackData track : matchedJet->GetTrackData()) 
    newJet->AddTrackData(&track);
  for (FatJet::SVData sv : matchedJet->GetSVData()) 
    newJet->AddSVData(&sv);
  for (FatJet::LeptonData lep : matchedJet->GetElectronData())
    newJet->AddElectronData(&lep);
  for (FatJet::LeptonData lep : matchedJet->GetMuonData())
    newJet->AddMuonData(&lep);
  
  for (float sjBtag : matchedJet->GetSubJetBtags())
    newJet->AddSubJetBtag(sjBtag);

  return;
}

//--------------------------------------------------------------------------------------------------
template<typename JETTYPE> void PuppiJetMod<JETTYPE>::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis.
    
  // Create the new output collection
  fJets = new Array<JETTYPE>(16,fJetsName);
  
  // Publish collection for further usage in the analysis
  if (fPublishOutput)
    PublishObj(fJets);
  
  if (fJetAlgorithm==kCambridgeAachen)
    fJetDef = new fastjet::JetDefinition(fastjet::cambridge_algorithm, fR0);
  else if (fJetAlgorithm==kKT)
    fJetDef = new fastjet::JetDefinition(fastjet::kt_algorithm, fR0);
  else
    fJetDef = new fastjet::JetDefinition(fastjet::antikt_algorithm,fR0);

  // Initialize area caculation (done with ghost particles)
  int activeAreaRepeats = 1;
  double ghostArea = 0.01;
  double ghostEtaMax = 5.0;
  if (fDoGhosts) {
    fActiveArea = new fastjet::GhostedAreaSpec(ghostEtaMax,activeAreaRepeats,ghostArea);
    fAreaDefinition = new fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts,*fActiveArea);
  }


  fStopwatch = new TStopwatch();

  return;
}

//--------------------------------------------------------------------------------------------------
template<typename JETTYPE> void PuppiJetMod<JETTYPE>::SlaveTerminate()
{
  Info("SlaveTerminate","Matched %i/%i jets",globalMatchedJetCounter,globalJetCounter);
  RetractObj(fJets->GetName());

  if (fJets)
    delete fJets;

  if (fJetDef)
    delete fJetDef;

  if (fDoGhosts) {
    if (fActiveArea)
      delete fActiveArea;
    if (fAreaDefinition)
      delete fAreaDefinition;
  }

  if (fStopwatch)
    delete fStopwatch;
}



}


#endif
