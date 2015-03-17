#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Mods/interface/VTagMod.h"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/contrib/Njettiness.hh"
#include "fastjet/contrib/Nsubjettiness.hh"

using namespace mithep;

ClassImp(mithep::VTagMod)

//--------------------------------------------------------------------------------------------------
VTagMod::VTagMod(const char *name, const char *title) : 
  BaseMod          (name,title),
  fPFCandidatesName(Names::gkPFCandidatesBrn),
  fGoodVTagsName   (ModNames::gkGoodVTagsName),
  fConeSize        (0.8)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void VTagMod::Process()
{
  // Load the branches we want to work with
  LoadBranch(fPFCandidatesName);

  // Do not even continue if particle flow candidates are not there
  assert(fPFCandidates);

  // Push all particle flow candidates into fastjet particle collection
  std::vector<fastjet::PseudoJet> lFJParts;
  for(UInt_t i=0; i<fPFCandidates->GetEntries(); ++i) {      
    const PFCandidate *pfCand = fPFCandidates->At(i);
    lFJParts.push_back(fastjet::PseudoJet(pfCand->Px(),pfCand->Py(),pfCand->Pz(),pfCand->E()));
    lFJParts.back().set_user_index(i);
  }

  // Setup the cluster for fastjet
  fastjet::ClusterSequenceArea *lClustering =
    new fastjet::ClusterSequenceArea(lFJParts,*fCAJetDef,*fAreaDefinition);

  // Produce a new set of jets based on the fastjet particle collection and the defined clustering
  std::vector<fastjet::PseudoJet> lOutJets = sorted_by_pt(lClustering->inclusive_jets(0.0));

  // You get a 4-vector (PseudoJet)
  fastjet::PseudoJet lJet;
  if (lOutJets.size() > 0) {
    lJet = (*fPruner)(lOutJets[0]);
  }

  delete lClustering;
}

//--------------------------------------------------------------------------------------------------
void VTagMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here, we just request the
  // particle flow collection branch.

  ReqEventObject(fPFCandidatesName, fPFCandidates, kTRUE);

  //Default pruning parameters
  fPruner          = new fastjet::Pruner(fastjet::cambridge_algorithm,0.1,0.5);       // CMS Default
  
  //CA constructor (fConeSize = 0.8 for CA8)
  fCAJetDef       = new fastjet::JetDefinition(fastjet::cambridge_algorithm, fConeSize);
  
  //Area caculation (done with ghost particles)
  int    activeAreaRepeats = 1;
  double ghostArea         = 0.01;
  double ghostEtaMax       = 7.0;
  
  fActiveArea     = new fastjet::ActiveAreaSpec(ghostEtaMax,activeAreaRepeats,         ghostArea);
  fAreaDefinition = new fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts, *fActiveArea);
}

//--------------------------------------------------------------------------------------------------
void VTagMod::Terminate()
{
  // Run finishing code on the computer (slave) that did the analysis
}

//--------------------------------------------------------------------------------------------------
float VTagMod::GetTau(fastjet::PseudoJet &iJet,int iN, float iKappa)
{
  // Calculate the tau variable for the given pseudojet
  fastjet::contrib::Nsubjettiness 
    nSubNKT(iN,fastjet::contrib::Njettiness::onepass_kt_axes,iKappa,fConeSize,fConeSize);
  
  return nSubNKT(iJet);
}
