#include "MitPhysics/Mods/interface/MVAMetMod.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataTree/interface/MetCol.h"
#include "MitAna/DataTree/interface/ParticleCol.h"
#include "MitAna/DataTree/interface/Names.h"
#include "MitPhysics/Init/interface/ModNames.h"

#include "TSystem.h"

using namespace mithep;

ClassImp(mithep::MVAMetMod)

//--------------------------------------------------------------------------------------------------
MVAMetMod::MVAMetMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fMVAMetName("MetMVA"),  
  fJetsName  ("correctPFJets"),
  fPFCandName(Names::gkPFCandidatesBrn),
  fVertexName(ModNames::gkGoodVertexesName),
  fPFMetName ("PFMet"),
  fRhoName   ("Rho"),
  fMVAMet(0)
{
  throw std::runtime_error("Broken Mod: MVAMetMod and Utils/MVAMet needs a major rework.");
}

//-------------------------------------------------------------------------------------------------
void MVAMetMod::Process()
{
  // Process entries of the tree. 

  auto* jets = GetObject<JetCol>(fJetsName); //corrected Jets
  //LoadBranch(fJetsName);
  auto* cands = GetObject<PFCandidateCol>(fPFCandName);
  auto* vertices = GetObject<VertexOArr>(fVertexName);
  auto* met = GetObject<PFMetCol>(fPFMetName);
  auto* rhoCol = GetObject<PileupEnergyDensityCol>(fRhoName);

  if (!jets || !cands || !vertices || !met) {
    SendError(kAbortModule, "Process", 
              "Pointer to input jet collection %s is null.",
              fJetsName.Data());
    return;
  }
  // lepton selection
  Float_t lPt0 = 0; Float_t lEta0 = 0; Float_t lPhi0 = 0;
  Float_t lPt1 = 0; Float_t lEta1 = 0; Float_t lPhi1 = 0;

  ParticleOArr *leptons = GetObjThisEvt<ParticleOArr>(ModNames::gkMergedLeptonsName);
  if(leptons->GetEntries() >= 1){
    lPt0 = leptons->At(0)->Pt();  lEta0 = leptons->At(0)->Eta(); lPhi0 = leptons->At(0)->Phi();
  }
  if(leptons->GetEntries() >= 2){
    lPt1 = leptons->At(1)->Pt();  lEta1 = leptons->At(1)->Eta(); lPhi1 = leptons->At(1)->Phi();
  }

  MetOArr *MVAMet = new MetOArr; 
  MVAMet->SetName(fMVAMetName);
  
  PFJetOArr *thePFJets = new PFJetOArr; 
  for(UInt_t i=0; i<jets->GetEntries(); i++) {
    const PFJet *ptJet = dynamic_cast<const PFJet*>(jets->At(i));
    thePFJets->Add(ptJet);
  }

  Met lMVAMet = fMVAMet->GetMet(  false,
                                  lPt0,lPhi0,lEta0,
                                  lPt1,lPhi1,lEta1,
                                  met->At(0),
                                  cands,vertices->At(0),vertices,rhoCol->At(0)->Rho(),
                                  thePFJets,
                                  int(vertices->GetEntries()));

  MVAMet->Add(&lMVAMet);

  // sort according to pt
  MVAMet->Sort();

  // add to event for other modules to use
  AddObjThisEvt(MVAMet);

  // delete temporal PFJet
  delete thePFJets;

}

//--------------------------------------------------------------------------------------------------
void MVAMetMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we typically initialize histograms and other analysis objects and request
  // branches. For this module, we request a branch of the MitTree.

  using std::string;

  TString dataDir(gSystem->Getenv("MIT_DATA"));
  if (dataDir.Length() == 0) {
    SendError(kAbortModule, "SlaveBegin", "MIT_DATA environment is not set.");
    return;
  }

  fMVAMet    = new MVAMet();
  fMVAMet->Initialize(dataDir + "/mva_JetID_lowpt.weights.xml",
                      dataDir + "/mva_JetID_highpt.weights.xml",
                      dataDir + "/JetIDMVA_JetIdParams.py",
                      dataDir + "/gbrmet_52.root",
                      dataDir + "/gbrmetphi_52.root",
                      dataDir + "/gbrmetu1cov_52.root",
                      dataDir + "/gbrmetu2cov_52.root"
                      );
}

//--------------------------------------------------------------------------------------------------
void MVAMetMod::SlaveTerminate()
{
  delete fMVAMet;
}
