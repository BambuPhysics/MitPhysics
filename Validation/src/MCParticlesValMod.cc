// $Id: MCParticlesValMod.cc,v 1.2 2009/06/15 15:00:23 loizides Exp $

#include "MitPhysics/Validation/interface/MCParticlesValMod.h"
#include "MitAna/DataTree/interface/MCParticleCol.h"
#include "MitAna/DataTree/interface/Names.h"
#include <TH1D.h>

using namespace mithep;

ClassImp(mithep::MCParticlesValMod)

//--------------------------------------------------------------------------------------------------
MCParticlesValMod::MCParticlesValMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fPartName(Names::gkMCPartBrn),
  fParticles(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void MCParticlesValMod::Process()
{
  // Process entries of the tree. For this module, we just load
  // the branch and fill the histograms.

  LoadBranch(GetPartName());

  Int_t ents=fParticles->GetEntries();
  for(Int_t i=0;i<ents;++i) {
     const MCParticle *p = fParticles->At(i);
     fHs[0]->Fill(p->Pt());
     fHs[1]->Fill(p->Eta());
     if (p->IsGenerated()) {
       fHs[2]->Fill(p->Pt());
       fHs[3]->Fill(p->Eta());
     }
     if (p->IsSimulated()) {
       fHs[4]->Fill(p->Pt());
       fHs[5]->Fill(p->Eta());
     }

     Double_t mdiff = (p->Mass() - p->PdgMass()) / (p->PdgMass()+0.01);
     fHs[6]->Fill(mdiff);
     if (p->IsGenerated())
       fHs[7]->Fill(mdiff);
     if (p->IsSimulated())
       fHs[8]->Fill(mdiff);
     if (mdiff>0.01) {
       std::cout << "----" << mdiff << "-> " << p->Mass() << " " << p->PdgMass() << ": " << std::endl;
       p->Print();
       fHs[9]->Fill(p->AbsPdgId());
     }
  }
}

//--------------------------------------------------------------------------------------------------
void MCParticlesValMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual
  // analysis. Here, we typically initialize histograms and 
  // other analysis objects and request branches. For this module,
  // we request a branch of the MitTree.

  ReqBranch(GetPartName(), fParticles);

  AddTH1(fHs[0],"hPtHist",";p_{t} [GeV]; #",100,0.,250.);
  AddTH1(fHs[1],"hEtaHist",";#eta;#",160,-8.,8.);
  AddTH1(fHs[2],"hGenPtHist",";p_{t} [GeV];#",100,0.,250.);
  AddTH1(fHs[3],"hGenEtaHist",";#eta;#",160,-8.,8.);
  AddTH1(fHs[4],"hSimPtHist",";p_{t} [GeV];#",100,0.,250.);
  AddTH1(fHs[5],"hSimEtaHist",";#eta;#",160,-8.,8.);
  AddTH1(fHs[6],"hMassDiff",";#eta;#",1000,-5.,5.);
  AddTH1(fHs[7],"hGenMassDiff",";#eta;#",1000,-5.,5.);
  AddTH1(fHs[8],"hSimMassDiff",";#eta;#",1000,-5.,5.);
  AddTH1(fHs[9],"hDiffPdg",";pdg;#",500,0.5,500.5);
}

//--------------------------------------------------------------------------------------------------
void MCParticlesValMod::SlaveTerminate()
{
  // Run finishing code on the computer (slave) that did the
  // analysis. For this module, we dont do anything here.
}
