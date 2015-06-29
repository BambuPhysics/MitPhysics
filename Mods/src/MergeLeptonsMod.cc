#include <TH1D.h>
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/ParticleCol.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Mods/interface/MergeLeptonsMod.h"

using namespace mithep;

ClassImp(mithep::MergeLeptonsMod)

//--------------------------------------------------------------------------------------------------
mithep::MergeLeptonsMod::MergeLeptonsMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fElName(ModNames::gkCleanElectronsName),
  fMuName(ModNames::gkCleanMuonsName),
  fMergedName(ModNames::gkMergedLeptonsName),
  fColOut(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void mithep::MergeLeptonsMod::BeginRun()
{
}

//--------------------------------------------------------------------------------------------------
void mithep::MergeLeptonsMod::Process()
{
  // Merge the two input collections and publish merged collection. 

  auto* elIn = GetObject<ElectronCol>(fElName);
  auto* muIn = GetObject<MuonCol>(fMuName);

  // determine how many there are in total 
  UInt_t nents = 0;
  if (elIn) 
    nents += elIn->GetEntries();
  if (muIn) 
    nents += muIn->GetEntries();

  // book collection with right length
  fColOut = new mithep::ParticleOArr(nents, GetMergedName());

  if (elIn)
    fColOut->Add(elIn);
  if (muIn)
    fColOut->Add(muIn);

  // sort according to pt
  fColOut->Sort();

  // add to event for other modules to use
  AddObjThisEvt(fColOut);

  // fill histograms
  if (elIn)
    fRecoWElectrons->Fill(elIn->GetEntries());

  if (muIn)
    fRecoWMuons->Fill(muIn->GetEntries());
}

//--------------------------------------------------------------------------------------------------
void mithep::MergeLeptonsMod::SlaveBegin()
{
  AddTH1(fRecoWMuons,    "fRecoWMuons",    "Number of Reconstructed Muons;N_{muons};#",
	 10,-0.5,9.5);
  AddTH1(fRecoWElectrons,"fRecoWElectrons","Number of Reconstructed Electrons;N_{electrons};#",
	 10,-0.5,9.5);
}
//--------------------------------------------------------------------------------------------------
void mithep::MergeLeptonsMod::SlaveTerminate()
{
}
