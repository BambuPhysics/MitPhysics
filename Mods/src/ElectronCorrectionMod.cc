#include <TH1D.h>
#include <TRandom3.h>
#include <TNtuple.h>
#include "MitAna/DataTree/interface/Names.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/EventHeader.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Mods/interface/ElectronCorrectionMod.h"

using namespace mithep;

ClassImp(mithep::ElectronCorrectionMod)

//--------------------------------------------------------------------------------------------------
ElectronCorrectionMod::ElectronCorrectionMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fElectronsFromBranch(kTRUE),
  fDoNotCorrect(kFALSE),
  fIsMC(kFALSE),
  fCorrectData(kTRUE),
  fSmearMC(kTRUE),
  fInElectronName(Names::gkElectronBrn),
  fCorrectionType("FromJosh"),
  fCorrectionFileName("scales.txt"),
  fCorrectedElectronsName(ModNames::gkCorrectedElectronsName),
  fElectrons(0),
  fScale(1.0),
  fResolution(0.0),
  fRand(0)
{
  // Constructor.
  fCorrections = new TNtuple("fCorrections","Scale factors and resolutions for MC -> Data",
			     "etaMin:etaMax:scale:resolution");

}
//--------------------------------------------------------------------------------------------------
ElectronCorrectionMod::~ElectronCorrectionMod()
{
  // Destructor
  delete fCorrections;
}
//--------------------------------------------------------------------------------------------------
void ElectronCorrectionMod::Begin()
{
  // Run startup code on the client machine. For this module, we dont do
  // anything here.
}

//--------------------------------------------------------------------------------------------------
void ElectronCorrectionMod::Process()
{
  LoadEventObject(fInElectronName,fElectrons);
  fIsMC = GetEventHeader()->IsMC();
  ElectronOArr *CorrectedElectrons = new ElectronOArr;
  CorrectedElectrons->SetOwner(kTRUE);
  CorrectedElectrons->SetName(fCorrectedElectronsName);

  Double_t rescaledSmearedPt; //Temporary variable, just for aesthetics

  // It is dumb to do this for every event, but I can't get it to work if I put it in slaveterminate
  // or in the constructor. huh?!
  if(fSmearMC && fIsMC) {
    fRand = new TRandom3();
    fRand->SetSeed();
    rescaledSmearedPt=0.0;
  }

  // Copy to new array, correct the Pt
  for(UInt_t i=0; i<fElectrons->GetEntries(); ++i) {
    Electron *el = new Electron(*(fElectrons->At(i)));
    CorrectedElectrons->AddOwned(el);

    //Set fScale and fResolution
    SetCorrectionValues(el);
    
    if(fDoNotCorrect) return;
    //Change the Pt
    if(fCorrectData && !fIsMC) {//If this is data *and* if we're supposed to correct data
      el->SetPtEtaPhi(fScale*el->Pt(),el->Eta(),el->Phi());
    }
    else {
      if(fSmearMC) {
	//rescale and smear the MC Pt
	rescaledSmearedPt  = fRand->Gaus((1/fScale)*el->Pt(), fResolution);
	el->SetPtEtaPhi(rescaledSmearedPt,el->Eta(),el->Phi());
      }
      else {
	//Just rescale the MC Pt
	el->SetPtEtaPhi((1/fScale)*el->Pt(),el->Eta(),el->Phi());
      }
    }
  }

  AddObjThisEvt(CorrectedElectrons);

  delete fRand;
}

//--------------------------------------------------------------------------------------------------
void ElectronCorrectionMod::SlaveBegin()
{
  ReqEventObject(fInElectronName,fElectrons,fElectronsFromBranch);
  fCorrections->ReadFile(fCorrectionFileName);
}

//--------------------------------------------------------------------------------------------------
void ElectronCorrectionMod::SlaveTerminate()
{
  // Run finishing code on the computer (slave) that did the analysis. 
}

//--------------------------------------------------------------------------------------------------
void ElectronCorrectionMod::Terminate()
{
  // Run finishing code on the client computer. For this module, we dont do
  // anything here.
}

//--------------------------------------------------------------------------------------------------
void ElectronCorrectionMod::SetCorrectionValues(Electron *el)
{
  Float_t etaMin=0.0, etaMax=0.0;

  fCorrections->SetBranchAddress("etaMin",&etaMin);
  fCorrections->SetBranchAddress("etaMax",&etaMax);
  fCorrections->SetBranchAddress("scale",&fScale);
  fCorrections->SetBranchAddress("resolution",&fResolution);
  Float_t etaAbs = fabs(el->Eta());
  UInt_t NBinsEta = fCorrections->GetEntries();
  for (UInt_t etaBin=0;etaBin<NBinsEta;etaBin++) {
    fCorrections->GetEntry(etaBin);
    if (etaAbs>=etaMin && etaAbs<etaMax)
      break;
    if (etaBin==NBinsEta-1) {
      fScale = 1.0;
      fResolution = 2.5; // Just making this number up . . .
      SendError(kWarning,"Process",
		Form("Electron not in any eta bin! Setting scale to %5f and resolution to %5f",
		     fScale,fResolution));
    }
  }
}
