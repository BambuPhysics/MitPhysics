#include "MitPhysics/Utils/interface/PhotonTools.h"
#include "MitPhysics/Utils/interface/MVAVBF.h"
#include "MitPhysics/Utils/interface/ElectronTools.h"
#include "MitPhysics/Utils/interface/IsolationTools.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/StableData.h"
#include <TMath.h> 
#include <TFile.h>
#include <TRandom3.h>
#include <TSystem.h>
#include "TMVA/Tools.h"//MVA
#include "TMVA/Reader.h"//MVA

ClassImp(mithep::MVAVBF)
  
  using namespace mithep;

//--------------------------------------------------------------------------------------------------
MVAVBF::MVAVBF():
  fReader(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void MVAVBF::InitializeMVA() {
  
  if (fReader) delete fReader;  
  
  fReader = new TMVA::Reader( "!Color:!Silent:Error" );       

  TString dataDir(gSystem->Getenv("MIT_DATA"));
  if (dataDir.Length() == 0)
    throw std::runtime_error("MIT_DATA environment is not set.");
  
  TString Weights;
  
  Weights = dataDir + "/TMVA_vbf_6var_mjj100_diphopt_phopt_BDTG.weights.xml";
  
  fReader->AddVariable("jet1pt",&_jet1pt);
  fReader->AddVariable("jet2pt",&_jet2pt);
  fReader->AddVariable("abs(jet1eta-jet2eta)",&_deltajeteta);
  fReader->AddVariable("mj1j2",&_dijetmass);
  fReader->AddVariable("zepp",&_zeppenfeld);
  fReader->AddVariable("dphi",&_dphidijetgg);
  fReader->AddVariable("diphopt/diphoM",&_diphoptOverdiphomass);
  fReader->AddVariable("pho1pt/diphoM",&_pho1ptOverdiphomass);
  fReader->AddVariable("pho2pt/diphoM",&_pho2ptOverdiphomass);

  fReader->BookMVA("BDT method",Weights);

  assert(fReader);
}

float MVAVBF::GetMVAbdtValue(float jet1pt, float jet2pt, float deltajeteta, float dijetmass, float zeppenfeld, float dphidijetgg, float diphoptOverdiphomass, float pho1ptOverdiphomass, float pho2ptOverdiphomass) {

  _jet1pt= jet1pt;
  _jet2pt= jet2pt;
  _deltajeteta= deltajeteta;
  _dijetmass= dijetmass;
  _zeppenfeld= zeppenfeld;
  _dphidijetgg= dphidijetgg;
  _diphoptOverdiphomass= diphoptOverdiphomass;
  _pho1ptOverdiphomass= pho1ptOverdiphomass;
  _pho2ptOverdiphomass= pho2ptOverdiphomass;

  TMVA::Reader* reader = NULL;
  reader = fReader;
  assert(reader); 
  
  return (reader->EvaluateMVA("BDT method"));
}

