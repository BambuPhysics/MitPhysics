#include <vector>
#include <TFile.h>
#include <TRandom3.h>
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "MitPhysics/Utils/interface/TauIsoMVA.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "CondFormats/EgammaObjects/interface/GBRForest.h"
// #include "Cintex/Cintex.h"
#include <utility>

ClassImp(mithep::TauIsoMVA)

using namespace mithep;

//--------------------------------------------------------------------------------------------------
TauIsoMVA::TauIsoMVA() {
  fReader = 0;
}
//--------------------------------------------------------------------------------------------------
TauIsoMVA::~TauIsoMVA() {
  delete fReader;
}
//--------------------------------------------------------------------------------------------------
void TauIsoMVA::Initialize( TString iWeightFile) { 
  fReader = new TMVA::Reader("!Color:!Silent");
  fReader->AddVariable("tau_nchiso",  (Float_t *)0);
  fReader->AddVariable("tau_ngiso",   (Float_t *)0);
  fReader->AddVariable("tau_nneuiso", (Float_t *)0);
  fReader->AddVariable("ring_ch_0*tau_pt", (Float_t *)0);
  fReader->AddVariable("ring_ch_1*tau_pt", (Float_t *)0);
  fReader->AddVariable("ring_ch_2*tau_pt", (Float_t *)0);
  fReader->AddVariable("ring_ch_3*tau_pt", (Float_t *)0);
  fReader->AddVariable("ring_ch_4*tau_pt", (Float_t *)0);
  fReader->AddVariable("ring_g_0*tau_pt", (Float_t *)0);
  fReader->AddVariable("ring_g_1*tau_pt", (Float_t *)0);
  fReader->AddVariable("ring_g_2*tau_pt", (Float_t *)0);
  fReader->AddVariable("ring_g_3*tau_pt", (Float_t *)0);
  fReader->AddVariable("ring_g_4*tau_pt", (Float_t *)0);
  fReader->AddVariable("ring_neu_0*tau_pt", (Float_t *)0);
  fReader->AddVariable("ring_neu_1*tau_pt", (Float_t *)0);
  fReader->AddVariable("ring_neu_2*tau_pt", (Float_t *)0);
  fReader->AddVariable("ring_neu_3*tau_pt", (Float_t *)0);
  fReader->AddVariable("ring_neu_4*tau_pt", (Float_t *)0);
  fReader->AddVariable("shape_ch_eta",    (Float_t *)0);
  fReader->AddVariable("shape_ch_phi",    (Float_t *)0);
  fReader->AddVariable("shape_ch_etaeta", (Float_t *)0);
  fReader->AddVariable("shape_ch_phiphi", (Float_t *)0);
  fReader->AddVariable("shape_ch_etaphi", (Float_t *)0);
  fReader->AddVariable("shape_g_eta",    (Float_t *)0);
  fReader->AddVariable("shape_g_phi",    (Float_t *)0);
  fReader->AddVariable("shape_g_etaeta", (Float_t *)0);
  fReader->AddVariable("shape_g_phiphi", (Float_t *)0);
  fReader->AddVariable("shape_g_etaphi", (Float_t *)0);
  fReader->AddVariable("shape_neu_eta",    (Float_t *)0);
  fReader->AddVariable("shape_neu_phi",    (Float_t *)0);
  fReader->AddVariable("shape_neu_etaeta", (Float_t *)0);
  fReader->AddVariable("shape_neu_phiphi", (Float_t *)0);
  fReader->AddVariable("shape_neu_etaphi", (Float_t *)0);
  fReader->AddVariable("rho", (Float_t *)0);
  fReader->AddSpectator("tau_pt",  (Float_t *)0);
  fReader->AddSpectator("tau_eta", (Float_t *)0);
  fReader->AddSpectator("tau_iso", (Float_t *)0);
  fReader->AddSpectator("gen_pt",  (Float_t *)0);
  fReader->AddSpectator("jet_pt",  (Float_t *)0);
  fReader->AddSpectator("pv",      (Float_t *)0);
  fReader->BookMVA     ("BDTG", iWeightFile);
  fGBR = false;
}
void TauIsoMVA::InitializeGBR( TString iWeightFile) {
  // Is this necessary? (Y.I. Apr 2 2015)
  // ROOT::Cintex::Cintex::Enable();   
  TFile *lForest = new TFile(iWeightFile,"READ");
  fGBRReader = (GBRForest*)lForest->Get("gbrfTauIso");
  lForest->Close();
  fGBR = true;
}
double   TauIsoMVA::MVAValue(const PFTau *iTau,double iRho) {
  IsoRings lRings        = computeIsoRings(iTau);
  std::vector<float> mvainput = lRings.getVector();
  mvainput.push_back( iRho);
  if(!fGBR) mvainput.insert(mvainput.end(), 6, 0); 
  if(!fGBR) return fReader->EvaluateMVA(mvainput, "BDTG");
  return fGBRReader->GetClassifier(&mvainput[0]);
}
TauIsoMVA::IsoRings TauIsoMVA::computeIsoRings(const PFTau *iTau) {
  std::vector<int>                 niso    (3);
  std::vector<std::vector<float> > rings   (3, std::vector<float>(5));
  std::vector<std::vector<float> > shapes  (3, std::vector<float>(5));
  std::vector<float>               isoptsum(3);

  for(UInt_t i0 = 0; i0 < iTau->NIsoPFCandS(); i0++) {
    const PFCandidate *pCand = iTau->IsoPFCand(i0);
    float deta = iTau->Eta() - pCand->Eta();
    float dphi =      MathUtils::DeltaPhi((double)iTau ->Phi(), (double)pCand->Phi());
    float dr   = fabs(MathUtils::DeltaR  ((double)iTau ->Phi(), (double)iTau ->Eta(),
					  (double)pCand->Phi(), (double)pCand->Eta()));
    int pftype = 0;
    
      
      if(pCand->Charge()       != 0)                   pftype = 0;
      else if(pCand->PFType()  == PFCandidate::eGamma) pftype = 1;
      else                                            pftype = 2;
      
      // Number of isolation candidates by type	
	
	niso[pftype]++;
	
	// Isolation Rings			
	  
	  if(dr < 0.1)      rings[pftype][0] += pCand->Pt();
	  else if(dr < 0.2) rings[pftype][1] += pCand->Pt();
	  else if(dr < 0.3) rings[pftype][2] += pCand->Pt();
	  else if(dr < 0.4) rings[pftype][3] += pCand->Pt();
	  else if(dr < 0.5) rings[pftype][4] += pCand->Pt();
	  
	  // Angle Shape Variables		
	    
	    shapes[pftype][0] += pCand->Pt() * deta;
	    shapes[pftype][1] += pCand->Pt() * dphi;
	    shapes[pftype][2] += pCand->Pt() * deta*deta;
	    shapes[pftype][3] += pCand->Pt() * dphi*dphi;
	    shapes[pftype][4] += pCand->Pt() * deta*dphi;
	    isoptsum[pftype]  += pCand->Pt();
  }
  
  // Mean and variance of angle variables are weighted by pT	
    
    for(uint i = 0; i < shapes.size(); i++)
      {
	for(uint j = 0; j < shapes[i].size(); j++)
	  {
	    shapes[i][j] = isoptsum[i] > 0 ? fabs(shapes[i][j]/isoptsum[i]) : 0;
	  }
      }
    
    // Fill IsoRing object			
      
      IsoRings isoRings;
      isoRings.niso = niso;
      isoRings.rings = rings;
      isoRings.shapes = shapes;
      
      return isoRings;
}
