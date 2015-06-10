#include "MitPhysics/Utils/interface/PhotonTools.h"
#include "MitPhysics/Utils/interface/MVATools.h"
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

using namespace mithep;

//--------------------------------------------------------------------------------------------------
MVATools::MVATools():
  fReaderEndcap(0),
  fReaderBarrel(0),
  
  fMVAType (MVATools::kNone)

  // ------------------------------------------------------------------------------
  // fab: These guys should all go away....
  //MVA Variables v4
//   HoE(0),
//   covIEtaIEta(0),
//   tIso1abs(0),
//   tIso3abs(0),
//   tIso2abs(0),
//   R9(0),
  
//   absIsoEcal(0),
//   absIsoHcal(0),
//   RelEMax(0),
//   RelETop(0),
//   RelEBottom(0),
//   RelELeft(0),
//   RelERight(0),
//   RelE2x5Max(0),
//   RelE2x5Top(0),
//   RelE2x5Bottom(0),
//   RelE2x5Left(0),
//   RelE2x5Right(0),
//   RelE5x5(0),
  
//   EtaWidth(0),
//   PhiWidth(0),
//   CoviEtaiPhi(0),
//   CoviPhiiPhi(0),
  
//   NVertexes(0),
//   RelPreshowerEnergy(0),

//   //MVA v2 and v1
//   RelIsoEcal(0),
//   RelIsoHcal(0),
//   tIso1(0),
//   tIso3(0),
//   tIso2(0),
//   ScEta(0.)
  
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void MVATools::InitializeMVA(MVATools::IdMVAType type) {
  
  fMVAType = type;

  if (fReaderEndcap) delete fReaderEndcap;  
  if (fReaderBarrel) delete fReaderBarrel; 

  // no MVA needed if none requested  
  if( type == kNone ) {    
    return;
  }
  
  fReaderEndcap = new TMVA::Reader( "!Color:!Silent:Error" );       
  fReaderBarrel = new TMVA::Reader( "!Color:!Silent:Error" );   

  TString BarrelWeights;
  TString EndcapWeights;

  varNames.resize(0);
  
  switch (type) {
    
  case k2011IdMVA_HZg:

    EndcapWeights =  (gSystem->Getenv("CMSSW_BASE")+
		      TString("/src/MitPhysics/data/")+
		      TString("PhotonId_lowPt_EE_BDTG.")+
		      TString("weights.xml"));
    BarrelWeights =  (gSystem->Getenv("CMSSW_BASE")+
		      TString("/src/MitPhysics/data/")+
		      TString("PhotonId_lowPt_EB_BDTG.")+
		      TString("weights.xml"));


    varNames.push_back("sigieie");
    varNames.push_back("covieip");
    varNames.push_back("s4r"    );
    varNames.push_back("r9"     );
    varNames.push_back("sigeta" );
    varNames.push_back("sigphi" );
    varNames.push_back("pfgios" );
    varNames.push_back("pfciso" );
    varNames.push_back("rho"    );
    varNames.push_back("sceta"  );
    varNames.push_back("rawe"   );

    mvaVars.resize(varNames.size());

    std::cout<<"  Adding stuff here.... "<<std::endl;

    for( unsigned int iV = 0; iV < mvaVars.size(); ++iV) {      
      mvaVarMapEB.insert(  std::pair<std::string,unsigned int>(varNames[iV], iV) );
      mvaVarMapEE.insert(  std::pair<std::string,unsigned int>(varNames[iV], iV) );
    }
    
    std::cout<<"  ... done "<<std::endl;

    break;

  case k2011IdMVA:

    EndcapWeights =  (gSystem->Getenv("CMSSW_BASE")+
		      TString("/src/MitPhysics/data/TMVAClassificationPhotonID_")+
		      TString("Endcap_PassPreSel_Variable_10_BDTnCuts2000_BDT.")+
		      TString("weights.xml"));
    BarrelWeights =  (gSystem->Getenv("CMSSW_BASE")+
		      TString("/src/MitPhysics/data/TMVAClassificationPhotonID_")+
		      TString("Barrel_PassPreSel_Variable_10_BDTnCuts2000_BDT.")+
		      TString("weights.xml"));

    // set up the variable names
    mvaVars.resize(12);
    varNames.push_back("HoE"        );
    varNames.push_back("covIEtaIEta");
    varNames.push_back("tIso1abs"   );
    varNames.push_back("tIso3abs"   );
    varNames.push_back("tIso2abs"   );
    varNames.push_back("R9"         );
    varNames.push_back("absIsoEcal" );
    varNames.push_back("absIsoHcal" );
    varNames.push_back("NVertexes"  );
    varNames.push_back("ScEta"      );
    varNames.push_back("EtaWidth"   );
    varNames.push_back("PhiWidth"   );

    for( unsigned int iV = 0; iV < mvaVars.size(); ++iV) {      
      mvaVarMapEB.insert(  std::pair<std::string,unsigned int>(varNames[iV], iV) );
      mvaVarMapEE.insert(  std::pair<std::string,unsigned int>(varNames[iV], iV) );
    }
    
    break;

  case k2012IdMVA:
  case k2012IdMVA_globe:
    
    EndcapWeights =      (gSystem->Getenv("CMSSW_BASE")+
			  TString("/src/MitPhysics/data/")+
			  TString("2012ICHEP_PhotonID_Endcap_BDT.")+
			  TString("weights_PSCorr.xml"));
    BarrelWeights =      (gSystem->Getenv("CMSSW_BASE")+
			  TString("/src/MitPhysics/data/")+
			  TString("2012ICHEP_PhotonID_Barrel_BDT.")+
			  TString("weights.xml"));

    mvaVars.resize(12);
    varNames.push_back("ph.r9"			      );
    varNames.push_back("ph.sigietaieta"		      );
    varNames.push_back("ph.scetawidth"		      );
    varNames.push_back("ph.scphiwidth"		      );
    varNames.push_back("ph.idmva_CoviEtaiPhi"	      );
    varNames.push_back("ph.idmva_s4ratio"	      );
    varNames.push_back("ph.idmva_GammaIso"	      );
    varNames.push_back("ph.idmva_ChargedIso_selvtx"   );
    varNames.push_back("ph.idmva_ChargedIso_worstvtx" );
    varNames.push_back("ph.sceta"		      );
    varNames.push_back("rho"			      );
    varNames.push_back("ph.idmva_PsEffWidthSigmaRR"   );

    for( unsigned int iV = 0; iV < mvaVars.size() - 1; ++iV) {
      mvaVarMapEB.insert(  std::pair<std::string,unsigned int>(varNames[iV], iV) );
      mvaVarMapEE.insert(  std::pair<std::string,unsigned int>(varNames[iV], iV) );
    }
    
    // pre-shower only used for Endcaps
    mvaVarMapEE.insert( std::pair<std::string,unsigned int> ( varNames[mvaVars.size() - 1], mvaVars.size() - 1) );

    break;

  case k2013FinalIdMVA_8TeV:
    
    EndcapWeights =      (gSystem->Getenv("CMSSW_BASE")+
			  TString("/src/MitPhysics/data/2013FinalPaper_PhotonID_Weight/8TeV/")+
			  TString("2013FinalPaper_PhotonID_Endcap_BDT_TrainRangePT15_8TeV.")+
			  TString("weights.xml"));
    BarrelWeights =      (gSystem->Getenv("CMSSW_BASE")+
			  TString("/src/MitPhysics/data/2013FinalPaper_PhotonID_Weight/8TeV/")+
			  TString("2013FinalPaper_PhotonID_Barrel_BDT_TrainRangePT15_8TeV.")+
			  TString("weights.xml"));

    mvaVars.resize(13);
    varNames.push_back("ph.scrawe"		      );
    varNames.push_back("ph.r9"			      );
    varNames.push_back("ph.sigietaieta"		      );
    varNames.push_back("ph.scetawidth"		      );
    varNames.push_back("ph.scphiwidth"		      );
    varNames.push_back("ph.idmva_CoviEtaiPhi"	      );
    varNames.push_back("ph.idmva_s4ratio"	      );
    varNames.push_back("ph.idmva_GammaIso"	      );
    varNames.push_back("ph.idmva_ChargedIso_selvtx"   );
    varNames.push_back("ph.idmva_ChargedIso_worstvtx" );
    varNames.push_back("ph.sceta"		      );
    varNames.push_back("rho"			      );
    varNames.push_back("ph.idmva_PsEffWidthSigmaRR"   );

    for( unsigned int iV = 0; iV < mvaVars.size() - 1; ++iV) {
      mvaVarMapEB.insert(  std::pair<std::string,unsigned int>(varNames[iV], iV) );
      mvaVarMapEE.insert(  std::pair<std::string,unsigned int>(varNames[iV], iV) );
    }
    
    // pre-shower only used for Endcaps
    mvaVarMapEE.insert( std::pair<std::string,unsigned int> ( varNames[mvaVars.size() - 1], mvaVars.size() - 1) );

    break;

  case k2013FinalIdMVA_7TeV:
    
    EndcapWeights =      (gSystem->Getenv("CMSSW_BASE")+
			  TString("/src/MitPhysics/data/2013FinalPaper_PhotonID_Weight/7TeV/")+
			  TString("2013FinalPaper_PhotonID_Endcap_BDT_TrainRangePT15_7TeV.")+
			  TString("weights.xml"));
    BarrelWeights =      (gSystem->Getenv("CMSSW_BASE")+
			  TString("/src/MitPhysics/data/2013FinalPaper_PhotonID_Weight/7TeV/")+
			  TString("2013FinalPaper_PhotonID_Barrel_BDT_TrainRangePT15_7TeV.")+
			  TString("weights.xml"));

    mvaVars.resize(13);
    varNames.push_back("ph.scrawe"		      );
    varNames.push_back("ph.r9"			      );
    varNames.push_back("ph.sigietaieta"		      );
    varNames.push_back("ph.scetawidth"		      );
    varNames.push_back("ph.scphiwidth"		      );
    varNames.push_back("ph.idmva_CoviEtaiPhi"	      );
    varNames.push_back("ph.idmva_s4ratio"	      );
    varNames.push_back("ph.idmva_GammaIso"	      );
    varNames.push_back("ph.idmva_ChargedIso_selvtx"   );
    varNames.push_back("ph.idmva_ChargedIso_worstvtx" );
    varNames.push_back("ph.sceta"		      );
    varNames.push_back("rho"			      );
    varNames.push_back("ph.idmva_PsEffWidthSigmaRR"   );

    for( unsigned int iV = 0; iV < mvaVars.size() - 1; ++iV) {
      mvaVarMapEB.insert(  std::pair<std::string,unsigned int>(varNames[iV], iV) );
      mvaVarMapEE.insert(  std::pair<std::string,unsigned int>(varNames[iV], iV) );
    }
    
    // pre-shower only used for Endcaps
    mvaVarMapEE.insert( std::pair<std::string,unsigned int> ( varNames[mvaVars.size() - 1], mvaVars.size() - 1) );

    break;
    
  default:
    // no variables... better never called..
    std::cerr<<" MVATools: Trying to initialize with unknown type."<<std::endl;
    break;
  }


  // looping over both maps and adding Vars to BDT readers
  for( unsigned int iV = 0; iV < varNames.size(); ++iV ){
    std::map<std::string,unsigned int>::const_iterator it = mvaVarMapEB.find(varNames[iV]);
    if ( it != mvaVarMapEB.end() )
      fReaderBarrel -> AddVariable( (it->first).c_str(), &(mvaVars[it->second]));
    it = mvaVarMapEE.find(varNames[iV]);
    if ( it != mvaVarMapEE.end() )
      fReaderEndcap -> AddVariable( (it->first).c_str(), &(mvaVars[it->second]));
  }

  fReaderEndcap->BookMVA("BDT method",EndcapWeights);
  fReaderBarrel->BookMVA("BDT method",BarrelWeights);
  
  assert(fReaderEndcap);
  assert(fReaderBarrel);
  
}

// //--------------------------------------------------------------------------------------------------
// void MVATools::InitializeMVA(int VariableType, TString EndcapWeights,TString BarrelWeights) {
  
//   if (fReaderEndcap) delete fReaderEndcap;  
//   if (fReaderBarrel) delete fReaderBarrel; 
  
//   fReaderEndcap = new TMVA::Reader( "!Color:!Silent:Error" );       
//   fReaderBarrel = new TMVA::Reader( "!Color:!Silent:Error" );   
  
//   TMVA::Reader *readers[2];
//   readers[0]  = fReaderEndcap;
//   readers[1]  = fReaderBarrel;  
  
//   for (UInt_t i=0; i<2; ++i) {

//     if(VariableType==0||VariableType==1||VariableType==2){
//       readers[i]->AddVariable( "HoE", &HoE );
//       readers[i]->AddVariable( "covIEtaIEta", &covIEtaIEta );
//       readers[i]->AddVariable( "tIso1", &tIso1 );
//       readers[i]->AddVariable( "tIso3", &tIso3 );
//       readers[i]->AddVariable( "tIso2", &tIso2 );
//       readers[i]->AddVariable( "R9", &R9 );
//     }
    
//     if(VariableType==3||VariableType==4){
//       readers[i]->AddVariable( "HoE", &HoE );
//       readers[i]->AddVariable( "covIEtaIEta", &covIEtaIEta );
//       readers[i]->AddVariable( "tIso1abs", &tIso1abs );
//       readers[i]->AddVariable( "tIso3abs", &tIso3abs );
//       readers[i]->AddVariable( "tIso2abs", &tIso2abs );
//       readers[i]->AddVariable( "R9", &R9 );
//     }
    
//     if(VariableType==1||VariableType==2){
//       readers[i]->AddVariable( "RelIsoEcal", &RelIsoEcal );
//       readers[i]->AddVariable( "RelIsoHcal", &RelIsoHcal );
//       readers[i]->AddVariable( "RelEMax", &RelEMax );
//       readers[i]->AddVariable( "RelETop", &RelETop );
//       readers[i]->AddVariable( "RelEBottom", &RelEBottom );
//       readers[i]->AddVariable( "RelELeft", &RelELeft );
//       readers[i]->AddVariable( "RelERight", &RelERight );
//       readers[i]->AddVariable( "RelE2x5Max", &RelE2x5Max );
//       readers[i]->AddVariable( "RelE2x5Top", &RelE2x5Top );
//       readers[i]->AddVariable( "RelE2x5Bottom", &RelE2x5Bottom );
//       readers[i]->AddVariable( "RelE2x5Left", &RelE2x5Left );
//       readers[i]->AddVariable( "RelE2x5Right", &RelE2x5Right );
//       readers[i]->AddVariable( "RelE5x5", &RelE5x5 );
//     }
    
//     if(VariableType==3||VariableType==4){
//       readers[i]->AddVariable( "absIsoEcal", &absIsoEcal );
//       readers[i]->AddVariable( "absIsoHcal", &absIsoHcal );
//       readers[i]->AddVariable( "RelEMax", &RelEMax );
//       readers[i]->AddVariable( "RelETop", &RelEBottom);
//       readers[i]->AddVariable( "RelEBottom", &RelEBottom );
//       readers[i]->AddVariable( "RelELeft", &RelELeft );
//       readers[i]->AddVariable( "RelERight", &RelERight );
//       readers[i]->AddVariable( "RelE2x5Max", &RelE2x5Max );
//       readers[i]->AddVariable( "RelE2x5Top", &RelE2x5Top );
//       readers[i]->AddVariable( "RelE2x5Bottom", &RelE2x5Bottom );
//       readers[i]->AddVariable( "RelE2x5Left", &RelE2x5Left );
//       readers[i]->AddVariable( "RelE2x5Right", &RelE2x5Right );
//       readers[i]->AddVariable( "RelE5x5", &RelE5x5 );
//     }
    
//     if(VariableType==2||VariableType==3||VariableType==4){
//       readers[i]->AddVariable( "EtaWidth", &EtaWidth );
//       readers[i]->AddVariable( "PhiWidth", &PhiWidth );
//       readers[i]->AddVariable( "CoviEtaiPhi", &CoviEtaiPhi );
//       readers[i]->AddVariable( "CoviPhiiPhi", &CoviPhiiPhi );
//       if(VariableType==4){
// 	readers[i]->AddVariable( "NVertexes", &NVertexes );
//       }
//       if(i==0){
// 	readers[i]->AddVariable( "RelPreshowerEnergy", &RelPreshowerEnergy );
//       }
//     }

//     if(VariableType==6){
//       readers[i]->AddVariable( "HoE", &HoE );
//       readers[i]->AddVariable( "covIEtaIEta", &covIEtaIEta );
//       readers[i]->AddVariable( "tIso1abs", &tIso1abs );
//       readers[i]->AddVariable( "tIso3abs", &tIso3abs );
//       readers[i]->AddVariable( "tIso2abs", &tIso2abs );
//       readers[i]->AddVariable( "R9", &R9 );
//       readers[i]->AddVariable( "absIsoEcal", &absIsoEcal );
//       readers[i]->AddVariable( "absIsoHcal", &absIsoHcal );
//       readers[i]->AddVariable( "RelE5x5", &RelE5x5 );
//       readers[i]->AddVariable( "EtaWidth", &EtaWidth );
//       readers[i]->AddVariable( "PhiWidth", &PhiWidth );
//       readers[i]->AddVariable( "CoviEtaiPhi", &CoviEtaiPhi );
//       readers[i]->AddVariable( "CoviPhiiPhi", &CoviPhiiPhi );
//       readers[i]->AddVariable( "NVertexes", &NVertexes );
//       if(i==0){
//       readers[i]->AddVariable( "RelPreshowerEnergy", &RelPreshowerEnergy );
//       }
//     }
    
//     if(VariableType==7){
//       readers[i]->AddVariable( "HoE", &HoE );
//       readers[i]->AddVariable( "covIEtaIEta", &covIEtaIEta );
//       readers[i]->AddVariable( "tIso1abs", &tIso1abs );
//       readers[i]->AddVariable( "tIso3abs", &tIso3abs );
//       readers[i]->AddVariable( "tIso2abs", &tIso2abs );
//       readers[i]->AddVariable( "R9", &R9 );
//       readers[i]->AddVariable( "absIsoEcal", &absIsoEcal );
//       readers[i]->AddVariable( "absIsoHcal", &absIsoHcal );
//       readers[i]->AddVariable( "NVertexes", &NVertexes );
//       readers[i]->AddVariable( "ScEta", &ScEta );
//     }    

//     if(VariableType==10){
//       readers[i]->AddVariable( "HoE", &HoE );
//       readers[i]->AddVariable( "covIEtaIEta", &covIEtaIEta );
//       readers[i]->AddVariable( "tIso1abs", &tIso1abs );
//       readers[i]->AddVariable( "tIso3abs", &tIso3abs );
//       readers[i]->AddVariable( "tIso2abs", &tIso2abs );
//       readers[i]->AddVariable( "R9", &R9 );
//       readers[i]->AddVariable( "absIsoEcal", &absIsoEcal );
//       readers[i]->AddVariable( "absIsoHcal", &absIsoHcal );
//       readers[i]->AddVariable( "NVertexes", &NVertexes );
//       readers[i]->AddVariable( "ScEta", &ScEta );
//       readers[i]->AddVariable( "EtaWidth", &EtaWidth );
//       readers[i]->AddVariable( "PhiWidth", &PhiWidth );      
//     }    

//     if(VariableType==1201){
//       /*readers[i]->AddVariable( "myphoton_pfchargedisogood03", &myphoton_pfchargedisogood03);
//       readers[i]->AddVariable( "myphoton_pfchargedisobad03", &myphoton_pfchargedisobad03); 
//       readers[i]->AddVariable( "myphoton_pfphotoniso03", &myphoton_pfphotoniso03 );
//       readers[i]->AddVariable( "myphoton_sieie", &myphoton_sieie ); 
//       readers[i]->AddVariable( "myphoton_sieip", &myphoton_sieip );
//       readers[i]->AddVariable( "myphoton_etawidth", &myphoton_etawidth ); 
//       readers[i]->AddVariable( "myphoton_phiwidth", &myphoton_phiwidth );
//       readers[i]->AddVariable( "myphoton_r9", &myphoton_r9 ); 
//       readers[i]->AddVariable( "myphoton_s4ratio", &myphoton_s4ratio );
//       readers[i]->AddVariable( "myphoton_SCeta", &myphoton_SCeta ); 
//       readers[i]->AddVariable( "event_rho", &event_rho );
//       if(i==0){
// 	readers[i]->AddVariable( "myphoton_ESEffSigmaRR", &myphoton_ESEffSigmaRR);
// 	}*/
//       readers[i]->AddVariable( "ph.r9", &myphoton_r9 ); 
//       readers[i]->AddVariable( "ph.sigietaieta", &myphoton_sieie ); 
//       readers[i]->AddVariable( "ph.scetawidth", &myphoton_etawidth ); 
//       readers[i]->AddVariable( "ph.scphiwidth", &myphoton_phiwidth );
//       readers[i]->AddVariable( "ph.idmva_CoviEtaiPhi", &myphoton_sieip );
//       readers[i]->AddVariable( "ph.idmva_s4ratio", &myphoton_s4ratio );
//       readers[i]->AddVariable( "ph.idmva_GammaIso", &myphoton_pfphotoniso03 );
//       readers[i]->AddVariable( "ph.idmva_ChargedIso_selvtx", &myphoton_pfchargedisogood03);
//       readers[i]->AddVariable( "ph.idmva_ChargedIso_worstvtx", &myphoton_pfchargedisobad03); 
//       readers[i]->AddVariable( "ph.sceta", &myphoton_SCeta ); 
//       readers[i]->AddVariable( "rho", &event_rho );
//       if(i==0){
// 	//readers[i]->AddVariable( "1.00023*ph.idmva_PsEffWidthSigmaRR + 0.0913", &myphoton_ESEffSigmaRR);
//         readers[i]->AddVariable( "ph.idmva_PsEffWidthSigmaRR", &myphoton_ESEffSigmaRR);
//       }
//     }
    
//   }
  
//   fReaderEndcap->BookMVA("BDT method",EndcapWeights);
//   fReaderBarrel->BookMVA("BDT method",BarrelWeights);
  
//   assert(fReaderEndcap);
//   assert(fReaderBarrel);
  
// }

// *** REMOVED THIS COMPLETELY. If a module wants to cut on BDT, it should 1,) compute the BDT value (using GetMVAbdtValue(...) ) and then make the cut itself...

// Bool_t MVATools::PassMVASelection(const Photon* p,const Vertex* vtx,const TrackCol* trackCol,const VertexCol* vtxCol,Double_t _tRho,Float_t bdtCutBarrel, Float_t bdtCutEndcap, const ElectronCol* els, Bool_t applyElectronVeto) {
  
//   //initilize the bool value
//   PassMVA=kFALSE;
  
//   Float_t photon_bdt =  MVATools::GetMVAbdtValue_2011(p,vtx,trackCol,vtxCol, _tRho, els, applyElectronVeto);
  
//   if (isbarrel) {
//     if(photon_bdt>bdtCutBarrel){
//       PassMVA=kTRUE;	
//     }
//   }
//   else {
//     if(photon_bdt>bdtCutEndcap){
//       PassMVA=kTRUE;	
//     }
//   } 
//   return PassMVA;
// }

// //---------------------------------------------------------------------------------
// Int_t MVATools::PassElectronVetoInt(const Photon* p, const ElectronCol* els) {
  
//   // these values are taken from the H2GGlobe code... (actually from Marco/s mail)
//   float cic4_allcuts_temp_sublead[] = { 
//     3.8,         2.2,         1.77,        1.29,
//     11.7,        3.4,         3.9,         1.84,
//     3.5,         2.2,         2.3,         1.45,
//     0.0106,      0.0097,      0.028,       0.027,
//     0.082,       0.062,       0.065,       0.048,
//     0.94,        0.36,        0.94,        0.32,
//     1.,          0.062,       0.97,        0.97,
//     1.5,         1.5,         1.5,         1.5 };  // the last line is PixelmatchVeto and un-used
  
//   //initilize the bool value
//   PassElecVetoInt=0;
  
//   dRTrack = PhotonTools::ElectronVetoCiC(p, els);
  
//   ScEta_MVA=p->SCluster()->Eta();
  
//   isbarrel = (fabs(ScEta_MVA)<1.4442);

//   R9 = p->R9();
//   //R9 = p->E33()/p->SCluster()->RawEnergy();
  
//   // check which category it is ...
//   _tCat = 1;
//   if ( !isbarrel ) _tCat = 3;
//   if ( R9 < 0.94 ) _tCat++;
  
//   //Electron Veto
//   if(dRTrack > cic4_allcuts_temp_sublead[_tCat-1+6*4]){
//     PassElecVetoInt=1;
//   }
  
//   return  PassElecVetoInt;
  
// }

//--------------------------------------------------------------------------------------------------


Double_t MVATools::GetMVAbdtValue(const Photon* p, const Vertex* vtx, const TrackCol* trackCol, const VertexCol* vtxCol, Double_t _tRho, const PFCandidateCol *fPFCands, const ElectronCol* els, Bool_t applyElectronVeto) {

  // if there's no reader, or the type is kNone, return the default values of -99.
  //if( ( !fReaderBarrel || !fReaderEndcap ) ) return -199.;
  if( ( !fReaderBarrel || !fReaderEndcap ) || fMVAType == kNone ) return -99.;

  // we compute the variable names... make sure no confilcts when adding new variables...
  
  // check if it's a Barrel or EE photon
  bool isBarrel = ( p->SCluster()->AbsEta() < 1.5 );

  std::map<std::string,unsigned int>* theVars = ( isBarrel ? &mvaVarMapEB : &mvaVarMapEE );  
  
  // loop over all the variables in the map... and keep count (to make sure we have filled all variables)
  unsigned int varCounter = 0;
  for( std::map<std::string,unsigned int>::const_iterator iV = theVars->begin(); iV != theVars->end(); ++iV ) {
    
    TString theVarName  = TString(iV->first);  
    float* theVarValue  = &(mvaVars[iV->second]);           // pointer to the variable...
    
    if( 
       !theVarName.CompareTo("HoE") 
       ) {
      (*theVarValue) = p->HadOverEm();
      varCounter++;
    } else if (
	       !theVarName.CompareTo("covIEtaIEta") || !theVarName.CompareTo("ph.sigietaieta") || !theVarName.CompareTo("sigieie") 
	       ) {
      (*theVarValue) = p->CoviEtaiEta();
      varCounter++;
    } else if (
	       !theVarName.CompareTo("R9") || !theVarName.CompareTo("ph.r9") || !theVarName.CompareTo("r9")
	       ) {
      (*theVarValue) = p->R9();
      varCounter++;
    } else if (
	       !theVarName.CompareTo("ScEta") || !theVarName.CompareTo("ph.sceta") || !theVarName.CompareTo("sceta")
	       ) {
      (*theVarValue) = p->SCluster()->Eta();
      varCounter++;
    } else if (
	       !theVarName.CompareTo("rho")
	       ) {
      (*theVarValue) = _tRho;
       varCounter++;
    } else if (
	       !theVarName.CompareTo("tIso1abs")
	       ) {      
      double _ecalIso3 = p->EcalRecHitIsoDr03();
      double _hcalIso4 = p->HcalTowerSumEtDr04();  
      double _trackIso1 = IsolationTools::CiCTrackIsolation(p,vtx, 0.3, 0.02, 0.0, 0.0, 0.1, 1.0, trackCol, NULL, NULL, (!applyElectronVeto ? els : NULL) );//Question Ming:whyfPV->At(0) instead of selected vertex using ranking method?
      (*theVarValue) = _ecalIso3+_hcalIso4+_trackIso1 - 0.17*_tRho;
      varCounter++;
    } else if (
	       !theVarName.CompareTo("tIso2abs")
	       ) {
      double _ecalIso4 = p->EcalRecHitIsoDr04();
      double _hcalIso4 = p->HcalTowerSumEtDr04();
      unsigned int wVtxInd = 0;
      double _trackIso2 = IsolationTools::CiCTrackIsolation(p,vtx, 0.4, 0.02, 0.0, 0.0, 0.1, 1.0, trackCol, &wVtxInd,vtxCol, (!applyElectronVeto ? els : NULL) );
      (*theVarValue) = _ecalIso4+_hcalIso4+_trackIso2 - 0.52*_tRho;
      varCounter++;
    } else if (
	       !theVarName.CompareTo("tIso3abs")
	       ) {
      (*theVarValue) = IsolationTools::CiCTrackIsolation(p,vtx, 0.3, 0.02, 0.0, 0.0, 0.1, 1.0, trackCol, NULL, NULL, (!applyElectronVeto ? els : NULL) );
      varCounter++;
    } else if (
	       !theVarName.CompareTo("absIsoEcal")
	       ) {
      double _ecalIso3 = p->EcalRecHitIsoDr03();
      (*theVarValue) =  (_ecalIso3-0.17*_tRho);
      varCounter++;
    } else if (
	       !theVarName.CompareTo("absIsoHcal")
	       ) {
      double _hcalIso4 = p->HcalTowerSumEtDr04();
      (*theVarValue) = (_hcalIso4-0.17*_tRho);
      varCounter++;
    } else if (
	       !theVarName.CompareTo("NVertexes")
	       ) {

      (*theVarValue) = vtxCol->GetEntries();
      varCounter++;
    } else if (
	       !theVarName.CompareTo("EtaWidth") || !theVarName.CompareTo("ph.scetawidth") || !theVarName.CompareTo("sigeta")
	       ) {
      (*theVarValue) = p->EtaWidth();
      varCounter++;
    } else if (
	       !theVarName.CompareTo("PhiWidth") || !theVarName.CompareTo("ph.scphiwidth") || !theVarName.CompareTo("sigphi")
	       ) {
      (*theVarValue) = p->PhiWidth();
      varCounter++;

    } else if (
	       !theVarName.CompareTo("ph.idmva_CoviEtaiPhi") || !theVarName.CompareTo("covieip")
	       ) {
      (*theVarValue) = p->SCluster()->Seed()->CoviEtaiPhi();
      varCounter++;
    } else if (
	       !theVarName.CompareTo("ph.idmva_s4ratio") || !theVarName.CompareTo("s4r") 
	       ) {
      (*theVarValue) = p->S4Ratio();
      //(*theVarValue) = p->SCluster()->Seed()->E2x2()/p->E55();
      varCounter++;
    } else if (
	       !theVarName.CompareTo("ph.idmva_GammaIso") || !theVarName.CompareTo("pfgiso") 
	       ) {
      (*theVarValue) = IsolationTools::PFGammaIsolation(p,0.3,0,fPFCands);
      varCounter++;
    } else if (
	       !theVarName.CompareTo("ph.idmva_ChargedIso_selvtx") || !theVarName.CompareTo("pfciso") 
	       ) {
      (*theVarValue) = IsolationTools::PFChargedIsolation(p,vtx,0.3,0,fPFCands);
      varCounter++;
    } else if (
	       !theVarName.CompareTo("ph.idmva_ChargedIso_worstvtx")
	       ) {
      unsigned int wVtxInd = 0;
      (*theVarValue) = IsolationTools::PFChargedIsolation(p,vtx,0.3,0,fPFCands,&wVtxInd,vtxCol); 
      varCounter++;
    } else if (
	       !theVarName.CompareTo("ph.idmva_PsEffWidthSigmaRR")
	       ) {
      (*theVarValue) = p->EffSigmaRR();
      varCounter++;
    } else if (
	       !theVarName.CompareTo("rawe")|| !theVarName.CompareTo("ph.scrawe") 
	       ) {
      (*theVarValue) = p->SCluster()->RawEnergy();
      varCounter++;
    } else {
      // a variable is not know... copmplain!
      std::cerr<<" ERROR: MVA Evaluation called with unknown variable name >"<<theVarName<<">."<<std::endl;
    }
  }
  
  // now all the variables should be filled... check!
  if( varCounter != theVars->size() )
    std::cerr<<" ERROR: MVA Evaludation called and not all variables are filled."<<std::endl;

  // we're ready to compute the MVA value
  TMVA::Reader* reader = NULL;
  if (isBarrel)
    reader = fReaderBarrel;
  else
    reader = fReaderEndcap;
  
  assert(reader); 
  
  return (reader->EvaluateMVA("BDT method"));
}

// Float_t MVATools::GetMVAbdtValue_2012_globe(const Photon* p,const Vertex* vtx,const TrackCol* trackCol,const VertexCol* vtxCol,Double_t _tRho, const PFCandidateCol *fPFCands,const ElectronCol* els,Bool_t applyElectronVeto) {
  
//   //get the variables used to compute MVA variables
//   ecalIso3 = p->EcalRecHitIsoDr03();
//   ecalIso4 = p->EcalRecHitIsoDr04();
//   hcalIso4 = p->HcalTowerSumEtDr04();
  
//   wVtxInd = 0;
  
//   trackIso1 = IsolationTools::CiCTrackIsolation(p,vtx, 0.3, 0.02, 0.0, 0.0, 0.1, 1.0, trackCol, NULL, NULL, (!applyElectronVeto ? els : NULL) );//Question Ming:whyfPV->At(0) instead of selected vertex using ranking method?
    
//   // track iso worst vtx
//   trackIso2 = IsolationTools::CiCTrackIsolation(p,vtx, 0.4, 0.02, 0.0, 0.0, 0.1, 1.0, trackCol, &wVtxInd,vtxCol, (!applyElectronVeto ? els : NULL) );
  
//   combIso1 = ecalIso3+hcalIso4+trackIso1 - 0.17*_tRho;
//   combIso2 = ecalIso4+hcalIso4+trackIso2 - 0.52*_tRho;
  
//   RawEnergy = p->SCluster()->RawEnergy();
  
//   ScEta = p->SCluster()->Eta();
  
//   //mva varialbes v1 and v2
//   tIso1 = (combIso1) *50./p->Et();
//   tIso3 = (trackIso1)*50./p->Et();
//   tIso2 = (combIso2) *50./(p->MomVtx(vtxCol->At(wVtxInd)->Position()).Pt());
//   RelIsoEcal=(ecalIso3-0.17*_tRho)/p->Et();
//   RelIsoHcal=(hcalIso4-0.17*_tRho)/p->Et();

//   //compute mva variables for v3
//   HoE = p->HadOverEm();
//   covIEtaIEta = p->CoviEtaiEta();
//   tIso1abs = combIso1;
//   tIso3abs = trackIso1;
//   tIso2abs = combIso2;
//   R9 = p->R9();

//   absIsoEcal=(ecalIso3-0.17*_tRho);
//   absIsoHcal=(hcalIso4-0.17*_tRho);
//   RelEMax=p->SCluster()->Seed()->EMax()/RawEnergy;
//   RelETop=p->SCluster()->Seed()->ETop()/RawEnergy;
//   RelEBottom=p->SCluster()->Seed()->EBottom()/RawEnergy;
//   RelELeft=p->SCluster()->Seed()->ELeft()/RawEnergy;
//   RelERight=p->SCluster()->Seed()->ERight()/RawEnergy;
//   RelE2x5Max=p->SCluster()->Seed()->E2x5Max()/RawEnergy;
//   RelE2x5Top=p->SCluster()->Seed()->E2x5Top()/RawEnergy;
//   RelE2x5Bottom=p->SCluster()->Seed()->E2x5Bottom()/RawEnergy;
//   RelE2x5Left=p->SCluster()->Seed()->E2x5Left()/RawEnergy;
//   RelE2x5Right=p->SCluster()->Seed()->E2x5Right()/RawEnergy;
//   RelE5x5=p->SCluster()->Seed()->E5x5()/RawEnergy;
  
//   EtaWidth=p->EtaWidth();
//   PhiWidth=p->PhiWidth();
//   CoviEtaiPhi=p->SCluster()->Seed()->CoviEtaiPhi(); 
//   CoviPhiiPhi=p->SCluster()->Seed()->CoviPhiiPhi();

//   RelPreshowerEnergy=p->SCluster()->PreshowerEnergy()/RawEnergy;
//   NVertexes=vtxCol->GetEntries();
  
//   //spectator variables
//   Pt_MVA=p->Pt();
//   ScEta_MVA=p->SCluster()->Eta();

//   //

//   isbarrel = (fabs(ScEta_MVA)<1.4442);
  
//   //variable 1201
//   myphoton_pfchargedisogood03=IsolationTools::PFChargedIsolation(p,vtx,0.3,0,fPFCands);
//   myphoton_pfchargedisobad03=IsolationTools::PFChargedIsolation(p,vtx,0.3,0,fPFCands,&wVtxInd,vtxCol); 
//   myphoton_pfphotoniso03=IsolationTools::PFGammaIsolation(p,0.3,0,fPFCands);
//   myphoton_sieie=covIEtaIEta; 
//   myphoton_sieip=CoviEtaiPhi;
//   myphoton_etawidth=EtaWidth; 
//   myphoton_phiwidth=PhiWidth;
//   myphoton_r9=R9; 
//   myphoton_s4ratio=p->S4Ratio();
//   myphoton_SCeta=ScEta_MVA; 
//   event_rho= _tRho;

//   myphoton_ESEffSigmaRR=-99;

//   if(!isbarrel){
//     myphoton_ESEffSigmaRR=p->EffSigmaRR();
//   }
  
//   if (isbarrel) {
//     reader = fReaderBarrel;
//   }
//   else {
//     reader = fReaderEndcap;
//   }
  
//   assert(reader); 

//   double bdt = reader->EvaluateMVA("BDT method");

//   /* printf("HoE: %f\n",HoE);
//   printf("covIEtaIEta: %f\n",covIEtaIEta);
//   printf("tIso1abs: %f\n",tIso1abs);
//   printf("tIso3abs: %f\n",tIso3abs);
//   printf("tIso2abs: %f\n",tIso2abs);
  
//   printf("absIsoEcal: %f\n",absIsoEcal);
//   printf("absIsoHcal: %f\n",absIsoHcal);
//   printf("RelEMax: %f\n",RelEMax);
//   printf("RelETop: %f\n",RelETop);
//   printf("RelEBottom: %f\n",RelEBottom);
//   printf("RelELeft: %f\n",RelELeft);
//   printf("RelERight: %f\n",RelERight);
//   printf("RelE2x5Max: %f\n",RelE2x5Max);
//   printf("RelE2x5Top: %f\n",RelE2x5Top);
//   printf("RelE2x5Bottom: %f\n",RelE2x5Bottom);
//   printf("RelE2x5Left: %f\n",RelE2x5Left);
//   printf("RelE2x5Right;: %f\n",RelE2x5Right);
//   printf("RelE5x5: %f\n",RelE5x5);
  
//   printf("EtaWidth: %f\n",EtaWidth);
//   printf("PhiWidth: %f\n",PhiWidth);
//   printf("CoviEtaiPhi: %f\n",CoviEtaiPhi);
//   printf("CoviPhiiPhi: %f\n",CoviPhiiPhi);
  
//   if (!isbarrel) {
//     printf("RelPreshowerEnergy: %f\n",RelPreshowerEnergy);
//     }*/
  
//   return bdt;
// }

// Float_t MVATools::GetMVAbdtValue_2011(const Photon* p,const Vertex* vtx,const TrackCol* trackCol,const VertexCol* vtxCol,Double_t _tRho,const ElectronCol* els,Bool_t applyElectronVeto) {
  
//   //get the variables used to compute MVA variables
//   ecalIso3 = p->EcalRecHitIsoDr03();
//   ecalIso4 = p->EcalRecHitIsoDr04();
//   hcalIso4 = p->HcalTowerSumEtDr04();
  
//   wVtxInd = 0;
  
//   trackIso1 = IsolationTools::CiCTrackIsolation(p,vtx, 0.3, 0.02, 0.0, 0.0, 0.1, 1.0, trackCol, NULL, NULL, (!applyElectronVeto ? els : NULL) );//Question Ming:whyfPV->At(0) instead of selected vertex using ranking method?
    
//   // track iso worst vtx
//   trackIso2 = IsolationTools::CiCTrackIsolation(p,vtx, 0.4, 0.02, 0.0, 0.0, 0.1, 1.0, trackCol, &wVtxInd,vtxCol, (!applyElectronVeto ? els : NULL) );
  
//   combIso1 = ecalIso3+hcalIso4+trackIso1 - 0.17*_tRho;
//   combIso2 = ecalIso4+hcalIso4+trackIso2 - 0.52*_tRho;
  
//   RawEnergy = p->SCluster()->RawEnergy();
  
//   ScEta = p->SCluster()->Eta();
  
//   //mva varialbes v1 and v2
//   tIso1 = (combIso1) *50./p->Et();
//   tIso3 = (trackIso1)*50./p->Et();
//   tIso2 = (combIso2) *50./(p->MomVtx(vtxCol->At(wVtxInd)->Position()).Pt());
//   RelIsoEcal=(ecalIso3-0.17*_tRho)/p->Et();
//   RelIsoHcal=(hcalIso4-0.17*_tRho)/p->Et();

//   //compute mva variables for v3
//   HoE = p->HadOverEm();
//   covIEtaIEta = p->CoviEtaiEta();
//   tIso1abs = combIso1;
//   tIso3abs = trackIso1;
//   tIso2abs = combIso2;
//   R9 = p->R9();

//   absIsoEcal=(ecalIso3-0.17*_tRho);
//   absIsoHcal=(hcalIso4-0.17*_tRho);
//   RelEMax=p->SCluster()->Seed()->EMax()/RawEnergy;
//   RelETop=p->SCluster()->Seed()->ETop()/RawEnergy;
//   RelEBottom=p->SCluster()->Seed()->EBottom()/RawEnergy;
//   RelELeft=p->SCluster()->Seed()->ELeft()/RawEnergy;
//   RelERight=p->SCluster()->Seed()->ERight()/RawEnergy;
//   RelE2x5Max=p->SCluster()->Seed()->E2x5Max()/RawEnergy;
//   RelE2x5Top=p->SCluster()->Seed()->E2x5Top()/RawEnergy;
//   RelE2x5Bottom=p->SCluster()->Seed()->E2x5Bottom()/RawEnergy;
//   RelE2x5Left=p->SCluster()->Seed()->E2x5Left()/RawEnergy;
//   RelE2x5Right=p->SCluster()->Seed()->E2x5Right()/RawEnergy;
//   RelE5x5=p->SCluster()->Seed()->E5x5()/RawEnergy;
  
//   EtaWidth=p->EtaWidth();
//   PhiWidth=p->PhiWidth();
//   CoviEtaiPhi=p->SCluster()->Seed()->CoviEtaiPhi(); 
//   CoviPhiiPhi=p->SCluster()->Seed()->CoviPhiiPhi();

//   RelPreshowerEnergy=p->SCluster()->PreshowerEnergy()/RawEnergy;
//   NVertexes=vtxCol->GetEntries();
  
//   //spectator variables
//   Pt_MVA=p->Pt();
//   ScEta_MVA=p->SCluster()->Eta();

//   //

//   isbarrel = (fabs(ScEta_MVA)<1.4442);
  
//   if (isbarrel) {
//     reader = fReaderBarrel;
//   }
//   else {
//     reader = fReaderEndcap;
//   }
  
//   assert(reader); 

//   double bdt = reader->EvaluateMVA("BDT method");

//   /* printf("HoE: %f\n",HoE);
//   printf("covIEtaIEta: %f\n",covIEtaIEta);
//   printf("tIso1abs: %f\n",tIso1abs);
//   printf("tIso3abs: %f\n",tIso3abs);
//   printf("tIso2abs: %f\n",tIso2abs);
  
//   printf("absIsoEcal: %f\n",absIsoEcal);
//   printf("absIsoHcal: %f\n",absIsoHcal);
//   printf("RelEMax: %f\n",RelEMax);
//   printf("RelETop: %f\n",RelETop);
//   printf("RelEBottom: %f\n",RelEBottom);
//   printf("RelELeft: %f\n",RelELeft);
//   printf("RelERight: %f\n",RelERight);
//   printf("RelE2x5Max: %f\n",RelE2x5Max);
//   printf("RelE2x5Top: %f\n",RelE2x5Top);
//   printf("RelE2x5Bottom: %f\n",RelE2x5Bottom);
//   printf("RelE2x5Left: %f\n",RelE2x5Left);
//   printf("RelE2x5Right;: %f\n",RelE2x5Right);
//   printf("RelE5x5: %f\n",RelE5x5);
  
//   printf("EtaWidth: %f\n",EtaWidth);
//   printf("PhiWidth: %f\n",PhiWidth);
//   printf("CoviEtaiPhi: %f\n",CoviEtaiPhi);
//   printf("CoviPhiiPhi: %f\n",CoviPhiiPhi);
  
//   if (!isbarrel) {
//     printf("RelPreshowerEnergy: %f\n",RelPreshowerEnergy);
//     }*/
  
//   return bdt;
// }
