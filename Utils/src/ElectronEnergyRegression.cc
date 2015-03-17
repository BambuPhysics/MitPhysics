#include "MitPhysics/Utils/interface/ElectronEnergyRegression.h"
#include "CondFormats/EgammaObjects/interface/GBRForest.h"
#include <cmath>
#include <cassert>
#include <iostream>

ClassImp(mithep::ElectronEnergyRegression)

using namespace mithep;


ElectronEnergyRegression::ElectronEnergyRegression() : 
  fIsInitialized(kFALSE),
  fVersionType(kNoTrkVar),
  forestCorrection_eb(0), 
  forestCorrection_ee(0), 
  forestUncertainty_eb(0), 
  forestUncertainty_ee(0) {
}

ElectronEnergyRegression::~ElectronEnergyRegression() {}
// Destructor does nothing


void ElectronEnergyRegression::initialize(std::string weightsFile, 
                                                  ElectronEnergyRegression::ElectronEnergyRegressionType type) {

  // Loading forest object according to different versions
  TFile file(weightsFile.c_str());

  if (type == kNoTrkVar || type == kWithTrkVarV1 || type == kWithTrkVarV2) {
    forestCorrection_eb = (GBRForest*) file.Get("EBCorrection");
    forestCorrection_ee = (GBRForest*) file.Get("EECorrection");
    forestUncertainty_eb = (GBRForest*) file.Get("EBUncertainty");
    forestUncertainty_ee = (GBRForest*) file.Get("EEUncertainty");

    // Just checking
    assert(forestCorrection_eb);
    assert(forestCorrection_ee);
    assert(forestUncertainty_eb);
    assert(forestUncertainty_ee);
  }

  // Updating type and marking as initialized
  fVersionType = type;
  fIsInitialized = kTRUE;
}


double ElectronEnergyRegression::regressionValueNoTrkVar(
                                                                 double SCRawEnergy,
                                                                 double scEta,
                                                                 double scPhi,
                                                                 double R9,
                                                                 double etawidth,
                                                                 double phiwidth,
                                                                 double NClusters,
                                                                 double HoE,
                                                                 double rho,
                                                                 double vertices,
                                                                 double EtaSeed,
                                                                 double PhiSeed,
                                                                 double ESeed,
                                                                 double E3x3Seed,
                                                                 double E5x5Seed,
                                                                 double see,
                                                                 double spp,
                                                                 double sep,
                                                                 double EMaxSeed,
                                                                 double E2ndSeed,
                                                                 double ETopSeed,
                                                                 double EBottomSeed,
                                                                 double ELeftSeed,
                                                                 double ERightSeed,
                                                                 double E2x5MaxSeed,
                                                                 double E2x5TopSeed,
                                                                 double E2x5BottomSeed,
                                                                 double E2x5LeftSeed,
                                                                 double E2x5RightSeed,
                                                                 double IEtaSeed,
                                                                 double IPhiSeed,
                                                                 double EtaCrySeed,
                                                                 double PhiCrySeed,
                                                                 double PreShowerOverRaw, 
                                                                 bool printDebug) 
{
  // Checking if instance has been initialized
  if (fIsInitialized == kFALSE) {
    printf("ElectronEnergyRegression instance not initialized !!!");
    return 0;
  }

  // Checking if type is correct
  if (!(fVersionType == kNoTrkVar)) {
    std::cout << "Error: Regression VersionType " << fVersionType << " is not supported to use function regressionValueNoTrkVar.\n";
    return 0;
  }

  // Now applying regression according to version and (endcap/barrel)
  float *vals = (fabs(scEta) <= 1.479) ? new float[38] : new float[31];
  if (fabs(scEta) <= 1.479) {		// Barrel
    vals[0]  = SCRawEnergy;
    vals[1]  = scEta;
    vals[2]  = scPhi;
    vals[3]  = R9;
    vals[4]  = E5x5Seed/SCRawEnergy;
    vals[5]  = etawidth;
    vals[6]  = phiwidth;
    vals[7]  = NClusters;
    vals[8]  = HoE;
    vals[9]  = rho;
    vals[10] = vertices;
    vals[11] = EtaSeed - scEta;
    vals[12] = atan2(sin(PhiSeed-scPhi),cos(PhiSeed-scPhi));
    vals[13] = ESeed/SCRawEnergy;
    vals[14] = E3x3Seed/ESeed;
    vals[15] = E5x5Seed/ESeed;
    vals[16] = see;
    vals[17] = spp;
    vals[18] = sep;
    vals[19] = EMaxSeed/ESeed;
    vals[20] = E2ndSeed/ESeed;
    vals[21] = ETopSeed/ESeed;
    vals[22] = EBottomSeed/ESeed;
    vals[23] = ELeftSeed/ESeed;
    vals[24] = ERightSeed/ESeed;
    vals[25] = E2x5MaxSeed/ESeed;
    vals[26] = E2x5TopSeed/ESeed;
    vals[27] = E2x5BottomSeed/ESeed;
    vals[28] = E2x5LeftSeed/ESeed;
    vals[29] = E2x5RightSeed/ESeed;
    vals[30] = IEtaSeed;
    vals[31] = IPhiSeed;
    vals[32] = ((int) IEtaSeed)%5;
    vals[33] = ((int) IPhiSeed)%2;
    vals[34] = (abs(IEtaSeed)<=25)*(((int)IEtaSeed)%25) + (abs(IEtaSeed)>25)*(((int) (IEtaSeed-25*abs(IEtaSeed)/IEtaSeed))%20);
    vals[35] = ((int) IPhiSeed)%20;
    vals[36] = EtaCrySeed;
    vals[37] = PhiCrySeed;
  }
  else {	// Endcap
    vals[0]  = SCRawEnergy;
    vals[1]  = scEta;
    vals[2]  = scPhi;
    vals[3]  = R9;
    vals[4]  = E5x5Seed/SCRawEnergy;
    vals[5]  = etawidth;
    vals[6]  = phiwidth;
    vals[7]  = NClusters;
    vals[8]  = HoE;
    vals[9]  = rho;
    vals[10] = vertices;
    vals[11] = EtaSeed - scEta;
    vals[12] = atan2(sin(PhiSeed-scPhi),cos(PhiSeed-scPhi));
    vals[13] = ESeed/SCRawEnergy;
    vals[14] = E3x3Seed/ESeed;
    vals[15] = E5x5Seed/ESeed;
    vals[16] = see;
    vals[17] = spp;
    vals[18] = sep;
    vals[19] = EMaxSeed/ESeed;
    vals[20] = E2ndSeed/ESeed;
    vals[21] = ETopSeed/ESeed;
    vals[22] = EBottomSeed/ESeed;
    vals[23] = ELeftSeed/ESeed;
    vals[24] = ERightSeed/ESeed;
    vals[25] = E2x5MaxSeed/ESeed;
    vals[26] = E2x5TopSeed/ESeed;
    vals[27] = E2x5BottomSeed/ESeed;
    vals[28] = E2x5LeftSeed/ESeed;
    vals[29] = E2x5RightSeed/ESeed;
    vals[30] = PreShowerOverRaw;
  }

  // Now evaluating the regression
  double regressionResult = 0;
  Int_t BinIndex = -1;

  if (fVersionType == kNoTrkVar) {
    if (fabs(scEta) <= 1.479) { 
      regressionResult = SCRawEnergy * forestCorrection_eb->GetResponse(vals); 
      BinIndex = 0;
    }
    else {
      regressionResult = (SCRawEnergy*(1+PreShowerOverRaw)) * forestCorrection_ee->GetResponse(vals);
      BinIndex = 1;
    }
  }

  //print debug
  if (printDebug) {    
    if ( fabs(scEta) <= 1.479) {
      std::cout << "Barrel :";
      for (uint v=0; v < 38; ++v) std::cout << vals[v] << ", ";
      std::cout << "\n";
    }
    else {
      std::cout << "Endcap :";
      for (uint v=0; v < 32; ++v) std::cout << vals[v] << ", ";
      std::cout << "\n";
    }
    std::cout << "BinIndex : " << BinIndex << "\n";
    std::cout << "SCRawEnergy = " << SCRawEnergy << " : PreShowerOverRaw = " << PreShowerOverRaw << std::endl;
    std::cout << "regression energy = " << regressionResult << std::endl;
  }
  

  // Cleaning up and returning
  delete[] vals;
  return regressionResult;
}

double ElectronEnergyRegression::regressionUncertaintyNoTrkVar(
                                                                       double SCRawEnergy,
                                                                       double scEta,
                                                                       double scPhi,
                                                                       double R9,
                                                                       double etawidth,
                                                                       double phiwidth,
                                                                       double NClusters,
                                                                       double HoE,
                                                                       double rho,
                                                                       double vertices,
                                                                       double EtaSeed,
                                                                       double PhiSeed,
                                                                       double ESeed,
                                                                       double E3x3Seed,
                                                                       double E5x5Seed,
                                                                       double see,
                                                                       double spp,
                                                                       double sep,
                                                                       double EMaxSeed,
                                                                       double E2ndSeed,
                                                                       double ETopSeed,
                                                                       double EBottomSeed,
                                                                       double ELeftSeed,
                                                                       double ERightSeed,
                                                                       double E2x5MaxSeed,
                                                                       double E2x5TopSeed,
                                                                       double E2x5BottomSeed,
                                                                       double E2x5LeftSeed,
                                                                       double E2x5RightSeed,
                                                                       double IEtaSeed,
                                                                       double IPhiSeed,
                                                                       double EtaCrySeed,
                                                                       double PhiCrySeed,
                                                                       double PreShowerOverRaw, 
                                                                       bool printDebug) 
{
  // Checking if instance has been initialized
  if (fIsInitialized == kFALSE) {
    printf("ElectronEnergyRegression instance not initialized !!!");
    return 0;
  }

  // Checking if type is correct
  if (!(fVersionType == kNoTrkVar)) {
    std::cout << "Error: Regression VersionType " << fVersionType << " is not supported to use function regressionValueNoTrkVar.\n";
    return 0;
  }

  // Now applying regression according to version and (endcap/barrel)
  float *vals = (fabs(scEta) <= 1.479) ? new float[38] : new float[31];
  if (fabs(scEta) <= 1.479) {		// Barrel
    vals[0]  = SCRawEnergy;
    vals[1]  = scEta;
    vals[2]  = scPhi;
    vals[3]  = R9;
    vals[4]  = E5x5Seed/SCRawEnergy;
    vals[5]  = etawidth;
    vals[6]  = phiwidth;
    vals[7]  = NClusters;
    vals[8]  = HoE;
    vals[9]  = rho;
    vals[10] = vertices;
    vals[11] = EtaSeed - scEta;
    vals[12] = atan2(sin(PhiSeed-scPhi),cos(PhiSeed-scPhi));
    vals[13] = ESeed/SCRawEnergy;
    vals[14] = E3x3Seed/ESeed;
    vals[15] = E5x5Seed/ESeed;
    vals[16] = see;
    vals[17] = spp;
    vals[18] = sep;
    vals[19] = EMaxSeed/ESeed;
    vals[20] = E2ndSeed/ESeed;
    vals[21] = ETopSeed/ESeed;
    vals[22] = EBottomSeed/ESeed;
    vals[23] = ELeftSeed/ESeed;
    vals[24] = ERightSeed/ESeed;
    vals[25] = E2x5MaxSeed/ESeed;
    vals[26] = E2x5TopSeed/ESeed;
    vals[27] = E2x5BottomSeed/ESeed;
    vals[28] = E2x5LeftSeed/ESeed;
    vals[29] = E2x5RightSeed/ESeed;
    vals[30] = IEtaSeed;
    vals[31] = IPhiSeed;
    vals[32] = ((int) IEtaSeed)%5;
    vals[33] = ((int) IPhiSeed)%2;
    vals[34] = (abs(IEtaSeed)<=25)*(((int)IEtaSeed)%25) + (abs(IEtaSeed)>25)*(((int) (IEtaSeed-25*abs(IEtaSeed)/IEtaSeed))%20);
    vals[35] = ((int) IPhiSeed)%20;
    vals[36] = EtaCrySeed;
    vals[37] = PhiCrySeed;
  }
  else {	// Endcap
    vals[0]  = SCRawEnergy;
    vals[1]  = scEta;
    vals[2]  = scPhi;
    vals[3]  = R9;
    vals[4]  = E5x5Seed/SCRawEnergy;
    vals[5]  = etawidth;
    vals[6]  = phiwidth;
    vals[7]  = NClusters;
    vals[8]  = HoE;
    vals[9]  = rho;
    vals[10] = vertices;
    vals[11] = EtaSeed - scEta;
    vals[12] = atan2(sin(PhiSeed-scPhi),cos(PhiSeed-scPhi));
    vals[13] = ESeed/SCRawEnergy;
    vals[14] = E3x3Seed/ESeed;
    vals[15] = E5x5Seed/ESeed;
    vals[16] = see;
    vals[17] = spp;
    vals[18] = sep;
    vals[19] = EMaxSeed/ESeed;
    vals[20] = E2ndSeed/ESeed;
    vals[21] = ETopSeed/ESeed;
    vals[22] = EBottomSeed/ESeed;
    vals[23] = ELeftSeed/ESeed;
    vals[24] = ERightSeed/ESeed;
    vals[25] = E2x5MaxSeed/ESeed;
    vals[26] = E2x5TopSeed/ESeed;
    vals[27] = E2x5BottomSeed/ESeed;
    vals[28] = E2x5LeftSeed/ESeed;
    vals[29] = E2x5RightSeed/ESeed;
    vals[30] = PreShowerOverRaw;
  }

  // Now evaluating the regression
  double regressionResult = 0;
  Int_t BinIndex = -1;

  if (fVersionType == kNoTrkVar) {
    if (fabs(scEta) <= 1.479) { 
      regressionResult = SCRawEnergy * forestUncertainty_eb->GetResponse(vals); 
      BinIndex = 0;
    }
    else {
      regressionResult = (SCRawEnergy*(1+PreShowerOverRaw)) * forestUncertainty_ee->GetResponse(vals);
      BinIndex = 1;
    }
  }

  //print debug
  if (printDebug) {    
    if (fabs(scEta) <= 1.479) {
      std::cout << "Barrel :";
      for (uint v=0; v < 38; ++v) std::cout << vals[v] << ", ";
      std::cout << "\n";
    }
    else {
      std::cout << "Endcap :";
      for (uint v=0; v < 32; ++v) std::cout << vals[v] << ", ";
      std::cout << "\n";
    }
    std::cout << "BinIndex : " << BinIndex << "\n";
    std::cout << "SCRawEnergy = " << SCRawEnergy << " : PreShowerOverRaw = " << PreShowerOverRaw << std::endl;
    std::cout << "regression energy uncertainty = " << regressionResult << std::endl;
  }
  

  // Cleaning up and returning
  delete[] vals;
  return regressionResult;
}

double ElectronEnergyRegression::regressionValueWithTrkVarV1( std::vector<double> &inputvars, 
                                                              bool printDebug) 
{
  // Checking if instance has been initialized
  if (fIsInitialized == kFALSE) {
    printf("ElectronEnergyRegression instance not initialized !!!");
    return 0;
  }

  // Checking if fVersionType is correct
  assert(fVersionType == kWithTrkVarV1);

  //Check that inputvars vector has the right number of inputs
  assert(inputvars.size() == 42);

  //assign variables from inputvars array to named variables
  double SCRawEnergy  = inputvars[0];
  double scEta  = inputvars[1];
  double scPhi  = inputvars[2];
  double R9  = inputvars[3];
  double etawidth  = inputvars[4];
  double phiwidth  = inputvars[5];
  double NClusters  = inputvars[6];
  double HoE  = inputvars[7];
  double rho  = inputvars[8];
  double vertices  = inputvars[9];
  double EtaSeed  = inputvars[10];
  double PhiSeed  = inputvars[11];
  double ESeed  = inputvars[12];
  double E3x3Seed  = inputvars[13];
  double E5x5Seed  = inputvars[14];
  double see  = inputvars[15];
  double spp  = inputvars[16];
  double sep  = inputvars[17];
  double EMaxSeed  = inputvars[18];
  double E2ndSeed  = inputvars[19];
  double ETopSeed  = inputvars[20];
  double EBottomSeed  = inputvars[21];
  double ELeftSeed  = inputvars[22];
  double ERightSeed  = inputvars[23];
  double E2x5MaxSeed  = inputvars[24];
  double E2x5TopSeed  = inputvars[25];
  double E2x5BottomSeed  = inputvars[26];
  double E2x5LeftSeed  = inputvars[27];
  double E2x5RightSeed  = inputvars[28];
  double IEtaSeed  = inputvars[29];
  double IPhiSeed  = inputvars[30];
  double EtaCrySeed  = inputvars[31];
  double PhiCrySeed  = inputvars[32];
  double PreShowerOverRaw  = inputvars[33];
  int    IsEcalDriven  = inputvars[34];
  double GsfTrackPIn  = inputvars[35];
  double fbrem  = inputvars[36];
  double Charge  = inputvars[37];
  double EoP  = inputvars[38];
  double TrackMomentumError  = inputvars[39];
  double EcalEnergyError  = inputvars[40];
  int    Classification  = inputvars[41]; 

  float *vals = (fabs(scEta) <= 1.479) ? new float[46] : new float[39];
  if (fabs(scEta) <= 1.479) {		// Barrel
    vals[0]  = SCRawEnergy;
    vals[1]  = scEta;
    vals[2]  = scPhi;
    vals[3]  = R9;
    vals[4]  = E5x5Seed/SCRawEnergy;
    vals[5]  = etawidth;
    vals[6]  = phiwidth;
    vals[7]  = NClusters;
    vals[8]  = HoE;
    vals[9]  = rho;
    vals[10] = vertices;
    vals[11] = EtaSeed - scEta;
    vals[12] = atan2(sin(PhiSeed-scPhi),cos(PhiSeed-scPhi));
    vals[13] = ESeed/SCRawEnergy;
    vals[14] = E3x3Seed/ESeed;
    vals[15] = E5x5Seed/ESeed;
    vals[16] = see;
    vals[17] = spp;
    vals[18] = sep;
    vals[19] = EMaxSeed/ESeed;
    vals[20] = E2ndSeed/ESeed;
    vals[21] = ETopSeed/ESeed;
    vals[22] = EBottomSeed/ESeed;
    vals[23] = ELeftSeed/ESeed;
    vals[24] = ERightSeed/ESeed;
    vals[25] = E2x5MaxSeed/ESeed;
    vals[26] = E2x5TopSeed/ESeed;
    vals[27] = E2x5BottomSeed/ESeed;
    vals[28] = E2x5LeftSeed/ESeed;
    vals[29] = E2x5RightSeed/ESeed;
    vals[30] = IsEcalDriven;
    vals[31] = GsfTrackPIn;
    vals[32] = fbrem;
    vals[33] = Charge;
    vals[34] = EoP;
    vals[35] = TrackMomentumError/GsfTrackPIn;
    vals[36] = EcalEnergyError/SCRawEnergy;
    vals[37] = Classification;
    vals[38] = IEtaSeed;
    vals[39] = IPhiSeed;
    vals[40] = ((int) IEtaSeed)%5;
    vals[41] = ((int) IPhiSeed)%2;
    vals[42] = (abs(IEtaSeed)<=25)*(((int)IEtaSeed)%25) + (abs(IEtaSeed)>25)*(((int) (IEtaSeed-25*abs(IEtaSeed)/IEtaSeed))%20);
    vals[43] = ((int) IPhiSeed)%20;
    vals[44] = EtaCrySeed;
    vals[45] = PhiCrySeed;
  }

  else {	// Endcap
    vals[0]  = SCRawEnergy;
    vals[1]  = scEta;
    vals[2]  = scPhi;
    vals[3]  = R9;
    vals[4]  = E5x5Seed/SCRawEnergy;
    vals[5]  = etawidth;
    vals[6]  = phiwidth;
    vals[7]  = NClusters;
    vals[8]  = HoE;
    vals[9]  = rho;
    vals[10] = vertices;
    vals[11] = EtaSeed - scEta;
    vals[12] = atan2(sin(PhiSeed-scPhi),cos(PhiSeed-scPhi));
    vals[13] = ESeed/SCRawEnergy;
    vals[14] = E3x3Seed/ESeed;
    vals[15] = E5x5Seed/ESeed;
    vals[16] = see;
    vals[17] = spp;
    vals[18] = sep;
    vals[19] = EMaxSeed/ESeed;
    vals[20] = E2ndSeed/ESeed;
    vals[21] = ETopSeed/ESeed;
    vals[22] = EBottomSeed/ESeed;
    vals[23] = ELeftSeed/ESeed;
    vals[24] = ERightSeed/ESeed;
    vals[25] = E2x5MaxSeed/ESeed;
    vals[26] = E2x5TopSeed/ESeed;
    vals[27] = E2x5BottomSeed/ESeed;
    vals[28] = E2x5LeftSeed/ESeed;
    vals[29] = E2x5RightSeed/ESeed;
    vals[30] = IsEcalDriven;
    vals[31] = GsfTrackPIn;
    vals[32] = fbrem;
    vals[33] = Charge;
    vals[34] = EoP;
    vals[35] = TrackMomentumError/GsfTrackPIn;
    vals[36] = EcalEnergyError/SCRawEnergy;
    vals[37] = Classification;
    vals[38] = PreShowerOverRaw;
  }

  // Now evaluating the regression
  double regressionResult = 0;

  if (fVersionType == kWithTrkVarV1) {
    if (fabs(scEta) <= 1.479) regressionResult = SCRawEnergy * forestCorrection_eb->GetResponse(vals);
    else regressionResult = (SCRawEnergy*(1+PreShowerOverRaw)) * forestCorrection_ee->GetResponse(vals);
  }


  //print debug
  if (printDebug) {
    if (fabs(scEta) <= 1.479) {
      std::cout << "Barrel :";
      for (uint v=0; v < 46; ++v) std::cout << vals[v] << ", ";
      std::cout << "\n";
    }
    else {
      std::cout << "Endcap :";
      for (uint v=0; v < 39; ++v) std::cout << vals[v] << ", ";
      std::cout << "\n";
    }
    std::cout << "SCRawEnergy = " << SCRawEnergy << " : PreShowerOverRaw = " << PreShowerOverRaw << std::endl;
    std::cout << "regression energy = " << regressionResult << std::endl;
  }

  // Cleaning up and returning
  delete[] vals;
  return regressionResult;
}





double ElectronEnergyRegression::regressionUncertaintyWithTrkVarV1( std::vector<double> &inputvars, 
                                                                    bool printDebug) 
{
  // Checking if instance has been initialized
  if (fIsInitialized == kFALSE) {
    printf("ElectronEnergyRegression instance not initialized !!!");
    return 0;
  }

  // Checking if fVersionType is correct
  assert(fVersionType == kWithTrkVarV1);

  // Checking if fVersionType is correct
  assert(inputvars.size() == 42);

  double SCRawEnergy  = inputvars[0];
  double scEta  = inputvars[1];
  double scPhi  = inputvars[2];
  double R9  = inputvars[3];
  double etawidth  = inputvars[4];
  double phiwidth  = inputvars[5];
  double NClusters  = inputvars[6];
  double HoE  = inputvars[7];
  double rho  = inputvars[8];
  double vertices  = inputvars[9];
  double EtaSeed  = inputvars[10];
  double PhiSeed  = inputvars[11];
  double ESeed  = inputvars[12];
  double E3x3Seed  = inputvars[13];
  double E5x5Seed  = inputvars[14];
  double see  = inputvars[15];
  double spp  = inputvars[16];
  double sep  = inputvars[17];
  double EMaxSeed  = inputvars[18];
  double E2ndSeed  = inputvars[19];
  double ETopSeed  = inputvars[20];
  double EBottomSeed  = inputvars[21];
  double ELeftSeed  = inputvars[22];
  double ERightSeed  = inputvars[23];
  double E2x5MaxSeed  = inputvars[24];
  double E2x5TopSeed  = inputvars[25];
  double E2x5BottomSeed  = inputvars[26];
  double E2x5LeftSeed  = inputvars[27];
  double E2x5RightSeed  = inputvars[28];
  double IEtaSeed  = inputvars[29];
  double IPhiSeed  = inputvars[30];
  double EtaCrySeed  = inputvars[31];
  double PhiCrySeed  = inputvars[32];
  double PreShowerOverRaw  = inputvars[33];
  int    IsEcalDriven  = inputvars[34];
  double GsfTrackPIn  = inputvars[35];
  double fbrem  = inputvars[36];
  double Charge  = inputvars[37];
  double EoP  = inputvars[38];
  double TrackMomentumError  = inputvars[39];
  double EcalEnergyError  = inputvars[40];
  int    Classification  = inputvars[41]; 


  float *vals = (fabs(scEta) <= 1.479) ? new float[46] : new float[39];
  if (fabs(scEta) <= 1.479) {		// Barrel
    vals[0]  = SCRawEnergy;
    vals[1]  = scEta;
    vals[2]  = scPhi;
    vals[3]  = R9;
    vals[4]  = E5x5Seed/SCRawEnergy;
    vals[5]  = etawidth;
    vals[6]  = phiwidth;
    vals[7]  = NClusters;
    vals[8]  = HoE;
    vals[9]  = rho;
    vals[10] = vertices;
    vals[11] = EtaSeed - scEta;
    vals[12] = atan2(sin(PhiSeed-scPhi),cos(PhiSeed-scPhi));
    vals[13] = ESeed/SCRawEnergy;
    vals[14] = E3x3Seed/ESeed;
    vals[15] = E5x5Seed/ESeed;
    vals[16] = see;
    vals[17] = spp;
    vals[18] = sep;
    vals[19] = EMaxSeed/ESeed;
    vals[20] = E2ndSeed/ESeed;
    vals[21] = ETopSeed/ESeed;
    vals[22] = EBottomSeed/ESeed;
    vals[23] = ELeftSeed/ESeed;
    vals[24] = ERightSeed/ESeed;
    vals[25] = E2x5MaxSeed/ESeed;
    vals[26] = E2x5TopSeed/ESeed;
    vals[27] = E2x5BottomSeed/ESeed;
    vals[28] = E2x5LeftSeed/ESeed;
    vals[29] = E2x5RightSeed/ESeed;
    vals[30] = IsEcalDriven;
    vals[31] = GsfTrackPIn;
    vals[32] = fbrem;
    vals[33] = Charge;
    vals[34] = EoP;
    vals[35] = TrackMomentumError/GsfTrackPIn;
    vals[36] = EcalEnergyError/SCRawEnergy;
    vals[37] = Classification;
    vals[38] = IEtaSeed;
    vals[39] = IPhiSeed;
    vals[40] = ((int) IEtaSeed)%5;
    vals[41] = ((int) IPhiSeed)%2;
    vals[42] = (abs(IEtaSeed)<=25)*(((int)IEtaSeed)%25) + (abs(IEtaSeed)>25)*(((int) (IEtaSeed-25*abs(IEtaSeed)/IEtaSeed))%20);
    vals[43] = ((int) IPhiSeed)%20;
    vals[44] = EtaCrySeed;
    vals[45] = PhiCrySeed;
  }

  else {	// Endcap
    vals[0]  = SCRawEnergy;
    vals[1]  = scEta;
    vals[2]  = scPhi;
    vals[3]  = R9;
    vals[4]  = E5x5Seed/SCRawEnergy;
    vals[5]  = etawidth;
    vals[6]  = phiwidth;
    vals[7]  = NClusters;
    vals[8]  = HoE;
    vals[9]  = rho;
    vals[10] = vertices;
    vals[11] = EtaSeed - scEta;
    vals[12] = atan2(sin(PhiSeed-scPhi),cos(PhiSeed-scPhi));
    vals[13] = ESeed/SCRawEnergy;
    vals[14] = E3x3Seed/ESeed;
    vals[15] = E5x5Seed/ESeed;
    vals[16] = see;
    vals[17] = spp;
    vals[18] = sep;
    vals[19] = EMaxSeed/ESeed;
    vals[20] = E2ndSeed/ESeed;
    vals[21] = ETopSeed/ESeed;
    vals[22] = EBottomSeed/ESeed;
    vals[23] = ELeftSeed/ESeed;
    vals[24] = ERightSeed/ESeed;
    vals[25] = E2x5MaxSeed/ESeed;
    vals[26] = E2x5TopSeed/ESeed;
    vals[27] = E2x5BottomSeed/ESeed;
    vals[28] = E2x5LeftSeed/ESeed;
    vals[29] = E2x5RightSeed/ESeed;
    vals[30] = IsEcalDriven;
    vals[31] = GsfTrackPIn;
    vals[32] = fbrem;
    vals[33] = Charge;
    vals[34] = EoP;
    vals[35] = TrackMomentumError/GsfTrackPIn;
    vals[36] = EcalEnergyError/SCRawEnergy;
    vals[37] = Classification;
    vals[38] = PreShowerOverRaw;
  }

  // Now evaluating the regression
  double regressionResult = 0;

  if (fVersionType == kWithTrkVarV1) {
    if (fabs(scEta) <= 1.479) regressionResult = SCRawEnergy * forestUncertainty_eb->GetResponse(vals);
    else regressionResult = (SCRawEnergy*(1+PreShowerOverRaw)) * forestUncertainty_ee->GetResponse(vals);
  }

  //print debug
  if (printDebug) {
    if (fabs(scEta) <= 1.479) {
      std::cout << "Barrel :";
      for (uint v=0; v < 46; ++v) std::cout << vals[v] << ", ";
      std::cout << "\n";
    }
    else {
      std::cout << "Endcap :";
      for (uint v=0; v < 39; ++v) std::cout << vals[v] << ", ";
      std::cout << "\n";
    }
    std::cout << " SCRawEnergy = " << SCRawEnergy << " : PreShowerOverRaw = " << PreShowerOverRaw << std::endl;
    std::cout << "regression energy uncertainty = " << regressionResult << std::endl;
  }


  // Cleaning up and returning
  delete[] vals;
  return regressionResult;
}






double ElectronEnergyRegression::regressionValueWithTrkVarV2( std::vector<double> &inputvars, 
                                                              bool printDebug) 
{
  // Checking if instance has been initialized
  if (fIsInitialized == kFALSE) {
    printf("ElectronEnergyRegression instance not initialized !!!");
    return 0;
  }

  // Checking if fVersionType is correct
  assert(fVersionType == kWithTrkVarV2);

  // Checking if fVersionType is correct
  assert(inputvars.size() == 49);

  double SCRawEnergy  = inputvars[0];
  double scEta  = inputvars[1];
  double scPhi  = inputvars[2];
  double R9  = inputvars[3];
  double etawidth  = inputvars[4];
  double phiwidth  = inputvars[5];
  double NClusters  = inputvars[6];
  double HoE  = inputvars[7];
  double rho  = inputvars[8];
  double vertices  = inputvars[9];
  double EtaSeed  = inputvars[10];
  double PhiSeed  = inputvars[11];
  double ESeed  = inputvars[12];
  double E3x3Seed  = inputvars[13];
  double E5x5Seed  = inputvars[14];
  double see  = inputvars[15];
  double spp  = inputvars[16];
  double sep  = inputvars[17];
  double EMaxSeed  = inputvars[18];
  double E2ndSeed  = inputvars[19];
  double ETopSeed  = inputvars[20];
  double EBottomSeed  = inputvars[21];
  double ELeftSeed  = inputvars[22];
  double ERightSeed  = inputvars[23];
  double E2x5MaxSeed  = inputvars[24];
  double E2x5TopSeed  = inputvars[25];
  double E2x5BottomSeed  = inputvars[26];
  double E2x5LeftSeed  = inputvars[27];
  double E2x5RightSeed  = inputvars[28];
  double IEtaSeed  = inputvars[29];
  double IPhiSeed  = inputvars[30];
  double EtaCrySeed  = inputvars[31];
  double PhiCrySeed  = inputvars[32];
  double PreShowerOverRaw  = inputvars[33];
  int    IsEcalDriven  = inputvars[34];
  double GsfTrackPIn  = inputvars[35];
  double fbrem  = inputvars[36];
  double Charge  = inputvars[37];
  double EoP  = inputvars[38];
  double TrackMomentumError  = inputvars[39];
  double EcalEnergyError  = inputvars[40];
  int    Classification  = inputvars[41]; 
  double detaIn  = inputvars[42];
  double dphiIn  = inputvars[43];
  double detaCalo  = inputvars[44];
  double dphiCalo  = inputvars[45];
  double GsfTrackChiSqr  = inputvars[46];
  double KFTrackNLayers  = inputvars[47];
  double ElectronEnergyOverPout  = inputvars[48];

  float *vals = (fabs(scEta) <= 1.479) ? new float[53] : new float[46];
  if (fabs(scEta) <= 1.479) {		// Barrel
    vals[0]  = SCRawEnergy;
    vals[1]  = scEta;
    vals[2]  = scPhi;
    vals[3]  = R9;
    vals[4]  = E5x5Seed/SCRawEnergy;
    vals[5]  = etawidth;
    vals[6]  = phiwidth;
    vals[7]  = NClusters;
    vals[8]  = HoE;
    vals[9]  = rho;
    vals[10] = vertices;
    vals[11] = EtaSeed - scEta;
    vals[12] = atan2(sin(PhiSeed-scPhi),cos(PhiSeed-scPhi));
    vals[13] = ESeed/SCRawEnergy;
    vals[14] = E3x3Seed/ESeed;
    vals[15] = E5x5Seed/ESeed;
    vals[16] = see;
    vals[17] = spp;
    vals[18] = sep;
    vals[19] = EMaxSeed/ESeed;
    vals[20] = E2ndSeed/ESeed;
    vals[21] = ETopSeed/ESeed;
    vals[22] = EBottomSeed/ESeed;
    vals[23] = ELeftSeed/ESeed;
    vals[24] = ERightSeed/ESeed;
    vals[25] = E2x5MaxSeed/ESeed;
    vals[26] = E2x5TopSeed/ESeed;
    vals[27] = E2x5BottomSeed/ESeed;
    vals[28] = E2x5LeftSeed/ESeed;
    vals[29] = E2x5RightSeed/ESeed;
    vals[30] = IsEcalDriven;
    vals[31] = GsfTrackPIn;
    vals[32] = fbrem;
    vals[33] = Charge;
    vals[34] = EoP;
    vals[35] = TrackMomentumError/GsfTrackPIn;
    vals[36] = EcalEnergyError/SCRawEnergy;
    vals[37] = Classification;
    vals[38] = detaIn;
    vals[39] = dphiIn;
    vals[40] = detaCalo;
    vals[41] = dphiCalo;
    vals[42] = GsfTrackChiSqr;
    vals[43] = KFTrackNLayers;
    vals[44] = ElectronEnergyOverPout;
    vals[45] = IEtaSeed;
    vals[46] = IPhiSeed;
    vals[47] = ((int) IEtaSeed)%5;
    vals[48] = ((int) IPhiSeed)%2;
    vals[49] = (abs(IEtaSeed)<=25)*(((int)IEtaSeed)%25) + (abs(IEtaSeed)>25)*(((int) (IEtaSeed-25*abs(IEtaSeed)/IEtaSeed))%20);
    vals[50] = ((int) IPhiSeed)%20;
    vals[51] = EtaCrySeed;
    vals[52] = PhiCrySeed;
  }

  else {	// Endcap
    vals[0]  = SCRawEnergy;
    vals[1]  = scEta;
    vals[2]  = scPhi;
    vals[3]  = R9;
    vals[4]  = E5x5Seed/SCRawEnergy;
    vals[5]  = etawidth;
    vals[6]  = phiwidth;
    vals[7]  = NClusters;
    vals[8]  = HoE;
    vals[9]  = rho;
    vals[10] = vertices;
    vals[11] = EtaSeed - scEta;
    vals[12] = atan2(sin(PhiSeed-scPhi),cos(PhiSeed-scPhi));
    vals[13] = ESeed/SCRawEnergy;
    vals[14] = E3x3Seed/ESeed;
    vals[15] = E5x5Seed/ESeed;
    vals[16] = see;
    vals[17] = spp;
    vals[18] = sep;
    vals[19] = EMaxSeed/ESeed;
    vals[20] = E2ndSeed/ESeed;
    vals[21] = ETopSeed/ESeed;
    vals[22] = EBottomSeed/ESeed;
    vals[23] = ELeftSeed/ESeed;
    vals[24] = ERightSeed/ESeed;
    vals[25] = E2x5MaxSeed/ESeed;
    vals[26] = E2x5TopSeed/ESeed;
    vals[27] = E2x5BottomSeed/ESeed;
    vals[28] = E2x5LeftSeed/ESeed;
    vals[29] = E2x5RightSeed/ESeed;
    vals[30] = IsEcalDriven;
    vals[31] = GsfTrackPIn;
    vals[32] = fbrem;
    vals[33] = Charge;
    vals[34] = EoP;
    vals[35] = TrackMomentumError/GsfTrackPIn;
    vals[36] = EcalEnergyError/SCRawEnergy;
    vals[37] = Classification;
    vals[38] = detaIn;
    vals[39] = dphiIn;
    vals[40] = detaCalo;
    vals[41] = dphiCalo;
    vals[42] = GsfTrackChiSqr;
    vals[43] = KFTrackNLayers;
    vals[44] = ElectronEnergyOverPout;
    vals[45] = PreShowerOverRaw;
  }

  // Now evaluating the regression
  double regressionResult = 0;

  if (fVersionType == kWithTrkVarV2) {
    if (fabs(scEta) <= 1.479) regressionResult = SCRawEnergy * forestCorrection_eb->GetResponse(vals);
    else regressionResult = (SCRawEnergy*(1+PreShowerOverRaw)) * forestCorrection_ee->GetResponse(vals);
  }


  //print debug
  if (printDebug) {
    if (fabs(scEta) <= 1.479) {
      std::cout << "Barrel :";
      for (uint v=0; v < 53; ++v) std::cout << vals[v] << ", ";
      std::cout << "\n";
    }
    else {
      std::cout << "Endcap :";
      for (uint v=0; v < 46; ++v) std::cout << vals[v] << ", ";
      std::cout << "\n";
    }
    std::cout << "SCRawEnergy = " << SCRawEnergy << " : PreShowerOverRaw = " << PreShowerOverRaw << std::endl;
    std::cout << "regression energy = " << regressionResult << std::endl;
  }

  // Cleaning up and returning
  delete[] vals;
  return regressionResult;
}






double ElectronEnergyRegression::regressionUncertaintyWithTrkVarV2( std::vector<double> &inputvars, 
                                                                    bool printDebug) 
{
  // Checking if instance has been initialized
  if (fIsInitialized == kFALSE) {
    printf("ElectronEnergyRegression instance not initialized !!!");
    return 0;
  }

  // Checking if fVersionType is correct
  assert(fVersionType == kWithTrkVarV2);

  // Checking if fVersionType is correct
  assert(inputvars.size() == 49);

  double SCRawEnergy  = inputvars[0];
  double scEta  = inputvars[1];
  double scPhi  = inputvars[2];
  double R9  = inputvars[3];
  double etawidth  = inputvars[4];
  double phiwidth  = inputvars[5];
  double NClusters  = inputvars[6];
  double HoE  = inputvars[7];
  double rho  = inputvars[8];
  double vertices  = inputvars[9];
  double EtaSeed  = inputvars[10];
  double PhiSeed  = inputvars[11];
  double ESeed  = inputvars[12];
  double E3x3Seed  = inputvars[13];
  double E5x5Seed  = inputvars[14];
  double see  = inputvars[15];
  double spp  = inputvars[16];
  double sep  = inputvars[17];
  double EMaxSeed  = inputvars[18];
  double E2ndSeed  = inputvars[19];
  double ETopSeed  = inputvars[20];
  double EBottomSeed  = inputvars[21];
  double ELeftSeed  = inputvars[22];
  double ERightSeed  = inputvars[23];
  double E2x5MaxSeed  = inputvars[24];
  double E2x5TopSeed  = inputvars[25];
  double E2x5BottomSeed  = inputvars[26];
  double E2x5LeftSeed  = inputvars[27];
  double E2x5RightSeed  = inputvars[28];
  double IEtaSeed  = inputvars[29];
  double IPhiSeed  = inputvars[30];
  double EtaCrySeed  = inputvars[31];
  double PhiCrySeed  = inputvars[32];
  double PreShowerOverRaw  = inputvars[33];
  int    IsEcalDriven  = inputvars[34];
  double GsfTrackPIn  = inputvars[35];
  double fbrem  = inputvars[36];
  double Charge  = inputvars[37];
  double EoP  = inputvars[38];
  double TrackMomentumError  = inputvars[39];
  double EcalEnergyError  = inputvars[40];
  int    Classification  = inputvars[41]; 
  double detaIn  = inputvars[42];
  double dphiIn  = inputvars[43];
  double detaCalo  = inputvars[44];
  double dphiCalo  = inputvars[45];
  double GsfTrackChiSqr  = inputvars[46];
  double KFTrackNLayers  = inputvars[47];
  double ElectronEnergyOverPout  = inputvars[48];

  float *vals = (fabs(scEta) <= 1.479) ? new float[53] : new float[46];
  if (fabs(scEta) <= 1.479) {		// Barrel
    vals[0]  = SCRawEnergy;
    vals[1]  = scEta;
    vals[2]  = scPhi;
    vals[3]  = R9;
    vals[4]  = E5x5Seed/SCRawEnergy;
    vals[5]  = etawidth;
    vals[6]  = phiwidth;
    vals[7]  = NClusters;
    vals[8]  = HoE;
    vals[9]  = rho;
    vals[10] = vertices;
    vals[11] = EtaSeed - scEta;
    vals[12] = atan2(sin(PhiSeed-scPhi),cos(PhiSeed-scPhi));
    vals[13] = ESeed/SCRawEnergy;
    vals[14] = E3x3Seed/ESeed;
    vals[15] = E5x5Seed/ESeed;
    vals[16] = see;
    vals[17] = spp;
    vals[18] = sep;
    vals[19] = EMaxSeed/ESeed;
    vals[20] = E2ndSeed/ESeed;
    vals[21] = ETopSeed/ESeed;
    vals[22] = EBottomSeed/ESeed;
    vals[23] = ELeftSeed/ESeed;
    vals[24] = ERightSeed/ESeed;
    vals[25] = E2x5MaxSeed/ESeed;
    vals[26] = E2x5TopSeed/ESeed;
    vals[27] = E2x5BottomSeed/ESeed;
    vals[28] = E2x5LeftSeed/ESeed;
    vals[29] = E2x5RightSeed/ESeed;
    vals[30] = IsEcalDriven;
    vals[31] = GsfTrackPIn;
    vals[32] = fbrem;
    vals[33] = Charge;
    vals[34] = EoP;
    vals[35] = TrackMomentumError/GsfTrackPIn;
    vals[36] = EcalEnergyError/SCRawEnergy;
    vals[37] = Classification;
    vals[38] = detaIn;
    vals[39] = dphiIn;
    vals[40] = detaCalo;
    vals[41] = dphiCalo;
    vals[42] = GsfTrackChiSqr;
    vals[43] = KFTrackNLayers;
    vals[44] = ElectronEnergyOverPout;
    vals[45] = IEtaSeed;
    vals[46] = IPhiSeed;
    vals[47] = ((int) IEtaSeed)%5;
    vals[48] = ((int) IPhiSeed)%2;
    vals[49] = (abs(IEtaSeed)<=25)*(((int)IEtaSeed)%25) + (abs(IEtaSeed)>25)*(((int) (IEtaSeed-25*abs(IEtaSeed)/IEtaSeed))%20);
    vals[50] = ((int) IPhiSeed)%20;
    vals[51] = EtaCrySeed;
    vals[52] = PhiCrySeed;
  }

  else {	// Endcap
    vals[0]  = SCRawEnergy;
    vals[1]  = scEta;
    vals[2]  = scPhi;
    vals[3]  = R9;
    vals[4]  = E5x5Seed/SCRawEnergy;
    vals[5]  = etawidth;
    vals[6]  = phiwidth;
    vals[7]  = NClusters;
    vals[8]  = HoE;
    vals[9]  = rho;
    vals[10] = vertices;
    vals[11] = EtaSeed - scEta;
    vals[12] = atan2(sin(PhiSeed-scPhi),cos(PhiSeed-scPhi));
    vals[13] = ESeed/SCRawEnergy;
    vals[14] = E3x3Seed/ESeed;
    vals[15] = E5x5Seed/ESeed;
    vals[16] = see;
    vals[17] = spp;
    vals[18] = sep;
    vals[19] = EMaxSeed/ESeed;
    vals[20] = E2ndSeed/ESeed;
    vals[21] = ETopSeed/ESeed;
    vals[22] = EBottomSeed/ESeed;
    vals[23] = ELeftSeed/ESeed;
    vals[24] = ERightSeed/ESeed;
    vals[25] = E2x5MaxSeed/ESeed;
    vals[26] = E2x5TopSeed/ESeed;
    vals[27] = E2x5BottomSeed/ESeed;
    vals[28] = E2x5LeftSeed/ESeed;
    vals[29] = E2x5RightSeed/ESeed;
    vals[30] = IsEcalDriven;
    vals[31] = GsfTrackPIn;
    vals[32] = fbrem;
    vals[33] = Charge;
    vals[34] = EoP;
    vals[35] = TrackMomentumError/GsfTrackPIn;
    vals[36] = EcalEnergyError/SCRawEnergy;
    vals[37] = Classification;
    vals[38] = detaIn;
    vals[39] = dphiIn;
    vals[40] = detaCalo;
    vals[41] = dphiCalo;
    vals[42] = GsfTrackChiSqr;
    vals[43] = KFTrackNLayers;
    vals[44] = ElectronEnergyOverPout;
    vals[45] = PreShowerOverRaw;
  }

  // Now evaluating the regression
  double regressionResult = 0;

  if (fVersionType == kWithTrkVarV2) {
    if (fabs(scEta) <= 1.479) regressionResult = SCRawEnergy * forestUncertainty_eb->GetResponse(vals);
    else regressionResult = (SCRawEnergy*(1+PreShowerOverRaw)) * forestUncertainty_ee->GetResponse(vals);
  }

  //print debug
  if (printDebug) {
    if (fabs(scEta) <= 1.479) {
      std::cout << "Barrel :";
      for (uint v=0; v < 53; ++v) std::cout << vals[v] << ", ";
      std::cout << "\n";
    }
    else {
      std::cout << "Endcap :";
      for (uint v=0; v < 46; ++v) std::cout << vals[v] << ", ";
      std::cout << "\n";
    }
    std::cout << "SCRawEnergy = " << SCRawEnergy << " : PreShowerOverRaw = " << PreShowerOverRaw << std::endl;
    std::cout << "regression energy uncertainty = " << regressionResult << std::endl;
  }


  // Cleaning up and returning
  delete[] vals;
  return regressionResult;
}

