//--------------------------------------------------------------------------------------------------
//
// ElectronEnergyRegression
//
// Helper Class for applying electron energy regression calculation
//
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_UTILS_ELECTRONENERGYREGRESSIONEVALUATE_H
#define MITPHYSICS_UTILS_ELECTRONENERGYREGRESSIONEVALUATE_H

// For applying regression
#include "TFile.h"

class GBRForest;

namespace mithep {
  class ElectronEnergyRegression{
    public:
      ElectronEnergyRegression();
      ~ElectronEnergyRegression();

      enum ElectronEnergyRegressionType {
        kNoTrkVar,
        kWithTrkVarV1,
        kWithTrkVarV2
      };

      void initialize(std::string weightsFile,
                      ElectronEnergyRegression::ElectronEnergyRegressionType type);

      bool isInitialized() const {return fIsInitialized;}
                
      // Evaluates regression without tracker variables
      double regressionValueNoTrkVar(
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
        bool printDebug = false);

      // Evaluates regression without tracker variables
      double regressionUncertaintyNoTrkVar(
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
        bool printDebug = false);

      // Evaluates regression using tracker variables
      double regressionValueWithTrkVarV1( std::vector<double> &inputvars,    
                                          bool printDebug = false );


      // Evaluates regression using tracker variables
      double regressionUncertaintyWithTrkVarV1(	std::vector<double> &inputvars,			
                                                bool printDebug = false );

      // Evaluates regression using tracker variables
      double regressionValueWithTrkVarV2( std::vector<double> &inputvars,    
                                          bool printDebug = false );


      // Evaluates regression using tracker variables
      double regressionUncertaintyWithTrkVarV2(	std::vector<double> &inputvars,			
                                                bool printDebug = false );

    private:
      bool fIsInitialized;
      ElectronEnergyRegression::ElectronEnergyRegressionType fVersionType;
      GBRForest *forestCorrection_eb;           //! Pointer to the GBRForest for barrel
      GBRForest *forestCorrection_ee;		//! Pointer to the GBRForest for endcap

      GBRForest *forestUncertainty_eb;	//!
      GBRForest *forestUncertainty_ee;	//!	
  };
}
#endif
