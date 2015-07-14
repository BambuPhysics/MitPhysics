from MitAna.TreeMod.bambu import mithep
import os

jetCorrectionMod = mithep.JetCorrectionMod(
    InputName = 'AKt4PFJetsCHS',
    CorrectedJetsName = 'CorrectedJets',
    RhoAlgo = mithep.PileupEnergyDensity.kFixedGridFastjetAll
)
jetCorrectionMod.AddCorrectionFromFile(os.environ['MIT_DATA'] + "/MCRUN2_74_V9_L1FastJet_AK4PFchs.txt")
jetCorrectionMod.AddCorrectionFromFile(os.environ['MIT_DATA'] + "/MCRUN2_74_V9_L2Relative_AK4PFchs.txt")
jetCorrectionMod.AddCorrectionFromFile(os.environ['MIT_DATA'] + "/MCRUN2_74_V9_L3Absolute_AK4PFchs.txt")
# currently correction file must be added in L1-L2-L3 order and cannot be removed
# -> will have to copy-paste this code to use any other jet collection / correction..
