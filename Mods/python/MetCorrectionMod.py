from MitAna.TreeMod.bambu import mithep, analysis
import os

metCorrectionMod = mithep.MetCorrectionMod(
    InputName = 'PFMet',
    OutputName = 'PFType1CorrectedMet',
    JetsName = 'AKt4PFJets',
    RhoAlgo = mithep.PileupEnergyDensity.kFixedGridFastjetAll,
    MaxEMFraction = 0.9,
    SkipMuons = True
)
metCorrectionMod.ApplyType0(False)
metCorrectionMod.ApplyType1(True)
metCorrectionMod.ApplyShift(False)
metCorrectionMod.IsData(analysis.isRealData)
if analysis.isRealData:
    metCorrectionMod.AddJetCorrectionFromFile(os.environ['MIT_DATA'] + "/74X_dataRun2_Prompt_v1_L1FastJet_AK4PFchs.txt")
    metCorrectionMod.AddJetCorrectionFromFile(os.environ['MIT_DATA'] + "/74X_dataRun2_Prompt_v1_L2Relative_AK4PFchs.txt")
    metCorrectionMod.AddJetCorrectionFromFile(os.environ['MIT_DATA'] + "/74X_dataRun2_Prompt_v1_L3Absolute_AK4PFchs.txt")
    metCorrectionMod.AddJetCorrectionFromFile(os.environ['MIT_DATA'] + "/74X_dataRun2_Prompt_v1_L2L3Residual_AK4PFchs.txt")
else:
    metCorrectionMod.AddJetCorrectionFromFile(os.environ['MIT_DATA'] + "/MCRUN2_74_V9_L1FastJet_AK4PF.txt")
    metCorrectionMod.AddJetCorrectionFromFile(os.environ['MIT_DATA'] + "/MCRUN2_74_V9_L2Relative_AK4PF.txt")
    metCorrectionMod.AddJetCorrectionFromFile(os.environ['MIT_DATA'] + "/MCRUN2_74_V9_L3Absolute_AK4PF.txt")
