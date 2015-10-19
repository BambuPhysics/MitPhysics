from MitAna.TreeMod.bambu import mithep
import os

jetIdMod = mithep.JetIdMod(
    InputName = mithep.ModNames.gkCorrectedJetsName,
    OutputName = 'GoodJets',
    PFId = mithep.JetTools.kPFLoose,
    MVATrainingSet = mithep.JetIDMVA.k74CHS, # set to mithep.JetIDMVA.nMVATrainingTypes to turn off MVA
    MVACutWP = mithep.JetIDMVA.kLoose,
    MVACutsFile = os.environ['MIT_DATA'] + '/JetId/jetIDCuts_150807.dat',
    UseClassicBetaForMVA = False,
    PtMin = 30.,
    EtaMax = 2.5,
    MaxChargedEMFraction = 0.99,
    MaxNeutralEMFraction = 0.99,
    MaxNeutralHadronFraction = 0.99,
    MaxMuonFraction = 0.80,
    MinNPFCandidates = 2,
    MinNChargedPFCandidates = 1,
)
jetIdMod.SetMVAWeightsFile(os.environ['MIT_DATA'] + '/JetId/TMVAClassificationCategory_BDTG.weights_jteta_0_2.xml', 0)
jetIdMod.SetMVAWeightsFile(os.environ['MIT_DATA'] + '/JetId/TMVAClassificationCategory_BDTG.weights_jteta_2_2p5.xml', 1)
jetIdMod.SetMVAWeightsFile(os.environ['MIT_DATA'] + '/JetId/TMVAClassificationCategory_BDTG.weights_jteta_2p5_3.xml', 2)
jetIdMod.SetMVAWeightsFile(os.environ['MIT_DATA'] + '/JetId/TMVAClassificationCategory_BDTG.weights_jteta_3_5.xml', 3)

# UseClassicBetaForMVA: set to true for MiniAOD input. Makes the MVA value agree
# with what CMSSW would give if Jet Id was calculated on MiniAOD. This is not the
# same as the pre-calculated MVA value stored in MiniAOD and accessible with
# jet.userFloat().
