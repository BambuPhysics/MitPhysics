from MitAna.TreeMod.bambu import mithep
import os

jetIdMod = mithep.JetIdMod(
    InputName = mithep.ModNames.gkCorrectedJetsName,
    OutputName = 'GoodJets',
    ApplyPFLooseId = True,
    MVATrainingSet = mithep.JetIDMVA.k53BDTCHSFullPlusRMS,
    MVACutWP = mithep.JetIDMVA.kLoose,
    MVAWeightsFile = os.environ['MIT_DATA'] + '/TMVAClassification_5x_BDT_chsFullPlusRMS.weights.xml',
    MVACutsFile = os.environ['MIT_DATA'] + '/jetIDCuts_121221.dat',
    UseClassicBetaForMVA = False,
    PtMin = 30.,
    EtaMax = 2.5
)

# UseClassicBetaForMVA: set to true for MiniAOD input. Makes the MVA value agree
# with what CMSSW would give if Jet Id was calculated on MiniAOD. This is not the
# same as the pre-calculated MVA value stored in MiniAOD and accessible with
# jet.userFloat().
