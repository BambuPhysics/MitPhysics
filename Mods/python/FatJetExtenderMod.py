from MitAna.TreeMod.bambu import mithep

fatJetExtenderMod = mithep.FatJetExtenderMod(
    InputName = "FatJets",
    OutputName = 'XlFatJets',
    ProcessNJets = 4,
    ConeSize = 0.8
)
