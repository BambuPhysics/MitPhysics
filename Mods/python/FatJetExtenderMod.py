from MitAna.TreeMod.bambu import mithep
import os

fatJetExtenderMod = mithep.fatJetExtenderMod(
    InputName = "FatJets",
    OutputName = 'XlFatJets',
    ProcessNJets = 4,
    ConeSize = 0.8,
    PublishOutput = True
)
