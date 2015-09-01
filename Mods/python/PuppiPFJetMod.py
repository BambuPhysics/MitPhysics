from MitAna.TreeMod.bambu import mithep
import os

puppiPFJetMod = mithep.PuppiPFJetMod(
    InputName = "PuppiParticles",
    OutputName = 'PuppiPFJets',
    ProcessNJets = 10,
    R0 = 0.4
)
