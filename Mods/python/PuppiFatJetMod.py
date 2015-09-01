from MitAna.TreeMod.bambu import mithep
import os

puppiFatJetMod = mithep.PuppiFatJetMod(
    InputName = "PuppiParticles",
    OutputName = 'PuppiFatJets',
    ProcessNJets = 4,
    R0 = 0.8
)
