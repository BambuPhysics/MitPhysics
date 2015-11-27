from MitAna.TreeMod.bambu import mithep

puppiPFJetMod = mithep.PuppiPFJetMod(
    InputName = "PuppiParticles",
    OutputName = 'PuppiPFJets',
    ProcessNJets = 10,
    R0 = 0.4,
    JetAlgorithm = mithep.PuppiPFJetMod.kAntiKT,
    DoMatching = True,
    MatchingJetsName = 'AKt4PFJetsCHS'
)
