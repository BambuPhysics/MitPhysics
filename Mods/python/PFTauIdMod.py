from MitAna.TreeMod.bambu import mithep

pfTauIdMod = mithep.PFTauIdMod(
    OutputName = 'LooseTaus',
    PtMin = 18.,
    EtaMax = 2.3
)
pfTauIdMod.AddDiscriminator(mithep.PFTau.kDiscriminationByDecayModeFindingNewDMs)
