from MitAna.TreeMod.bambu import mithep

photonIdMod = mithep.PhotonIdMod(
    OutputName = 'LoosePhotons',
    IdType = mithep.PhotonTools.kPhys14Loose,
    IsoType = mithep.PhotonTools.kPhys14Loose,
    PtMin = 15.,
    EtaMax = 2.5
)
