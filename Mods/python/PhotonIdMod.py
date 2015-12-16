from MitAna.TreeMod.bambu import mithep

photonIdMod = mithep.PhotonIdMod(
    OutputName = 'LoosePhotons',
    IdType = mithep.PhotonTools.kSpring15Loose,
    IsoType = mithep.PhotonTools.kSpring15LooseIso,
    PtMin = 15.,
    EtaMax = 2.5
)
