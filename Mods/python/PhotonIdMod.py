from MitAna.TreeMod.bambu import mithep

photonIdMod = mithep.PhotonIdMod(
    OutputName = 'LoosePhotons',
    IdType = mithep.PhotonTools.kSummer15Loose,
    IsoType = mithep.PhotonTools.kSummer15LooseIso,
    PtMin = 15.,
    EtaMax = 2.5
)
