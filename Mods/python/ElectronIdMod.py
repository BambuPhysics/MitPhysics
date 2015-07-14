from MitAna.TreeMod.bambu import mithep

electronIdMod = mithep.ElectronIdMod(
    OutputName = 'VetoElectrons',
    IdType = mithep.ElectronTools.kPhys14Veto,
    IsoType = mithep.ElectronTools.kPhys14VetoIso,
    ApplyEcalFiducial = True,
    WhichVertex = 0,
    PtMin = 10.,
    EtaMax = 2.5,
    ConversionsName = 'Conversions'
)
