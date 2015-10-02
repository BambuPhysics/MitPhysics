from MitAna.TreeMod.bambu import mithep

electronIdMod = mithep.ElectronIdMod(
    OutputName = 'VetoElectrons',
    IdType = mithep.ElectronTools.kNoId,
    IsoType = mithep.ElectronTools.kNoIso,
    ApplyEcalFiducial = True,
    ApplyD0Cut = False,
    ApplyDZCut = False,
    WhichVertex = 0,
    PtMin = 10.,
    EtaMax = 2.5,
    ApplyConvFilterType1 = False,
    ApplyConvFilterType2 = False,
    ApplyNExpectedHitsInnerCut = False,
    InvertNExpectedHitsInnerCut = False,
    ApplySpikeRemoval = False,
    ChargeFilter = False,
    ApplyTriggerMatching = False,
    ApplyEcalSeeded = False,
    RhoAlgo = mithep.PileupEnergyDensity.kFixedGridFastjetAll,
    EtMin = 0.,
    IdLikelihoodCut = -999.,
    ConversionsName = 'Conversions'
)
