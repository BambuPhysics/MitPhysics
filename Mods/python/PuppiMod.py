import os
from MitAna.TreeMod.bambu import mithep

puppiMod = mithep.PuppiMod(
    EtaConfigName = os.getenv("CMSSW_BASE")+'/src/MitPhysics/data/PuppiEta_150701.cfg',
    OutputName = 'PuppiParticles',
    RMin = 0.02,
    R0 = 0.3,
    Alpha = 2.0,
    Beta = 2.0,
    D0Cut = 0.1,
    DZCut = 0.3,
    MinWeightCut = 0.01,
    RMSScaleFactor = 1.0,
    TrackUncertainty = 1.0,
    KeepPileup = False,
    Invert = False,
    ApplyCHS = True,
    ApplyLowPUCorr = True,
    EtaForAlgo = 2.5,
    Dump = False
)
