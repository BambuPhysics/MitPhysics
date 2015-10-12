from MitAna.TreeMod.bambu import mithep, analysis
import os

metCorrectionMod = mithep.MetCorrectionMod(
    InputName = 'PFMet',
    OutputName = 'PFType1CorrectedMet',
    JetsName = 'AKt4PFJets',
    RhoAlgo = mithep.PileupEnergyDensity.kFixedGridFastjetAll,
    MaxEMFraction = 0.9,
    SkipMuons = True,
    MinJetPt = 15.
)
metCorrectionMod.ApplyType0(False)
metCorrectionMod.ApplyType1(True)
metCorrectionMod.ApplyShift(False)
metCorrectionMod.IsData(analysis.isRealData)
