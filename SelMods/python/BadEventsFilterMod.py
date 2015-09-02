from MitAna.TreeMod.bambu import mithep, analysis

badEventsFilterMod = mithep.BadEventsFilterMod('BadEventsFilterMod',
    EEBadScFilter = True,
    CSCTightHaloFilter = True,
    HBHENoiseFilter = analysis.isRealData,
    FillHist = True
)
