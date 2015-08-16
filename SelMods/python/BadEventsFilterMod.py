from MitAna.TreeMod.bambu import mithep

badEventsFilterMod = mithep.BadEventsFilterMod('BadEventsFilterMod',
    EEBadScFilter = True,
    CSCTightHaloFilter = True,
    HBHENoiseFilter = False
)
