from MitAna.TreeMod.bambu import mithep, analysis

puppiMod = mithep.PuppiMod()
#puppiMod.SetKeepPileup(True)
puppiMod.SetR0(0.4)
puppiMod.SetAlpha(2)
puppiMod.SetDump(True)
puppiMod.SetD0Cut(0.5)

analysis.SetUseHLT(False)
analysis.setSequence(puppiMod)
