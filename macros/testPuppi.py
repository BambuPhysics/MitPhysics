from MitAna.TreeMod.bambu import mithep, analysis

puppiMod = mithep.PuppiMod()
#puppiMod.SetKeepPileup(True)
puppiMod.SetR0(0.4)
puppiMod.SetAlpha(2)
puppiMod.SetBeta(2)
puppiMod.SetDump(True)
puppiMod.SetDZCut(0.3)
puppiMod.SetD0Cut(2)

analysis.SetUseHLT(False)
analysis.setSequence(puppiMod)
