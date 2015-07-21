from MitAna.TreeMod.bambu import mithep, analysis
from MitPhysics.Mods.PuppiMod import puppiMod

puppiMod.SetEtaConfigName("/home/dabercro/cms/cmssw/041/CMSSW_7_4_6/src/MitPhysics/data/PuppiEta_150701.cfg")
puppiMod.SetR0(0.4)
puppiMod.SetAlpha(2)
puppiMod.SetBeta(2)
puppiMod.SetDump(True)
puppiMod.SetDZCut(0.3)
puppiMod.SetD0Cut(1)

analysis.SetUseHLT(False)
analysis.setSequence(puppiMod)
