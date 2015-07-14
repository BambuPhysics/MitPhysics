from MitAna.TreeMod.bambu import mithep

goodPVFilterMod = mithep.GoodPVFilterMod(
    MinVertexNTracks = 0,
    MinNDof = 4,
    MaxAbsZ = 24.,
    MaxRho = 2.,
    IsMC = True,
    VertexesName = mithep.Names.gkPVBrn
)
