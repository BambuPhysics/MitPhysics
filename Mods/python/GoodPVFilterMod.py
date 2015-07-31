from MitAna.TreeMod.bambu import mithep, analysis

goodPVFilterMod = mithep.GoodPVFilterMod(
    MinVertexNTracks = 0,
    MinNDof = 4,
    MaxAbsZ = 24.,
    MaxRho = 2.,
    VertexesName = mithep.Names.gkPVBrn
)
