module PlanarFields

export
    Grid,
    GridAxis,
    PlanarField

using TypeUtils
using OffsetArrays

using Unitful
using Unitful:
    Length,
    °, rad,
    km, m, cm, mm, µm, nm

using TwoDimensional
using TwoDimensional: ShapeElement, dot, ⋅

using Base: axes1
using Base.Broadcast: broadcasted

include("types.jl")
include("grids.jl")
include("fields.jl")

end
