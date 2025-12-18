module ClassIRFier
using Reexport

include(joinpath(@__DIR__, "VARUtils", "VARUtils.jl"))
@reexport using .VARUtils


include(joinpath(@__DIR__, "Heads", "Heads.jl"))
@reexport using .Heads
end
