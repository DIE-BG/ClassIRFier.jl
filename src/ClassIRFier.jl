module ClassIRFier
using Reexport

include(joinpath(@__DIR__, "VARUtils", "VARUtils.jl"))
@reexport using .VARUtils


include(joinpath(@__DIR__, "Heads", "Heads.jl"))
@reexport using .Heads

include(joinpath(@__DIR__, "SmoothIRF", "SmoothIRF.jl"))
@reexport using .SmoothIRF

end
