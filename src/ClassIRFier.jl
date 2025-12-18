using DrWatson
@quickactivate "ClassIRFier"

module ClassIRFier
    using Reexport

    include(srcdir("VARUtils", "VARUtils.jl"))
    @reexport using .VARUtils


    include(srcdir("Heads", "Heads.jl"))
    @reexport using .Heads
end
