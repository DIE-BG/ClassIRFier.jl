using DrWatson
@quickactivate "ClassIRFier.jl"

using Distributions
using CairoMakie

includet(srcdir("VARUtils.jl"))
using Main.VARUtils

IRF_DIR = mkpath(datadir("irf"))

# simulation grid for the
grid = Dict(
    :K => 2,
    :p => [1, 2],
    :dist => [Uniform(-1, 1), Normal(0, 1), Uniform(0, 1)],
    :h => [40, 50],
) |> dict_list

# definition of the simulation given the parameters
function makesim(nsim::Int = 100; kwargs...)
    h = kwargs[:h]
    K = kwargs[:K]
    p = kwargs[:p]
    dist = kwargs[:dist]

    @assert h > 0
    @assert K >= 2
    @assert p >= 1

    irf = []

    @info "K=$(K), p=$(p) h=$(h), dist=$(dist)"
    for i in 1:nsim
        v = simulate_irf(VAR(K, p, dist), h)
        irf = [irf..., v...]
    end

    return irf
end


irf = []
for (i, g) in enumerate(grid)
    @info "Config: $(i)/$(length(grid))"
    irf = [irf..., makesim(; g...)...]
end

safesave(joinpath(IRF_DIR, "irf.jld2"), Dict(:sim => irf))

# exploring the irf realizations
using CairoMakie
f = Figure()
ax = Axis(f[1, 1])
for g in irf
    lines!(ax, 1:length(g), g)
end
f
