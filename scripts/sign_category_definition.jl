using DrWatson
@quickactivate "ClassIRFier.jl"
using MacroModelling
import StatsPlots
using CairoMakie

# Definition of an AR(1) model
@model AR begin
    y[0] = x[0]
    x[0] = a + b * x[-1] + c * e[x]
end

@parameters AR begin
    a = 0
    b = 0.5
    c = 1
end


# this functions checks the "sign" of the IRF by the criteria in notebooks/sign.qmd
function s(ϕ; δ = 1e-4 , η = 1e-4)
    A_pos = sum(maximum.([ϕ, 0]))
    A_neg = sum(maximum.([-1 .* ϕ, 0]))
    A = A_pos + A_neg
    if A < δ
        return 1
    elseif (A ≥ δ) & (A_pos/A ≥ 1 - η)
        return 2
    elseif (A ≥ δ) & (A_neg/A ≥ 1 - η)
        return 3
    else
        return 4
    end
end

# simulate IRF for a grid of parameters for b
grid = Dict(
    :periods => collect(30:10:60),
    :b => collect(-1:0.013:1),
    :c => collect(0.1:0.01:10),
    :negative_shock => [true, false],
) |> dict_list

function single_sim(simconf)
    
    r = get_irf(
        AR, 
        parameters = (
            :b => simconf[:b],
            :c => simconf[:c],
        ); 
        periods = simconf[:periods],
        negative_shock = simconf[:negative_shock],

    )
    
    return r(:x, :, :e) |> collect
end

single_sim(grid[1])

function make_sim(config_list)
    max_t = maximum(
        [conf[:periods] for conf in config_list]
    )
    responses = Matrix{Float64}(undef, max_t, length(config_list))
    
    for (i, conf) in enumerate(config_list)
        responses[1:conf[:periods], i] = single_sim(conf)
    end
    return responses
end

# generating the simulated irf
responses = make_sim(grid)


# respondes categories
responses_cat = mapslices(s, responses; dims=1)

# plotting the positive ones
mask_responses = responses_cat .== 2

f = Figure()
ax = Axis(f[1, 1])
for col in 1:size(responses, 2)
    if mask_responses[col]
        lines!(ax, 1:t_sim, responses[:, col])
    end
end
f