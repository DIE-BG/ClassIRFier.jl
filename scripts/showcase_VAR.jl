using DrWatson
@quickactivate "ClassIRFier"

using ClassIRFier
using UnicodePlots

TEST_DATA_DIR = datadir("test")

# loading some test IRFs whit different levels of smoothness
irf0 = zeros(Float32, 40)
irf1 = Float32.(load(joinpath(TEST_DATA_DIR, "test_4.jld2"))["irf"])
irf2 = load(joinpath(TEST_DATA_DIR, "test_7.jld2"))["irf"]
irf3 = load(joinpath(TEST_DATA_DIR, "test_8.jld2"))["irf"]

function head_smoothness1(ϕ::Vector{<:AbstractFloat})
    # Total Variation TV
    TV = sum(abs.(diff(ϕ)))
    # Normalized Total Variation TV_n
    TV_n = TV / sum(abs.(ϕ))
    # Smoothness measure (the higher, the smoother)
    return 1 / (1 + TV_n)
end

function head_smoothness2(ϕ::Vector{<:AbstractFloat})
    # Roughness
    R = sum((ϕ |> diff |> diff) .^ 2)
    # Normalized Roughness R_n
    R_n = R / sum(ϕ .^ 2)
    # Smoothness measure (the higher, the smoother)
    return 1 / (1 + R_n)
end

# the distributions helps to control the shape of the IRFs, determining the
# oscillation and sign properties

function test_heads(ϕ::Vector{Float32})

    plt1 = lineplot(ϕ, name = "ϕ")
    plt2 = lineplot(abs.(diff(ϕ)), name = "|Δϕ|")
    plt3 = lineplot((ϕ |> diff |> diff) .^ 2, name = "(Δ²ϕ)²")

    UnicodePlots.show(plt1)
    println("\n")
    UnicodePlots.show(plt2)
    println("\n")
    UnicodePlots.show(plt3)
    println("\n")

    logs = """Test Heads Results:
        - Sign: $(head_sign(ϕ))
        - Oscillation: $(head_oscillation(ϕ))
        - Smoothness 1 (TV-based): $(head_smoothness1(ϕ))
        - Smoothness 2 (Roughness-based): $(head_smoothness2(ϕ))
            ⋅ Smoothness avg: $((head_smoothness1(ϕ) + head_smoothness2(ϕ)) / 2)
    """
    @info logs

    return nothing
end

test_heads(irf0)
test_heads(irf1)
test_heads(irf2)
test_heads(irf3)
