"""
    head_sign(ϕ; δ = 1.0e-4, η = 1.0e-4)

Classifies an impulse response function (IRF) vector `ϕ` according to the dominance of its positive and negative areas.

# Arguments
- `ϕ`: A vector representing the IRF.
- `δ`: Threshold below which the total area is considered negligible (default: `1.0e-4`).
- `η`: Tolerance for dominance; if one side covers at least `1-η` of the total area, it is considered dominant (default: `1.0e-4`).

# Returns
- `:negligible` if the total area is less than `δ`.
- `:positive` if the positive area covers at least `1-η` of the total area.
- `:negative` if the negative area covers at least `1-η` of the total area.
- `:mixed` if neither side clearly dominates.
"""
function head_sign(ϕ::Vector{Float32}; δ = 1.0f-4, η = 1.0f-4)

    # Compute the area of the positive and negative parts
    A_pos = sum(maximum.([ϕ, 0.0f0]))
    A_neg = sum(maximum.([-1.0f0 .* ϕ, 0.0f0]))
    # Total area
    A = A_pos + A_neg

    # Classify based on the computed areas

    # If the total area is less than δ, the signal is considered negligible
    if A < δ
        return :negligible
        # If the positive area is not negligible and covers at least (1-η) of the total area, classify as positive
    elseif (A ≥ δ) & (A_pos / A ≥ 1.0f0 - η)
        return :positive
        # If the negative area is not negligible and covers at least (1-η) of the total area, classify as negative
    elseif (A ≥ δ) & (A_neg / A ≥ 1.0f0 - η)
        return :negative
        # If neither side dominates, classify as mixed
    else
        return :mixed
    end
end
