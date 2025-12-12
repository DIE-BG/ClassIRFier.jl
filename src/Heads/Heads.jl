module Heads

export head_sign, head_sign_training_target


"""
    head_sign(ϕ; δ = 1.0e-4, η = 1.0e-4)

Classify an IRF `ϕ` based on the balance between its positive and negative components.

# Arguments
- `ϕ`: An IRF as a vector.
- `δ`: magnitude threshold below which the signal is considered negligible (default `1.0e-4`).
- `η`: dominance tolerance for considering either the positive or negative part as nearly the entire signal (default `1.0e-4`).

# Returns
- `1` when the total magnitude is below `δ`.
- `2` when the positive part dominates.
- `3` when the negative part dominates.
- `4` when neither side clearly dominates.
"""
function head_sign(ϕ::Vector{Float32}; δ = 1.0f-4, η = 1.0f-4)

    # computing the the area of the positive and negative parts
    A_pos = sum(maximum.([ϕ, 0.0f0]))
    A_neg = sum(maximum.([-1.0f0 .* ϕ, 0.0f0]))
    # total area
    A = A_pos + A_neg

    # classification based on the areas

    # if the total area is less that a δ threshold we consider the signal negligible
    if A < δ
        return :negligible
        # if the positive is not negligible and it covers at least (1-η) of the
        # total area we consider it positive
    elseif (A ≥ δ) & (A_pos / A ≥ 1.0f0 - η)
        return :positive
        # if the negative is not negligible and it covers at least (1-η) of the
        # total area we consider it negative
    elseif (A ≥ δ) & (A_neg / A ≥ 1.0f0 - η)
        return :negative
        # if is not negligible but neither side dominates we consider it mixed
    else
        return :mixed
    end
end

end
