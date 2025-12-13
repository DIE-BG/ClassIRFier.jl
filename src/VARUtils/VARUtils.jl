"""
    module VARUtils

Utility functions and types for simulating and working with Vector Autoregressive (VAR) models and their impulse response functions (IRFs).

# Description
The `VARUtils` module provides tools to generate stable VAR models, simulate their impulse response functions (IRFs), and perform related operations. These utilities are designed to support the simulation and analysis of VAR processes, particularly for generating data to train and test models in the broader `ClassIRFier.jl` package.

# Exports
- `VAR`: Type representing a VAR model.
- `simulate_irf`: Functions to simulate IRFs from a VAR model.
"""
module VARUtils

using LinearAlgebra
using Distributions

export VAR, simulate_irf

abstract type AbstractVAR end
"""
    VAR <: AbstractVAR

Type representing a Vector Autoregressive (VAR) model.

# Fields
- `K::Int`: Number of variables in the VAR.
- `p::Int`: Lag order.
- `A::Matrix{Float32}`: Companion matrix of the VAR process.
"""
struct VAR <: AbstractVAR
    K::Int
    p::Int
    A::Matrix{Float32}

    """
        VAR(K, p, A)

    Construct a VAR object by specifying the number of variables `K`, lag order `p`, and a companion matrix `A` directly.

    # Arguments
    - `K::Int`: Number of variables.
    - `p::Int`: Lag order.
    - `A::Matrix{Float32}`: Companion matrix.
    """
    function VAR(K::Int, p::Int, A::Matrix{Float32})
        @assert p >= 1 "The `p` parameter should be positive"
        @assert p <= 10 "Lag order `p` should not exceed 10."

        @assert K >= 2 "At least 2 variables are required."
        @assert K <= 10 "Number of variables `K` should not exceed 10."
        @assert size(A, 1) == K * p && size(A, 2) == K * p "Companion matrix dimensions must be (K*p, K*p)"
        @assert check_stability(A) "The provided companion matrix A is not stable."
        return new(K, p, A)
    end

    """
        VAR(K, p, [dist])

    Construct a stable VAR object by specifying the number of variables `K`, lag order `p`, and optionally a distribution `dist` for the random coefficients (default: `Normal()`).

    The function ensures the generated VAR is stable (all eigenvalues of the companion matrix are inside the unit circle).

    # Arguments
    - `K::Int`: Number of variables (must be between 2 and 10).
    - `p::Int`: Lag order (must be between 1 and 10).
    - `dist::UnivariateDistribution`: Distribution for random coefficients (default: `Normal()`).
    """
    function VAR(K::Int, p::Int, dist::UnivariateDistribution = Normal())
        @assert p >= 1 "The `p` parameter should be positive"
        @assert p <= 10 "Lag order `p` should not exceed 10."

        @assert K >= 2 "At least 2 variables are required."
        @assert K <= 10 "Number of variables `K` should not exceed 10."

        # Generate a stable companion matrix A
        A = generate_A(K, p, dist)
        flag = !check_stability(A)

        while flag
            A = generate_A(K, p, dist)
            flag = !check_stability(A)
        end

        return new(K, p, A)
    end
end

"""
    VAR(coefs::Vector{Matrix{Float32}})

Construct a VAR object from a list of K×K coefficient matrices, assembling the companion matrix automatically.

# Arguments
- `coefs::Vector{Matrix{Float32}}`: List of p coefficient matrices, each of size K×K (Float32).

# Returns
- `VAR`: Instantiated VAR object with the corresponding companion matrix.
"""
function VAR(coefs::Vector{Matrix{Float32}})
    p = length(coefs)
    @assert p >= 1 "At least one coefficient matrix is required."
    K = size(coefs[1], 1)
    @assert all(size(A, 1) == K && size(A, 2) == K for A in coefs) "All matrices must be K×K and of type Float32."

    # Build the companion matrix A (size Kp × Kp)
    A = zeros(Float32, K * p, K * p)
    # Fill the first K rows with the coefficient matrices
    for i in 1:p
        A[1:K, ((i - 1) * K + 1):(i * K)] .= coefs[i]
    end
    # Fill the lower block with identity matrices
    if p > 1
        A[(K + 1):(K * p), 1:(K * (p - 1))] .= Matrix{Float32}(I, K * (p - 1), K * (p - 1))
    end

    @assert check_stability(A) "The companion matrix constructed from the coefficient matrices is not stable."
    return VAR(K, p, A)
end

"""
    check_stability(A)

Check if the companion matrix `A` is stable (all eigenvalues inside the unit circle).

# Arguments
- `A::Matrix{<:AbstractFloat}`: Companion matrix.

# Returns
- `Bool`: `true` if stable, `false` otherwise.
"""
function check_stability(A::Matrix{F}) where {F <: AbstractFloat}
    e = eigvals(A)
    n = norm.(e)
    mask = n .< 1
    return all(mask)
end

"""
    generate_A(K, p, dist)

Generate a random companion matrix for a VAR(K, p) process using the specified distribution.

# Arguments
- `K::Int`: Number of variables.
- `p::Int`: Lag order.
- `dist::UnivariateDistribution`: Distribution for random coefficients.

# Returns
- `Matrix{Float32}`: Generated companion matrix.
"""
function generate_A(K::Int, p::Int, dist::UnivariateDistribution)::Matrix{Float32}
    # The majority of the companion matrix is zeros
    A = zeros(Float32, K * p, K * p)

    # Simulate random coefficients for the first K rows
    A[1:K, 1:(K * p)] = reshape(rand(dist, K * K * p), K, K * p)

    # Below the coefficients, fill with identity matrices to complete the companion form
    A[(K + 1):(K * p), 1:((K * p) - K)] = Matrix{Float32}(I, (K * p) - K, (K * p) - K)

    return A
end

"""
    simulate_irf(v, h, std)

Simulate the impulse response functions (IRFs) of a VAR model to a random shock.

# Arguments
- `v::VAR`: VAR model.
- `h::Int`: Horizon (number of periods to simulate).
- `std::Vector{Float32}`: Shock sizes for each variable (only one variable will receive the shock, chosen randomly).

# Returns
- `Vector{Vector{Float32}}`: List of IRFs, one per variable.
"""
function simulate_irf(v::VAR, h::Int, std::Vector{F}) where {F <: Float32}
    @assert length(std) == v.K "You need as many std as variables are in the model"

    response = Matrix{F}(undef, h, v.K * v.p)

    # The shock will occur in a single variable, chosen randomly
    shock_idx = sample(1:v.K)

    u = zeros(F, v.K * v.p, 1)
    u[shock_idx] = std[shock_idx]

    for t in 1:h
        response[t, :] = (v.A^t) * u
    end

    # Only keep the responses for the original K variables
    response = response[:, 1:v.K]

    return [vec(response[:, i]) for i in 1:size(response, 2)]
end

"""
    simulate_irf(v, h)

Simulate the IRF of a VAR model with a shock of size 1 in a randomly chosen variable.

# Arguments
- `v::VAR`: VAR model.
- `h::Int`: Horizon (number of periods to simulate).

# Returns
- `Vector{Vector{Float32}}`: List of IRFs, one per variable.
"""
function simulate_irf(v::VAR, h::Int)
    std = ones(Float32, v.K)
    return simulate_irf(v, h, std)
end

end # VARUtils
