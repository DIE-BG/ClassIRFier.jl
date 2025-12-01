using LinearAlgebra
using Distributions

module VARUtils

    export VAR, simulate_irf

    abstract type AbstractVAR end

    """
        VAR <: AbstractVAR
    Includes the number of variables `K`, lag order `p` and companion matrix `A`
    """
    struct VAR <: AbstractVAR
        K::Int
        p::Int
        A::Matrix{Float64}

        """
            VAR(K, p, [dist])
        Instantiate a VAR object just giving the number of variables `K` and lag order `p`.
        """
        function VAR(K::Int, p::Int, dist::UnivariateDistribution = Normal())
            @assert p >= 1 "The `p` parameter should be positive"
            @assert p <= 10 "Why do you need a `p` grater than 10?"

            @assert K >= 2 "You need at least 2 variables"
            @assert K <= 10 "Why do you need more than 10 variables?"

            A = generate_A(K, p, dist)
            flag = !check_stability(A)

            while flag
                A = generate_A(K, p, dist)
                flag = !check_stability(A)
            end

            return new(K, p, A)
        end
    end

    function check_stability(A::Matrix{Float64})
        e = eigvals(A)
        n = norm.(e)
        mask = n .< 1
        return all(mask)
    end

    function generate_A(K::Int, p::Int, dist::UnivariateDistribution)
        # the majority of the companion matrix is zeros
        A = zeros(K * p, K * p)

        # we simulate a random draw for the coefficients
        A[1:K, 1:(K * p)] = reshape(rand(dist, K * K * p), K, K * p)

        # below the coefficients we have identities matrices
        A[(K + 1):(K * p), 1:((K * p) - K)] = Matrix{Float64}(I, (K * p) - K, (K * p) - K)

        return A
    end

    """
        simulate_irf(v, h, std)
    Simulate the IRF of a VAR model. The shocks size is given in std.

    Only one variable will receive the shock, but the variable is selected randomly.
    """
    function simulate_irf(v::VAR, h::Int, std::Vector{Float64})
        @assert length(std) == v.K "You need as many std as variables are in the model"
        response = Matrix{Float64}(undef, h, v.K * v.p)

        # the shock will happen in a single variable randomly
        shock_idx = sample(1:v.K)

        u = zeros(v.K * v.p, 1)
        u[shock_idx] = std[shock_idx]

        for t in 1:h
            response[t, :] = (v.A^t) * u
        end

        return response[:, 1:v.K]
    end

    """
        simulate_irf(v, h)
    Simulate the IRF with a shock of size 1.
    """
    function simulate_irf(v::VAR, h::Int)
        std = ones(v.K)
        return simulate_irf(v, h, std)
    end
end
