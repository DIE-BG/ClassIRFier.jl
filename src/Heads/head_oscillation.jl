"""
        head_oscillation(ϕ::Vector{Float32};
                                         atol::Float32 = 1f-4,
                                         tail_len::Int = 3) -> Symbol

Returns `:oscillatory` if there is at least 1 sign change in the IRF (after:
- trimming the tail that has converged to 0, defined as a final streak of `tail_len`
    consecutive values with `abs(ϕ[t]) ≤ atol`;
- ignoring values with `abs(ϕ[t]) ≤ atol` when evaluating signs).

Otherwise, returns `:no_oscillatory`.
"""
function head_oscillation(
        ϕ::Vector{Float32};
        atol::Float32 = 1.0f-4,
        tail_len::Int = 3
    )::Symbol
    n = length(ϕ)
    n < 2 && return :no_oscillatory

    # 1) Trim tail that has converged to 0:
    #    we look for the final streak (from the end backwards) with abs ≤ atol.
    #    If its length >= tail_len, we trim from the start of that streak.
    last_keep = n
    run = 0
    @inbounds for i in n:-1:1
        if abs(ϕ[i]) ≤ atol
            run += 1
        else
            break
        end
    end
    if run ≥ tail_len
        last_keep = n - run
    end

    # If we trim everything (or it's too short), there is no oscillation.
    last_keep < 2 && return :no_oscillatory

    # 2) Count effective sign changes ignoring |x| ≤ atol
    prev_sign = Int8(0)
    @inbounds for i in 1:last_keep
        x = ϕ[i]
        if abs(x) ≤ atol
            continue
        end
        s = x > 0.0f0 ? Int8(1) : Int8(-1)
        if prev_sign != 0 && s != prev_sign
            return :oscillatory
        end
        prev_sign = s
    end

    return :no_oscillatory
end
