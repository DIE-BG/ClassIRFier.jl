"""
Module `SmoothIRF`

Provides tools for generating and scaling smooth impulse response functions (IRFs) for time series models.
Includes methods for constructing smooth IRFs using the normal distribution and scaling them for further analysis or visualization.

# Exports
- `_scale_irf`: Scale and normalize an IRF vector.
- `_gen_smooth_irf`: Generate a smooth IRF using a normal distribution.

# Usage
```julia
using .SmoothIRF
irf = _gen_smooth_irf(20, 10.0f0, 5.0f0)
scaled_irf = _scale_irf(irf, 2.0f0)
```
"""
module SmoothIRF

using Distributions: Normal, pdf

export _scale_irf, _gen_smooth_irf

"""
_scale_irf(ϕ::Vector{F}, scale::F = 1.0f0) where {F <: Float32}

Scale and normalize an impulse response function (IRF) vector.

# Arguments
- `ϕ::Vector{F}`: IRF values to be scaled (vector of Float32).
- `scale::F`: Scaling factor (default: 1.0f0).

# Returns
- `Vector{F}`: Scaled and normalized IRF vector.
"""

function _scale_irf(ϕ::Vector{F}, scale::F = 1.0f0) where {F <: Float32}
    return scale .* ϕ ./ maximum(ϕ)
end


"""
_gen_smooth_irf(h::Int, max_position::F = 1.0f0, smoothness::F = 8.0f0, scale::F = 1.0f0) where {F <: Float32}

Generate a smooth impulse response function (IRF) of length `h` using a normal distribution.

# Arguments
- `h::Int`: Length of the IRF (must be > 0).
- `max_position::F`: Mean (center) of the normal distribution (default: 1.0f0).
- `smoothness::F`: Standard deviation (spread/smoothness, default: 8.0f0).
- `scale::F`: Scaling factor (default: 1.0f0).

# Returns
- `Vector{F}`: Scaled, smooth IRF vector.

# Example
```julia
irf = _gen_smooth_irf(20, 10.0f0, 5.0f0)
```
"""
function _gen_smooth_irf(h::Int, max_position::F = 1.0f0, smoothness::F = 8.0f0, scale::F = 1.0f0) where {F <: Float32}
    @assert h > 0
    dist = Normal(max_position, smoothness)
    irf = pdf.(dist, 1:h)
    return _scale_irf(irf, scale)
end

end
