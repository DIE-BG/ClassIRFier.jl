using DrWatson
@quickactivate "ClassIRFier"

using ClassIRFier
using Distributions
using CairoMakie

# the distributions helps to control the shape of the IRFs, determining the
# oscillation and sign properties
var = VAR(2, 2, Uniform(-0.2, 3))
# with a stable VAR, we can simulate IRFs. the shock size i s set to 1 by default
sim_irf = simulate_irf(var, 40)


@info "Sign 1:" head_sign(sim_irf[1])
@info "Sign 2:" head_sign(sim_irf[2])

@info "Oscillation 1:" head_oscillation(sim_irf[1])
@info "Oscillation 2:" head_oscillation(sim_irf[2])

fig = Figure()
ax1 = Axis(fig[1, 1], title = "IRF: 1")
ax2 = Axis(fig[2, 1], title = "IRF: 2")
lines!(ax1, sim_irf[1])
lines!(ax2, sim_irf[2])
fig
