using Base.Test, OrdinaryDiffEq, HyperbolicDiffEq

u₀func(x) = sinpi(x)

balance_law = Burgers()
@inferred UniformPeriodicMesh1D(0., 2., 50)
meshx = UniformPeriodicMesh1D(0., 2., 50)
@inferred LobattoLegendre(4)
basis = LobattoLegendre(4)
tspan = (0.0, 0.3)

@inferred UniformPeriodicFluxDiffDisc1D(balance_law, meshx, basis, EnergyConservativeFlux(), GodunovFlux())
semidisc = UniformPeriodicFluxDiffDisc1D(balance_law, meshx, basis, EnergyConservativeFlux(), GodunovFlux())
semidiscretise(semidisc, u₀func, tspan)
ode = semidiscretise(semidisc, u₀func, tspan)
sol = solve(ode, SSPRK104(), dt=max_dt(tspan[1], ode.u0, semidisc), save_everystep=false)
@inferred evaluate_coefficients(sol.u[end], semidisc)
x, u = evaluate_coefficients(sol.u[end], semidisc)

@inferred UniformPeriodicFluxDiffDisc1D(balance_law, meshx, basis, L2L4ConservativeFlux(), GodunovFlux())
semidisc = UniformPeriodicFluxDiffDisc1D(balance_law, meshx, basis, L2L4ConservativeFlux(), GodunovFlux())
ode = semidiscretise(semidisc, u₀func, tspan)
sol = solve(ode, SSPRK104(), dt=max_dt(tspan[1], ode.u0, semidisc), save_everystep=false)
@inferred evaluate_coefficients(sol.u[end], semidisc)
x, u = evaluate_coefficients(sol.u[end], semidisc)
