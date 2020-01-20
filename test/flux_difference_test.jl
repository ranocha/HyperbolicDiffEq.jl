using Test, OrdinaryDiffEq, HyperbolicDiffEq

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

@inferred UniformPeriodicFluxDiffDisc1D(balance_law, meshx, basis, L2L2sConservativeFlux(Val{3}()), GodunovFlux())
semidisc = UniformPeriodicFluxDiffDisc1D(balance_law, meshx, basis, L2L2sConservativeFlux(Val{3}()), GodunovFlux())
ode = semidiscretise(semidisc, u₀func, tspan)
sol = solve(ode, SSPRK104(), dt=max_dt(tspan[1], ode.u0, semidisc), save_everystep=false)
@inferred evaluate_coefficients(sol.u[end], semidisc)
x, u = evaluate_coefficients(sol.u[end], semidisc)



balance_law = HyperbolicDiffEq.Euler{Float64,2}(1.4)
ϱ₀(x,y)  = 1.
vx₀(x,y) = sin(x)*cos(y)
vy₀(x,y) = -cos(x)*sin(y)
p₀(x,y)  = 100/1.4 + (cos(2x)+cos(2y))*3/16
u₀(x,y) = conserved_variables(ϱ₀(x,y), vx₀(x,y), vy₀(x,y), p₀(x,y), balance_law)
xmin=0.; xmax=2π; ymin=0.; ymax=2π
fvol = KennedyGruberFlux()
fnum = fvol
p = 9
Nx = 5
Ny = 5
meshx = UniformPeriodicMesh1D(xmin, xmax, Nx)
meshy = UniformPeriodicMesh1D(ymin, ymax, Ny)
basis = LobattoLegendre(p)
tspan = (0., 14.)
semidisc = UniformPeriodicFluxDiffDisc2D(balance_law, meshx, meshy, basis, fvol, fnum, true)
ode = semidiscretise(semidisc, u₀, tspan)
tstops = range(tspan[1], tspan[end], length=100)
sol = solve(ode, SSPRK33(), dt=0.1/((2p+1)*max(Nx,Ny)), saveat=tstops, tstops=tstops)
