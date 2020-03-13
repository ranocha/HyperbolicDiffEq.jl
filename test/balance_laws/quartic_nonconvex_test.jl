using Test, OrdinaryDiffEq, HyperbolicDiffEq, DiffEqCallbacks

function runtest()
    xmin = -1.0
    xmax = +1.0
    p = 1
    N = 2^6
    tspan = (0., 0.02)
    ode_solver = OrdinaryDiffEq.SSPRK33()

    mesh = UniformPeriodicMesh1D(xmin, xmax, N)
    basis = LobattoLegendre(p)
    balance_law = QuarticNonconvex()
    fvol = EnergyConservativeFlux()
    fnum = FluxPlusDissipation(EnergyConservativeFlux(), LocalLaxFriedrichsDissipation(max_abs_speed))
    semidisc = UniformPeriodicFluxDiffDisc1D(balance_law, mesh, basis, fvol, fnum)
    ode = semidiscretise(semidisc, sinpi, tspan)

    save_func = (u,t,integrator) -> integrate(u->IntegralQuantitiesBurgers(u,balance_law), u, integrator.f.f)
    saved_values = SavedValues(Float64, IntegralQuantitiesBurgers{Float64})
    saving = SavingCallback(save_func, saved_values, saveat=range(tspan..., length=10^2))
    stepsize = StepsizeLimiter((u,p,t)->max_dt(t,u,semidisc), safety_factor=0.5)

    sol = solve(ode, ode_solver, dt=0.5*max_dt(tspan[1], ode.u0, semidisc), save_everystep=false, callback=CallbackSet(stepsize,saving))

    mass = first.(saving.affect!.saved_values.saveval)
    entropy = last.(saving.affect!.saved_values.saveval)
    @test mass[1] ≈ mass[end] atol=10*eps()
    @test entropy[1] >= entropy[end]

    nothing
end

runtest()
