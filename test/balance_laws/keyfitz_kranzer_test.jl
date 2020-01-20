using Test, OrdinaryDiffEq, HyperbolicDiffEq, DiffEqCallbacks

function runtest()
    xmin = -1.0
    xmax = +0.5
    p = 1
    N = 2^8
    tspan = (0., 1)
    ode_solver = OrdinaryDiffEq.Euler()
    dt(p,N) = 1/((p^2+1)*N)
    u01func(x) = x < 0 ? 1.5 : -1.725862; u02func(x) = x < 0 ? 0.0 : 1.276293
    u0func(x) = KeyfitzKranzerVar(u01func(x), u02func(x))

    mesh = UniformPeriodicMesh1D(xmin, xmax, N)
    basis = LobattoLegendre(p)
    balance_law = KeyfitzKranzer()
    fvol = EnergyConservativeFlux()
    fnum = FluxPlusDissipation(EnergyConservativeFlux(), LocalLaxFriedrichsDissipation())
    semidisc = UniformPeriodicFluxDiffDisc1D(balance_law, mesh, basis, fvol, fnum)
    ode = semidiscretise(semidisc, u0func, tspan)

    save_func = (u,t,integrator) -> integrate(u->IntegralQuantitiesKeyfitzKranzer(u,balance_law), u, integrator.f.f)
    saved_values = SavedValues(Float64, IntegralQuantitiesKeyfitzKranzer{Float64})
    saving = SavingCallback(save_func, saved_values, saveat=range(tspan..., length=10^3))
    stepsize = StepsizeLimiter((u,p,t)->max_dt(t,u,semidisc), safety_factor=1.0)

    sol = solve(ode, ode_solver, dt=dt(p,N), save_everystep=false, callback=CallbackSet(stepsize,saving))

    u1 = map(q->q.u1, saving.affect!.saved_values.saveval)
    u2 = map(q->q.u2, saving.affect!.saved_values.saveval)
    entropy = map(q->q.entropy, saving.affect!.saved_values.saveval)
    @test u1[1] ≈ u1[end]
    @test u2[1] ≈ u2[end]
    @test entropy[1] >= entropy[end]

    nothing
end

runtest()
