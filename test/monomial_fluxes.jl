using Base.Test, OrdinaryDiffEq, HyperbolicDiffEq, DiffEqCallbacks

function save_func(t, u, integrator)
    meshx = integrator.f.meshx
    basis = integrator.f.basis
    balance_law = integrator.f.balance_law
    integrate(u->IntegralQuantitiesBurgers(u,balance_law), u, meshx, basis)
end

function compute_solution(balance_law, p, N, basis_type, fvol_type, fnum_type, tspan, cfl)
    xmin, xmax = -1., 1.
    meshx = UniformPeriodicMesh1D(xmin, xmax, N)
    if basis_type == :LobattoLegendre
        basis = LobattoLegendre(p)
    elseif basis_type == :GaussLegendre
        basis = GaussLegendre(p)
    else
        error("basis_type $basis_type unknown.")
    end
    parallel = Val{:serial}()

    if fvol_type == :Central
        fvol = CentralFlux()
    elseif fvol_type == :EC
        fvol = EnergyConservativeFlux()
    else
        error("fvol_type $fvol_type unknown.")
    end
    if fnum_type == :Central
        fnumint = CentralFlux()
    elseif fnum_type == :EC
        fnumint = EnergyConservativeFlux()
    elseif fnum_type == :Godunov
        fnumint = GodunovFlux()
    elseif fnum_type == :LLF
        fnumint = LocalLaxFriedrichsFlux()
    elseif fnum_type == :HLL
        fnumint = HartenLaxVanLeerFlux()
    else
        error("fnum_type $fnum_type unknown.")
    end
    semidisc = UniformPeriodicFluxDiffDisc1D(balance_law, meshx, basis, fvol, fnumint)

    u₀= sinpi

    ode = semidiscretise(semidisc, u₀, tspan)
    tstops = linspace(tspan[1], tspan[end], 2)
    maxdt = (t,u) -> max_dt(t, u, semidisc, cfl)
    saved_values = SavedValues(Vector{Float64}(), Vector{IntegralQuantitiesBurgers{Float64}}())
    cb = CallbackSet(StepsizeLimiter(maxdt), SavingCallback(save_func, saved_values))
    sol = solve(ode, SSPRK104(), dt=maxdt(tspan[1],ode.u0), saveat=tstops, tstops=tstops, dense=false,
                  callback=cb)
    sol, saved_values
end

function mass_and_energy_differences(balance_law, basis_type, cfl)
    p = 4
    N = 20
    fvol_type = :EC
    fnum_type = :EC

    if typeof(balance_law) <: Burgers
        tspan = (0., 0.3)
    elseif typeof(balance_law) <: Cubic
        tspan = (0., 0.1)
    elseif typeof(balance_law) <: Quartic
        tspan = (0., 0.065)
    elseif typeof(balance_law) <: Quintic
        tspan = (0., 0.045)
    elseif typeof(balance_law) <: Sextic
        tspan = (0., 0.035)
    elseif typeof(balance_law) <: Septic
        tspan = (0., 0.025)
    elseif typeof(balance_law) <: Octic
        tspan = (0., 0.020)
    else
        error("Balance law $balance_law not supported.")
    end

    sol, saved_values = compute_solution(balance_law, p, N, basis_type, fvol_type, fnum_type, tspan, cfl)
    mass_diff = saved_values.saveval[1].mass-saved_values.saveval[end].mass
    ener_diff = saved_values.saveval[1].energy-saved_values.saveval[end].energy

    mass_diff, ener_diff
end

mass_diff, ener_diff = mass_and_energy_differences(Burgers(), :LobattoLegendre, 0.9)
@test abs(mass_diff) < 10*eps()
@test ener_diff < 5.e-8
mass_diff, ener_diff = mass_and_energy_differences(Burgers(), :LobattoLegendre, 0.09)
@test abs(mass_diff) < 10*eps()
@test ener_diff < 5.e-12
mass_diff, ener_diff = mass_and_energy_differences(Burgers(), :GaussLegendre, 0.9)
@test abs(mass_diff) < 10*eps()
@test ener_diff < 5.e-8
mass_diff, ener_diff = mass_and_energy_differences(Burgers(), :GaussLegendre, 0.09)
@test abs(mass_diff) < 10*eps()
@test ener_diff < 5.e-12

mass_diff, ener_diff = mass_and_energy_differences(Cubic(), :LobattoLegendre, 0.9)
@test abs(mass_diff) < 10*eps()
@test ener_diff < 5.e-8
mass_diff, ener_diff = mass_and_energy_differences(Cubic(), :LobattoLegendre, 0.09)
@test abs(mass_diff) < 10*eps()
@test ener_diff < 5.e-12
mass_diff, ener_diff = mass_and_energy_differences(Cubic(), :GaussLegendre, 0.9)
@test abs(mass_diff) < 10*eps()
@test ener_diff < 5.e-8
mass_diff, ener_diff = mass_and_energy_differences(Cubic(), :GaussLegendre, 0.09)
@test abs(mass_diff) < 10*eps()
@test ener_diff < 5.e-12

mass_diff, ener_diff = mass_and_energy_differences(Quartic(), :LobattoLegendre, 0.9)
@test abs(mass_diff) < 10*eps()
@test ener_diff < 5.e-8
mass_diff, ener_diff = mass_and_energy_differences(Quartic(), :LobattoLegendre, 0.09)
@test abs(mass_diff) < 10*eps()
@test ener_diff < 5.e-12
mass_diff, ener_diff = mass_and_energy_differences(Quartic(), :GaussLegendre, 0.9)
@test abs(mass_diff) < 10*eps()
@test ener_diff < 5.e-8
mass_diff, ener_diff = mass_and_energy_differences(Quartic(), :GaussLegendre, 0.09)
@test abs(mass_diff) < 10*eps()
@test ener_diff < 5.e-12

mass_diff, ener_diff = mass_and_energy_differences(Quintic(), :LobattoLegendre, 0.9)
@test abs(mass_diff) < 10*eps()
@test ener_diff < 5.e-8
mass_diff, ener_diff = mass_and_energy_differences(Quintic(), :LobattoLegendre, 0.09)
@test abs(mass_diff) < 10*eps()
@test ener_diff < 5.e-12
mass_diff, ener_diff = mass_and_energy_differences(Quintic(), :GaussLegendre, 0.9)
@test abs(mass_diff) < 10*eps()
@test ener_diff < 5.e-8
mass_diff, ener_diff = mass_and_energy_differences(Quintic(), :GaussLegendre, 0.09)
@test abs(mass_diff) < 10*eps()
@test ener_diff < 5.e-12

mass_diff, ener_diff = mass_and_energy_differences(Sextic(), :LobattoLegendre, 0.9)
@test abs(mass_diff) < 10*eps()
@test ener_diff < 5.e-8
mass_diff, ener_diff = mass_and_energy_differences(Sextic(), :LobattoLegendre, 0.09)
@test abs(mass_diff) < 10*eps()
@test ener_diff < 5.e-12
mass_diff, ener_diff = mass_and_energy_differences(Sextic(), :GaussLegendre, 0.9)
@test abs(mass_diff) < 10*eps()
@test ener_diff < 5.e-8
mass_diff, ener_diff = mass_and_energy_differences(Sextic(), :GaussLegendre, 0.09)
@test abs(mass_diff) < 10*eps()
@test ener_diff < 5.e-12

mass_diff, ener_diff = mass_and_energy_differences(Septic(), :LobattoLegendre, 0.9)
@test abs(mass_diff) < 10*eps()
@test ener_diff < 5.e-8
mass_diff, ener_diff = mass_and_energy_differences(Septic(), :LobattoLegendre, 0.09)
@test abs(mass_diff) < 10*eps()
@test ener_diff < 5.e-12
mass_diff, ener_diff = mass_and_energy_differences(Septic(), :GaussLegendre, 0.9)
@test abs(mass_diff) < 10*eps()
@test ener_diff < 5.e-8
mass_diff, ener_diff = mass_and_energy_differences(Septic(), :GaussLegendre, 0.09)
@test abs(mass_diff) < 10*eps()
@test ener_diff < 5.e-12

mass_diff, ener_diff = mass_and_energy_differences(Octic(), :LobattoLegendre, 0.9)
@test abs(mass_diff) < 10*eps()
@test ener_diff < 5.e-8
mass_diff, ener_diff = mass_and_energy_differences(Octic(), :LobattoLegendre, 0.09)
@test abs(mass_diff) < 10*eps()
@test ener_diff < 5.e-12
mass_diff, ener_diff = mass_and_energy_differences(Octic(), :GaussLegendre, 0.9)
@test abs(mass_diff) < 10*eps()
@test ener_diff < 5.e-8
mass_diff, ener_diff = mass_and_energy_differences(Octic(), :GaussLegendre, 0.09)
@test abs(mass_diff) < 10*eps()
@test ener_diff < 5.e-12
