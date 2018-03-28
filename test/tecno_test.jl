using OrdinaryDiffEq, HyperbolicDiffEq

function compute(K, N, fnum, compute_entropy_var, diss)
    xmin = -1.0
    xmax = +1.0
    tspan = (0., 1.0)
    ode_solver = SSPRK104()
    
    mesh = UniformPeriodicMesh1D(xmin, xmax, N)
    balance_law = Cubic()
    reconstruction = ModifiedENO(K)
    semidisc = ScalarUniformPeriodicTecno1D(balance_law, mesh, fnum, diss, reconstruction, compute_entropy_var)
    ode = semidiscretise(semidisc, sinpi, tspan)

    solve(ode, ode_solver, dt=1/N, save_everystep=false);
end

let K = 1
    N = 2^10
    fnum = CentralFlux()
    compute_entropy_var = (u,model) -> u^3
    diss = (wl, wr, model) -> -0.5*max_abs_speed(wl, wr, model) * (wr-wl) # TeCNO
    compute(K, N, fnum, compute_entropy_var, diss)

    fnum = CentralFlux()
    compute_entropy_var = (u,model) -> u^3
    diss = (wl, wr, model) -> -0.5*max_abs_speed(wl, wr, model) * abs(wr-wl)^(K-1) * (wr-wl) #ELW
    compute(K, N, fnum, compute_entropy_var, diss)

    fnum = EnergyConservativeFlux()
    compute_entropy_var = (u,model) -> u
    diss = (wl, wr, model) -> -0.5*max_abs_speed(wl, wr, model) * (wr-wl) # TeCNO
    compute(K, N, fnum, compute_entropy_var, diss)

    fnum = EnergyConservativeFlux()
    compute_entropy_var = (u,model) -> u
    diss = (wl, wr, model) -> -0.5*max_abs_speed(wl, wr, model) * abs(wr-wl)^(K-1) * (wr-wl) #ELW
    compute(K, N, fnum, compute_entropy_var, diss)
end

let K = 2
    N = 2^10
    fnum = CentralFlux()
    compute_entropy_var = (u,model) -> u^3
    diss = (wl, wr, model) -> -0.5*max_abs_speed(wl, wr, model) * (wr-wl) # TeCNO
    compute(K, N, fnum, compute_entropy_var, diss)

    fnum = CentralFlux()
    compute_entropy_var = (u,model) -> u^3
    diss = (wl, wr, model) -> -0.5*max_abs_speed(wl, wr, model) * abs(wr-wl)^(K-1) * (wr-wl) #ELW
    compute(K, N, fnum, compute_entropy_var, diss)

    fnum = EnergyConservativeFlux()
    compute_entropy_var = (u,model) -> u
    diss = (wl, wr, model) -> -0.5*max_abs_speed(wl, wr, model) * (wr-wl) # TeCNO
    compute(K, N, fnum, compute_entropy_var, diss)

    fnum = EnergyConservativeFlux()
    compute_entropy_var = (u,model) -> u
    diss = (wl, wr, model) -> -0.5*max_abs_speed(wl, wr, model) * abs(wr-wl)^(K-1) * (wr-wl) #ELW
    compute(K, N, fnum, compute_entropy_var, diss)
end

let K = 3
    N = 2^10
    fnum = CentralFlux()
    compute_entropy_var = (u,model) -> u^3
    diss = (wl, wr, model) -> -0.5*max_abs_speed(wl, wr, model) * (wr-wl) # TeCNO
    compute(K, N, fnum, compute_entropy_var, diss)

    fnum = CentralFlux()
    compute_entropy_var = (u,model) -> u^3
    diss = (wl, wr, model) -> -0.5*max_abs_speed(wl, wr, model) * abs(wr-wl)^(K-1) * (wr-wl) #ELW
    compute(K, N, fnum, compute_entropy_var, diss)

    fnum = EnergyConservativeFlux()
    compute_entropy_var = (u,model) -> u
    diss = (wl, wr, model) -> -0.5*max_abs_speed(wl, wr, model) * (wr-wl) # TeCNO
    compute(K, N, fnum, compute_entropy_var, diss)

    fnum = EnergyConservativeFlux()
    compute_entropy_var = (u,model) -> u
    diss = (wl, wr, model) -> -0.5*max_abs_speed(wl, wr, model) * abs(wr-wl)^(K-1) * (wr-wl) #ELW
    compute(K, N, fnum, compute_entropy_var, diss)
end
