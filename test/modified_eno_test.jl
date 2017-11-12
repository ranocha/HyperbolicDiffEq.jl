using Base.Test, OrdinaryDiffEq, HyperbolicDiffEq, Roots


function calc_error(balance_law, uₐₙₐ, tmin, tmax, fnum, reconstruction, N, parallel)
    u₀ = x -> uₐₙₐ(tspan[1], x)

    tspan = (tmin, tmax)
    mesh = UniformPeriodicMesh1D(0., 2., N)
    fv = UniformPeriodicReconstructedFV1D(balance_law, mesh, fnum, reconstruction, parallel)
    ode = semidiscretise(fv, u₀, tspan)
    if order(reconstruction) <= 4
        sol = solve(ode, SSPRK104(), dt=0.5*max_dt(tspan[1], ode.u0, fv),
                    save_everystep=false, dense=false)
    else
        sol = solve(ode, SSPRK104(), dt=0.01*max_dt(tspan[1], ode.u0, fv),
                    save_everystep=false, dense=false)
    end
    uana = compute_coefficients(x->uₐₙₐ(tspan[end],x), mesh, order(reconstruction))

    N, norm(sol[end] - uana, 1) / N
end

function calc_order_estimate(Ns_and_errors)
    orders = zeros(length(Ns_and_errors)-1)
    for i in 1:length(orders)
        orders[i] = -log(Ns_and_errors[i+1][2]/Ns_and_errors[i][2]) /
                     log(Ns_and_errors[i+1][1]/Ns_and_errors[i][1])
    end
    mean(orders)
end

function calc_order_estimate(balance_law, uₐₙₐ, tspan, fnum, reconstruction, Ns, parallel=Val{:serial}())
    calc_order_estimate(calc_error.(balance_law, uₐₙₐ, tspan..., fnum, reconstruction, Ns, parallel))
end

# test convergence orders
Ns_lo = 10 .* 2 .^ (2:6)
Ns_hi =  4 .* 2 .^ (2:5)


balance_law = ConstantLinearAdvection()
tspan = (0., 1.)
fnum = LocalLaxFriedrichsFlux()
uₐₙₐ = (t,x) -> sinpi(x-t)

@test calc_order_estimate(balance_law, uₐₙₐ, tspan, fnum, ModifiedENO(1), Ns_lo, Val{:serial}()) > 0.8
@test calc_order_estimate(balance_law, uₐₙₐ, tspan, fnum, ModifiedENO(2), Ns_lo, Val{:serial}()) > 1.6
@test calc_order_estimate(balance_law, uₐₙₐ, tspan, fnum, ModifiedENO(3), Ns_lo, Val{:serial}()) > 2.8
@test calc_order_estimate(balance_law, uₐₙₐ, tspan, fnum, ModifiedENO(4), Ns_lo, Val{:serial}()) > 3.6
@test calc_order_estimate(balance_law, uₐₙₐ, tspan, fnum, ModifiedENO(5), Ns_lo, Val{:serial}()) > 4.8
@test calc_order_estimate(balance_law, uₐₙₐ, tspan, fnum, ModifiedENO(6), Ns_hi, Val{:serial}()) > 5.6
@test calc_order_estimate(balance_law, uₐₙₐ, tspan, fnum, ModifiedENO(7), Ns_hi, Val{:serial}()) > 6.8
@test calc_order_estimate(balance_law, uₐₙₐ, tspan, fnum, ModifiedENO(8), Ns_hi, Val{:serial}()) > 7.5

@test calc_order_estimate(balance_law, uₐₙₐ, tspan, fnum, ModifiedENO(1), Ns_lo, Val{:threads}()) > 0.8
@test calc_order_estimate(balance_law, uₐₙₐ, tspan, fnum, ModifiedENO(2), Ns_lo, Val{:threads}()) > 1.6
@test calc_order_estimate(balance_law, uₐₙₐ, tspan, fnum, ModifiedENO(3), Ns_lo, Val{:threads}()) > 2.8
@test calc_order_estimate(balance_law, uₐₙₐ, tspan, fnum, ModifiedENO(4), Ns_lo, Val{:threads}()) > 3.6
@test calc_order_estimate(balance_law, uₐₙₐ, tspan, fnum, ModifiedENO(5), Ns_lo, Val{:threads}()) > 4.8
@test calc_order_estimate(balance_law, uₐₙₐ, tspan, fnum, ModifiedENO(6), Ns_hi, Val{:threads}()) > 5.6
@test calc_order_estimate(balance_law, uₐₙₐ, tspan, fnum, ModifiedENO(7), Ns_hi, Val{:threads}()) > 6.8
@test calc_order_estimate(balance_law, uₐₙₐ, tspan, fnum, ModifiedENO(8), Ns_hi, Val{:threads}()) > 7.5
