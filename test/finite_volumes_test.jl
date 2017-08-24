using Base.Test
using OrdinaryDiffEq
using HyperbolicDiffEq


function calc_error(balance_law, numflux, N, usethreads=true)
    uₐₙₐ = solve(RiemannProblem(balance_law, 0., 1., -0.5) *
                    RiemannProblem(balance_law, 1., 0., 0.5))
    tspan = (0., 0.5)
    u₀ = x -> uₐₙₐ(tspan[1], x)

    mesh = UniformPeriodicMesh1D(-2., 2., N)
    fv = FirstOrderFV(balance_law, mesh, numflux, usethreads)
    ode = semidiscretise(fv, u₀, tspan)
    sol = solve(ode, Euler(), dt=max_dt(tspan[1], ode.u0, fv), save_everystep=false)
    uana = compute_coefficients(x->uₐₙₐ(tspan[end],x), mesh)

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

function calc_order_estimate(balance_law, numflux, Ns, usethreads=true)
    calc_order_estimate(calc_error.(balance_law, numflux, Ns, usethreads))
end


# test convergence orders
Ns = 100 .* 2 .^ (1:7)
@test calc_order_estimate(Burgers(), godunov, Ns) > 0.8
@test calc_order_estimate(Burgers(), local_lax_friedrichs, Ns) > 0.8
@test calc_order_estimate(BuckleyLeverette(), godunov, Ns) > 0.8
@test calc_order_estimate(BuckleyLeverette(), local_lax_friedrichs, Ns) > 0.8
