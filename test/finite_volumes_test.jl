using Base.Test, OrdinaryDiffEq, HyperbolicDiffEq


function calc_error(balance_law, uₐₙₐ, tmin, tmax, numflux, N, usethreads=true)
    u₀ = x -> uₐₙₐ(tspan[1], x)

    tspan = (tmin, tmax)
    mesh = UniformPeriodicMesh1D(-2., 2., N)
    fv = FirstOrderFV(balance_law, mesh, numflux, usethreads)
    ode = semidiscretise(fv, u₀, tspan)
    sol = solve(ode, OrdinaryDiffEq.Euler(), dt=max_dt(tspan[1], ode.u0, fv), save_everystep=false)
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

function calc_order_estimate(balance_law, uₐₙₐ, tspan, numflux, Ns, usethreads=true)
    calc_order_estimate(calc_error.(balance_law, uₐₙₐ, tspan..., numflux, Ns, usethreads))
end


# test convergence orders
Ns = 100 .* 2 .^ (1:7)

balance_law = Burgers()
tspan = (0., 0.5)
uₐₙₐ = solve(RiemannProblem(balance_law, 0., 1., -0.5) *
                RiemannProblem(balance_law, 1., 0., 0.5))
@test calc_order_estimate(balance_law, uₐₙₐ, tspan, GodunovFlux(), Ns) > 0.8
@test calc_order_estimate(balance_law, uₐₙₐ, tspan, LocalLaxFriedrichsFlux(), Ns) > 0.8

balance_law = BuckleyLeverette()
tspan = (0., 0.5)
uₐₙₐ = solve(RiemannProblem(balance_law, 0., 1., -0.5) *
                RiemannProblem(balance_law, 1., 0., 0.5))
@test calc_order_estimate(balance_law, uₐₙₐ, tspan, GodunovFlux(), Ns) > 0.8
@test calc_order_estimate(balance_law, uₐₙₐ, tspan, LocalLaxFriedrichsFlux(), Ns) > 0.8

balance_law = ShallowWater()
tspan = (0., 0.1)
u1 = variables(balance_law)(1., 2.)
u2 = variables(balance_law)(2., 0.)
uₐₙₐ = solve(RiemannProblem(balance_law, u1, u2, -1.) *
                RiemannProblem(balance_law, u2, u1, 0.5))
@test calc_order_estimate(balance_law, uₐₙₐ, tspan, LocalLaxFriedrichsFlux(), Ns) > 0.75
@test calc_order_estimate(balance_law, uₐₙₐ, tspan, SuliciuFlux(), Ns) > 0.75
