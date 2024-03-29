using Test, OrdinaryDiffEq, HyperbolicDiffEq, Roots
using Statistics: mean
using LinearAlgebra: norm


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

    N, integrate(u->sum(abs.(u)), sol[end]-uana, mesh)
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

println("CentralReconstruction")
@test calc_order_estimate(balance_law, uₐₙₐ, tspan, fnum, CentralReconstruction(1), Ns_lo, Val{:serial}()) > 0.8
@test calc_order_estimate(balance_law, uₐₙₐ, tspan, fnum, CentralReconstruction(3), Ns_lo, Val{:serial}()) > 2.8
@test calc_order_estimate(balance_law, uₐₙₐ, tspan, fnum, CentralReconstruction(5), Ns_lo, Val{:serial}()) > 4.8
@test calc_order_estimate(balance_law, uₐₙₐ, tspan, fnum, CentralReconstruction(7), Ns_hi, Val{:serial}()) > 6.8

@test calc_order_estimate(balance_law, uₐₙₐ, tspan, fnum, CentralReconstruction(1), Ns_lo, Val{:threads}()) > 0.8
@test calc_order_estimate(balance_law, uₐₙₐ, tspan, fnum, CentralReconstruction(3), Ns_lo, Val{:threads}()) > 2.8
@test calc_order_estimate(balance_law, uₐₙₐ, tspan, fnum, CentralReconstruction(5), Ns_lo, Val{:threads}()) > 4.8
@test calc_order_estimate(balance_law, uₐₙₐ, tspan, fnum, CentralReconstruction(7), Ns_hi, Val{:threads}()) > 6.8

println("ModifiedENO")
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

println("MinL2Choice")
@test calc_order_estimate(balance_law, uₐₙₐ, tspan, fnum, ModifiedENO(Val{1}(), MinL2Choice()), Ns_lo, Val{:serial}()) > 0.8
@test calc_order_estimate(balance_law, uₐₙₐ, tspan, fnum, ModifiedENO(Val{2}(), MinL2Choice()), Ns_lo, Val{:serial}()) > 1.6
@test calc_order_estimate(balance_law, uₐₙₐ, tspan, fnum, ModifiedENO(Val{3}(), MinL2Choice()), Ns_lo, Val{:serial}()) > 2.8
@test 0.2 < calc_order_estimate(balance_law, uₐₙₐ, tspan, fnum, ModifiedENO(Val{4}(), MinL2Choice()), Ns_lo, Val{:serial}()) < 0.4
@test 0.2 < calc_order_estimate(balance_law, uₐₙₐ, tspan, fnum, ModifiedENO(Val{5}(), MinL2Choice()), Ns_lo, Val{:serial}()) < 0.4
@test 0.5 < calc_order_estimate(balance_law, uₐₙₐ, tspan, fnum, ModifiedENO(Val{6}(), MinL2Choice()), Ns_hi, Val{:serial}()) < 1.1
@test 0.5 < calc_order_estimate(balance_law, uₐₙₐ, tspan, fnum, ModifiedENO(Val{7}(), MinL2Choice()), Ns_hi, Val{:serial}()) < 1.1
@test 0.5 < calc_order_estimate(balance_law, uₐₙₐ, tspan, fnum, ModifiedENO(Val{8}(), MinL2Choice()), Ns_hi, Val{:serial}()) < 1.1

println("BiasedMinL2Choice")
@test calc_order_estimate(balance_law, uₐₙₐ, tspan, fnum, ModifiedENO(Val{1}(), BiasedMinL2Choice()), Ns_lo, Val{:serial}()) > 0.8
@test calc_order_estimate(balance_law, uₐₙₐ, tspan, fnum, ModifiedENO(Val{2}(), BiasedMinL2Choice()), Ns_lo, Val{:serial}()) > 1.6
@test calc_order_estimate(balance_law, uₐₙₐ, tspan, fnum, ModifiedENO(Val{3}(), BiasedMinL2Choice()), Ns_lo, Val{:serial}()) > 2.8
@test calc_order_estimate(balance_law, uₐₙₐ, tspan, fnum, ModifiedENO(Val{4}(), BiasedMinL2Choice()), Ns_lo, Val{:serial}()) > 3.6
@test calc_order_estimate(balance_law, uₐₙₐ, tspan, fnum, ModifiedENO(Val{5}(), BiasedMinL2Choice()), Ns_lo, Val{:serial}()) > 4.8

println("WENOJiangShu")
@test calc_order_estimate(balance_law, uₐₙₐ, tspan, fnum, WENOJiangShu{5}(), Ns_lo, Val{:serial}()) > 4.8


# Order reduction observed by Rogerson & Meiburg and Shu (1990).
balance_law = ConstantLinearAdvection()
tspan = (0., 1.)
fnum = LocalLaxFriedrichsFlux()
uₐₙₐ = (t,x) -> (sinpi(x-t)^2)^2
println("ClassicalChoiceENO")
@test calc_order_estimate(balance_law, uₐₙₐ, tspan, fnum, ModifiedENO(Val{2}(), ClassicalChoiceENO()), Ns_lo, Val{:serial}()) > 1.6
@test calc_order_estimate(balance_law, uₐₙₐ, tspan, fnum, ModifiedENO(Val{3}(), ClassicalChoiceENO()), Ns_lo, Val{:serial}()) > 2.2
@test calc_order_estimate(balance_law, uₐₙₐ, tspan, fnum, ModifiedENO(Val{4}(), ClassicalChoiceENO()), Ns_lo, Val{:serial}()) < 2.2
@test calc_order_estimate(balance_law, uₐₙₐ, tspan, fnum, ModifiedENO(Val{5}(), ClassicalChoiceENO()), Ns_lo, Val{:serial}()) < 1.9

println("BiasedENOChoice")
@test calc_order_estimate(balance_law, uₐₙₐ, tspan, fnum, ModifiedENO(Val{2}(), BiasedENOChoice()), Ns_lo, Val{:serial}()) > 1.6
@test calc_order_estimate(balance_law, uₐₙₐ, tspan, fnum, ModifiedENO(Val{3}(), BiasedENOChoice()), Ns_lo, Val{:serial}()) > 2.6
@test calc_order_estimate(balance_law, uₐₙₐ, tspan, fnum, ModifiedENO(Val{4}(), BiasedENOChoice()), Ns_lo, Val{:serial}()) > 3.4
@test calc_order_estimate(balance_law, uₐₙₐ, tspan, fnum, ModifiedENO(Val{5}(), BiasedENOChoice()), Ns_lo, Val{:serial}()) > 4.6



# test evaluation
mesh = UniformPeriodicMesh1D(0., 2., 50)

K = 1
fv = UniformPeriodicReconstructedFV1D(balance_law, mesh, fnum, ModifiedENO(K), Val{:serial})
u = compute_coefficients(sinpi, mesh, K)
xplot, uplot = evaluate_coefficients(u, fv)
@test norm(sinpi.(xplot) - uplot) < 5.e-1
K = 2
fv = UniformPeriodicReconstructedFV1D(balance_law, mesh, fnum, ModifiedENO(K), Val{:serial})
u = compute_coefficients(sinpi, mesh, K)
xplot, uplot = evaluate_coefficients(u, fv)
@test norm(sinpi.(xplot) - uplot) < 4.e-2
K = 3
fv = UniformPeriodicReconstructedFV1D(balance_law, mesh, fnum, ModifiedENO(K), Val{:serial})
u = compute_coefficients(sinpi, mesh, K)
xplot, uplot = evaluate_coefficients(u, fv)
@test norm(sinpi.(xplot) - uplot) < 2.e-3
K = 4
fv = UniformPeriodicReconstructedFV1D(balance_law, mesh, fnum, ModifiedENO(K), Val{:serial})
u = compute_coefficients(sinpi, mesh, K)
xplot, uplot = evaluate_coefficients(u, fv)
@test norm(sinpi.(xplot) - uplot) < 2.e-4
K = 5
fv = UniformPeriodicReconstructedFV1D(balance_law, mesh, fnum, ModifiedENO(K), Val{:serial})
u = compute_coefficients(sinpi, mesh, K)
xplot, uplot = evaluate_coefficients(u, fv)
@test norm(sinpi.(xplot) - uplot) < 7.e-6
K = 6
fv = UniformPeriodicReconstructedFV1D(balance_law, mesh, fnum, ModifiedENO(K), Val{:serial})
u = compute_coefficients(sinpi, mesh, K)
xplot, uplot = evaluate_coefficients(u, fv)
@test norm(sinpi.(xplot) - uplot) < 5.e-7
K = 7
fv = UniformPeriodicReconstructedFV1D(balance_law, mesh, fnum, ModifiedENO(K), Val{:serial})
u = compute_coefficients(sinpi, mesh, K)
xplot, uplot = evaluate_coefficients(u, fv)
@test norm(sinpi.(xplot) - uplot) < 3.e-8
K = 8
fv = UniformPeriodicReconstructedFV1D(balance_law, mesh, fnum, ModifiedENO(K), Val{:serial})
u = compute_coefficients(sinpi, mesh, K)
xplot, uplot = evaluate_coefficients(u, fv)
@test norm(sinpi.(xplot) - uplot) < 2.e-9
