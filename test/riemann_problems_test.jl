using Base.Test
using HyperbolicDiffEq
using Plots


# Burgers' equation

# Rarefaction wave
prob1 = HyperbolicDiffEq.RiemannProblem(HyperbolicDiffEq.Burgers(), 0., 1.)
print(DevNull, prob1)
@test prob1(0-eps()) ≈ 0
@test prob1(0+eps()) ≈ 1
sol = HyperbolicDiffEq.solve(prob1)
print(DevNull, sol)
@test sol(-1) ≈ 0
@test sol( 2) ≈ 1
ξ = linspace(0, 1, 10)
@test all(sol.(ξ) .== ξ)
@test sol(2., -0.5) ≈ 0
@test sol(2., 1.) ≈ 0.5
@test sol(2., 2.) ≈ 1
plot(sol)

# Stationary shock wave
prob2 = HyperbolicDiffEq.RiemannProblem(HyperbolicDiffEq.Burgers(), 1., -1., 1.)
print(DevNull, prob2)
@test prob2(1-eps()) ≈ 1
@test prob2(1+eps()) ≈ -1
sol = HyperbolicDiffEq.solve(prob2)
print(DevNull, sol)
@test sol(-1) ≈  1
@test sol( 1) ≈ -1
plot(sol)

# Tuples of Riemann problems
prob3 = RiemannProblem(Burgers(), -1., -2., 2.f0)
prob4 = RiemannProblem(Burgers(), -2., 5., 7//2)

@test_throws ErrorException prob3 * prob1
prob1 * prob2
(prob1 * prob2) * prob3
prob1 * (prob2 * prob3)
(prob1 * prob2) * (prob3 * prob4)
prob = prob1 * prob2
print(DevNull, prob)

sol = solve(prob)
print(DevNull, sol)
@test sol(1., -1.) ≈ 0
@test sol(1., 0.) ≈ 0
@test sol(1., 0.5) ≈ 0.5
@test sol(1., 1-eps()) ≈ 1
@test sol(1., 1+eps()) ≈ -1
@test sol(1., 2.) ≈ -1
plot(sol)
