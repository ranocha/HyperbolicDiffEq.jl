using Base.Test
using HyperbolicDiffEq
#using Plots # Plots fails on win32 [appveyor]


# Burgers' equation
model = Burgers()

# Rarefaction wave
prob1 = RiemannProblem(model, 0., 1.)
print(DevNull, prob1)
@test prob1(0-eps()) ≈ 0
@test prob1(0+eps()) ≈ 1
sol = solve(prob1)
print(DevNull, sol)
@test sol(-1) ≈ 0
@test sol( 2) ≈ 1
ξ = linspace(0, 1, 10)
@test all(sol.(ξ) .== ξ)
@test sol(2., -0.5) ≈ 0
@test sol(2., 1.) ≈ 0.5
@test sol(2., 2.) ≈ 1
#plot(sol) # Plots fails on win32 [appveyor]
@test flux(sol(0), sol.prob.model) ≈ godunov(sol.prob.uₗ, sol.prob.uᵣ, sol.prob.model)

# Stationary shock wave
prob2 = RiemannProblem(model, 1., -1., 1.)
print(DevNull, prob2)
@test prob2(1-eps()) ≈ 1
@test prob2(1+eps()) ≈ -1
sol = solve(prob2)
print(DevNull, sol)
@test sol(-1) ≈  1
@test sol( 1) ≈ -1
#plot(sol) # Plots fails on win32 [appveyor]
@test flux(sol(0), sol.prob.model) ≈ godunov(sol.prob.uₗ, sol.prob.uᵣ, sol.prob.model)

# Tuples of Riemann problems
prob3 = RiemannProblem(model, -1., -2., 2.f0)
prob4 = RiemannProblem(model, -2., 5., 7//2)

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
#plot(sol) # Plots fails on win32 [appveyor]

us = linspace(-1, 1, 10)
for (uₗ,uᵣ) in Iterators.product(us,us)
  prob = RiemannProblem(model, uₗ, uᵣ)
  sol = solve(prob)
  @test isapprox(godunov(uₗ, uᵣ, model), flux(sol(0), model), atol=1.e-30)
end


################################################################################

# Buckley-Leverette equation
model = BuckleyLeverette()
us = linspace(0, 1, 10)
for (uₗ,uᵣ) in Iterators.product(us,us)
  prob = RiemannProblem(model, uₗ, uᵣ)
  sol = solve(prob)
  @test godunov(uₗ, uᵣ, model) ≈ flux(sol(0), model)
end
