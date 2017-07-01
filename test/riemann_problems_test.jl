using Base.Test
using MappedArrays
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


################################################################################

# Shallow water equations
model = ShallowWater()

# Example 5.20 of Holden & Risebro, 2002
uₗ = ShallowWaterVar1D(1., 0.)
uᵣ = ShallowWaterVar1D(0., 0.)
prob = RiemannProblem(model, uₗ, uᵣ)
sol = solve(prob)
x = linspace(-3, 3)
u = sol.(1., x)
h = mappedarray(u->u.h, u)
v = mappedarray(u->(u.h ≈ 0 ? zero(u.hv) : u.hv/u.h), u)
@test all(h[x .< -1] .≈ 1)
@test all(h[-1 .< x .< 2] ≈ (2 .- x[-1 .< x .< 2]).^2 ./ 9)
@test all(h[x .> 2] .≈ 0)
@test all(v[x .< -1] .≈ 0)
@test all(v[-1 .< x .< 2] ≈ 2 .* (1 .+ x[-1 .< x .< 2]) ./ 3)
@test all(v[x .> 2] .≈ 0)

uₗ, uᵣ = uᵣ, uₗ
prob = RiemannProblem(model, uₗ, uᵣ)
sol = solve(prob)
x = linspace(-3, 3)
u = sol.(1., x)
h = mappedarray(u->u.h, u)
v = mappedarray(u->(u.h ≈ 0 ? zero(u.hv) : u.hv/u.h), u)
@test all(h[x .< -2] .≈ 0)
@test all(h[-2 .< x .< 1] ≈ (2 .+ x[-2 .< x .< 1]).^2 ./ 9)
@test all(h[x .> 1] .≈ 1)
@test all(v[x .< -2] .≈ 0)
@test all(v[-2 .< x .< 1] ≈ -2 .* (1 .- x[-2 .< x .< 1]) ./ 3)
@test all(v[x .> 1] .≈ 0)

# Example 5.21 of Holden & Risebro, 2002
h₀ = 1.
v₀ = 2.5
uₗ = ShallowWaterVar1D(h₀, -v₀)
uᵣ = ShallowWaterVar1D(h₀,  v₀)
prob = RiemannProblem(model, uₗ, uᵣ)
sol = solve(prob)
x = linspace(-1-(v₀+sqrt(h₀)), 1+v₀+sqrt(h₀))
u = sol.(1., x)
h = mappedarray(u->u.h, u)
v = mappedarray(u->(u.h ≈ 0 ? zero(u.hv) : u.hv/u.h), u)
idx1 = x .< -(v₀+sqrt(h₀))
idx2 = -(v₀+sqrt(h₀)) .< x .< (2*sqrt(h₀)-v₀)
idx3 = (2*sqrt(h₀)-v₀) .< x .< (v₀-2*sqrt(h₀))
idx4 = (v₀-2*sqrt(h₀)) .< x .< (v₀+sqrt(h₀))
idx5 = (v₀+sqrt(h₀)) .< x
@test all(h[idx1] .≈ h₀)
@test all(h[idx2] .≈ (-v₀ .+ 2 .* sqrt(h₀) .- x[idx2]).^2 ./ 9)
@test all(h[idx3] .≈ 0)
@test all(h[idx4] .≈ (v₀ .- 2 .* sqrt(h₀) .- x[idx4]).^2 ./ 9)
@test all(h[idx5] .≈ h₀)
@test all(v[idx1] .≈ -v₀)
@test all(v[idx2] .≈ (-v₀ .+ 2 .* sqrt(h₀) .+ 2 .* x[idx2]) ./ 3)
@test all(v[idx3] .≈ 0)
@test all(v[idx4] .≈ (v₀ .- 2 .* sqrt(h₀) .+ 2 .* x[idx4]) ./ 3)
@test all(v[idx5] .≈ v₀)
