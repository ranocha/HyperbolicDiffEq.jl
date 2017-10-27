using Base.Test
using HyperbolicDiffEq

tic()

@time @testset "Riemann Problems" begin include("riemann_problems_test.jl") end
@time @testset "Finite Volumes" begin include("finite_volumes_test.jl") end
@time @testset "Monomial Fluxes" begin include("finite_volumes_test.jl") end

toc()
