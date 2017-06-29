using Base.Test
using HyperbolicDiffEq

tic()

@time @testset "Riemann Problems" begin include("riemann_problems_test.jl") end

toc()
