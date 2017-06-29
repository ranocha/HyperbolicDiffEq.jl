using HyperbolicDiffEq
using Base.Test

tic()

@time @testset "Riemann Problems" begin include("riemann_problems_test.jl") end

toc()
