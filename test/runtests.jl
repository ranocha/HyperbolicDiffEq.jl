using Test
using HyperbolicDiffEq

@elapsed begin

  @time @testset "Riemann Problems" begin include("riemann_problems_test.jl") end
  @time @testset "Finite Volumes" begin include("finite_volumes_test.jl") end
  @time @testset "Monomial Fluxes" begin include("monomial_fluxes_test.jl") end
  @time @testset "Flux Difference" begin include("flux_difference_test.jl") end
  @time @testset "ENO Reconstruction" begin include("modified_eno_test.jl") end
  @time @testset "TeCNO" begin include("tecno_test.jl") end
  @time @testset "Quartic Nonconvex Flux" begin include("balance_laws/quartic_nonconvex_test.jl") end
  @time @testset "Keyfitz-Kranzer System" begin include("balance_laws/keyfitz_kranzer_test.jl") end

end
