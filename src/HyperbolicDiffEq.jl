__precompile__()

module HyperbolicDiffEq

using Reexport

@reexport using PolynomialBases
import PolynomialBases: integrate, interpolate, interpolate!,
                        compute_coefficients, compute_coefficients!,
                        evaluate_coefficients, evaluate_coefficients!

using Roots
using ArgCheck
using Parameters
using StaticArrays
using MappedArrays

using RecipesBase
using LaTeXStrings

@reexport using DiffEqBase
@reexport using DiffEqPDEBase

# interfaces
import DiffEqBase: solve

import Base: *

import StaticArrays: similar_type


# types
"""
An abstract type for a computational mesh.
"""
abstract type AbstractMesh end

"""
An abstract type representing a balance law in `Dim` space dimensions.
"""
abstract type AbstractBalanceLaw{Dim} end

"""
An abstract type representing a scalar balance law using `T` as real type in
`Dim` space dimensions.
"""
abstract type ScalarBalanceLaw{T,Dim} <: AbstractBalanceLaw{Dim} end
variables{T,Dim}(model::ScalarBalanceLaw{T,Dim}) = T

"""
An abstract type representing a Riemann problem.
"""
abstract type AbstractRiemannProblem <: PDEProblem end

"""
An abstract type representing the solution of a Riemann problem.
"""
abstract type AbstractRiemannSolution end

"""
An abstract type for a semidiscretisation.
"""
abstract type AbstractSemidiscretisation end

"""
An abstract type for a reconstruction.
"""
abstract type AbstractReconstruction end

"""
An abstract type for a numerical flux.
"""
abstract type NumericalFlux end

"""
An abstract type for a dissipation operator.
"""
abstract type DissipationOperator end


include("meshes/meshes1d.jl")
include("meshes/compute_coefficients.jl")

include("finite_volumes/finite_volumes.jl")
include("finite_volumes/numerical_fluxes.jl")
include("finite_volumes/periodic_fv_1d.jl")
include("finite_volumes/reconstructions.jl")
include("finite_volumes/modified_eno.jl")
include("finite_volumes/central_reconstruction.jl")

include("flux_difference/uniform_flux_difference.jl")
include("flux_difference/uniform_flux_difference_1d.jl")

include("riemann_problems.jl")
include("balance_laws/constant_linear_advection.jl")
include("balance_laws/burgers.jl")
include("balance_laws/cubic.jl")
include("balance_laws/quartic.jl")
include("balance_laws/quintic.jl")
include("balance_laws/sextic.jl")
include("balance_laws/septic.jl")
include("balance_laws/octic.jl")
include("balance_laws/buckley_leverette.jl")
include("balance_laws/shallow_water.jl")
include("balance_laws/euler.jl")

include("utils.jl")


# models
export ConstantLinearAdvection, Burgers, IntegralQuantitiesBurgers,
       Cubic, Quartic, Quintic, Sextic, Septic, Octic,
       BuckleyLeverette
export ShallowWater, ShallowWaterVar1D
export Euler, EulerVar1D, EulerVar2D, EulerVar3D, IntegralQuantitiesEuler

export flux, max_abs_speed, naive_max_abs_speed, variables, kinetic_energy,
       entropy, conserved_variables, primitive_variables, entropy_variables,
       satisfies_physical_constraints, flux_potential

export RiemannProblem

# meshes
export UniformMesh1D, UniformPeriodicMesh1D
export compute_coefficients, compute_coefficients!
export evaluate_coefficients, evaluate_coefficients!

# general semidiscretisations
export semidiscretise, max_dt

# finite volume methods
export FirstOrderFV, UniformPeriodicReconstructedFV1D
export CentralReconstruction
export ModifiedENO, ClassicalChoiceENO, BiasedENOChoice,
       MinL2Choice, BiasedMinL2Choice, LexMinLegendreChoice
export GodunovFlux, LocalLaxFriedrichsFlux, HartenLaxVanLeerFlux, HLL, SuliciuFlux,
       KineticFlux
export CentralFlux, MorinishiFlux, DucrosEtAlFlux, KennedyGruberFlux, PirozzoliFlux
export EnergyConservativeFlux, EnergyConservativeFlux1Param, EnergyConservativeFlux2Param,
       IsmailRoeFluxEC, ChandrashekarFluxEC, RanochaFluxECandKEP,
       L2L4ConservativeFlux, L2L2sConservativeFlux
export FluxPlusDissipation, LocalLaxFriedrichsDissipation, ScalarDissipation,
       MatrixDissipation

# flux difference methods
export UniformFluxDiffDisc1D, UniformPeriodicFluxDiffDisc1D,
        UniformPeriodicFluxDiffDisc2D, UniformPeriodicFluxDiffDisc3D

# utilities
export logmean, integrate, order


end # module
