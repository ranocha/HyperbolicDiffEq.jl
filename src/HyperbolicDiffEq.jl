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

# interfaces
import DiffEqBase: solve

import Base: *

import StaticArrays: similar_type


# types
"""
An abstract type for a computational mesh.
"""
abstract type AbstractMesh end

# for broadcasting; treat as scalar
Base.broadcastable(mesh::AbstractMesh) = Ref(mesh)

"""
An abstract type representing a balance law in `Dim` space dimensions.
"""
abstract type AbstractBalanceLaw{Dim} end

# for broadcasting; treat as scalar
Base.broadcastable(model::AbstractBalanceLaw) = Ref(model)

"""
An abstract type representing a scalar balance law using `T` as real type in
`Dim` space dimensions.
"""
abstract type ScalarBalanceLaw{T,Dim} <: AbstractBalanceLaw{Dim} end
variables(model::ScalarBalanceLaw{T}) where {T} = T

"""
An abstract type representing a Riemann problem.
"""
abstract type AbstractRiemannProblem end

# for broadcasting; treat as scalar
Base.broadcastable(prob::AbstractRiemannProblem) = Ref(prob)

"""
An abstract type representing the solution of a Riemann problem.
"""
abstract type AbstractRiemannSolution end

# for broadcasting; treat as scalar
Base.broadcastable(sol::AbstractRiemannSolution) = Ref(sol)

"""
An abstract type for a semidiscretisation.
"""
abstract type AbstractSemidiscretisation end

# for broadcasting; treat as scalar
Base.broadcastable(semidisc::AbstractSemidiscretisation) = Ref(semidisc)

"""
An abstract type for a reconstruction.
"""
abstract type AbstractReconstruction end

# for broadcasting; treat as scalar
Base.broadcastable(recons::AbstractReconstruction) = Ref(recons)

"""
An abstract type for a numerical flux.
"""
abstract type NumericalFlux end

# for broadcasting; treat as scalar
Base.broadcastable(flux::NumericalFlux) = Ref(flux)

"""
An abstract type for a dissipation operator.
"""
abstract type DissipationOperator end

# for broadcasting; treat as scalar
Base.broadcastable(op::DissipationOperator) = Ref(op)


include("meshes/meshes1d.jl")
include("meshes/compute_coefficients.jl")

include("finite_volumes/finite_volumes.jl")
include("finite_volumes/numerical_fluxes.jl")
include("finite_volumes/periodic_fv_1d.jl")
include("finite_volumes/reconstructions.jl")
include("finite_volumes/modified_eno.jl")
include("finite_volumes/weno_jiang_shu.jl")
include("finite_volumes/central_reconstruction.jl")

include("tecno/tecno.jl")

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
include("balance_laws/quartic_nonconvex.jl")
include("balance_laws/shallow_water.jl")
include("balance_laws/euler.jl")
include("balance_laws/keyfitz_kranzer.jl")

include("utils.jl")


# models
export ConstantLinearAdvection, Burgers, IntegralQuantitiesBurgers,
       Cubic, Quartic, Quintic, Sextic, Septic, Octic,
       BuckleyLeverette,
       QuarticNonconvex
export ShallowWater, ShallowWaterVar1D
export Euler, EulerVar1D, EulerVar2D, EulerVar3D, IntegralQuantitiesEuler
export KeyfitzKranzer, KeyfitzKranzerVar, IntegralQuantitiesKeyfitzKranzer

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
export FirstOrderFV, UniformPeriodicReconstructedFV1D, ScalarUniformPeriodicTecno1D
export CentralReconstruction
export ModifiedENO, ClassicalChoiceENO, BiasedENOChoice,
       MinL2Choice, BiasedMinL2Choice, LexMinLegendreChoice
export WENOJiangShu
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
