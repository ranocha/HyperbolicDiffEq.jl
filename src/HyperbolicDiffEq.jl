__precompile__()

module HyperbolicDiffEq

using Roots
using Parameters
using StaticArrays
using MappedArrays

using RecipesBase
using LaTeXStrings

using Reexport
@reexport using DiffEqBase
@reexport using DiffEqPDEBase

# interfaces
import DiffEqBase: solve

import Base: show, *, start, done, next, promote_rule, convert

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
variables{T}(model::ScalarBalanceLaw{T,1}) = T

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


include("meshes/meshes1d.jl")
include("meshes/compute_coefficients.jl")

include("finite_volumes/finite_volumes.jl")
include("finite_volumes/numerical_fluxes.jl")

include("riemann_problems.jl")
include("balance_laws/burgers.jl")
include("balance_laws/buckley_leverette.jl")
include("balance_laws/shallow_water.jl")



# models
export Burgers, BuckleyLeverette
export ShallowWater, ShallowWaterVar1D

export flux, max_abs_speed, variables

export RiemannProblem

# meshes
export UniformPeriodicMesh1D
export compute_coefficients, compute_coefficients!
export evaluate_coefficients, evaluate_coefficients!

# finite volume methods
export FirstOrderFV
export semidiscretise, max_dt
export local_lax_friedrichs, godunov, suliciu


end # module
