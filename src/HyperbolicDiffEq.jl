__precompile__()

module HyperbolicDiffEq

import Jacobi

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
An abstract type for a numerical flux.
"""
abstract type NumericalFlux end

"""
An abstract type for a basis, e.g. a polynomial basis.
"""
abstract type AbstractBasis end
abstract type NodalBasis <: AbstractBasis end

include("bases/bases1d.jl")

include("meshes/meshes1d.jl")
include("meshes/compute_coefficients.jl")

include("finite_volumes/finite_volumes.jl")
include("finite_volumes/numerical_fluxes.jl")

include("flux_difference/uniform_flux_difference.jl")

include("riemann_problems.jl")
include("balance_laws/burgers.jl")
include("balance_laws/buckley_leverette.jl")
include("balance_laws/shallow_water.jl")
include("balance_laws/euler.jl")



# models
export Burgers, BuckleyLeverette
export ShallowWater, ShallowWaterVar1D
export Euler, EulerVar2D

export flux, max_abs_speed, variables

export RiemannProblem

# bases
export LobattoLegendre

# meshes
export UniformPeriodicMesh1D
export compute_coefficients, compute_coefficients!
export evaluate_coefficients, evaluate_coefficients!

# finite volume methods
export FirstOrderFV
export semidiscretise, max_dt
export local_lax_friedrichs, godunov, suliciu, EnergyConservativeFlux, CentralFlux,
        SuliciuFlux

# flux difference methods
export UniformPeriodicFluxDiffDisc2D


end # module
