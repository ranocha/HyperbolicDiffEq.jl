module HyperbolicDiffEq

using Roots
using Parameters

using RecipesBase
using LaTeXStrings

using Reexport
@reexport using DiffEqBase
@reexport using DiffEqPDEBase

# interfaces
import DiffEqBase: solve

import Base: show, *, start, done, next, promote_rule, convert


# types
"""
An abstract type representing a balance law in `Dim` space dimensions.
"""
abstract type AbstractBalanceLaw{Dim} end

"""
An abstract type representing a scalar balance law using `T` as real type in
`Dim` space dimensions.
"""
abstract type ScalarBalanceLaw{T,Dim} <: AbstractBalanceLaw{Dim} end

"""
An abstract type representing a Riemann problem.
"""
abstract type AbstractRiemannProblem <: PDEProblem end

"""
An abstract type representing the solution of a Riemann problem.
"""
abstract type AbstractRiemannSolution end


include("riemann_problems.jl")
include("balance_laws/burgers.jl")


# models
export Burgers, BuckleyLeverette

export godunov
export flux, max_abs_speed

export RiemannProblem, *

end # module
