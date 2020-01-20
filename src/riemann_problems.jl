
"""
    RiemannProblem{Model<:AbstractBalanceLaw{1}, U, T<:Real}

A Riemann problem in one space dimension for `model` with left and right values
`uₗ`,`uᵣ` at `x₀`,`t₀`.
"""
struct RiemannProblem{Model<:AbstractBalanceLaw{1}, U, T<:Real} <: AbstractRiemannProblem
  model::Model
  uₗ::U
  uᵣ::U
  x₀::T
  t₀::T
end

"""
    RiemannProblem(model::AbstractBalanceLaw{1}, uₗ, uᵣ, x₀::Real=0, t₀::Real=0)

Create the `RiemannProblem` for `model` with left and right values  `uₗ`,`uᵣ` at
`x₀`,`t₀`.
"""
function RiemannProblem(model::AbstractBalanceLaw{1}, uₗ, uᵣ, x₀::Real=0, t₀::Real=0)
  @assert(typeof(uₗ) == typeof(uᵣ))
  x₀, t₀ == promote(x₀, t₀)
  RiemannProblem{typeof(model),typeof(uₗ),typeof(x₀)}(model, uₗ, uᵣ, x₀, t₀)
end


function Base.convert(::Type{RiemannProblem{Model,U,T1}}, prob::RiemannProblem{Model,U,T2}) where {Model,U,T1,T2}
  RiemannProblem{Model,U,T1}(prob.model, prob.uₗ, prob.uᵣ, convert(T1,prob.x₀), convert(T1,prob.t₀))
end

function Base.promote_rule(::Type{RiemannProblem{Model,U,T1}}, ::Type{RiemannProblem{Model,U,T2}}) where {Model,U,T1,T2}
  RiemannProblem{Model,U,promote_type(T1,T2)}
end


function (prob::RiemannProblem)(x::Real)
  @unpack x₀, uₗ, uᵣ = prob

  x < x₀ ? uₗ : uᵣ
end


function Base.show(io::IO, prob::RiemannProblem)
  print(io,
    "Riemann problem for ", prob.model, "\n  ",
    "with uₗ = ", prob.uₗ, " and uᵣ = ", prob.uᵣ, "\n  ",
    "at (x₀,t₀) = (", prob.x₀, ", ", prob.t₀, ")")
end


################################################################################


"""
An abstract type representing the solution of a scalar Riemann problem.
"""
abstract type ScalarRiemannSolution <: AbstractRiemannSolution end


function Base.show(io::IO, sol::AbstractRiemannSolution)
  σ⁻, σ⁺ =  minmax_speeds(sol)
  print(io, "Solution  of the ", sol.prob,
        "\n  speeds: σ⁻ = ", σ⁻, ", σ⁺ = ", σ⁺ )
end


@recipe function f(sol::AbstractRiemannSolution)
  xguide --> L"\xi"
  label  --> L"u"
  legend --> false

  σ⁻, σ⁺ =  minmax_speeds(sol)
  if σ⁻ ≈ σ⁺
    σ = one(σ⁻)
  else
    σ = min(abs(σ⁻ - σ⁺), max(abs(σ⁻), abs(σ⁺))/10)
  end
  σ⁻ = σ⁻ - σ
  σ⁺ = σ⁺ + σ

  ξ = range(σ⁻, σ⁺, length=10^3)
  u = sol.(ξ)

  ((ξ, u, sol.prob.model),)
end


@recipe function f(ξumodel::Tuple{Ξ,U,AbstractBalanceLaw{1}}) where {Ξ,U}
    ξ, u, model = ξumodel
    ((ξ, u), )
end


################################################################################

"""
    RiemannProblemTuple{N, Model<:AbstractBalanceLaw{1}, U, T<:Real}

An NTuple of consecutive Riemann problems.
"""
struct RiemannProblemTuple{N, Model<:AbstractBalanceLaw{1}, U, T<:Real}
  tup::NTuple{N, RiemannProblem{Model,U,T}}

  function RiemannProblemTuple{N,Model,U,T}(tup::NTuple{N, RiemannProblem{Model,U,T}}) where{N,Model,U,T}
    for i in 1:N-1
      tup[i].t₀ == tup[i+1].t₀ ||
        error("The Riemann problems must have the same initial time `t₀`.")
    end
    for i in 1:N-1
      tup[i].x₀ < tup[i+1].x₀ ||
        error("The Riemann problems must be sorted in increasing order.")
    end
    for i in 1:N-1
      tup[i].uᵣ == tup[i+1].uₗ ||
        error("The Riemann problems must be compatible.")
    end
    new(tup)
  end
end

function RiemannProblemTuple(tup::NTuple{N,RiemannProblem{Model,U,T}}) where {N,Model,U,T}
  RiemannProblemTuple{N,Model,U,T}(tup)
end

Base.length(tup::RiemannProblemTuple) = length(tup.tup)
Base.iterate(tup::RiemannProblemTuple) = iterate(tup.tup)
Base.iterate(tup::RiemannProblemTuple, state) = iterate(tup.tup, state)


function *(prob1::RiemannProblem{Model,U}, prob2::RiemannProblem{Model,U}) where {Model,U}
  RiemannProblemTuple(promote(prob1, prob2))
end

function *(tup::RiemannProblemTuple{N,Model,U}, prob::RiemannProblem{Model,U}) where {N,Model,U}
  RiemannProblemTuple(promote(tup..., prob))
end

function *(prob::RiemannProblem{Model,U}, tup::RiemannProblemTuple{N,Model,U}) where {N,Model,U}
  RiemannProblemTuple(promote(prob, tup...))
end

function *(tup1::RiemannProblemTuple{N,Model,U}, tup2::RiemannProblemTuple{M,Model,U}) where {N,M,Model,U}
  RiemannProblemTuple(promote(tup1..., tup2...))
end


function Base.show(io::IO, prob::RiemannProblemTuple{N,Model,U,T}) where {N,Model,U,T}
  println(io, "Tuple of ", N, " consecutive Riemann problems:")
  for i in 1:N
    print(io, prob.tup[i], "\n")
  end
end


"""
An NTuple of consecutive Riemann solutions.
"""
struct RiemannSolutionTuple{N,Sol<:AbstractRiemannSolution}
  tup::NTuple{N, Sol}

  function RiemannSolutionTuple{N,Sol}(tup::NTuple{N, Sol}) where {N,Sol}
    for i in 1:N-1
      tup[i].prob.x₀ < tup[i+1].prob.x₀ ||
        error("The Riemann problems must have the same initial time `t₀`.")
    end
    for i in 1:N-1
      tup[i].prob.x₀ < tup[i+1].prob.x₀ ||
        error("The Riemann problems must be sorted in increasing order.")
    end
    new(tup)
  end
end

function RiemannSolutionTuple(tup::NTuple{N,Sol}) where {N,Sol<:AbstractRiemannSolution}
  RiemannSolutionTuple{N,Sol}(tup)
end


function Base.show(io::IO, prob::RiemannSolutionTuple{N,Sol}) where {N,Sol}
  println(io, "Tuple of ", N, " solutions to consecutive Riemann problems:")
  for i in 1:N
    println(io, prob.tup[i])
  end
end

# for broadcasting; treat as scalar
Base.broadcastable(sol::RiemannSolutionTuple) = Ref(sol)


function solve(prob::RiemannProblemTuple)
  RiemannSolutionTuple(map(solve, prob.tup))
end


@inline get_t₀(sol::RiemannSolutionTuple) = sol.tup[1].prob.t₀


"""
Tha maximal time before any of the solutions interact.
"""
function compute_tmax(sol::RiemannSolutionTuple{N,Sol}) where {N,Sol}
  tup = sol.tup
  t₀ = get_t₀(sol)

  tmax = typemax(tup[1].prob.t₀)
  for i in 1:N-1
    _, σₗ⁺ = minmax_speeds(tup[i])
    σᵣ⁻, _ = minmax_speeds(tup[i+1])
    if σₗ⁺ > σᵣ⁻
      tmax = min(tmax, t₀+(tup[i+1].prob.x₀-tup[i].prob.x₀)/(σₗ⁺-σᵣ⁻))
    end
  end
  tmax
end


function (sol::RiemannSolutionTuple{N,Sol})(t::Real, x::Real) where {N,Sol}
  tup = sol.tup
  t₀ = sol.tup[1].prob.t₀

  t > compute_tmax(sol) && error("The solutions of the Riemann problems interact.")

  for i in 1:N
    σ⁻, σ⁺ = minmax_speeds(tup[i])
    if (x-tup[i].prob.x₀)/(t-t₀) <= σ⁺
      return tup[i](t,x)
    end
  end
  tup[end](t,x)
end


@recipe function f(sol::RiemannSolutionTuple)
  xguide --> L"x"
  yguide --> L"u"

  t₀ = get_t₀(sol)
  tmax = compute_tmax(sol)
  t = t₀ + (tmax-t₀)/2

  label --> "\$ u($t) \$"
  legend --> false

  σ⁻, _ = minmax_speeds(sol.tup[1])
  xmin = sol.tup[1].prob.x₀ + (t-t₀)*σ⁻
  _, σ⁺ = minmax_speeds(sol.tup[end])
  xmax = sol.tup[end].prob.x₀ + (t-t₀)*σ⁺

  if xmin ≈ xmax
    Δx = one(xmin)
  else
    Δx = min(abs(xmin - xmax), max(abs(xmin), abs(xmax))/10)
  end
  xmin = xmin - Δx
  xmax = xmax + Δx

  x = range(xmin, xmax, length=10^3)
  u = sol.(t,x)

  ((x, u, sol.tup[1].prob.model),)
end
