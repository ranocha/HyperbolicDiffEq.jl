
"""
    Burgers{T,Dim}

Burgers' equation
``
  \\partial_t u + \\partial_x \\frac{u^2}{2} = 0
``
in `Dim` space dimensions using `T` as scalar type.
"""
struct Burgers{T,Dim} <: ScalarBalanceLaw{T,Dim} end

function Burgers(T=Float64, Dim=1)
  Burgers{T,Dim}()
end

function Base.show(io::IO, model::Burgers{T,Dim}) where {T,Dim}
  print(io, "Burgers' equation {T=", T, ", Dim=", Dim, "}",
            " with flux f(u) = u^2 / 2 * (1,...,1)/sqrt(Dim)")
end


"""
    IntegralQuantitiesBurgers{T<:Real}

Some integrated quantities of interest for Burgers' equation. Can be used in
callbacks of DifferentialEquations.jl.
"""
struct IntegralQuantitiesBurgers{T<:Real} <: FieldVector{2,T}
    mass::T
    energy::T
end

function IntegralQuantitiesBurgers(u, model::ScalarBalanceLaw)
    IntegralQuantitiesBurgers(u, u^2)
end


"""
    flux(u, model::Burgers)

Compute the flux of `u` for `model`.
"""
@inline flux(u, model::Burgers) = u^2/2

"""
    speed(u::Real, model::Burgers)

Compute the speed f'(`u`) for `model`.
"""
@inline speed(u, model::Burgers) = u


@inline flux(u, model::Burgers{T,2}, direction) where {T} = u^2 / (2*sqrt(2))


################################################################################

"""
    (::GodunovFlux)(uₗ::T, uᵣ::T, model::Burgers{T,1}) where {T}

Compute Godunov's flux between `uₗ` and `uᵣ` for `model`.
"""
function (::GodunovFlux)(uₗ::T, uᵣ::T, model::Burgers{T,1}) where {T}
  if uₗ < uᵣ
    if uₗ < 0 && 0 < uᵣ
      zero(T)
    else
      min(flux(uₗ, model), flux(uᵣ, model))
    end
  else
    max(flux(uₗ, model), flux(uᵣ, model))
  end
end


function (::GodunovFlux)(uₗ::T, uᵣ::T, model::Burgers{T,2}, direction) where {T}
  if uₗ < uᵣ
    if uₗ < 0 && 0 < uᵣ
      zero(T)
    else
      min(flux(uₗ, model, direction), flux(uᵣ, model, direction))
    end
  else
    max(flux(uₗ, model, direction), flux(uᵣ, model, direction))
  end
end


function (flux::EnergyConservativeFlux)(uₗ::T, uᵣ::T, model::Burgers{T,1}) where T<:Real
    (uₗ^2 + uₗ*uᵣ + uᵣ^2) / 6
end

function (flux::EnergyConservativeFlux)(uₗ::T, uᵣ::T, model::Burgers{T,2}, direction) where T<:Real
    (uₗ^2 + uₗ*uᵣ + uᵣ^2) / (6*sqrt(6))
end


@inline function (flux::L2L4ConservativeFlux)(uₗ::T, uᵣ::T, model::Burgers{T,1}) where T<:Real
    ( 18 * jump_pow_u_r_over_jump_u(uₗ, uᵣ, Val{5}()) + 5 * jump_pow_u_r_over_jump_u(uₗ, uᵣ, Val{3}()) ) /
        ( 60 * jump_pow_u_r_over_jump_u(uₗ, uᵣ, Val{3}()) + 30 )
end

@generated function (flux::L2L2sConservativeFlux{s})(uₗ, uᵣ, model::Burgers{T,1}) where {s,T<:Real}
    ex = :(
        ( $(6s^2-3s) * jump_pow_u_r_over_jump_u(uₗ, uᵣ, Val{$(2s+1)}())
        + $(2s+1) * jump_pow_u_r_over_jump_u(uₗ, uᵣ, Val{3}()) ) /
            ( $(12s^2+6s) * jump_pow_u_r_over_jump_u(uₗ, uᵣ, Val{$(2s-1)}()) + $(12s+6) )
    )
    return :(Base.@_inline_meta; $ex)
end


################################################################################

@inline function add_numerical_fluxes_inner_loop2!(du, fluxes, u, balance_law::Burgers,
                                                    Nx, basis, boundaries_included::Val{false},
                                                    jacx, parallel)
    Pp1 = length(basis.nodes)
    Rl = basis.interp_left
    Rr = basis.interp_right
    one_3 = 1 / 3
    one_6 = 1 / 6

    # add numerical fluxes
    @inbounds for ix in Base.OneTo(Nx)
        Rul = zero(eltype(u))
        Rur = zero(eltype(u))
        for nx in Base.OneTo(Pp1)
            Rul += Rl[nx]*u[nx,ix]
            Rur += Rr[nx]*u[nx,ix]
        end
        Ru2l = zero(eltype(u))
        Ru2r = zero(eltype(u))
        for nx in Base.OneTo(Pp1)
            tmp = u[nx,ix]^2
            Ru2l += Rl[nx]*tmp
            Ru2r += Rr[nx]*tmp
        end

        #@. du[:,ix] += ((fluxes[ix] - one_3 * Ru2l - one_6 * Rul^2) * Rl
        #                - (fluxes[ix+1] - one_3 * Ru2r - one_6 * Rur^2) * Rr
        #                ) * jacx / basis.weights
        for nx in Base.OneTo(Pp1)
            du[nx,ix] += ((fluxes[ix] - one_3 * Ru2l - one_6 * Rul^2) * Rl[nx]
                            - (fluxes[ix+1] - one_3 * Ru2r - one_6 * Rur^2) * Rr[nx]
                            ) * jacx / basis.weights[nx]
        end
    end

    nothing
end

################################################################################

"""
    BurgersRiemannSolution{T,T1}

The solution of a Riemann problem `prob` for Burgers' equation.
"""
struct BurgersRiemannSolution{T,T1} <: ScalarRiemannSolution
  prob::RiemannProblem{Burgers{T,1},T,T1}
  σ⁻::T
  σ⁺::T
end


"""
    minmax_speeds(sol::BurgersRiemannSolution)

Return the minimal and maximal speeds `σ⁻, σ⁺` that appear in the solution `sol`.
"""
function minmax_speeds(sol::BurgersRiemannSolution)
  sol.σ⁻, sol.σ⁺
end


"""
    (sol::BurgersRiemannSolution)(ξ::Real)

Evaluate the solution `sol` at the value `ξ` of the self-similarity variable
`ξ = (x - x₀) / (t - t₀)`.
"""
function (sol::BurgersRiemannSolution)(ξ::Real)
  @unpack σ⁻, σ⁺ = sol
  @unpack uₗ, uᵣ = sol.prob

  if ξ < σ⁻
    uₗ
  elseif ξ < σ⁺
    uₗ + (ξ-σ⁻)/(σ⁺-σ⁻) * (uᵣ-uₗ)
  else
    uᵣ
  end
end


"""
    (sol::BurgersRiemannSolution)(t::Real, x::Real)

Evaluate the solution `sol` at the time and space coordinates `t` and `x`.
"""
function (sol::BurgersRiemannSolution)(t::Real, x::Real)
  @unpack x₀, t₀ = sol.prob

  sol((x-x₀)/(t-t₀))
end


"""
    solve(prob::RiemannProblem{Burgers{T,1},T,T1}) where {T,T1}

Compute the solution of the Riemann prolem `prob`.
"""
function solve(prob::RiemannProblem{Burgers{T,1},T,T1}) where {T,T1}
  @unpack uₗ, uᵣ = prob
  if uₗ > uᵣ
    σ⁻ = σ⁺ = (uₗ + uᵣ) / 2
  else
    σ⁻ = uₗ
    σ⁺ = uᵣ
  end
  BurgersRiemannSolution(prob, σ⁻, σ⁺)
end
