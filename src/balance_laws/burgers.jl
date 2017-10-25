
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

function show{T,Dim}(io::IO, model::Burgers{T,Dim})
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

function IntegralQuantitiesBurgers(u, model::Burgers)
    IntegralQuantitiesBurgers(u, u^2)
end


"""
    flux{T}(u, model::Burgers{T,1})

Compute the flux of `u` for `model`.
"""
@inline flux{T}(u, model::Burgers{T,1}) = u^2/2

"""
    max_abs_speed(u, model::Burgers)

Compute the maximal absolute value of speed at `u` for `model`.
"""
@inline max_abs_speed(u, model::Burgers) = abs(u)



@inline flux{T}(u, model::Burgers{T,2}, direction) = u^2 / (2*sqrt(2))


################################################################################

"""
    (::GodunovFlux){T}(uₗ::T, uᵣ::T, model::Burgers{T,1})

Compute Godunov's flux between `uₗ` and `uᵣ` for `model`.
"""
function (::GodunovFlux){T}(uₗ::T, uᵣ::T, model::Burgers{T,1})
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


function (::GodunovFlux){T}(uₗ::T, uᵣ::T, model::Burgers{T,2}, direction)
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


################################################################################

@inline function add_numerical_fluxes_inner_loop2!(du, fluxes, u, balance_law::Burgers,
                                                    Nx, basis::GaussLegendre, jacx, parallel)
    Rl = reshape(interpolation_matrix(-1, basis), length(basis.nodes))
    Rr = reshape(interpolation_matrix(+1, basis), length(basis.nodes))
    MinvRl = Rl ./ basis.weights
    MinvRr = Rr ./ basis.weights
    utmp = zeros(eltype(u), size(u,1))
    one_3 = 1 / 3
    one_6 = 1 / 6

    # add numerical fluxes
    @inbounds for ix in Base.OneTo(Nx)
        @views @. utmp = u[:,ix]
        Rul = dot(Rl, utmp)
        Rur = dot(Rr, utmp)
        @views @. utmp = u[:,ix]^2
        Ru2l = dot(Rl, utmp)
        Ru2r = dot(Rr, utmp)

        @. du[:,ix] += ((fluxes[ix] - one_3 * Ru2l - one_6 * Rul^2) * MinvRl
                        - (fluxes[ix+1] - one_3 * Ru2r - one_6 * Rur^2) * MinvRr
                        ) * jacx
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
    solve{T,T1}(prob::RiemannProblem{Burgers{T,1},T,T1})

Compute the solution of the Riemann prolem `prob`.
"""
function solve{T,T1}(prob::RiemannProblem{Burgers{T,1},T,T1})
  @unpack uₗ, uᵣ = prob
  if uₗ > uᵣ
    σ⁻ = σ⁺ = (uₗ + uᵣ) / 2
  else
    σ⁻ = uₗ
    σ⁺ = uᵣ
  end
  BurgersRiemannSolution(prob, σ⁻, σ⁺)
end
