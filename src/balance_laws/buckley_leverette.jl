
"""
    BuckleyLeverette{T,Dim}

A struct representing the Buckley-Leverette equation in `Dim` dimensions using
`T` as scalar type.
"""
struct BuckleyLeverette{T} <: ScalarBalanceLaw{T,1} end

function BuckleyLeverette(T=Float64)
  BuckleyLeverette{T}()
end

function show{T}(io::IO, model::BuckleyLeverette{T})
  print(io, "Buckley-Leverette equation {T=", T, "}",
            " with flux f(u) = u^2 / (u^2 + (1-u)^2)")
end


"""
    flux(u, model::BuckleyLeverette)

Compute the flux of `u` for `model`.
"""
@inline flux(u, model::BuckleyLeverette) = u^2 / (u^2 + (1-u)^2)

"""
    speed(u::Real, model::BuckleyLeverette)

Compute the speed f'(`u`) for `model`.
"""
@inline speed(u, model::BuckleyLeverette) = 2*u*(1-u) / (u^2 + (1-u)^2)^2

"""
    max_abs_speed(u, model::BuckleyLeverette)

Compute the maximal absolute value of speed at `u` for `model`.
"""
@inline max_abs_speed(u, model::BuckleyLeverette) = abs(speed(u, model))


################################################################################

"""
    local_lax_friedrichs(uₗ, uᵣ, model::BuckleyLeverette)

Compute the local Lax-Friedrichs flux between `uₗ` and `uᵣ` for `model`.
"""
function local_lax_friedrichs(u, uᵣ, model::BuckleyLeverette)
  # the flux is not convex; the maximal speed is 0.5 at u=0.5
  critical_u = one(T)/2
  if min(uₗ,uᵣ) < critical_u && max(uₗ,uᵣ) > critical_u
    λ = max_abs_speed(critical_u, model)
  else
    λ = max(max_abs_speed(uₗ,model), max_abs_speed(uᵣ,model))
  end
  ( flux(uₗ, model) + flux(uᵣ, model) ) / 2 - λ/2 * (uᵣ-uₗ)
end

"""
    godunov(uₗ, uᵣ, model::BuckleyLeverette)

Compute Godunov's flux between `uₗ` and `uᵣ` for `model`.
"""
function godunov(uₗ, uᵣ, model::BuckleyLeverette)
  flux(uₗ, model)
end


################################################################################

"""
    BuckleyLeveretteRiemannSolution{T,T1}

The solution of a Riemann problem `prob` for the Buckley-Leverette equation.
"""
struct BuckleyLeveretteRiemannSolution{T,T1} <: ScalarRiemannSolution
  prob::RiemannProblem{BuckleyLeverette{T},T,T1}
  σ⁻::T
  σ⁺::T
end


"""
    minmax_speeds(sol::BuckleyLeveretteRiemannSolution)

Return the minimal and maximal speeds `σ⁻, σ⁺` that appear in the solution `sol`.
"""
function minmax_speeds(sol::BuckleyLeveretteRiemannSolution)
  sol.σ⁻, sol.σ⁺
end


"""
    (sol::BuckleyLeveretteRiemannSolution)(ξ::Real)

Evaluate the solution `sol` at the value `ξ` of the self-similarity variable
`ξ = (x - x₀) / (t - t₀)`.
"""
function (sol::BuckleyLeveretteRiemannSolution)(ξ::Real)
  @unpack σ⁻, σ⁺ = sol
  @unpack uₗ, uᵣ = sol.prob

  if ξ < σ⁻
    uₗ
  elseif ξ < σ⁺
    # rarefaction wave
    if uₗ < uᵣ
      ( 1 - sqrt(-1 + 4/(1+sqrt(1+4ξ))) ) / 2
    else
      ( 1 + sqrt(-1 + 4/(1+sqrt(1+4ξ))) ) / 2
    end
  else
    uᵣ
  end
end


"""
    (sol::BuckleyLeveretteRiemannSolution)(t::Real, x::Real)

Evaluate the solution `sol` at the time and space coordinates `t` and `x`.
"""
function (sol::BuckleyLeveretteRiemannSolution)(t::Real, x::Real)
  @unpack x₀, t₀ = sol.prob

  sol((x-x₀)/(t-t₀))
end


"""
    shockspeed(uₗ, uᵣ, model::BuckleyLeverette)

Compute the speed f a shock with left and right state `uₗ`, `uᵣ` for `model`.
"""
function shockspeed(uₗ, uᵣ, model::BuckleyLeverette)
  (-2*uₗ*uᵣ + uₗ + uᵣ) /
    (4*uₗ^2*uᵣ^2 - 4*uₗ^2*uᵣ + 2*uₗ^2 - 4*uₗ*uᵣ^2 + 4*uₗ*uᵣ - 2*uₗ + 2*uᵣ^2 - 2*uᵣ + 1)
end


"""
    solve{T,T1}(prob::RiemannProblem{BuckleyLeverette{T,1},T,T1})

Compute the solution of the Riemann prolem `prob`.
"""
function solve{T,T1}(prob::RiemannProblem{BuckleyLeverette{T,1},T,T1})
  @unpack uₗ, uᵣ, model = prob
  u_crit = one(T) / 2

  if uₗ < 0 || uₗ > 1 || uᵣ < 0 || uᵣ > 1
    error("The values uₗ, uᵣ are not in the valid interval [0,1].")
  end

  # f(u) = u^2 / (u^2 + (1-u)^2) is convex on [0, 0.5] and concave on [0.5, 1]
  if uₗ < uᵣ
    # find the lower convex envelope
    if uₗ < u_crit && uᵣ < u_crit
      # f is convex
      # the solutions contains a single rarefaction wave
      σ⁻ = speed(uₗ, model)
      σ⁺ = speed(uᵣ, model)
    elseif uₗ > u_crit && uᵣ > u_crit
      # f is concave
      # the solutions contains a single shock
      σ⁻ = σ⁺ = shockspeed(uₗ, uᵣ, model)
    else
      # f has both convex and concave parts
      σ⁻ = speed(uₗ, model)
      umin = min(uₗ, uᵣ)
      umax = max(uₗ, uᵣ)
      σ⁺ = speed(
        Roots.fzero(u -> speed(u, model) - shockspeed(u, uᵣ, model), [umin, umax]),
        model)
    end
  else
    # find the upper concave envelope
    if uₗ < u_crit && uᵣ < u_crit
      # f is convex
      # the solutions contains a single shock
      σ⁻ = σ⁺ = shockspeed(uₗ, uᵣ, model)
    elseif uₗ > u_crit && uᵣ > u_crit
      # f is concave
      # the solutions contains a single rarefaction wave
      σ⁻ = speed(uₗ, model)
      σ⁺ = speed(uᵣ, model)
    else
      # f has both convex and concave parts
      σ⁻ = speed(uₗ, model)
      umin = min(uₗ, uᵣ); umin += eps(umin) # otherwise fzero takes umin if umin=0
      umax = max(uₗ, uᵣ)
      σ⁺ = speed(
        Roots.fzero(u -> speed(u, model) - shockspeed(u, uᵣ, model), [umin, umax]),
        model)
    end
  end

  BuckleyLeveretteRiemannSolution(prob, σ⁻, σ⁺)
end
