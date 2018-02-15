
"""
    Cubic{T}

The scalar conservation law with flux
``
  \\partial_t u + \\partial_x u^3 = 0
``
in one space dimensions using `T` as scalar type.
"""
struct Cubic{T} <: ScalarBalanceLaw{T,1} end

function Cubic(T=Float64)
  Cubic{T}()
end

function Base.show{T}(io::IO, model::Cubic{T})
  print(io, "Scalar conservation law {T=", T, "} with flux f(u) = u^3")
end


"""
    flux{T}(u, model::Cubic{T,1})

Compute the flux of `u` for `model`.
"""
@inline flux{T}(u, model::Cubic{T}) = u^3

"""
    speed(u::Real, model::Cubic)

Compute the speed f'(`u`) for `model`.
"""
@inline speed(u, model::Cubic) = 3u^2


################################################################################

"""
    (::GodunovFlux){T}(uₗ::T, uᵣ::T, model::Cubic{T})

Compute Godunov's flux between `uₗ` and `uᵣ` for `model`.
"""
function (::GodunovFlux){T}(uₗ::T, uᵣ::T, model::Cubic{T})
  flux(uₗ, model)
end


"""
    (::EnergyConservativeFlux){T}(uₗ::T, uᵣ::T, model::Cubic{T})

Compute the energy (L₂ entropy) conservative flux between `uₗ` and `uᵣ` for `model`.
"""
function (flux::EnergyConservativeFlux)(uₗ::T, uᵣ::T, model::Cubic{T}) where T<:Real
    (uₗ^3 + uₗ^2*uᵣ + uₗ*uᵣ^2 + uᵣ^3) / 4
end


################################################################################

@inline function add_numerical_fluxes_inner_loop2!(du, fluxes, u, balance_law::Cubic,
                                                    Nx, basis, boundaries_included::Val{false},
                                                    jacx, parallel)
    Pp1 = length(basis.nodes)
    Rl = basis.interp_left
    Rr = basis.interp_right
    one_2 = 1 / 2
    one_4 = 1 / 4

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
        Ru3l = zero(eltype(u))
        Ru3r = zero(eltype(u))
        for nx in Base.OneTo(Pp1)
            tmp = u[nx,ix]^3
            Ru3l += Rl[nx]*tmp
            Ru3r += Rr[nx]*tmp
        end

        for nx in Base.OneTo(Pp1)
            du[nx,ix] += (  (fluxes[ix  ]   - one_2 * Ru3l
                                            - one_4 * u[nx,ix] * Ru2l
                                            - one_4 * Rul^3
                                            + one_4 * u[nx,ix] * Rul^2
                                            - one_4 * Ru2l * Rul) * Rl[nx]
                          - (fluxes[ix+1]   - one_2 * Ru3r
                                            - one_4 * u[nx,ix] * Ru2r
                                            - one_4 * Rur^3
                                            + one_4 * u[nx,ix] * Rur^2
                                            - one_4 * Ru2r * Rur) * Rr[nx]
                            ) * jacx / basis.weights[nx]
        end
    end

    nothing
end


################################################################################

"""
    CubicRiemannSolution{T,T1}

The solution of a Riemann problem `prob` for the cubic conservation law.
"""
immutable CubicRiemannSolution{T,T1} <: ScalarRiemannSolution
  prob::RiemannProblem{Cubic{T},T,T1}
  σ⁻::T
  σ⁺::T
end


"""
    minmax_speeds(sol::CubicRiemannSolution)

Return the minimal and maximal speeds `σ⁻, σ⁺` that appear in the solution `sol`.
"""
function minmax_speeds(sol::CubicRiemannSolution)
  sol.σ⁻, sol.σ⁺
end


"""
    (sol::CubicRiemannSolution)(ξ::Real)

Evaluate the solution `sol` at the value `ξ` of the self-similarity variable
`ξ = (x - x₀) / (t - t₀)`.
"""
function (sol::CubicRiemannSolution)(ξ::Real)
  @unpack σ⁻, σ⁺ = sol
  @unpack uₗ, uᵣ = sol.prob

  if ξ < σ⁻
    uₗ
  elseif ξ < σ⁺
    # rarefaction wave
    sign(uᵣ)*sqrt(ξ/3)
  else
    uᵣ
  end
end


"""
    (sol::CubicRiemannSolution)(t::Real, x::Real)

Evaluate the solution `sol` at the time and space coordinates `t` and `x`.
"""
function (sol::CubicRiemannSolution)(t::Real, x::Real)
  @unpack x₀, t₀ = sol.prob

  sol((x-x₀)/(t-t₀))
end


function shockspeed(uₗ, uᵣ, model::Cubic)
  uₗ^2 + uₗ*uᵣ + uᵣ^2
end


"""
    solve{T,T1}(prob::RiemannProblem{Cubic{T},T,T1})

Compute the solution of the Riemann prolem `prob`.
"""
function solve{T,T1}(prob::RiemannProblem{Cubic{T},T,T1})
  @unpack uₗ, uᵣ, model = prob
  u_crit = zero(T)

  # f(u) = u^3 is concave on (-Inf,0) and convex on (0,Inf)
  if uₗ < uᵣ
    # find the lower convex envelope
    if uₗ < u_crit && uᵣ < u_crit
      # f is concave
      # the solutions contains a single shock
      σ⁻ = σ⁺ = shockspeed(uₗ, uᵣ, model)
    elseif uₗ > u_crit && uᵣ > u_crit
      # f is convex
      # the solutions contains a single rarefaction wave
      σ⁻ = speed(uₗ, model)
      σ⁺ = speed(uᵣ, model)
    else
      # f has both concave and convex parts
      σ⁻ = speed(-uₗ/2, model)
      σ⁺ = speed(uᵣ, model)
    end
  else
    # find the upper concave envelope
    if uₗ < u_crit && uᵣ < u_crit
      # f is concave
      # the solutions contains a single rarefaction wave
      σ⁻ = speed(uₗ, model)
      σ⁺ = speed(uᵣ, model)
    elseif uₗ > u_crit && uᵣ > u_crit
      # f is convex
      # the solutions contains a single shock
      σ⁻ = σ⁺ = shockspeed(uₗ, uᵣ, model)
    else
      # f has both convex and concave parts
      σ⁻ = speed(-uₗ/2, model)
      σ⁺ = speed(uᵣ, model)
    end
  end

  CubicRiemannSolution(prob, σ⁻, σ⁺)
end
