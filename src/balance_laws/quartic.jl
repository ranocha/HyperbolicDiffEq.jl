
"""
    Quartic{T}

The scalar conservation law with flux
``
  \\partial_t u + \\partial_x u^4 = 0
``
in one space dimensions using `T` as scalar type.
"""
struct Quartic{T} <: ScalarBalanceLaw{T,1} end

function Quartic(T=Float64)
  Quartic{T}()
end

function Base.show(io::IO, model::Quartic{T}) where {T}
  print(io, "Scalar conservation law {T=", T, "} with flux f(u) = u^4")
end


"""
    flux(u, model::Quartic)

Compute the flux of `u` for `model`.
"""
@inline flux(u, model::Quartic) = (u^2)^2

"""
    speed(u::Real, model::Quartic)

Compute the speed f'(`u`) for `model`.
"""
@inline speed(u, model::Quartic) = 4u^3


################################################################################

"""
    (::GodunovFlux)(uₗ::T, uᵣ::T, model::Quartic{T}) where {T}

Compute Godunov's flux between `uₗ` and `uᵣ` for `model`.
"""
function (::GodunovFlux)(uₗ::T, uᵣ::T, model::Quartic{T}) where {T}
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


"""
    (::EnergyConservativeFlux)(uₗ::T, uᵣ::T, model::Quartic{T}) where {T}

Compute the energy (L₂ entropy) conservative flux between `uₗ` and `uᵣ` for `model`.
"""
function (flux::EnergyConservativeFlux)(uₗ::T, uᵣ::T, model::Quartic{T}) where T<:Real
    ((uₗ^2)^2 + uₗ^3*uᵣ + uₗ^2*uᵣ^2 + uₗ*uᵣ^3 + (uᵣ^2)^2) / 5
end


################################################################################

@inline function add_numerical_fluxes_inner_loop2!(du, fluxes, u, balance_law::Quartic,
                                                    Nx, basis, boundaries_included::Val{false},
                                                    jacx, parallel)
    Pp1 = length(basis.nodes)
    Rl = basis.interp_left
    Rr = basis.interp_right
    one_5 = 1 / 5
    two_5 = 2 / 5

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
        Ru4l = zero(eltype(u))
        Ru4r = zero(eltype(u))
        for nx in Base.OneTo(Pp1)
            tmp = (u[nx,ix]^2)^2
            Ru4l += Rl[nx]*tmp
            Ru4r += Rr[nx]*tmp
        end

        for nx in Base.OneTo(Pp1)
            du[nx,ix] += (  (fluxes[ix  ]   - two_5 * Ru4l
                                            - one_5 * Ru2l^2
                                            - one_5 * Ru2l * Rul^2
                                            - one_5 * (Rul^2)^2
                                            - two_5 * u[nx,ix] * Ru3l
                                            + one_5 * u[nx,ix] * Ru2l * Rul
                                            + one_5 * u[nx,ix] * Rul^3) * Rl[nx]
                          - (fluxes[ix+1]   - two_5 * Ru4r
                                            - one_5 * Ru2r^2
                                            - one_5 * Ru2r * Rur^2
                                            - one_5 * (Rur^2)^2
                                            - two_5 * u[nx,ix] * Ru3r
                                            + one_5 * u[nx,ix] * Ru2r * Rur
                                            + one_5 * u[nx,ix] * Rur^3) * Rr[nx]
                            ) * jacx / basis.weights[nx]
        end
    end

    nothing
end
