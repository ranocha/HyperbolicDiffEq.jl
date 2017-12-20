
"""
    Octic{T}

The scalar conservation law with flux
``
  \\partial_t u + \\partial_x u^8 = 0
``
in one space dimensions using `T` as scalar type.
"""
struct Octic{T} <: ScalarBalanceLaw{T,1} end

function Octic(T=Float64)
  Octic{T}()
end

function Base.show{T}(io::IO, model::Octic{T})
  print(io, "Scalar conservation law {T=", T, "} with flux f(u) = u^8")
end


"""
    flux{T}(u, model::Octic{T,1})

Compute the flux of `u` for `model`.
"""
@inline flux{T}(u, model::Octic{T}) = ((u^2)^2)^2

"""
    speed(u::Real, model::Octic)

Compute the speed f'(`u`) for `model`.
"""
@inline speed(u, model::Octic) = 8*(u^2)^2*u^3


################################################################################

"""
    (::GodunovFlux){T}(uₗ::T, uᵣ::T, model::Octic{T})

Compute Godunov's flux between `uₗ` and `uᵣ` for `model`.
"""
function (::GodunovFlux){T}(uₗ::T, uᵣ::T, model::Octic{T})
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
    (::EnergyConservativeFlux){T}(uₗ::T, uᵣ::T, model::Octic{T})

Compute the energy (L₂ entropy) conservative flux between `uₗ` and `uᵣ` for `model`.
"""
function (flux::EnergyConservativeFlux)(uₗ::T, uᵣ::T, model::Octic{T}) where T<:Real
    ( ((uₗ^2)^2)^2 + (uₗ^2)^2*uₗ^3*uᵣ + (uₗ^3)^2*uᵣ^2 + uₗ^2*uₗ^3*uᵣ^3 + (uₗ^2)^2*(uᵣ^2)^2
    + uₗ^3*uᵣ^2*uᵣ^3 + uₗ^2*(uᵣ^3)^2 + uₗ*(uᵣ^2)^2*uᵣ^3 + ((uᵣ^2)^2)^2 ) / 9
end


################################################################################

@inline function add_numerical_fluxes_inner_loop2!(du, fluxes, u, balance_law::Octic,
                                                    Nx, basis::GaussLegendre,
                                                    jacx, parallel)
    Pp1 = length(basis.nodes)
    Rl = basis.interp_left
    Rr = basis.interp_right
    one_9 = 1 / 9
    two_9 = 2 / 9

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
        Ru5l = zero(eltype(u))
        Ru5r = zero(eltype(u))
        for nx in Base.OneTo(Pp1)
            tmp = u[nx,ix]^2*u[nx,ix]^3
            Ru5l += Rl[nx]*tmp
            Ru5r += Rr[nx]*tmp
        end
        Ru6l = zero(eltype(u))
        Ru6r = zero(eltype(u))
        for nx in Base.OneTo(Pp1)
            tmp = (u[nx,ix]^3)^2
            Ru6l += Rl[nx]*tmp
            Ru6r += Rr[nx]*tmp
        end
        Ru7l = zero(eltype(u))
        Ru7r = zero(eltype(u))
        for nx in Base.OneTo(Pp1)
            tmp = (u[nx,ix]^3)^2*u[nx,ix]
            Ru7l += Rl[nx]*tmp
            Ru7r += Rr[nx]*tmp
        end
        Ru8l = zero(eltype(u))
        Ru8r = zero(eltype(u))
        for nx in Base.OneTo(Pp1)
            tmp = ((u[nx,ix]^2)^2)^2
            Ru8l += Rl[nx]*tmp
            Ru8r += Rr[nx]*tmp
        end

        for nx in Base.OneTo(Pp1)
            du[nx,ix] += (  (fluxes[ix  ]   - two_9 * Ru8l
                                            - two_9 * u[nx,ix] * Ru7l
                                            - two_9 * u[nx,ix]^2 * Ru6l
                                            - two_9 * u[nx,ix]^3 * Ru5l
                                            - one_9 * ((Rul^2)^2)^2
                                            + one_9 * u[nx,ix] * Rul^3*(Rul^2)^2
                                            - one_9 * Ru2l * (Rul^3)^2
                                            + one_9 * u[nx,ix]^2 * (Rul^3)^2
                                            - one_9 * Ru3l * Rul^2*Rul^3
                                            + one_9 * u[nx,ix]^3 * Rul^2*Rul^3
                                            - one_9 * Ru4l * (Rul^2)^2
                                            + one_9 * u[nx,ix] * Ru4l * Rul^3
                                            - one_9 * Ru4l * Ru2l * Rul^2
                                            + one_9 * u[nx,ix]^2 * Ru4l * Rul^2
                                            - one_9 * Ru4l * Ru3l * Rul
                                            + one_9 * u[nx,ix]^3 * Ru4l * Rul
                                            - one_9 * Ru4l^2) * Rl[nx]
                          - (fluxes[ix+1]   - two_9 * Ru8r
                                            - two_9 * u[nx,ix] * Ru7r
                                            - two_9 * u[nx,ix]^2 * Ru6r
                                            - two_9 * u[nx,ix]^3 * Ru5r
                                            - one_9 * ((Rur^2)^2)^2
                                            + one_9 * u[nx,ix] * Rur^3*(Rur^2)^2
                                            - one_9 * Ru2r * (Rur^3)^2
                                            + one_9 * u[nx,ix]^2 * (Rur^3)^2
                                            - one_9 * Ru3r * Rur^2*Rur^3
                                            + one_9 * u[nx,ix]^3 * Rur^2*Rur^3
                                            - one_9 * Ru4r * (Rur^2)^2
                                            + one_9 * u[nx,ix] * Ru4r * Rur^3
                                            - one_9 * Ru4r * Ru2r * Rur^2
                                            + one_9 * u[nx,ix]^2 * Ru4r * Rur^2
                                            - one_9 * Ru4r * Ru3r * Rur
                                            + one_9 * u[nx,ix]^3 * Ru4r * Rur
                                            - one_9 * Ru4r^2) * Rr[nx]
                            ) * jacx / basis.weights[nx]
        end
    end

    nothing
end
