
"""
    Sextic{T}

The scalar conservation law with flux
``
  \\partial_t u + \\partial_x u^6 = 0
``
in one space dimensions using `T` as scalar type.
"""
struct Sextic{T} <: ScalarBalanceLaw{T,1} end

function Sextic(T=Float64)
  Sextic{T}()
end

function Base.show(io::IO, model::Sextic{T}) where {T}
  print(io, "Scalar conservation law {T=", T, "} with flux f(u) = u^6")
end


"""
    flux(u, model::Sextic)

Compute the flux of `u` for `model`.
"""
@inline flux(u, model::Sextic) = (u^3)^2

"""
    speed(u::Real, model::Sextic)

Compute the speed f'(`u`) for `model`.
"""
@inline speed(u, model::Sextic) = 6*u^2*u^3


################################################################################

"""
    (::GodunovFlux)(uₗ::T, uᵣ::T, model::Sextic{T}) where {T}

Compute Godunov's flux between `uₗ` and `uᵣ` for `model`.
"""
function (::GodunovFlux)(uₗ::T, uᵣ::T, model::Sextic{T}) where {T}
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
    (::EnergyConservativeFlux)(uₗ::T, uᵣ::T, model::Sextic{T}) where {T}

Compute the energy (L₂ entropy) conservative flux between `uₗ` and `uᵣ` for `model`.
"""
function (flux::EnergyConservativeFlux)(uₗ::T, uᵣ::T, model::Sextic{T}) where T<:Real
    ((uₗ^3)^2 + uₗ^2*uₗ^3*uᵣ + (uₗ^2)^2*uᵣ^2 + uₗ^3*uᵣ^3 + uₗ^2*(uᵣ^2)^2 + uₗ*uᵣ^2*uᵣ^3 + (uᵣ^3)^2) / 7
end


################################################################################

@inline function add_numerical_fluxes_inner_loop2!(du, fluxes, u, balance_law::Sextic,
                                                    Nx, basis, boundaries_included::Val{false},
                                                    jacx, parallel)
    Pp1 = length(basis.nodes)
    Rl = basis.interp_left
    Rr = basis.interp_right
    one_7 = 1 / 7
    two_7 = 2 / 7

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

        for nx in Base.OneTo(Pp1)
            du[nx,ix] += (  (fluxes[ix  ]   - two_7 * Ru6l
                                            - two_7 * u[nx,ix] * Ru5l
                                            - two_7 * u[nx,ix]^2 * Ru4l
                                            - one_7 * (Rul^3)^2
                                            + one_7 * u[nx,ix] * Rul^3*Rul^2
                                            - one_7 * Ru2l * (Rul^2)^2
                                            + one_7 * u[nx,ix]^2 * (Rul^2)^2
                                            - one_7 * Ru3l * Rul^3
                                            + one_7 * u[nx,ix] * Ru3l * Rul^2
                                            - one_7 * Ru3l * Ru2l * Rul
                                            + one_7 * u[nx,ix]^2 * Ru3l * Rul
                                            - one_7 * Ru3l^2) * Rl[nx]
                          - (fluxes[ix+1]   - two_7 * Ru6r
                                            - two_7 * u[nx,ix] * Ru5r
                                            - two_7 * u[nx,ix]^2 * Ru4r
                                            - one_7 * (Rur^3)^2
                                            + one_7 * u[nx,ix] * Rur^3*Rur^2
                                            - one_7 * Ru2r * (Rur^2)^2
                                            + one_7 * u[nx,ix]^2 * (Rur^2)^2
                                            - one_7 * Ru3r * Rur^3
                                            + one_7 * u[nx,ix] * Ru3r * Rur^2
                                            - one_7 * Ru3r * Ru2r * Rur
                                            + one_7 * u[nx,ix]^2 * Ru3r * Rur
                                            - one_7 * Ru3r^2) * Rr[nx]
                            ) * jacx / basis.weights[nx]
        end
    end

    nothing
end
