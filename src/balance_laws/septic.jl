
"""
    Septic{T}

The scalar conservation law with flux
``
  \\partial_t u + \\partial_x u^7 = 0
``
in one space dimensions using `T` as scalar type.
"""
struct Septic{T} <: ScalarBalanceLaw{T,1} end

function Septic(T=Float64)
  Septic{T}()
end

function Base.show(io::IO, model::Septic{T}) where {T}
  print(io, "Scalar conservation law {T=", T, "} with flux f(u) = u^7")
end


"""
    flux(u, model::Septic)

Compute the flux of `u` for `model`.
"""
@inline flux(u, model::Septic) = (u^2)^2*u^3

"""
    speed(u::Real, model::Septic)

Compute the speed f'(`u`) for `model`.
"""
@inline speed(u, model::Septic) = 7*(u^3)^2


################################################################################

"""
    (::GodunovFlux){T}(uₗ::T, uᵣ::T, model::Septic{T})

Compute Godunov's flux between `uₗ` and `uᵣ` for `model`.
"""
function (::GodunovFlux)(uₗ::T, uᵣ::T, model::Septic{T}) where {T}
  flux(uₗ, model)
end


"""
    (::EnergyConservativeFlux)(uₗ::T, uᵣ::T, model::Septic{T}) where {T}

Compute the energy (L₂ entropy) conservative flux between `uₗ` and `uᵣ` for `model`.
"""
function (flux::EnergyConservativeFlux)(uₗ::T, uᵣ::T, model::Septic{T}) where T<:Real
    ((uₗ^2)^2*uₗ^3 + (uₗ^3)^2*uᵣ + uₗ^2*uₗ^3*uᵣ^2 + (uₗ^2)^2*uᵣ^3 + uₗ^3*(uᵣ^2)^2
    + uₗ^2*uᵣ^2*uᵣ^3 + uₗ*(uᵣ^3)^2 + (uᵣ^2)^2*uᵣ^3) / 8
end


################################################################################

@inline function add_numerical_fluxes_inner_loop2!(du, fluxes, u, balance_law::Septic,
                                                    Nx, basis, boundaries_included::Val{false},
                                                    jacx, parallel)
    Pp1 = length(basis.nodes)
    Rl = basis.interp_left
    Rr = basis.interp_right
    one_4 = 1 / 4
    one_8 = 1 / 8

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

        for nx in Base.OneTo(Pp1)
            du[nx,ix] += (  (fluxes[ix  ]   - one_4 * Ru7l
                                            - one_4 * u[nx,ix] * Ru6l
                                            - one_4 * u[nx,ix]^2 * Ru5l
                                            - one_8 * u[nx,ix]^3 * Ru4l
                                            - one_8 * Rul^3*(Rul^2)^2
                                            + one_8 * u[nx,ix] * (Rul^3)^2
                                            - one_8 * Ru2l * Rul^2*Rul^3
                                            + one_8 * u[nx,ix]^2 * Rul^2*Rul^3
                                            - one_8 * Ru3l * (Rul^2)^2
                                            + one_8 * u[nx,ix]^3 * (Rul^2)^2
                                            - one_8 * Ru4l * Rul^3
                                            + one_8 * u[nx,ix] * Ru4l * Rul^2
                                            - one_8 * Ru4l * Ru2l * Rul
                                            + one_8 * u[nx,ix]^2 * Ru4l * Rul
                                            - one_8 * Ru4l * Ru3l) * Rl[nx]
                          - (fluxes[ix+1]   - one_4 * Ru7r
                                            - one_4 * u[nx,ix] * Ru6r
                                            - one_4 * u[nx,ix]^2 * Ru5r
                                            - one_8 * u[nx,ix]^3 * Ru4r
                                            - one_8 * Rur^3*(Rur^2)^2
                                            + one_8 * u[nx,ix] * (Rur^3)^2
                                            - one_8 * Ru2r * Rur^2*Rur^3
                                            + one_8 * u[nx,ix]^2 * Rur^2*Rur^3
                                            - one_8 * Ru3r * (Rur^2)^2
                                            + one_8 * u[nx,ix]^3 * (Rur^2)^2
                                            - one_8 * Ru4r * Rur^3
                                            + one_8 * u[nx,ix] * Ru4r * Rur^2
                                            - one_8 * Ru4r * Ru2r * Rur
                                            + one_8 * u[nx,ix]^2 * Ru4r * Rur
                                            - one_8 * Ru4r * Ru3r) * Rr[nx]
                            ) * jacx / basis.weights[nx]
        end
    end

    nothing
end
