
"""
    Quintic{T}

The scalar conservation law with flux
``
  \\partial_t u + \\partial_x u^5 = 0
``
in one space dimensions using `T` as scalar type.
"""
struct Quintic{T} <: ScalarBalanceLaw{T,1} end

function Quintic(T=Float64)
  Quintic{T}()
end

function Base.show{T}(io::IO, model::Quintic{T})
  print(io, "Scalar conservation law {T=", T, "} with flux f(u) = u^5")
end


"""
    flux{T}(u, model::Quintic{T,1})

Compute the flux of `u` for `model`.
"""
@inline flux{T}(u, model::Quintic{T}) = u^2*u^3

"""
    speed(u::Real, model::Quintic)

Compute the speed f'(`u`) for `model`.
"""
@inline speed(u, model::Quintic) = 5*(u^2)^2


################################################################################

"""
    (::GodunovFlux){T}(uₗ::T, uᵣ::T, model::Quintic{T})

Compute Godunov's flux between `uₗ` and `uᵣ` for `model`.
"""
function (::GodunovFlux){T}(uₗ::T, uᵣ::T, model::Quintic{T})
  flux(uₗ, model)
end


"""
    (::EnergyConservativeFlux){T}(uₗ::T, uᵣ::T, model::Quintic{T})

Compute the energy (L₂ entropy) conservative flux between `uₗ` and `uᵣ` for `model`.
"""
function (flux::EnergyConservativeFlux)(uₗ::T, uᵣ::T, model::Quintic{T}) where T<:Real
    (uₗ^2*uₗ^3 + (uₗ^2)^2*uᵣ + uₗ^3*uᵣ^2 + uₗ^2*uᵣ^3 + uₗ*(uᵣ^2)^2 + uᵣ^2*uᵣ^3) / 6
end


################################################################################

@inline function add_numerical_fluxes_inner_loop2!(du, fluxes, u, balance_law::Quintic,
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

        for nx in Base.OneTo(Pp1)
            du[nx,ix] += (  (fluxes[ix  ]   - one_3 * Ru5l
                                            - one_3 * u[nx,ix] * Ru4l
                                            - one_6 * u[nx,ix]^2 * Ru3l
                                            - one_6 * Rul^3*Rul^2
                                            + one_6 * u[nx,ix] * (Rul^2)^2
                                            - one_6 * Ru2l * Rul^3
                                            + one_6 * u[nx,ix]^2 * Rul^3
                                            - one_6 * Ru3l * Rul^2
                                            + one_6 * u[nx,ix] * Ru3l * Rul
                                            - one_6 * Ru3l * Ru2l) * Rl[nx]
                          - (fluxes[ix+1]   - one_3 * Ru5r
                                            - one_3 * u[nx,ix] * Ru4r
                                            - one_6 * u[nx,ix]^2 * Ru3r
                                            - one_6 * Rur^3*Rur^2
                                            + one_6 * u[nx,ix] * (Rur^2)^2
                                            - one_6 * Ru2r * Rur^3
                                            + one_6 * u[nx,ix]^2 * Rur^3
                                            - one_6 * Ru3r * Rur^2
                                            + one_6 * u[nx,ix] * Ru3r * Rur
                                            - one_6 * Ru3r * Ru2r) * Rr[nx]
                            ) * jacx / basis.weights[nx]
        end
    end

    nothing
end
