
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
                                                    Nx, basis::GaussLegendre,
                                                    jacx, parallel)
    Pp1 = length(basis.nodes)
    Rl = reshape(interpolation_matrix(-1, basis), Pp1)
    Rr = reshape(interpolation_matrix(+1, basis), Pp1)
    utmp = zeros(eltype(u), size(u,1))
    one_3 = 1 / 3
    one_6 = 1 / 6

    # add numerical fluxes
    @inbounds for ix in Base.OneTo(Nx)
        for nx in 1:Pp1
            utmp[nx] = u[nx,ix]
        end
        Rul = dot(Rl, utmp)
        Rur = dot(Rr, utmp)
        for nx in 1:Pp1
            utmp[nx] = u[nx,ix]^2
        end
        Ru2l = dot(Rl, utmp)
        Ru2r = dot(Rr, utmp)
        for nx in 1:Pp1
            utmp[nx] = u[nx,ix]^3
        end
        Ru3l = dot(Rl, utmp)
        Ru3r = dot(Rr, utmp)
        for nx in 1:Pp1
            utmp[nx] = (u[nx,ix]^2)^2
        end
        Ru4l = dot(Rl, utmp)
        Ru4r = dot(Rr, utmp)
        for nx in 1:Pp1
            utmp[nx] = u[nx,ix]^3*u[nx,ix]^2
        end
        Ru5l = dot(Rl, utmp)
        Ru5r = dot(Rr, utmp)

        for nx in 1:Pp1
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
