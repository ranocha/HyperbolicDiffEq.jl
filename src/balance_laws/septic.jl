
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

function Base.show{T}(io::IO, model::Septic{T})
  print(io, "Scalar conservation law {T=", T, "} with flux f(u) = u^7")
end


"""
    flux{T}(u, model::Septic{T,1})

Compute the flux of `u` for `model`.
"""
@inline flux{T}(u, model::Septic{T}) = (u^2)^2*u^3

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
function (::GodunovFlux){T}(uₗ::T, uᵣ::T, model::Septic{T})
  flux(uₗ, model)
end


"""
    (::EnergyConservativeFlux){T}(uₗ::T, uᵣ::T, model::Septic{T})

Compute the energy (L₂ entropy) conservative flux between `uₗ` and `uᵣ` for `model`.
"""
function (flux::EnergyConservativeFlux)(uₗ::T, uᵣ::T, model::Septic{T}) where T<:Real
    ((uₗ^2)^2*uₗ^3 + (uₗ^3)^2*uᵣ + uₗ^2*uₗ^3*uᵣ^2 + (uₗ^2)^2*uᵣ^3 + uₗ^3*(uᵣ^2)^2
    + uₗ^2*uᵣ^2*uᵣ^3 + uₗ*(uᵣ^3)^2 + (uᵣ^2)^2*uᵣ^3) / 8
end


################################################################################

@inline function add_numerical_fluxes_inner_loop2!(du, fluxes, u, balance_law::Septic,
                                                    Nx, basis::GaussLegendre,
                                                    jacx, parallel)
    Pp1 = length(basis.nodes)
    Rl = reshape(interpolation_matrix(-1, basis), Pp1)
    Rr = reshape(interpolation_matrix(+1, basis), Pp1)
    utmp = zeros(eltype(u), size(u,1))
    one_4 = 1 / 4
    one_8 = 1 / 8

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
            utmp[nx] = (u[nx,ix]^3)^2
        end
        Ru6l = dot(Rl, utmp)
        Ru6r = dot(Rr, utmp)
        for nx in 1:Pp1
            utmp[nx] = u[nx,ix]^3*(u[nx,ix]^2)^2
        end
        Ru7l = dot(Rl, utmp)
        Ru7r = dot(Rr, utmp)

        for nx in 1:Pp1
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
