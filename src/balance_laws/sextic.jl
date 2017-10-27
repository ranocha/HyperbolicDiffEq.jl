
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

function show{T}(io::IO, model::Sextic{T})
  print(io, "Scalar conservation law {T=", T, "} with flux f(u) = u^6")
end


"""
    flux{T}(u, model::Sextic{T,1})

Compute the flux of `u` for `model`.
"""
@inline flux{T}(u, model::Sextic{T}) = (u^3)^2

"""
    speed(u::Real, model::Sextic)

Compute the speed f'(`u`) for `model`.
"""
@inline speed(u, model::Sextic) = 6*u^2*u^3


################################################################################

"""
    (::GodunovFlux){T}(uₗ::T, uᵣ::T, model::Sextic{T})

Compute Godunov's flux between `uₗ` and `uᵣ` for `model`.
"""
function (::GodunovFlux){T}(uₗ::T, uᵣ::T, model::Sextic{T})
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
    (::EnergyConservativeFlux){T}(uₗ::T, uᵣ::T, model::Sextic{T})

Compute the energy (L₂ entropy) conservative flux between `uₗ` and `uᵣ` for `model`.
"""
function (flux::EnergyConservativeFlux)(uₗ::T, uᵣ::T, model::Sextic{T}) where T<:Real
    ((uₗ^3)^2 + uₗ^2*uₗ^3*uᵣ + (uₗ^2)^2*uᵣ^2 + uₗ^3*uᵣ^3 + uₗ^2*(uᵣ^2)^2 + uₗ*uᵣ^2*uᵣ^3 + (uᵣ^3)^2) / 7
end


################################################################################

@inline function add_numerical_fluxes_inner_loop2!(du, fluxes, u, balance_law::Sextic,
                                                    Nx, basis::GaussLegendre,
                                                    jacx, parallel)
    Pp1 = length(basis.nodes)
    Rl = reshape(interpolation_matrix(-1, basis), Pp1)
    Rr = reshape(interpolation_matrix(+1, basis), Pp1)
    utmp = zeros(eltype(u), size(u,1))
    one_7 = 1 / 7
    two_7 = 2 / 7

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
