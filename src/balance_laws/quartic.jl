
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

function Base.show{T}(io::IO, model::Quartic{T})
  print(io, "Scalar conservation law {T=", T, "} with flux f(u) = u^4")
end


"""
    flux{T}(u, model::Quartic{T,1})

Compute the flux of `u` for `model`.
"""
@inline flux{T}(u, model::Quartic{T}) = (u^2)^2

"""
    speed(u::Real, model::Quartic)

Compute the speed f'(`u`) for `model`.
"""
@inline speed(u, model::Quartic) = 4u^3


################################################################################

"""
    (::GodunovFlux){T}(uₗ::T, uᵣ::T, model::Quartic{T})

Compute Godunov's flux between `uₗ` and `uᵣ` for `model`.
"""
function (::GodunovFlux){T}(uₗ::T, uᵣ::T, model::Quartic{T})
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
    (::EnergyConservativeFlux){T}(uₗ::T, uᵣ::T, model::Quartic{T})

Compute the energy (L₂ entropy) conservative flux between `uₗ` and `uᵣ` for `model`.
"""
function (flux::EnergyConservativeFlux)(uₗ::T, uᵣ::T, model::Quartic{T}) where T<:Real
    ((uₗ^2)^2 + uₗ^3*uᵣ + uₗ^2*uᵣ^2 + uₗ*uᵣ^3 + (uᵣ^2)^2) / 5
end


################################################################################

@inline function add_numerical_fluxes_inner_loop2!(du, fluxes, u, balance_law::Quartic,
                                                    Nx, basis::GaussLegendre,
                                                    jacx, parallel)
    Pp1 = length(basis.nodes)
    Rl = reshape(interpolation_matrix(-1, basis), Pp1)
    Rr = reshape(interpolation_matrix(+1, basis), Pp1)
    utmp = zeros(eltype(u), size(u,1))
    one_5 = 1 / 5
    two_5 = 2 / 5

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
