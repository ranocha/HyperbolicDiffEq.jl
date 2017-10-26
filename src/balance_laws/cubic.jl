
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

function show{T}(io::IO, model::Cubic{T})
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
                                                    Nx, basis::GaussLegendre,
                                                    jacx, parallel)
    Pp1 = length(basis.nodes)
    Rl = reshape(interpolation_matrix(-1, basis), Pp1)
    Rr = reshape(interpolation_matrix(+1, basis), Pp1)
    utmp = zeros(eltype(u), size(u,1))
    one_2 = 1 / 2
    one_4 = 1 / 4

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
