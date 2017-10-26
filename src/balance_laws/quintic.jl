
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

function show{T}(io::IO, model::Quintic{T})
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
