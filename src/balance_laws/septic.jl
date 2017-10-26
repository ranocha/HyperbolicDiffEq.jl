
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

function show{T}(io::IO, model::Septic{T})
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
