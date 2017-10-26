
"""
    Octic{T}

The scalar conservation law with flux
``
  \\partial_t u + \\partial_x u^8 = 0
``
in one space dimensions using `T` as scalar type.
"""
struct Octic{T} <: ScalarBalanceLaw{T,1} end

function Octic(T=Float64)
  Octic{T}()
end

function show{T}(io::IO, model::Octic{T})
  print(io, "Scalar conservation law {T=", T, "} with flux f(u) = u^8")
end


"""
    flux{T}(u, model::Octic{T,1})

Compute the flux of `u` for `model`.
"""
@inline flux{T}(u, model::Octic{T}) = ((u^2)^2)^2

"""
    speed(u::Real, model::Octic)

Compute the speed f'(`u`) for `model`.
"""
@inline speed(u, model::Octic) = 8*(u^2)^2*u^3


################################################################################

"""
    (::GodunovFlux){T}(uₗ::T, uᵣ::T, model::Octic{T})

Compute Godunov's flux between `uₗ` and `uᵣ` for `model`.
"""
function (::GodunovFlux){T}(uₗ::T, uᵣ::T, model::Octic{T})
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
    (::EnergyConservativeFlux){T}(uₗ::T, uᵣ::T, model::Octic{T})

Compute the energy (L₂ entropy) conservative flux between `uₗ` and `uᵣ` for `model`.
"""
function (flux::EnergyConservativeFlux)(uₗ::T, uᵣ::T, model::Octic{T}) where T<:Real
    ( ((uₗ^2)^2)^2 + (uₗ^2)^2*uₗ^3*uᵣ + (uₗ^3)^2*uᵣ^2 + uₗ^2*uₗ^3*uᵣ^3 + (uₗ^2)^2*(uᵣ^2)^2
    + uₗ^3*uᵣ^2*uᵣ^3 + uₗ^2*(uᵣ^3)^2 + uₗ*(uᵣ^2)^2*uᵣ^3 + ((uᵣ^2)^2)^2 ) / 9
end
