
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

function show{T}(io::IO, model::Quartic{T})
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
