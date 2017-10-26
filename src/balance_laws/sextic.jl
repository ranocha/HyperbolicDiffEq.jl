
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
