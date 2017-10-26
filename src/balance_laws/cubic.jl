
"""
    Cubic{T}

The scalar conservation law with cubic flux
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
  print(io, "Scalar conservation law {T=", T, "} with cubic flux f(u) = u^3")
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
