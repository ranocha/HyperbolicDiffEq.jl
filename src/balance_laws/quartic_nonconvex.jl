
"""
    QuarticNonconvex{T}

The scalar conservation law
``
  \\partial_t u + \\partial_x ( u(t,x)^4 - 10 u(t,x)^2 + 3 u(t,x) ) = 0
``
in one space dimensions using `T` as scalar type.
"""
struct QuarticNonconvex{T} <: ScalarBalanceLaw{T,1} end

function QuarticNonconvex(T=Float64)
  QuarticNonconvex{T}()
end

function Base.show(io::IO, model::QuarticNonconvex{T}) where {T}
  print(io, "Scalar conservation law {T=", T, "} with flux f(u) = u^4 - 10 u^2 + 3 u")
end


"""
    flux(u, model::QuarticNonconvex)

Compute the flux of `u` for `model`.
"""
@inline flux(u, model::QuarticNonconvex) = u^2 * (u^2 - 10) + 3u

"""
    speed(u::Real, model::QuarticNonconvex)

Compute the speed f'(`u`) for `model`.
"""
@inline speed(u, model::QuarticNonconvex) = 4*u^3 - 20*u + 3


function min_max_speed(ul::T, ur::T, model::QuarticNonconvex{T}) where T<:Real
  min_speed = zero(T)
  max_speed = zero(T)

  if ur > 1.2909944487358056284 && -2.5819888974716112568 < ul < 1.2909944487358056284
    min_speed = T(-14.213259316477408379)
  elseif (ur > 1.2909944487358056284 && (ul >= 1.2909944487358056284 || ul <= -2.5819888974716112568)) || (ur <= 1.2909944487358056284 && 2 * ul + ur + sqrt(20 - 3 * ur^2) <= 0) || 1.2909944487358056284 + ur <= 0
    min_speed = 3 + 4 * ul * (-5 + ul^2)
  else
    min_speed = 3 + 4 * ur * (-5 + ur^2)
  end

  if (1.2909944487358056284 + ul <= 0 && 1.2909944487358056284 + ur > 0 && ur < 1.2909944487358056284) || (1.2909944487358056284 + ul < 0 && 1.2909944487358056284 <= ur < 2.5819888974716112568)
    max_speed = T(20.213259316477408379)
  elseif (1.2909944487358056284 + ul == 0 && ur == 2.5819888974716112568) || (1.2909944487358056284 + ul >= 0 && (ur ≈ 1.2909944487358056284 || (sqrt(20 - 3 * ur^2) >= 2 * ul + ur && 1.2909944487358056284 < ur < 2.5819888974716112568))) || (ur < 1.2909944487358056284 && 1.2909944487358056284 + ul > 0)
    max_speed = 3 + 4 * ul * (-5 + ul^2)
  else
    max_speed = 3 + 4 * ur * (-5 + ur^2)
  end

  min_speed, max_speed
end

function max_abs_speed(ul::T, ur::T, model::QuarticNonconvex{T}) where T<:Real
  min_speed, max_speed = min_max_speed(ul, ur, model)
  max(abs(min_speed), abs(max_speed))
end


################################################################################


"""
    (::EnergyConservativeFlux)(uₗ::T, uᵣ::T, model::QuarticNonconvex{T}) where {T}

Compute the energy (L₂ entropy) conservative flux between `uₗ` and `uᵣ` for `model`.
"""
function (flux::EnergyConservativeFlux)(ul::T, ur::T, model::QuarticNonconvex{T}) where T<:Real
  (ul^4 + ur*ul^3 + ur^2*ul^2 + ur^3*ul + ur^4) / 5 - 10 * (ul^2 + ul*ur + ur^2) / 3 + 3 * (ul + ur) / 2
end
