
"""
    ConstantLinearAdvection{T,Dim}

Linear advection equation with constant coefficients
``
  \\partial_t u + \\sum_{i=1}^{\\mathrm{Dim}} \\partial_{x_i} u = 0
``
in `Dim` space dimensions using `T` as scalar type.
"""
struct ConstantLinearAdvection{T,Dim} <: ScalarBalanceLaw{T,Dim} end

function ConstantLinearAdvection(T=Float64, Dim=1)
  ConstantLinearAdvection{T,Dim}()
end

function Base.show{T,Dim}(io::IO, model::ConstantLinearAdvection{T,Dim})
  print(io, "Linear advection equation with constant coefficient {T=", T, ", Dim=", Dim, "}",
            " with flux f(u) = u * (1,...,1)")
end


"""
    flux{T}(u, model::ConstantLinearAdvection{T,1})

Compute the flux of `u` for `model`.
"""
@inline flux{T}(u, model::ConstantLinearAdvection{T,1}) = u

"""
    speed(u::Real, model::ConstantLinearAdvection)

Compute the speed f'(`u`) for `model`.
"""
@inline speed(u, model::ConstantLinearAdvection) = one(u)


@inline flux{T}(u, model::ConstantLinearAdvection{T,2}, direction) = u


################################################################################

"""
    (::GodunovFlux){T}(uₗ::T, uᵣ::T, model::ConstantLinearAdvection{T,1})

Compute Godunov's flux between `uₗ` and `uᵣ` for `model`.
"""
function (::GodunovFlux){T}(uₗ::T, uᵣ::T, model::ConstantLinearAdvection{T,1})
    uₗ
end
