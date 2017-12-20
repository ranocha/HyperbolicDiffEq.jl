
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

function (flux::EnergyConservativeFlux)(uₗ::T, uᵣ::T, model::ConstantLinearAdvection{T,1}) where T<:Real
    (uₗ + uᵣ) / 2
end


################################################################################

@inline function add_numerical_fluxes_inner_loop2!(du, fluxes, u, balance_law::ConstantLinearAdvection,
                                                    Nx, basis::GaussLegendre,
                                                    jacx, parallel)
    Pp1 = length(basis.nodes)
    Rl = basis.interp_left
    Rr = basis.interp_right

    # add numerical fluxes
    @inbounds for ix in Base.OneTo(Nx)
        Rul = zero(eltype(u))
        Rur = zero(eltype(u))
        for nx in Base.OneTo(Pp1)
            Rul += Rl[nx]*u[nx,ix]
            Rur += Rr[nx]*u[nx,ix]
        end

        for nx in Base.OneTo(Pp1)
            du[nx,ix] += ((fluxes[ix] - Rul) * Rl[nx] - (fluxes[ix+1] - Rur) * Rr[nx]
                            ) * jacx / basis.weights[nx]
        end
    end

    nothing
end


################################################################################

"""
    ConstantLinearAdvectionRiemannSolution{T,T1}

The solution of a Riemann problem `prob` for the linear advection equation with
constant coefficients.
"""
struct ConstantLinearAdvectionRiemannSolution{T,T1} <: ScalarRiemannSolution
    prob::RiemannProblem{ConstantLinearAdvection{T,1},T,T1}
    σ::T
end


"""
    minmax_speeds(sol::ConstantLinearAdvectionRiemannSolution)

Return the minimal and maximal speeds `σ⁻, σ⁺` that appear in the solution `sol`.
"""
function minmax_speeds(sol::ConstantLinearAdvectionRiemannSolution)
    sol.σ, sol.σ
end


"""
    (sol::ConstantLinearAdvectionRiemannSolution)(ξ::Real)

Evaluate the solution `sol` at the value `ξ` of the self-similarity variable
`ξ = (x - x₀) / (t - t₀)`.
"""
function (sol::ConstantLinearAdvectionRiemannSolution)(ξ::Real)
    @unpack σ = sol
    @unpack uₗ, uᵣ = sol.prob

    if ξ < σ
        uₗ
    else
        uᵣ
    end
end


"""
    (sol::ConstantLinearAdvectionRiemannSolution)(t::Real, x::Real)

Evaluate the solution `sol` at the time and space coordinates `t` and `x`.
"""
function (sol::ConstantLinearAdvectionRiemannSolution)(t::Real, x::Real)
    @unpack x₀, t₀ = sol.prob

    sol((x-x₀)/(t-t₀))
end


"""
    solve{T,T1}(prob::RiemannProblem{ConstantLinearAdvection{T,1},T,T1})

Compute the solution of the Riemann prolem `prob`.
"""
function solve{T,T1}(prob::RiemannProblem{ConstantLinearAdvection{T,1},T,T1})
    ConstantLinearAdvectionRiemannSolution(prob, one(T))
end
