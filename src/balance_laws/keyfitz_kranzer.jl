
doc"
    KeyfitzKranzer{T<:Real}

The Keyfitz-Kranzer system
$$
    \partial_t u_1 + \partial_x (u_1^3 - u_2) = 0,
    \partial_t u_2 + \partial_x (\frac{1}{3}u_2^3 - u_1) = 0,
$$
using `T` as scalar type.
"
struct KeyfitzKranzer{T} <: AbstractBalanceLaw{1} end

KeyfitzKranzer(::Type{T}=Float64) where {T} = KeyfitzKranzer{T}()

function Base.show(io::IO, model::KeyfitzKranzer{T}) where {T}
  print(io, "Keyfitz-Kranzer system {T=", T, "}")
end


"""
Conserved variables (u1, u2) of the Keyfitz-Kranzer system.
"""
struct KeyfitzKranzerVar{T} <: FieldVector{2,T}
  u1::T
  u2::T
end

function (::Type{KeyfitzKranzerVar{T}}){T}(val::Real)
  KeyfitzKranzerVar{T}(val, val)
end

function similar_type{T}(::KeyfitzKranzerVar{T})
  KeyfitzKranzerVar{T}
end

@inline variables(model::KeyfitzKranzer{T}) where {T} = KeyfitzKranzerVar{T}

@inline function entropy(u::KeyfitzKranzerVar, model::KeyfitzKranzer)
    @unpack u1, u2 = u

    exp(u1^2/2 - u2)
end

@inline function entropy_variables(u::KeyfitzKranzerVar, model::KeyfitzKranzer)
    @unpack u1, u2 = u
    U = entropy(u, model)

    SVector(u1*U, -U)
end

@inline function flux_potential(u::KeyfitzKranzerVar, model::KeyfitzKranzer)
    @unpack u1, u2 = u
    U = entropy(u, model)

    (2*u1^3/3 - u1*u2) * U
end

@inline function flux(u::KeyfitzKranzerVar, model::KeyfitzKranzer)
  @unpack u1, u2 = u

  SVector(u1^2 - u2, u1^3/3 - u1)
end

@inline function max_abs_speed(u::KeyfitzKranzerVar, model::KeyfitzKranzer)
  @unpack u1, u2 = u

  abs(u1) + 1
end

#= TODO: Estimate the maximal speed in the solution of the Riemann problem with
left and right states `uₗ`, `uᵣ` for `model`.
@inline function max_abs_speed(uₗ::KeyfitzKranzerVar, uᵣ::KeyfitzKranzerVar, model::KeyfitzKranzer)
    λ₋, λ₊ = min_max_speed(uₗ, uᵣ, model)

    max(abs(λ₋), abs(λ₊))
end

@inline function min_max_speed(uₗ::KeyfitzKranzerVar, uᵣ::KeyfitzKranzerVar, model::KeyfitzKranzer)


    λ₋, λ₊
end
=#


"""
    IntegralQuantitiesKeyfitzKranzer{T<:Real}

Some integrated quantities of interest for the Keyfitz Kranzer system. Can be
used in callbacks of DifferentialEquations.jl.
"""
struct IntegralQuantitiesKeyfitzKranzer{T<:Real} <: FieldVector{3,T}
    u1::T
    u2::T
    entropy::T
end

function IntegralQuantitiesKeyfitzKranzer(u::KeyfitzKranzerVar{T}, model::KeyfitzKranzer) where {T}
    @unpack u1, u2 = u

    IntegralQuantitiesKeyfitzKranzer{T}(u1, u2, entropy(u, model))
end


function (fvol::EnergyConservativeFlux)(uₗ::KeyfitzKranzerVar, uᵣ::KeyfitzKranzerVar,
                                        model::KeyfitzKranzer)
    u1₋, u2₋ = uₗ
    u1₊, u2₊ = uᵣ
    U₋ = entropy(uₗ, model)
    U₊ = entropy(uᵣ, model)

    u1 = (u1₋ + u1₊) / 2
    U = (U₋ + U₊) / 2
    Ulog = logmean(U₋, U₊)

    fnum_u1 = ( 5*u1₊^2 + 2*u1₋*u1₊ + 5*u1₋^2) / 12 - (u2₋ + u2₊) / 2
    fnum_u2 = u1 * fnum_u1 - (u1₋^3 + u1₊^3) / 3 + (u1₋*u2₋ + u1₊*u2₊) / 2 - u1 * U / Ulog

    SVector(fnum_u1, fnum_u2)
end



@recipe function f{Ξ,T}(ξu::Tuple{Ξ,Vector{KeyfitzKranzerVar{T}}})
  ξ, u = ξu

  h  = mappedarray(u->u.u1, u)
  hv = mappedarray(u->u.u2, u)

  size --> (1000, 400)
  layout --> (1,2)
  legend --> false

  @series begin
    subplot := 1
    ylabel --> L"u_1"
    label  --> L"u_1"
    ξ, h
  end

  @series begin
    subplot := 2
    ylabel --> L"u_2"
    label  --> L"u_2"
    ξ, hv
  end
end
