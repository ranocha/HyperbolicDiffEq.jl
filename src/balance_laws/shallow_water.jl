
"""
    ShallowWater{T<:Real,Dim}

The shallow water equations with gravitational constant `g`
in `Dim` dimensions using `T` as scalar tpe.
"""
struct ShallowWater{T<:Real,Dim} <: AbstractBalanceLaw{1}
  g::T
end

function ShallowWater(g=1., Dim=1)
  ShallowWater{typeof(g), Dim}(g)
end

function show{T,Dim}(io::IO, model::ShallowWater{T,Dim})
  print(io, "Shallow water equations {T=", T, ", Dim=", Dim, "}")
end


"""
Conserved variables (h, hv) of the shallow water system in one space dimension.
"""
struct ShallowWaterVar1D{T} <: FieldVector{2,T}
  h ::T
  hv::T
end

function (::Type{ShallowWaterVar1D{T}}){T}(val::Real)
  ShallowWaterVar1D{T}(val, val)
end

function similar_type{T}(::ShallowWaterVar1D{T})
  ShallowWaterVar1D{T}
end

@inline variables{T}(model::ShallowWater{T,1}) = ShallowWaterVar1D{T}

@inline function primitive_variables(u::ShallowWaterVar1D)
  @unpack h, hv = u
  if h ≈ 0
    v = zero(h)
  else
    v = hv/h
  end
  h, v
end

@inline function primitive_variables(u::ShallowWaterVar1D, model::ShallowWater)
  primitive_variables(u)
end

@inline function flux{T}(u::ShallowWaterVar1D{T}, model::ShallowWater{T,1})
  h, v = primitive_variables(u)
  @unpack g = model

  SVector{2,T}(h*v, h*v^2 + g*h^2/2)
end

@inline function max_speed{T}(u::ShallowWaterVar1D, model::ShallowWater{T,1})
  h, v = primitive_variables(u)
  @unpack g = model

  abs(v) + sqrt(g*h)
end


################################################################################

struct ShallowWaterRiemannSolution{T,T1} <: AbstractRiemannSolution
  prob::RiemannProblem{ShallowWater{T,1},ShallowWaterVar1D{T},T1}
  uₘ::ShallowWaterVar1D{T} # middle state
  σ₁⁻::T # slow speed of first family
  σ₁⁺::T # fast speed of first family
  σ₂⁻::T # slow speed of second family
  σ₂⁺::T # fast speed of second family

  function ShallowWaterRiemannSolution{T,T1}(prob::RiemannProblem{ShallowWater{T,1},ShallowWaterVar1D{T},T1}) where {T,T1}
    uₘ, σ₁⁻, σ₁⁺, σ₂⁻, σ₂⁺ = compute_state_and_speeds(prob.uₗ, prob.uᵣ, prob.model)
    new(prob, uₘ, σ₁⁻, σ₁⁺, σ₂⁻, σ₂⁺)
  end
end

function solve{T,T1}(prob::RiemannProblem{ShallowWater{T,1},ShallowWaterVar1D{T},T1})
  ShallowWaterRiemannSolution{T,T1}(prob)
end

function minmax_speeds(sol::ShallowWaterRiemannSolution)
  sol.σ₁⁻, sol.σ₂⁺
end

function (sol::ShallowWaterRiemannSolution)(ξ::Real)
  @unpack uₘ, σ₁⁻, σ₁⁺, σ₂⁻, σ₂⁺ = sol
  @unpack uₗ, uᵣ = sol.prob

  if ξ <= σ₁⁻
    uₗ
  elseif ξ <= σ₁⁺
    rarefaction_wave_1(ξ, uₗ)
  elseif ξ <= σ₂⁻
    uₘ
  elseif ξ <= σ₂⁺
    rarefaction_wave_2(ξ, uᵣ)
  else # σ₂⁺ < ξ
    uᵣ
  end
end

function rarefaction_wave_1{T}(ξ, uₗ::ShallowWaterVar1D{T})
  hₗ, vₗ = primitive_variables(uₗ)
  h  = (vₗ + 2*sqrt(hₗ) -  ξ)^2 / 9
  hv = (vₗ + 2*sqrt(hₗ) + 2ξ)*h / 3
  ShallowWaterVar1D{T}(h, hv)
end

function rarefaction_wave_2{T}(ξ, uₗ::ShallowWaterVar1D{T})
  hₗ, vₗ = primitive_variables(uₗ)
  h  = (vₗ - 2*sqrt(hₗ) -  ξ)^2 / 9
  hv = (vₗ - 2*sqrt(hₗ) + 2ξ)*h / 3
  ShallowWaterVar1D{T}(h, hv)
end


function (sol::ShallowWaterRiemannSolution)(t::Real, x::Real)
  @unpack x₀, t₀ = sol.prob

  sol((x-x₀)/(t-t₀))
end


function compute_state_and_speeds{T}(uₗ::ShallowWaterVar1D{T}, uᵣ::ShallowWaterVar1D{T}, model::ShallowWater)
  hₗ, vₗ = primitive_variables(uₗ)
  hᵣ, vᵣ = primitive_variables(uᵣ)

  if v_S1_SWE(hᵣ,hₗ,vₗ) <= vᵣ && vᵣ <= v_R2_SWE(hᵣ,hₗ,vₗ)
    # Case I: wave1 = shock, wave2 = rarefaction
    hₘ = Roots.fzero(h -> vᵣ-v_R2_SWE(hᵣ,h,v_S1_SWE(h,hₗ,vₗ)), hₗ, hᵣ)
     # this is better than vₘ = v_S1_SWE(hₘ,hₗ,vₗ) for hₗ≈0
    vₘ = v_R2_SWE(hₘ,hᵣ,vᵣ)
    uₘ = ShallowWaterVar1D{T}(hₘ, hₘ*vₘ)
    # this is better than σ₁⁻ = σ₁⁺ = vₗ - sqrt(hₘ+hₘ^2/hₗ) / sqrt(2) for hₗ≈0
    σ₁⁻ = σ₁⁺ = vₘ - sqrt(hₗ+hₗ^2/hₘ) / sqrt(2)
    σ₂⁻ = vₘ + sqrt(hₘ)
    σ₂⁺ = vᵣ + sqrt(hᵣ)
  elseif v_R1_SWE(hᵣ,hₗ,vₗ) <= vᵣ && v_R2_SWE(hᵣ,hₗ,vₗ) <= vᵣ
    # Case II: wave1 = rarefaction, wave2 = rarefaction
    if 2*(sqrt(hᵣ)+sqrt(hₗ)) >= vᵣ-vₗ
      # no vacuum
      hₘ = (2*(sqrt(hᵣ)+sqrt(hₗ)) -(vᵣ-vₗ))^2 / 16
      vₘ = v_R1_SWE(hₘ,hₗ,vₗ)
      uₘ = ShallowWaterVar1D{T}(hₘ,hₘ*vₘ)
      σ₁⁻ = vₗ - sqrt(hₗ)
      σ₁⁺ = vₘ - sqrt(hₘ)
      σ₂⁻ = vₘ + sqrt(hₘ)
      σ₂⁺ = vᵣ + sqrt(hᵣ)
    else
      # vacuum
      uₘ = ShallowWaterVar1D{T}(0,0)
      σ₁⁻ = vₗ - sqrt(hₗ)
      σ₁⁺ = vₗ + 2*sqrt(hₗ)
      σ₂⁻ = vᵣ - 2*sqrt(hᵣ)
      σ₂⁺ = vᵣ + sqrt(hᵣ)
    end
  elseif v_S2_SWE(hᵣ,hₗ,vₗ) <= vᵣ && vᵣ <= v_R1_SWE(hᵣ,hₗ,vₗ)
    # Case III: wave1 = rarefaction, wave2 = shock
    hₘ = Roots.fzero(h->vᵣ-v_S2_SWE(hᵣ,h,v_R1_SWE(h,hₗ,vₗ)), hᵣ, hₗ)
    vₘ = v_R1_SWE(hₘ,hₗ,vₗ)
    uₘ = ShallowWaterVar1D{T}(hₘ,hₘ*vₘ)
    σ₁⁻ = vₗ - sqrt(hₗ)
    σ₁⁺ = vₘ - sqrt(hₘ)
    σ₂⁻ = σ₂⁺ = vₘ + sqrt(hᵣ+hᵣ^2/hₘ) / sqrt(2)
  else #vᵣ <= v_S1_SWE(hᵣ,hₗ,vₗ) && vᵣ <= v_S2_SWE(hᵣ,hₗ,vₗ)
    # Case IV: wave1 = shock, wave2 = shock
    # TODO: Chose some interval?
    hₘ = Roots.fzero(h->vᵣ-v_S2_SWE(hᵣ,h,v_S1_SWE(h,hₗ,vₗ)), max(hₗ,hᵣ))
    vₘ = v_S1_SWE(hₘ,hₗ,vₗ)
    uₘ = ShallowWaterVar1D{T}(hₘ,hₘ*vₘ)
    σ₁⁻ = σ₁⁺ = vₗ - hₘ * sqrt(1/hₘ+1/hₗ) / sqrt(2)
    σ₂⁻ = σ₂⁺ = vₘ + hᵣ * sqrt(1/hᵣ+1/hₘ) / sqrt(2)
  end

  uₘ, σ₁⁻, σ₁⁺, σ₂⁻, σ₂⁺
end


# wave curves parameterised by the height h
v_R1_SWE(h, hₗ, vₗ) = vₗ - 2*(sqrt(h) - sqrt(hₗ))
v_R2_SWE(h, hₗ, vₗ) = vₗ + 2*(sqrt(h) - sqrt(hₗ))
v_S1_SWE(h, hₗ, vₗ) = vₗ - (h-hₗ) * sqrt(1/h+1/hₗ) / sqrt(2)
v_S2_SWE(h, hₗ, vₗ) = vₗ + (h-hₗ) * sqrt(1/h+1/hₗ) / sqrt(2)



@recipe function f(sol::ShallowWaterRiemannSolution)
  xguide --> L"\xi"

  σ⁻, σ⁺ =  minmax_speeds(sol)
  if σ⁻ ≈ σ⁺
    σ = one(σ⁻)
  else
    σ = min(abs(σ⁻ - σ⁺), max(abs(σ⁻), abs(σ⁺))/10)
  end
  σ⁻ = σ⁻ - σ
  σ⁺ = σ⁺ + σ

  ξ = linspace(σ⁻, σ⁺, 10^3)
  u = sol.(ξ)

  ((ξ, u),)
end


@recipe function f{Ξ,T}(ξu::Tuple{Ξ,Vector{ShallowWaterVar1D{T}}})
  ξ = ξu[1]
  u = ξu[2]

  h  = mappedarray(u->u.h, u)
  hv = mappedarray(u->u.hv, u)

  size --> (1000, 400)
  layout --> (1,2)
  legend --> false

  @series begin
    subplot := 1
    ylabel --> L"h"
    label  --> L"h"
    ξ, h
  end

  @series begin
    subplot := 2
    ylabel --> L"hv"
    label  --> L"hv"
    ξ, hv
  end
end
