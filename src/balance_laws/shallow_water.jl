
"""
    ShallowWater{T<:Real,Dim}

The shallow water equations with gravitational constant `g`
in `Dim` dimensions using `T` as scalar type.
"""
struct ShallowWater{T,Dim} <: AbstractBalanceLaw{1}
  g::T
end

function ShallowWater(g=1., Dim=1)
  ShallowWater{typeof(g), Dim}(g)
end

function show{T,Dim}(io::IO, model::ShallowWater{T,Dim})
  print(io, "Shallow water equations with g=", model.g, " {T=", T, ", Dim=", Dim, "}")
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

@inline function entropy_variables(u::ShallowWaterVar1D, model::ShallowWater)
    entropy_variables(primitive_variables(u, model)..., model)
end

@inline function entropy_variables(h, v, model::ShallowWater)
    @unpack g = model
    SVector(g*h - v^2/2, v)
end

@inline function flux{T}(u::ShallowWaterVar1D{T}, model::ShallowWater{T,1})
  h, v = primitive_variables(u)
  @unpack g = model

  SVector{2,T}(h*v, h*v^2 + g*h^2/2)
end

@inline function max_abs_speed{T}(u::ShallowWaterVar1D, model::ShallowWater{T,1})
  h, v = primitive_variables(u)
  @unpack g = model

  abs(v) + sqrt(g*h)
end

"""
    max_abs_speed(uₗ::ShallowWaterVar1D, uᵣ::ShallowWaterVar1D, model::ShallowWater)

Estimate the maximal speed in the solution of the Riemann problem with left and
right states `uₗ`, `uᵣ` for `model` via the two rarefaction approximation, see
Guermond and Popov (2016) Fast Estimation from above of the maximum wave speed in
the Riemann roblem for the Euler equations, Journal of Computational Physics 321,
pp. 908-926, and
Chen and Shu (2017) Entropy stable high order discontinuous Galerkin methods with
suitable quadrature rules for hyperbolic conservation laws, Journal of Computational
Physics 345, pp. 427-461.
"""
@inline function max_abs_speed(uₗ::ShallowWaterVar1D, uᵣ::ShallowWaterVar1D, model::ShallowWater)
    λ₋, λ₊ = min_max_speed(uₗ, uᵣ, model)

    max(abs(λ₋), abs(λ₊))
end

"""
    min_max_speed(uₗ::ShallowWaterVar1D, uᵣ::ShallowWaterVar1D, model::ShallowWater)

Estimate the maximal speed in the solution of the Riemann problem with left and
right states `uₗ`, `uᵣ` for `model` via the two rarefaction approximation, see
Guermond and Popov (2016) Fast Estimation from above of the maximum wave speed in
the Riemann roblem for the Euler equations, Journal of Computational Physics 321,
pp. 908-926, and
Chen and Shu (2017) Entropy stable high order discontinuous Galerkin methods with
suitable quadrature rules for hyperbolic conservation laws, Journal of Computational
Physics 345, pp. 427-461.
"""
@inline function min_max_speed(uₗ::ShallowWaterVar1D, uᵣ::ShallowWaterVar1D, model::ShallowWater)
    @unpack g = model
    hₗ, vₗ = primitive_variables(uₗ, model)
    aₗ = sqrt(g * hₗ)
    hᵣ, vᵣ = primitive_variables(uᵣ, model)
    aᵣ = sqrt(g * hᵣ)

    hₘ = (2*(aᵣ+aₗ) - (vᵣ-vₗ))^2 / 16g
    if hₘ <= hₗ
        # approximated by rarefaction wave
        λ₋ = vₗ - aₗ
    else
        # approximated by shock wave
        λ₋ = vₗ - hₘ * sqrt( g*(hₘ+hₗ) / (2*hₘ*hₗ) )
    end
    if hₘ <= hᵣ
        # approximated by rarefaction wave
        λ₊ = vᵣ + aᵣ
    else
        # approximated by shock wave
        λ₊ = vᵣ + hₘ * sqrt( g*(hₘ+hᵣ) / (2*hₘ*hᵣ) )
    end

    λ₋, λ₊
end


"""
    (::EnergyConservativeFlux1Param)(uₗ::ShallowWaterVar1D, uᵣ::ShallowWaterVar1D, model::ShallowWater)

The one-parameter family of energy conservative numerical fluxes for the shallow
water equations, see
Ranocha (2017) Shallow water equations: Split-form, entropy stable, well-balanced,
and positivity-preserving numerical methods, GEM - International Journal on
Geomathematics 8(1), pp. 85-133.
"""
function (fvol::EnergyConservativeFlux1Param)(uₗ::ShallowWaterVar1D, uᵣ::ShallowWaterVar1D,
                                              model::ShallowWater)
    h₋, v₋ = primitive_variables(uₗ, model)
    h₊, v₊ = primitive_variables(uᵣ, model)
    @unpack g = model
    @unpack a₁ = fvol

    fnum_h  = ( (3-a₁)*(h₋*v₋ + h₊*v₊) + (1+a₁)*(h₊*v₋ + h₋*v₊) ) / 8
    fnum_hv = g * ( (1+a₁)*(h₊^2 + h₋^2) + (2-2a₁)*(h₊*h₋) ) / 8 +
                ( (3-a₁)*(h₊*v₊^2 + h₋*v₋^2) + (1+a₁)*(h₊*v₋^2 + h₋*v₊^2) ) / 16 +
                (h₊*v₊*v₋ + h₋*v₊*v₋) / 4

    SVector(fnum_h, fnum_hv)
end


doc"
    (diss::ScalarDissipation)(uₗ::ShallowWaterVar1D, uᵣ::ShallowWaterVar1D, model::ShallowWater)

The scalar dissipation operator $- \frac{\lambda}{2} R \cdot R^T \cdot  (w_r - w_l)$
for the shallow water equations, where $\partial_w u = R \cdot R^T$ is evaluated
at the arithmetic mean values $h = (h_- + h_+) / 2$, $v = (v_- + v_+) / 2$.
"
function (diss::ScalarDissipation)(uₗ::ShallowWaterVar1D, uᵣ::ShallowWaterVar1D, model::ShallowWater)
    h₋, v₋ = primitive_variables(uₗ, model)
    h₊, v₊ = primitive_variables(uᵣ, model)
    @unpack g = model

    w₋ = entropy_variables(h₋, v₋, model)
    w₊ = entropy_variables(h₊, v₊, model)

    h = (h₋ + h₊) / 2
    v = (v₋ + v₊) / 2
    du_dw = SArray{Tuple{2,2}}(1, v, v, g*h+v^2) / g

    λ = diss.max_abs_speed(uₗ, uᵣ, model)

    -λ/2 * du_dw * (w₊ - w₋)
end

doc"
    (diss::MatrixDissipation)(uₗ::ShallowWaterVar1D, uᵣ::ShallowWaterVar1D, model::ShallowWater)

The scalar dissipation operator $- \frac{1}{2} R \cdot |\Lambda| \cdot R^T \cdot  (w_r - w_l)$
for the shallow water equations, where $f'(u) = R \cdot \Lambda \cdot R^{-1}$ is
evaluated at the arithmetic mean values $h = (h_- + h_+) / 2$, $v = (v_- + v_+) / 2$.
"
function (diss::MatrixDissipation)(uₗ::ShallowWaterVar1D, uᵣ::ShallowWaterVar1D, model::ShallowWater)
    h₋, v₋ = primitive_variables(uₗ, model)
    h₊, v₊ = primitive_variables(uᵣ, model)
    @unpack g = model

    w₋ = entropy_variables(h₋, v₋, model)
    w₊ = entropy_variables(h₊, v₊, model)

    h = (h₋ + h₊) / 2
    v = (v₋ + v₊) / 2
    sqrt_gh = sqrt(g*h)
    R = SArray{Tuple{2,2}}(1, v-sqrt_gh, 1, v+sqrt_gh) / sqrt(2g)
    Λ = SArray{Tuple{2,2}}(abs(v-sqrt_gh), 0, 0, abs(v+sqrt_gh))

    -R * Λ * R' * (w₊ - w₋) / 2
end


"""
    (::SuliciuFlux)(uₗ::ShallowWaterVar1D, uᵣ::ShallowWaterVar1D, model::ShallowWater{T,1})

The Suliciu relaxation solver / numerical flux, see
Bouchut (2004) Nonlinear Stability of Finite Volume Methods for Hyperbolic
Conservation Laws and Well-Balanced Schemes for Sources.
"""
function (::SuliciuFlux)(uₗ::ShallowWaterVar1D, uᵣ::ShallowWaterVar1D, model::ShallowWater{T,1}) where T
    hl, vl = primitive_variables(uₗ, model)
    hr, vr = primitive_variables(uᵣ, model)
    @unpack g = model

    sqrt_ghl = sqrt(g*hl)
    sqrt_ghr = sqrt(g*hr)

    # compute speeds
    cl_hl = zero(T)
    cr_hr = zero(T)
    if hl <= hr && 0 < hr
      cl_hl = sqrt_ghl + 1.5*max(0, 0.5*g*(hr*hr-hl*hl)/(hr*sqrt_ghr) + vl - vr )
      cr_hr = sqrt_ghr + 1.5*max(0, 0.5*g*(hl*hl-hr*hr)/(cl_hl*hl) + vl - vr )
    elseif hr <= hl && 0 < hl
      cr_hr = sqrt_ghr + 1.5*max(0, 0.5*g*(hl*hl-hr*hr)/(hl*sqrt_ghl) + vl - vr )
      cl_hl = sqrt_ghl + 1.5*max(0, 0.5*g*(hr*hr-hl*hl)/(cr_hr*hr) + vl - vr )
    else
      cl_hl = cr_hr = zero(T)
    end

    # compute intermediate values
    cl = cl_hl*hl
    cr = cr_hr*hr
    vls = vrs = ( cl*vl+cr*vr+0.5*g*(hl*hl-hr*hr) ) / ( cl+cr )
    vls = vrs = ifelse(isnan(vls), zero(T), vls)
    πls = πrs = ( 0.5*g*(cr*hl*hl+cl*hr*hr)-cl*cr*(vr-vl) ) / (cl+cr)
    πls = πrs = ifelse(isnan(πls), zero(T), πls)
    hls = 1 / ( 1/hl + (cr*(vr-vl)+0.5*g*(hl*hl-hr*hr))/(cl*(cl+cr)) )
    hls = ifelse(isnan(hls), zero(T), hls)
    hrs = 1 / ( 1/hr + (cl*(vr-vl)+0.5*g*(hr*hr-hl*hl))/(cr*(cl+cr)) )
    hrs = ifelse(isnan(hrs), zero(T), hrs)

    # compute fluxes
    fnum_h  = zero(T)
    fnum_hv = zero(T)
    if 0 <= vl-cl_hl
      fnum_h  = hl*vl
      fnum_hv = hl*vl*vl + 0.5*g*hl*hl
    elseif 0 <= vls
      fnum_h  = hls*vls
      fnum_hv = hls*vls*vls + πls
    elseif 0 <= vr+cr_hr
      fnum_h  = hrs*vrs
      fnum_hv = hrs*vrs*vrs + πrs
    else
      fnum_h  = hr*vr
      fnum_hv = hr*vr*vr + 0.5*g*hr*hr
    end

    SVector(fnum_h, fnum_hv)
end


"""
    (::KineticFlux)(uₗ::ShallowWaterVar1D, uᵣ::ShallowWaterVar1D, model::ShallowWater{T,1})

The kinetic relaxation solver / numerical flux, see
Perthame and Simeoni (2001) A kinetic scheme for the Saint-Venant system with a
source term, Calcolo 38(4), pp. 201-231.
"""
function (::KineticFlux)(uₗ::ShallowWaterVar1D, uᵣ::ShallowWaterVar1D, model::ShallowWater{T,1}) where T
    hₗ, vₗ = primitive_variables(uₗ, model)
    hᵣ, vᵣ = primitive_variables(uᵣ, model)
    @unpack g = model

    # integrals for \xi \geq 0
    fnum_h, fnum_hv = kinetic_integrals_ξ_geq_0(hₗ, vₗ, g)
    # use symmetry for \xi \leq 0
    val_h, val_hv   = kinetic_integrals_ξ_geq_0(hᵣ, -vᵣ, g)
    fnum_h  -= val_h
    fnum_hv += val_hv

    SVector(fnum_h, fnum_hv)
end

# \int_{\xi \geq 0} (\xi, xi^2)  M_l(\xi) \dif \xi  for (fnum_h, fnum_hv)
function kinetic_integrals_ξ_geq_0(h, v, g)
    val_h  = zero(h)
    val_hv = zero(h)
    sqrt_2gh = sqrt(2g*h)

    if h > 0 && v >= sqrt_2gh
        val_h += h*v

        val_hv += h*v*v + g*h*h/2
    elseif h > 0 && sqrt_2gh+v > 0 && sqrt_2gh-v > 0
        sqrt_2ghmv2 = sqrt(2g*h-v^2)
        atan_v_sqrt_2ghmv2 = atan(v/sqrt_2ghmv2)

        val_h += h * v / 2
        val_h += h * sqrt_2ghmv2 * 2 / (3π)
        val_h += v^2 * sqrt_2ghmv2 / (6g*π)
        val_h += h*v * atan_v_sqrt_2ghmv2 / π

        val_hv += ( h*v*v + g*h*h/2 ) / 2
        val_hv += h*v*sqrt_2ghmv2 * 13 / (12π)
        val_hv += v^3*sqrt_2ghmv2 / (12g*π)
        val_hv += ( h*v^2 + g*h^2/2 ) * atan_v_sqrt_2ghmv2 / π
    end

    val_h, val_hv
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
  @unpack g = sol.prob.model

  if ξ <= σ₁⁻
    uₗ
  elseif ξ <= σ₁⁺
    rarefaction_wave_1(ξ, uₗ, g)
  elseif ξ <= σ₂⁻
    uₘ
  elseif ξ <= σ₂⁺
    rarefaction_wave_2(ξ, uᵣ, g)
  else # σ₂⁺ < ξ
    uᵣ
  end
end

function rarefaction_wave_1{T}(ξ, uₗ::ShallowWaterVar1D{T}, g)
  hₗ, vₗ = primitive_variables(uₗ)
  h  = (vₗ + 2*sqrt(g*hₗ) -  ξ)^2 / 9g
  hv = (vₗ + 2*sqrt(g*hₗ) + 2ξ)*h / 3
  ShallowWaterVar1D{T}(h, hv)
end

function rarefaction_wave_2{T}(ξ, uₗ::ShallowWaterVar1D{T}, g)
  hₗ, vₗ = primitive_variables(uₗ)
  h  = (vₗ - 2*sqrt(g*hₗ) -  ξ)^2 / 9g
  hv = (vₗ - 2*sqrt(g*hₗ) + 2ξ)*h / 3
  ShallowWaterVar1D{T}(h, hv)
end


function (sol::ShallowWaterRiemannSolution)(t::Real, x::Real)
  @unpack x₀, t₀ = sol.prob

  sol((x-x₀)/(t-t₀))
end


function compute_state_and_speeds{T}(uₗ::ShallowWaterVar1D{T}, uᵣ::ShallowWaterVar1D{T}, model::ShallowWater)
  @unpack g = model
  hₗ, vₗ = primitive_variables(uₗ)
  hᵣ, vᵣ = primitive_variables(uᵣ)

  if v_S1_SWE(hᵣ,hₗ,vₗ,g) <= vᵣ && vᵣ <= v_R2_SWE(hᵣ,hₗ,vₗ,g)
    # Case I: wave1 = shock, wave2 = rarefaction
    if hₗ ≈ 0
      hₘ = hₗ
      vₘ = v_R2_SWE(hₘ,hᵣ,vᵣ,g)
      uₘ = ShallowWaterVar1D{T}(hₘ, hₘ*vₘ)
      σ₁⁻ = σ₁⁺ = vₘ - sqrt(g*hₗ)
      σ₂⁻ = vₘ + sqrt(g*hₘ)
      σ₂⁺ = vᵣ + sqrt(g*hᵣ)
    else
      hₘ = Roots.fzero(h -> vᵣ-v_R2_SWE(hᵣ,h,v_S1_SWE(h,hₗ,vₗ,g),g), hₗ, hᵣ)
      # this is better than vₘ = v_S1_SWE(hₘ,hₗ,vₗ) for hₗ≈0
      vₘ = v_R2_SWE(hₘ,hᵣ,vᵣ,g)
      uₘ = ShallowWaterVar1D{T}(hₘ, hₘ*vₘ)
      # this is better than σ₁⁻ = σ₁⁺ = vₗ - sqrt(hₘ+hₘ^2/hₗ) / sqrt(2) for hₗ≈0
      σ₁⁻ = σ₁⁺ = vₘ - sqrt(hₗ+hₗ^2/hₘ) * sqrt(g/2)
      σ₂⁻ = vₘ + sqrt(g*hₘ)
      σ₂⁺ = vᵣ + sqrt(g*hᵣ)
    end
  elseif v_R1_SWE(hᵣ,hₗ,vₗ,g) <= vᵣ && v_R2_SWE(hᵣ,hₗ,vₗ,g) <= vᵣ
    # Case II: wave1 = rarefaction, wave2 = rarefaction
    if 2*(sqrt(g*hᵣ)+sqrt(g*hₗ)) >= vᵣ-vₗ
      # no vacuum
      hₘ = (2*(sqrt(g*hᵣ)+sqrt(g*hₗ)) -(vᵣ-vₗ))^2 / 16g
      vₘ = v_R1_SWE(hₘ,hₗ,vₗ,g)
      uₘ = ShallowWaterVar1D{T}(hₘ,hₘ*vₘ)
      σ₁⁻ = vₗ - sqrt(g*hₗ)
      σ₁⁺ = vₘ - sqrt(g*hₘ)
      σ₂⁻ = vₘ + sqrt(g*hₘ)
      σ₂⁺ = vᵣ + sqrt(g*hᵣ)
    else
      # vacuum
      uₘ = ShallowWaterVar1D{T}(0,0)
      σ₁⁻ = vₗ - sqrt(g*hₗ)
      σ₁⁺ = vₗ + 2*sqrt(g*hₗ)
      σ₂⁻ = vᵣ - 2*sqrt(g*hᵣ)
      σ₂⁺ = vᵣ + sqrt(g*hᵣ)
    end
  elseif v_S2_SWE(hᵣ,hₗ,vₗ,g) <= vᵣ && vᵣ <= v_R1_SWE(hᵣ,hₗ,vₗ,g)
    # Case III: wave1 = rarefaction, wave2 = shock
    if hᵣ ≈ 0
      hₘ = hᵣ
      vₘ = v_R1_SWE(hₘ,hₗ,vₗ,g)
      uₘ = ShallowWaterVar1D{T}(hₘ,hₘ*vₘ)
      σ₁⁻ = vₗ - sqrt(g*hₗ)
      σ₁⁺ = vₘ - sqrt(g*hₘ)
      σ₂⁻ = σ₂⁺ = vₘ + sqrt(g*hᵣ)
    else
      hₘ = Roots.fzero(h->vᵣ-v_S2_SWE(hᵣ,h,v_R1_SWE(h,hₗ,vₗ,g),g), hᵣ, hₗ)
      vₘ = v_R1_SWE(hₘ,hₗ,vₗ,g)
      uₘ = ShallowWaterVar1D{T}(hₘ,hₘ*vₘ)
      σ₁⁻ = vₗ - sqrt(g*hₗ)
      σ₁⁺ = vₘ - sqrt(g*hₘ)
      σ₂⁻ = σ₂⁺ = vₘ + sqrt(hᵣ+hᵣ^2/hₘ) * sqrt(g/2)
    end
  else #vᵣ <= v_S1_SWE(hᵣ,hₗ,vₗ) && vᵣ <= v_S2_SWE(hᵣ,hₗ,vₗ)
    # Case IV: wave1 = shock, wave2 = shock
    # TODO: Chose some interval?
    hₘ = Roots.fzero(h->vᵣ-v_S2_SWE(hᵣ,h,v_S1_SWE(h,hₗ,vₗ,g),g), max(hₗ,hᵣ))
    vₘ = v_S1_SWE(hₘ,hₗ,vₗ,g)
    uₘ = ShallowWaterVar1D{T}(hₘ,hₘ*vₘ)
    σ₁⁻ = σ₁⁺ = vₗ - hₘ * sqrt(1/hₘ+1/hₗ) * sqrt(g/2)
    σ₂⁻ = σ₂⁺ = vₘ + hᵣ * sqrt(1/hᵣ+1/hₘ) * sqrt(g/2)
  end

  uₘ, σ₁⁻, σ₁⁺, σ₂⁻, σ₂⁺
end


# wave curves parameterised by the height h
v_R1_SWE(h, hₗ, vₗ, g) = vₗ - 2*(sqrt(g*h) - sqrt(g*hₗ))
v_R2_SWE(h, hₗ, vₗ, g) = vₗ + 2*(sqrt(g*h) - sqrt(g*hₗ))
v_S1_SWE(h, hₗ, vₗ, g) = vₗ - (h-hₗ) * sqrt( (1/h+1/hₗ)*g/2 )
v_S2_SWE(h, hₗ, vₗ, g) = vₗ + (h-hₗ) * sqrt( (1/h+1/hₗ)*g/2 )



@recipe function f{Ξ,T}(ξu::Tuple{Ξ,Vector{ShallowWaterVar1D{T}}})
  ξ, u = ξu

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
