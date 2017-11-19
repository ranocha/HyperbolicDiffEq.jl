
doc"
    GodunovFlux

Godunovs flux $f(u(x/t = 0))$, i.e. the flux of the value at zero of the solution
of the Riemann problem.
"
struct GodunovFlux <: NumericalFlux end


doc"
    LocalLaxFriedrichsFlux{MaxAbsSpeed}

The local Lax-Friedrichs flux $\frac{f(u_r) + f(u_l)}{2} - \frac{\lambda}{2} (u_r - u_l)$.
$\lambda$ is the maximal absolute value of the speed in the solution of a Riemann
problem, computed by `max_abs_speed::MaxAbsSpeed`.
"
struct LocalLaxFriedrichsFlux{MaxAbsSpeed} <: NumericalFlux
    max_abs_speed::MaxAbsSpeed
end

LocalLaxFriedrichsFlux() = LocalLaxFriedrichsFlux(max_abs_speed)

Base.@pure function (fnum::LocalLaxFriedrichsFlux)(uₗ, uᵣ, model::AbstractBalanceLaw{1})
    λ = fnum.max_abs_speed(uₗ, uᵣ, model)

    (flux(uₗ,model) + flux(uᵣ,model))/2 - λ/2 * (uᵣ - uₗ)
end

Base.@pure function (fnum::LocalLaxFriedrichsFlux)(uₗ, uᵣ, model::AbstractBalanceLaw, direction)
    λ = fnum.max_abs_speed(uₗ, uᵣ, model, direction)

    (flux(uₗ, model, direction) + flux(uᵣ, model, direction))/2 - λ/2 * (uᵣ - uₗ)
end


"""
    max_abs_speed(uₗ, uᵣ, model::AbstractBalanceLaw)

Compute the maximal absolute value of the speed in the solution of the Riemann
problem with states `uₗ`, `uᵣ` for `model`.
"""
@inline function max_abs_speed(uₗ, uᵣ, model::AbstractBalanceLaw)
    naive_max_abs_speed(uₗ, uᵣ, model)
end

@inline function max_abs_speed(uₗ, uᵣ, model::AbstractBalanceLaw, direction)
    naive_max_abs_speed(uₗ, uᵣ, model, direction)
end

@inline function naive_max_abs_speed(uₗ, uᵣ, model::AbstractBalanceLaw)
    max(max_abs_speed(uₗ,model), max_abs_speed(uᵣ,model))
end

@inline function naive_max_abs_speed(uₗ, uᵣ, model::AbstractBalanceLaw, direction)
    max(max_abs_speed(uₗ,model,direction), max_abs_speed(uᵣ,model,direction))
end

"""
    max_abs_speed{T}(u::T, model::ScalarBalanceLaw{T,1})

Compute the maximal absolute value of speed at `u` for `model`.
"""
@inline max_abs_speed{T}(u::T, model::ScalarBalanceLaw{T,1}) = abs(speed(u, model))

"""
    min_max_speed{T}(u::T, model::ScalarBalanceLaw{T,1})

Compute the minimal and maximal speed at `u` for `model`.
"""
@inline function min_max_speed{T}(u::T, model::ScalarBalanceLaw{T,1})
    λ = speed(u, model)
    λ, λ
end


doc"
    HartenLaxVanLeerFlux{MinMaxSpeed}

The HLL (Harten, Lax, van Leer) flux using `min_max_speed::MinMaxSpeed` to compute
 the (approximate) minimal and maximal speed in the solution of a Riemann problem.
"
struct HartenLaxVanLeerFlux{MinMaxSpeed} <: NumericalFlux
    min_max_speed::MinMaxSpeed
end

const HLL = HartenLaxVanLeerFlux

HartenLaxVanLeerFlux() = HartenLaxVanLeerFlux(min_max_speed)

@inline function (fnum::HartenLaxVanLeerFlux)(uₗ, uᵣ, model::AbstractBalanceLaw{1})
    λ₋, λ₊ = fnum.min_max_speed(uₗ, uᵣ, model)

    if 0 <= λ₋
        flux(uₗ, model)
    elseif 0 < λ₊
        ( λ₊*flux(uₗ, model) - λ₋*flux(uᵣ, model) + λ₋*λ₊*(uᵣ-uₗ) ) / (λ₊ - λ₋)
    else
        flux(uᵣ, model)
    end
end

@inline function (fnum::HartenLaxVanLeerFlux)(uₗ, uᵣ, model::AbstractBalanceLaw, direction)
    λ₋, λ₊ = fnum.min_max_speed(uₗ, uᵣ, model, direction)

    if 0 <= λ₋
        flux(uₗ, model, direction)
    elseif 0 < λ₊
        ( λ₊*flux(uₗ, model, direction) - λ₋*flux(uᵣ, model, direction) + λ₋*λ₊*(uᵣ-uₗ) ) / (λ₊ - λ₋)
    else
        flux(uᵣ, model, direction)
    end
end


"""
    min_max_speed(uₗ, uᵣ, model::AbstractBalanceLaw{1})

Compute the maximal absolute value of the speed in the solution of the Riemann
problem with states `uₗ`, `uᵣ` for `model`.
"""
@inline function min_max_speed(uₗ, uᵣ, model::AbstractBalanceLaw{1})
    λₗ₋, λₗ₊ = min_max_speed(uₗ, model)
    λᵣ₋, λᵣ₊ = min_max_speed(uᵣ, model)
    min(λₗ₋, λᵣ₋), max(λₗ₊, λᵣ₊)
end

@inline function min_max_speed(uₗ, uᵣ, model::AbstractBalanceLaw, direction)
    λₗ₋, λₗ₊ = min_max_speed(uₗ, model, direction)
    λᵣ₋, λᵣ₊ = min_max_speed(uᵣ, model, direction)
    min(λₗ₋, λᵣ₋), max(λₗ₊, λᵣ₊)
end


"""
    EnergyConservativeFlux

The "standard" energy conservative flux.
"""
struct EnergyConservativeFlux <: NumericalFlux end

"""
    EnergyConservativeFlux1Param

A one-parameter family of energy conservative fluxes.
"""
struct EnergyConservativeFlux1Param{T} <: NumericalFlux
    a₁::T
end

"""
    EnergyConservativeFlux2Param

A two-parameter family of energy conservative fluxes.
"""
struct EnergyConservativeFlux2Param{T} <: NumericalFlux
    a₁::T
    a₂::T
end

doc"
    L2L4ConservativeFlux

A numerical flux conserving the $L^2 \cap L^4$ entropy U(U) = u^2 + u^4.
"
struct L2L4ConservativeFlux <: NumericalFlux end

doc"
    L2L2sConservativeFlux

A numerical flux conserving the $L^2 \cap L^{2s}$ entropy U(U) = u^2 + u^(2s).
"
struct L2L2sConservativeFlux{s} <: NumericalFlux
    half_pow::Val{s}

    function L2L2sConservativeFlux(half_pow::Val{s}) where s
        @assert typeof(s) <: Integer && s >= 1
        new{s}(half_pow)
    end
end


"""
    jump_pow_u_r_over_jump_u(uₗ, uᵣ, ::Val{r})

Compute `(uᵣ^r - uₗ^r) / (uᵣ - uₗ)` for integer `r >= 2` via the expanded
expression `uᵣ^(r-1) + rᵣ^(r-2)*uₗ + ... + uᵣ*uₗ^(r-2) + uₗ^(r-1)`.
"""
@generated function jump_pow_u_r_over_jump_u(uₗ, uᵣ, ::Val{r}) where r
    @assert typeof(r) <: Integer && r >= 2
    ex = :( Base.FastMath.pow_fast(uₗ, $r-1) )
    ex = :( $ex + Base.FastMath.pow_fast(uᵣ, $r-1) )
    for pₗ in Base.OneTo(r-2)
        pᵣ = r-1 - pₗ
        ex = :( $ex + Base.FastMath.pow_fast(uₗ, $pₗ) * Base.FastMath.pow_fast(uᵣ, $pᵣ) )
    end
    return :(Base.@_inline_meta; $ex)
end


doc"
    CentralFlux

The central flux $\frac{f(u_r) + f(u_l)}{2}$.
"
struct CentralFlux <: NumericalFlux end

Base.@pure function (fnum::CentralFlux)(uₗ, uᵣ, model)
    (flux(uₗ, model) + flux(uᵣ, model)) / 2
end

Base.@pure function (fnum::CentralFlux)(uₗ, uᵣ, model, direction)
    (flux(uₗ, model, direction) + flux(uᵣ, model, direction)) / 2
end


"""
    MorinishiFlux

The numerical flux corresponding to the splitting of Morinishi (2010),
see Gassner, Winter, Kopriva (2016).
"""
struct MorinishiFlux <: NumericalFlux end

"""
    DucrosEtAlFlux

The numerical flux corresponding to the splitting of Ducros, Laporte, Soulères,
Guinot, Moinat, Caruelle (2000), see Gassner, Winter, Kopriva (2016).
"""
struct DucrosEtAlFlux <: NumericalFlux end

"""
    KennedyGruberFlux

The numerical flux corresponding to the splitting of Kennedy and Gruber (2008),
see Gassner, Winter, Kopriva (2016).
"""
struct KennedyGruberFlux <: NumericalFlux end

"""
    PirozzoliFlux

The numerical flux corresponding to the splitting of Pirozzoli (2011),
see Gassner, Winter, Kopriva (2016).
"""
struct PirozzoliFlux <: NumericalFlux end

"""
    SuliciuFlux

The Suliciu relaxation solver.
"""
struct SuliciuFlux <: NumericalFlux end

"""
    KineticFlux

A kinetic relaxation solver.
"""
struct KineticFlux <: NumericalFlux end

"""
    ChandrashekarFluxEC

The entropy conservative flux of Chandrashekar (2013).
"""
struct ChandrashekarFluxEC <: NumericalFlux end

"""
    IsmailRoeFluxEC

The entropy conservative flux of Ismail and Roe (2009).
"""
struct IsmailRoeFluxEC <: NumericalFlux end

"""
    RanochaFluxECandKEP

The entropy conservative and kinetic energy preserving flux of Ranocha (2017).
"""
struct RanochaFluxECandKEP <: NumericalFlux end



"""
    FluxPlusDissipation{Fnum::NumericalFlux, Diss} <: NumericalFlux

A numerical flux given by the (possibly central / symmetric) numerical flux
`fnum` and the dissipation operator `diss`.
"""
struct FluxPlusDissipation{Fnum<:NumericalFlux, Diss<:DissipationOperator} <: NumericalFlux
    fnum::Fnum
    diss::Diss
end

@inline function (fnumdiss::FluxPlusDissipation)(uₗ, uᵣ, model::AbstractBalanceLaw)
    @unpack fnum, diss = fnumdiss

    fnum(uₗ, uᵣ, model) + diss(uₗ, uᵣ, model)
end


doc"
    LocalLaxFriedrichsDissipation{MaxAbsSpeed}

The local Lax-Friedrichs dissipation $- \frac{\lambda}{2} (u_r - u_l)$.
$\lambda$ is the approximate maximal absolute value of the speed in the solution
of a Riemann problem, computed by `max_abs_speed::MaxAbsSpeed`.
"
struct LocalLaxFriedrichsDissipation{MaxAbsSpeed} <: DissipationOperator
    max_abs_speed::MaxAbsSpeed
end

LocalLaxFriedrichsDissipation() = LocalLaxFriedrichsDissipation(naive_max_abs_speed)

@inline function (diss::LocalLaxFriedrichsDissipation)(uₗ, uᵣ, model::AbstractBalanceLaw)
    λ = diss.max_abs_speed(uₗ, uᵣ, model)

    -λ * (uᵣ - uₗ) / 2
end



doc"
    ScalarDissipation{MaxAbsSpeed}

The scalar dissipation $- \frac{\lambda}{2} R \cdot R^T \cdot  (w_r - w_l)$.
$\lambda$ is the approximate maximal absolute value of the speed, computed by
`max_abs_speed::MaxAbsSpeed`.
"
struct ScalarDissipation{MaxAbsSpeed} <: DissipationOperator
    max_abs_speed::MaxAbsSpeed
end

ScalarDissipation() = ScalarDissipation(naive_max_abs_speed)


doc"
    MatrixDissipation

The matrix dissipation $- \frac{1}{2} R \cdot |\Lambda| \cdot R^T \cdot  (w_r - w_l)$.
"
struct MatrixDissipation <: DissipationOperator end
