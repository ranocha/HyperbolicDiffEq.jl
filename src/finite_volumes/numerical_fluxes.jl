
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
    max_abs_speed(uₗ, uᵣ, model::AbstractBalanceLaw{1})

Compute the maximal absolute value of the speed in the solution of the Riemann
problem with states `uₗ`, `uᵣ` for `model`.
"""
Base.@pure function max_abs_speed(uₗ, uᵣ, model::AbstractBalanceLaw{1})
    max(max_abs_speed(uₗ,model), max_abs_speed(uᵣ,model))
end

Base.@pure function max_abs_speed(uₗ, uᵣ, model::AbstractBalanceLaw, direction)
    max(max_abs_speed(uₗ,model,direction), max_abs_speed(uᵣ,model,direction))
end


"""
    EnergyConservativeFlux

The "standard" energy conservative flux.
"""
struct EnergyConservativeFlux <: NumericalFlux end


doc"
    CentralFlux

The central flux $\frac{f(u_r) + f(u_l)}{2}$.
"
struct CentralFlux <: NumericalFlux end

Base.@pure function (fnum::CentralFlux)(uₗ, uᵣ, model, direction)
    (flux(uₗ, model, direction) + flux(uᵣ, model, direction)) / 2
end


"""
    SuliciuFlux

The Suliciu relaxation solver.
"""
struct SuliciuFlux <: NumericalFlux end



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
