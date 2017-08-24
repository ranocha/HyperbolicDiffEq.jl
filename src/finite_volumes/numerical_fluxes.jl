
"""
    local_lax_friedrichs(uₗ, uᵣ, model::BuckleyLeverette)

Compute the local Lax-Friedrichs flux between `uₗ` and `uᵣ` for `model`.
"""
function local_lax_friedrichs(uₗ, uᵣ, model::AbstractBalanceLaw{1})
    λ = max_abs_speed(uₗ, uᵣ, model)

    (flux(uₗ,model) + flux(uᵣ,model))/2 - λ/2 * (uᵣ - uₗ)
end


"""
    max_abs_speed(uₗ, uᵣ, model::AbstractBalanceLaw{1})

Compute the maximal absolute value of the speed in the solution of the Riemann
problem with states `uₗ`, `uᵣ` for `model`.
"""
function max_abs_speed(uₗ, uᵣ, model::AbstractBalanceLaw{1})
    max(max_abs_speed(uₗ,model), max_abs_speed(uᵣ,model))
end
