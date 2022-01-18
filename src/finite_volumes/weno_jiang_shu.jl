
"""
    WENOJiangShu{K}

The classical weighted essentially non-oscillatory (WENO reconstruction of
Jiang and Shu on `K` cells.
"""
struct WENOJiangShu{K, T<:Real} <: AbstractReconstruction
  ε::T
end

WENOJiangShu{K}(ε) where {K} = WENOJiangShu{K, typeof(ε)}(ε)

WENOJiangShu{K}() where {K} = WENOJiangShu{K}(1.0e-6)

@inline stencil_width(::WENOJiangShu{K}) where {K} = K
@inline stencil_width_val(::WENOJiangShu{K}) where {K} = Val{K}()

@inline order(::WENOJiangShu{K}) where {K} = K


function (weno::WENOJiangShu{5})(edge_u, cell, u_m2, u_m1, u_0, u_p1, u_p2, balance_law, meshx)
    ε = weno.ε
    T = typeof(ε)

    # smoothness indicators of the three three-point stencils
    β0 = 13 * (u_0  - 2 * u_p1 + u_p2)^2 / 12 + (3 * u_0 - 4 * u_p1 + u_p2)^2 / 4
    β1 = 13 * (u_m1 - 2 * u_0  + u_p1)^2 / 12 + (u_m1 - u_p1)^2 / 4
    β2 = 13 * (u_m2 - 2 * u_m1 + u_0 )^2 / 12 + (u_m2 - 4 * u_m1 + 3 * u_0)^2 / 4

    # left edge
    edge_u_l0 = ( 11*u_0 - 7*u_p1 + 2*u_p2 ) / 6
    edge_u_l1 = ( 5*u_0 - u_p1 + 2*u_m1 ) / 6
    edge_u_l2 = ( 2*u_0 + 5*u_m1 - u_m2 ) / 6

    d_l0 = 1 * one(T) / 10
    d_l1 = 3 * one(T) / 5
    d_l2 = 3 * one(T) / 10

    α_l0 = d_l0 / (ε + β0)^2
    α_l1 = d_l1 / (ε + β1)^2
    α_l2 = d_l2 / (ε + β2)^2

    sum_α_l = α_l0 + α_l1 + α_l2

    ω_l0 = α_l0 / sum_α_l
    ω_l1 = α_l1 / sum_α_l
    ω_l2 = α_l2 / sum_α_l

    # right edge
    edge_u_r0 = ( 2*u_0 + 5*u_p1 - u_p2 ) / 6
    edge_u_r1 = ( 5*u_0 + 2*u_p1 - u_m1 ) / 6
    edge_u_r2 = ( 11*u_0 - 7*u_m1 + 2*u_m2 ) / 6

    d_r0 = 3 * one(T) / 10
    d_r1 = 3 * one(T) / 5
    d_r2 = 1 * one(T) / 10

    α_r0 = d_r0 / (ε + β0)^2
    α_r1 = d_r1 / (ε + β1)^2
    α_r2 = d_r2 / (ε + β2)^2

    sum_α_r = α_r0 + α_r1 + α_r2

    ω_r0 = α_r0 / sum_α_r
    ω_r1 = α_r1 / sum_α_r
    ω_r2 = α_r2 / sum_α_r

    # assign values
    @inbounds begin
        edge_u[1,cell] = ω_l0 * edge_u_l0 + ω_l1 * edge_u_l1 + ω_l2 * edge_u_l2
        edge_u[2,cell] = ω_r0 * edge_u_r0 + ω_r1 * edge_u_r1 + ω_r2 * edge_u_r2
    end

    nothing
end
