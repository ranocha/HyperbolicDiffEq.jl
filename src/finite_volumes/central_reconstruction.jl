
"""
    CentralReconstruction{K}

A simple, central, polynomial reconstruction (linear procedure) without adaptive
techniques using polynomials of degee `K-1`.
"""
struct CentralReconstruction{K} <: AbstractReconstruction
    function CentralReconstruction{K}() where {K}
        @argcheck (K>0 && isodd(K)) ArgumentError("`K` must be odd (K = $K).")
        new{K}()
    end
end

function CentralReconstruction(::Val{K}=Val{3}()) where {K}
    CentralReconstruction{K}()
end

function CentralReconstruction(K::Integer)
    CentralReconstruction{K}()
end


@inline stencil_width(::CentralReconstruction{K}) where {K} = K
@inline stencil_width_val(::CentralReconstruction{K}) where {K} = Val{K}()

@inline order(::CentralReconstruction{K}) where {K} = K


function (recons::CentralReconstruction{3})(edge_u, cell, u_m1, u_0, u_p1, balance_law, meshx)
    edge_u[1,cell] = ( 5*u_0 - u_p1 + 2*u_m1 ) / 6
    edge_u[2,cell] = ( 5*u_0 + 2*u_p1 - u_m1 ) / 6

    nothing
end

function (recons::CentralReconstruction{5})(edge_u, cell, u_m2, u_m1, u_0, u_p1, u_p2, balance_law, meshx)
    edge_u[1,cell] = ( 47*u_0 - 13*u_p1 + 2*u_p2 + 27*u_m1 - 3*u_m2 ) / 60
    edge_u[2,cell] = ( 47*u_0 + 27*u_p1 - 3*u_p2 - 13*u_m1 + 2*u_m2 ) / 60
   
    nothing
end

function (recons::CentralReconstruction{7})(edge_u, cell, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, balance_law, meshx)
    edge_u[1,cell] = ( 319*u_0 - 101*u_p1 + 25*u_p2 - 3*u_p3 + 214*u_m1 - 38*u_m2 + 4*u_m3 ) / 420
    edge_u[2,cell] = ( 319*u_0 + 214*u_p1 - 38*u_p2 + 4*u_p3 - 101*u_m1 + 25*u_m2 - 3*u_m3 ) / 420
    
    nothing
end


function interpolate!(uval, ξ, u_m1, u_0, u_p1, balance_law, mesh, recons::CentralReconstruction{3})
    for i in eachindex(uval)
        uval[i] = @evalpoly(ξ[i], 5*u_0/6 - u_p1/6 + u_m1/3, u_0 - u_m1, -u_0 + u_p1/2 + u_m1/2)
    end
    
    nothing
end

function interpolate!(uval, ξ, u_m2, u_m1, u_0, u_p1, u_p2, balance_law, mesh, recons::CentralReconstruction{5})
    for i in eachindex(uval)
        uval[i] = @evalpoly(ξ[i], 47*u_0/60 - 13*u_p1/60 + u_p2/30 + 9*u_m1/20 - u_m2/20, 5*u_0/4 - u_p1/12 - 5*u_m1/4 + u_m2/12, -u_0 + 3*u_p1/4 - u_p2/8 + u_m1/4 + u_m2/8, -u_0/2 + u_p1/6 + u_m1/2 - u_m2/6, u_0/4 - u_p1/6 + u_p2/24 - u_m1/6 + u_m2/24)
    end
    
    nothing
end

function interpolate!(uval, ξ, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, balance_law, mesh, recons::CentralReconstruction{7})
    for i in eachindex(uval)
        uval[i] = @evalpoly(ξ[i], 319*u_0/420 - 101*u_p1/420 + 5*u_p2/84 - u_p3/140 + 107*u_m1/210 - 19*u_m2/210 + u_m3/105, 49*u_0/36 - 5*u_p1/36 + u_p2/90 - 49*u_m1/36 + 5*u_m2/36 - u_m3/90, -23*u_0/24 + 7*u_p1/8 - 19*u_p2/80 + 7*u_p3/240 + u_m1/16 + 21*u_m2/80 - u_m3/30, -7*u_0/9 + 11*u_p1/36 - u_p2/36 + 7*u_m1/9 - 11*u_m2/36 + u_m3/36, 23*u_0/72 - 13*u_p1/48 + 5*u_p2/48 - u_p3/72 - u_m1/6 + u_m2/48 + u_m3/144, u_0/12 - u_p1/24 + u_p2/120 - u_m1/12 + u_m2/24 - u_m3/120, -u_0/36 + u_p1/48 - u_p2/120 + u_p3/720 + u_m1/48 - u_m2/120 + u_m3/720)
    end

    nothing
end
