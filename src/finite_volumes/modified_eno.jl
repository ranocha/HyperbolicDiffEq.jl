
"""
    ModifiedENO{K, ChooseStencil}

A reconstruction using the same polynomials as the essentially non-oscillatory
(ENO) reconstruction on `K` cells. The stencil is chosen by `choose_stencil`.
"""
struct ModifiedENO{K, ChooseStencil} <: AbstractReconstruction
    choose_stencil::ChooseStencil
end

function ModifiedENO{K}(::Val{K}=Val{2}(), choose_stencil=ClassicalChoiceENO())
    ModifiedENO{K,typeof(choose_stencil)}(choose_stencil)
end

function ModifiedENO(K::Integer, choose_stencil=ClassicalChoiceENO())
    ModifiedENO{K,typeof(choose_stencil)}(choose_stencil)
end


# ENO{K} uses K different polynomials of degree K-1
# -> the total stencil width is 2K-1
@inline stencil_width{K}(::ModifiedENO{K}) = 2K-1
@inline stencil_width_val(::ModifiedENO{1}) = Val{1}()
@inline stencil_width_val(::ModifiedENO{2}) = Val{3}()
@inline stencil_width_val(::ModifiedENO{3}) = Val{5}()
@inline stencil_width_val(::ModifiedENO{4}) = Val{7}()
@inline stencil_width_val(::ModifiedENO{5}) = Val{9}()
@inline stencil_width_val(::ModifiedENO{6}) = Val{11}()
@inline stencil_width_val(::ModifiedENO{7}) = Val{13}()
@inline stencil_width_val(::ModifiedENO{8}) = Val{15}()

@inline order{K}(::ModifiedENO{K}) = K


function (eno::ModifiedENO{2})(edge_u, cell, u_m1, u_0, u_p1, balance_law, meshx)
    idx = eno.choose_stencil(u_m1, u_0, u_p1)

    @inbounds if idx == -1
        edge_u[1,cell] = ( u_0 + u_m1 ) / 2
        edge_u[2,cell] = ( 3*u_0 - u_m1 ) / 2
    else # idx == 0
        edge_u[1,cell] = ( 3*u_0 - u_p1 ) / 2
        edge_u[2,cell] = ( u_0 + u_p1 ) / 2
    end
    nothing
end

function (eno::ModifiedENO{3})(edge_u, cell, u_m2, u_m1, u_0, u_p1, u_p2, balance_law, meshx)
    idx = eno.choose_stencil(u_m2, u_m1, u_0, u_p1, u_p2)

    @inbounds if idx == -2
        edge_u[1,cell] = ( 2*u_0 + 5*u_m1 - u_m2 ) / 6
        edge_u[2,cell] = ( 11*u_0 - 7*u_m1 + 2*u_m2 ) / 6
    elseif idx == -1
        edge_u[1,cell] = ( 5*u_0 - u_p1 + 2*u_m1 ) / 6
        edge_u[2,cell] = ( 5*u_0 + 2*u_p1 - u_m1 ) / 6
    else # idx == 0
        edge_u[1,cell] = ( 11*u_0 - 7*u_p1 + 2*u_p2 ) / 6
        edge_u[2,cell] = ( 2*u_0 + 5*u_p1 - u_p2 ) / 6
    end
    nothing
end

function (eno::ModifiedENO{4})(edge_u, cell, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, balance_law, meshx)
    idx = eno.choose_stencil(u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3)

    @inbounds if idx == -3
        edge_u[1,cell] = ( 3*u_0 + 13*u_m1 - 5*u_m2 + u_m3 ) / 12
        edge_u[2,cell] = ( 25*u_0 - 23*u_m1 + 13*u_m2 - 3*u_m3 ) / 12
    elseif idx == -2
        edge_u[1,cell] = ( 7*u_0 - u_p1 + 7*u_m1 - u_m2 ) / 12
        edge_u[2,cell] = ( 13*u_0 + 3*u_p1 - 5*u_m1 + u_m2 ) / 12
    elseif idx == -1
        edge_u[1,cell] = ( 13*u_0 - 5*u_p1 + u_p2 + 3*u_m1 ) / 12
        edge_u[2,cell] = ( 7*u_0 + 7*u_p1 - u_p2 - u_m1 ) / 12
    else # idx == 0
        edge_u[1,cell] = ( 25*u_0 - 23*u_p1 + 13*u_p2 - 3*u_p3 ) / 12
        edge_u[2,cell] = ( 3*u_0 + 13*u_p1 - 5*u_p2 + u_p3 ) / 12
    end
    nothing
end

function (eno::ModifiedENO{5})(edge_u, cell, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, balance_law, meshx)
    idx = eno.choose_stencil(u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4)

    @inbounds if idx == -4
        edge_u[1,cell] = ( 12*u_0 + 77*u_m1 - 43*u_m2 + 17*u_m3 - 3*u_m4 ) / 60
        edge_u[2,cell] = ( 137*u_0 - 163*u_m1 + 137*u_m2 - 63*u_m3 + 12*u_m4 ) / 60
    elseif idx == -3
        edge_u[1,cell] = ( 27*u_0 - 3*u_p1 + 47*u_m1 - 13*u_m2 + 2*u_m3 ) / 60
        edge_u[2,cell] = ( 77*u_0 + 12*u_p1 - 43*u_m1 + 17*u_m2 - 3*u_m3 ) / 60
    elseif idx == -2
        edge_u[1,cell] = ( 47*u_0 - 13*u_p1 + 2*u_p2 + 27*u_m1 - 3*u_m2 ) / 60
        edge_u[2,cell] = ( 47*u_0 + 27*u_p1 - 3*u_p2 - 13*u_m1 + 2*u_m2 ) / 60
    elseif idx == -1
        edge_u[1,cell] = ( 77*u_0 - 43*u_p1 + 17*u_p2 - 3*u_p3 + 12*u_m1 ) / 60
        edge_u[2,cell] = ( 27*u_0 + 47*u_p1 - 13*u_p2 + 2*u_p3 - 3*u_m1 ) / 60
    else # idx == 0
        edge_u[1,cell] = ( 137*u_0 - 163*u_p1 + 137*u_p2 - 63*u_p3 + 12*u_p4 ) / 60
        edge_u[2,cell] = ( 12*u_0 + 77*u_p1 - 43*u_p2 + 17*u_p3 - 3*u_p4 ) / 60
    end
    nothing
end

function (eno::ModifiedENO{6})(edge_u, cell, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, balance_law, meshx)
    idx = eno.choose_stencil(u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5)

    @inbounds if idx == -5
        edge_u[1,cell] = ( 10*u_0 + 87*u_m1 - 63*u_m2 + 37*u_m3 - 13*u_m4 + 2*u_m5 ) / 60
        edge_u[2,cell] = ( 147*u_0 - 213*u_m1 + 237*u_m2 - 163*u_m3 + 62*u_m4 - 10*u_m5 ) / 60
    elseif idx == -4
        edge_u[1,cell] = ( 22*u_0 - 2*u_p1 + 57*u_m1 - 23*u_m2 + 7*u_m3 - u_m4 ) / 60
        edge_u[2,cell] = ( 87*u_0 + 10*u_p1 - 63*u_m1 + 37*u_m2 - 13*u_m3 + 2*u_m4 ) / 60
    elseif idx == -3
        edge_u[1,cell] = ( 37*u_0 - 8*u_p1 + u_p2 + 37*u_m1 - 8*u_m2 + u_m3 ) / 60
        edge_u[2,cell] = ( 57*u_0 + 22*u_p1 - 2*u_p2 - 23*u_m1 + 7*u_m2 - u_m3 ) / 60
    elseif idx == -2
        edge_u[1,cell] = ( 57*u_0 - 23*u_p1 + 7*u_p2 - u_p3 + 22*u_m1 - 2*u_m2 ) / 60
        edge_u[2,cell] = ( 37*u_0 + 37*u_p1 - 8*u_p2 + u_p3 - 8*u_m1 + u_m2 ) / 60
    elseif idx == -1
        edge_u[1,cell] = ( 87*u_0 - 63*u_p1 + 37*u_p2 - 13*u_p3 + 2*u_p4 + 10*u_m1 ) / 60
        edge_u[2,cell] = ( 22*u_0 + 57*u_p1 - 23*u_p2 + 7*u_p3 - u_p4 - 2*u_m1 ) / 60
    else # idx == 0
        edge_u[1,cell] = ( 147*u_0 - 213*u_p1 + 237*u_p2 - 163*u_p3 + 62*u_p4 - 10*u_p5 ) / 60
        edge_u[2,cell] = ( 10*u_0 + 87*u_p1 - 63*u_p2 + 37*u_p3 - 13*u_p4 + 2*u_p5 ) / 60
    end
    nothing
end

function (eno::ModifiedENO{7})(edge_u, cell, u_m6, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, u_p6, balance_law, meshx)
    idx = eno.choose_stencil(u_m6, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, u_p6)

    @inbounds if idx == -6
        edge_u[1,cell] = ( 60*u_0 + 669*u_m1 - 591*u_m2 + 459*u_m3 - 241*u_m4 + 74*u_m5 - 10*u_m6 ) / 420
        edge_u[2,cell] = ( 1089*u_0 - 1851*u_m1 + 2559*u_m2 - 2341*u_m3 + 1334*u_m4 - 430*u_m5 + 60*u_m6 ) / 420
    elseif idx == -5
        edge_u[1,cell] = ( 130*u_0 - 10*u_p1 + 459*u_m1 - 241*u_m2 + 109*u_m3 - 31*u_m4 + 4*u_m5 ) / 420
        edge_u[2,cell] = ( 669*u_0 + 60*u_p1 - 591*u_m1 + 459*u_m2 - 241*u_m3 + 74*u_m4 - 10*u_m5 ) / 420
    elseif idx == -4
        edge_u[1,cell] = ( 214*u_0 - 38*u_p1 + 4*u_p2 + 319*u_m1 - 101*u_m2 + 25*u_m3 - 3*u_m4 ) / 420
        edge_u[2,cell] = ( 459*u_0 + 130*u_p1 - 10*u_p2 - 241*u_m1 + 109*u_m2 - 31*u_m3 + 4*u_m4 ) / 420
    elseif idx == -3
        edge_u[1,cell] = ( 319*u_0 - 101*u_p1 + 25*u_p2 - 3*u_p3 + 214*u_m1 - 38*u_m2 + 4*u_m3 ) / 420
        edge_u[2,cell] = ( 319*u_0 + 214*u_p1 - 38*u_p2 + 4*u_p3 - 101*u_m1 + 25*u_m2 - 3*u_m3 ) / 420
    elseif idx == -2
        edge_u[1,cell] = ( 459*u_0 - 241*u_p1 + 109*u_p2 - 31*u_p3 + 4*u_p4 + 130*u_m1 - 10*u_m2 ) / 420
        edge_u[2,cell] = ( 214*u_0 + 319*u_p1 - 101*u_p2 + 25*u_p3 - 3*u_p4 - 38*u_m1 + 4*u_m2 ) / 420
    elseif idx == -1
        edge_u[1,cell] = ( 669*u_0 - 591*u_p1 + 459*u_p2 - 241*u_p3 + 74*u_p4 - 10*u_p5 + 60*u_m1 ) / 420
        edge_u[2,cell] = ( 130*u_0 + 459*u_p1 - 241*u_p2 + 109*u_p3 - 31*u_p4 + 4*u_p5 - 10*u_m1 ) / 420
    else # idx == 0
        edge_u[1,cell] = ( 1089*u_0 - 1851*u_p1 + 2559*u_p2 - 2341*u_p3 + 1334*u_p4 - 430*u_p5 + 60*u_p6 ) / 420
        edge_u[2,cell] = ( 60*u_0 + 669*u_p1 - 591*u_p2 + 459*u_p3 - 241*u_p4 + 74*u_p5 - 10*u_p6 ) / 420
    end
    nothing
end

function (eno::ModifiedENO{8})(edge_u, cell, u_m7, u_m6, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, u_p6, u_p7, balance_law, meshx)
    idx = eno.choose_stencil(u_m7, u_m6, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, u_p6, u_p7)

    @inbounds if idx == -7
        edge_u[1,cell] = ( 105*u_0 + 1443*u_m1 - 1497*u_m2 + 1443*u_m3 - 1007*u_m4 + 463*u_m5 - 125*u_m6 + 15*u_m7 ) / 840
        edge_u[2,cell] = ( 2283*u_0 - 4437*u_m1 + 7323*u_m2 - 8357*u_m3 + 6343*u_m4 - 3065*u_m5 + 855*u_m6 - 105*u_m7 ) / 840
    elseif idx == -6
        edge_u[1,cell] = ( 225*u_0 - 15*u_p1 + 1023*u_m1 - 657*u_m2 + 393*u_m3 - 167*u_m4 + 43*u_m5 - 5*u_m6 ) / 840
        edge_u[2,cell] = ( 1443*u_0 + 105*u_p1 - 1497*u_m1 + 1443*u_m2 - 1007*u_m3 + 463*u_m4 - 125*u_m5 + 15*u_m6 ) / 840
    elseif idx == -5
        edge_u[1,cell] = ( 365*u_0 - 55*u_p1 + 5*u_p2 + 743*u_m1 - 307*u_m2 + 113*u_m3 - 27*u_m4 + 3*u_m5 ) / 840
        edge_u[2,cell] = ( 1023*u_0 + 225*u_p1 - 15*u_p2 - 657*u_m1 + 393*u_m2 - 167*u_m3 + 43*u_m4 - 5*u_m5 ) / 840
    elseif idx == -4
        edge_u[1,cell] = ( 533*u_0 - 139*u_p1 + 29*u_p2 - 3*u_p3 + 533*u_m1 - 139*u_m2 + 29*u_m3 - 3*u_m4 ) / 840
        edge_u[2,cell] = ( 743*u_0 + 365*u_p1 - 55*u_p2 + 5*u_p3 - 307*u_m1 + 113*u_m2 - 27*u_m3 + 3*u_m4 ) / 840
    elseif idx == -3
        edge_u[1,cell] = ( 743*u_0 - 307*u_p1 + 113*u_p2 - 27*u_p3 + 3*u_p4 + 365*u_m1 - 55*u_m2 + 5*u_m3 ) / 840
        edge_u[2,cell] = ( 533*u_0 + 533*u_p1 - 139*u_p2 + 29*u_p3 - 3*u_p4 - 139*u_m1 + 29*u_m2 - 3*u_m3 ) / 840
    elseif idx == -2
        edge_u[1,cell] = ( 1023*u_0 - 657*u_p1 + 393*u_p2 - 167*u_p3 + 43*u_p4 - 5*u_p5 + 225*u_m1 - 15*u_m2 ) / 840
        edge_u[2,cell] = ( 365*u_0 + 743*u_p1 - 307*u_p2 + 113*u_p3 - 27*u_p4 + 3*u_p5 - 55*u_m1 + 5*u_m2 ) / 840
    elseif idx == -1
        edge_u[1,cell] = ( 1443*u_0 - 1497*u_p1 + 1443*u_p2 - 1007*u_p3 + 463*u_p4 - 125*u_p5 + 15*u_p6 + 105*u_m1 ) / 840
        edge_u[2,cell] = ( 225*u_0 + 1023*u_p1 - 657*u_p2 + 393*u_p3 - 167*u_p4 + 43*u_p5 - 5*u_p6 - 15*u_m1 ) / 840
    else # idx == 0
        edge_u[1,cell] = ( 2283*u_0 - 4437*u_p1 + 7323*u_p2 - 8357*u_p3 + 6343*u_p4 - 3065*u_p5 + 855*u_p6 - 105*u_p7 ) / 840
        edge_u[2,cell] = ( 105*u_0 + 1443*u_p1 - 1497*u_p2 + 1443*u_p3 - 1007*u_p4 + 463*u_p5 - 125*u_p6 + 15*u_p7 ) / 840
    end
    nothing
end


function interpolate!(uval, ξ, u_m1, u_0, u_p1, balance_law, mesh, eno::ModifiedENO{2})
    idx = eno.choose_stencil(u_m1, u_0, u_p1)

    if idx == -1
        for i in eachindex(uval)
            uval[i] = @evalpoly(ξ[i], u_0/2 + u_m1/2, u_0 - u_m1)
        end
    else # idx == 0
        for i in eachindex(uval)
            uval[i] = @evalpoly(ξ[i], 3*u_0/2 - u_p1/2, -u_0 + u_p1)
        end
    end
    nothing
end

function interpolate!(uval, ξ, u_m2, u_m1, u_0, u_p1, u_p2, balance_law, mesh, eno::ModifiedENO{3})
    idx = eno.choose_stencil(u_m2, u_m1, u_0, u_p1, u_p2)

    if idx == -2
        for i in eachindex(uval)
            uval[i] = @evalpoly(ξ[i], u_0/3 + 5*u_m1/6 - u_m2/6, u_0 - u_m1, u_0/2 - u_m1 + u_m2/2)
        end
    elseif idx == -1
        for i in eachindex(uval)
            uval[i] = @evalpoly(ξ[i], 5*u_0/6 - u_p1/6 + u_m1/3, u_0 - u_m1, -u_0 + u_p1/2 + u_m1/2)
        end
    else # idx == 0
        for i in eachindex(uval)
            uval[i] = @evalpoly(ξ[i], 11*u_0/6 - 7*u_p1/6 + u_p2/3, -2*u_0 + 3*u_p1 - u_p2, u_0/2 - u_p1 + u_p2/2)
        end
    end
    nothing
end

function interpolate!(uval, ξ, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, balance_law, mesh, eno::ModifiedENO{4})
    idx = eno.choose_stencil(u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3)

    if idx == -3
        for i in eachindex(uval)
            uval[i] = @evalpoly(ξ[i], u_0/4 + 13*u_m1/12 - 5*u_m2/12 + u_m3/12, 11*u_0/12 - 3*u_m1/4 - u_m2/4 + u_m3/12, 3*u_0/4 - 7*u_m1/4 + 5*u_m2/4 - u_m3/4, u_0/6 - u_m1/2 + u_m2/2 - u_m3/6)
        end
    elseif idx == -2
        for i in eachindex(uval)
            uval[i] = @evalpoly(ξ[i], 7*u_0/12 - u_p1/12 + 7*u_m1/12 - u_m2/12, 5*u_0/4 - u_p1/12 - 5*u_m1/4 + u_m2/12, -u_0/4 + u_p1/4 - u_m1/4 + u_m2/4, -u_0/2 + u_p1/6 + u_m1/2 - u_m2/6)
        end
    elseif idx == -1
        for i in eachindex(uval)
            uval[i] = @evalpoly(ξ[i], 13*u_0/12 - 5*u_p1/12 + u_p2/12 + u_m1/4, 3*u_0/4 + u_p1/4 - u_p2/12 - 11*u_m1/12, -7*u_0/4 + 5*u_p1/4 - u_p2/4 + 3*u_m1/4, u_0/2 - u_p1/2 + u_p2/6 - u_m1/6)
        end
    else # idx == 0
        for i in eachindex(uval)
            uval[i] = @evalpoly(ξ[i], 25*u_0/12 - 23*u_p1/12 + 13*u_p2/12 - u_p3/4, -35*u_0/12 + 23*u_p1/4 - 15*u_p2/4 + 11*u_p3/12, 5*u_0/4 - 13*u_p1/4 + 11*u_p2/4 - 3*u_p3/4, -u_0/6 + u_p1/2 - u_p2/2 + u_p3/6)
        end
    end
    nothing
end

function interpolate!(uval, ξ, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, balance_law, mesh, eno::ModifiedENO{5})
    idx = eno.choose_stencil(u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4)

    if idx == -4
        for i in eachindex(uval)
            uval[i] = @evalpoly(ξ[i], u_0/5 + 77*u_m1/60 - 43*u_m2/60 + 17*u_m3/60 - u_m4/20, 5*u_0/6 - 5*u_m1/12 - 3*u_m2/4 + 5*u_m3/12 - u_m4/12, 7*u_0/8 - 9*u_m1/4 + 2*u_m2 - 3*u_m3/4 + u_m4/8, u_0/3 - 7*u_m1/6 + 3*u_m2/2 - 5*u_m3/6 + u_m4/6, u_0/24 - u_m1/6 + u_m2/4 - u_m3/6 + u_m4/24)
        end
    elseif idx == -3
        for i in eachindex(uval)
            uval[i] = @evalpoly(ξ[i], 9*u_0/20 - u_p1/20 + 47*u_m1/60 - 13*u_m2/60 + u_m3/30, 5*u_0/4 - u_p1/12 - 5*u_m1/4 + u_m2/12, u_0/4 + u_p1/8 - u_m1 + 3*u_m2/4 - u_m3/8, -u_0/2 + u_p1/6 + u_m1/2 - u_m2/6, -u_0/6 + u_p1/24 + u_m1/4 - u_m2/6 + u_m3/24)
        end
    elseif idx == -2
        for i in eachindex(uval)
            uval[i] = @evalpoly(ξ[i], 47*u_0/60 - 13*u_p1/60 + u_p2/30 + 9*u_m1/20 - u_m2/20, 5*u_0/4 - u_p1/12 - 5*u_m1/4 + u_m2/12, -u_0 + 3*u_p1/4 - u_p2/8 + u_m1/4 + u_m2/8, -u_0/2 + u_p1/6 + u_m1/2 - u_m2/6, u_0/4 - u_p1/6 + u_p2/24 - u_m1/6 + u_m2/24)
        end
    elseif idx == -1
        for i in eachindex(uval)
            uval[i] = @evalpoly(ξ[i], 77*u_0/60 - 43*u_p1/60 + 17*u_p2/60 - u_p3/20 + u_m1/5, 5*u_0/12 + 3*u_p1/4 - 5*u_p2/12 + u_p3/12 - 5*u_m1/6, -9*u_0/4 + 2*u_p1 - 3*u_p2/4 + u_p3/8 + 7*u_m1/8, 7*u_0/6 - 3*u_p1/2 + 5*u_p2/6 - u_p3/6 - u_m1/3, -u_0/6 + u_p1/4 - u_p2/6 + u_p3/24 + u_m1/24)
        end
    else # idx == 0
        for i in eachindex(uval)
            uval[i] = @evalpoly(ξ[i], 137*u_0/60 - 163*u_p1/60 + 137*u_p2/60 - 21*u_p3/20 + u_p4/5, -15*u_0/4 + 109*u_p1/12 - 35*u_p2/4 + 17*u_p3/4 - 5*u_p4/6, 17*u_0/8 - 27*u_p1/4 + 8*u_p2 - 17*u_p3/4 + 7*u_p4/8, -u_0/2 + 11*u_p1/6 - 5*u_p2/2 + 3*u_p3/2 - u_p4/3, u_0/24 - u_p1/6 + u_p2/4 - u_p3/6 + u_p4/24)
        end
    end
    nothing
end

function interpolate!(uval, ξ, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, balance_law, mesh, eno::ModifiedENO{6})
    idx = eno.choose_stencil(u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5)

    if idx == -5
        for i in eachindex(uval)
            uval[i] = @evalpoly(ξ[i], u_0/6 + 29*u_m1/20 - 21*u_m2/20 + 37*u_m3/60 - 13*u_m4/60 + u_m5/30, 137*u_0/180 - u_m1/18 - 53*u_m2/36 + 41*u_m3/36 - 4*u_m4/9 + 13*u_m5/180, 15*u_0/16 - 41*u_m1/16 + 21*u_m2/8 - 11*u_m3/8 + 7*u_m4/16 - u_m5/16, 17*u_0/36 - 67*u_m1/36 + 26*u_m2/9 - 20*u_m3/9 + 31*u_m4/36 - 5*u_m5/36, 5*u_0/48 - 23*u_m1/48 + 7*u_m2/8 - 19*u_m3/24 + 17*u_m4/48 - u_m5/16, u_0/120 - u_m1/24 + u_m2/12 - u_m3/12 + u_m4/24 - u_m5/120)
        end
    elseif idx == -4
        for i in eachindex(uval)
            uval[i] = @evalpoly(ξ[i], 11*u_0/30 - u_p1/30 + 19*u_m1/20 - 23*u_m2/60 + 7*u_m3/60 - u_m4/60, 43*u_0/36 - 13*u_p1/180 - 41*u_m1/36 - u_m2/36 + u_m3/18 - u_m4/90, 9*u_0/16 + u_p1/16 - 13*u_m1/8 + 11*u_m2/8 - 7*u_m3/16 + u_m4/16, -13*u_0/36 + 5*u_p1/36 + 2*u_m1/9 + u_m2/9 - 5*u_m3/36 + u_m4/36, -13*u_0/48 + u_p1/16 + 11*u_m1/24 - 3*u_m2/8 + 7*u_m3/48 - u_m4/48, -u_0/24 + u_p1/120 + u_m1/12 - u_m2/12 + u_m3/24 - u_m4/120)
        end
    elseif idx == -3
        for i in eachindex(uval)
            uval[i] = @evalpoly(ξ[i], 37*u_0/60 - 2*u_p1/15 + u_p2/60 + 37*u_m1/60 - 2*u_m2/15 + u_m3/60, 49*u_0/36 - 5*u_p1/36 + u_p2/90 - 49*u_m1/36 + 5*u_m2/36 - u_m3/90, -3*u_0/8 + 7*u_p1/16 - u_p2/16 - 3*u_m1/8 + 7*u_m2/16 - u_m3/16, -7*u_0/9 + 11*u_p1/36 - u_p2/36 + 7*u_m1/9 - 11*u_m2/36 + u_m3/36, u_0/24 - u_p1/16 + u_p2/48 + u_m1/24 - u_m2/16 + u_m3/48, u_0/12 - u_p1/24 + u_p2/120 - u_m1/12 + u_m2/24 - u_m3/120)
        end
    elseif idx == -2
        for i in eachindex(uval)
            uval[i] = @evalpoly(ξ[i], 19*u_0/20 - 23*u_p1/60 + 7*u_p2/60 - u_p3/60 + 11*u_m1/30 - u_m2/30, 41*u_0/36 + u_p1/36 - u_p2/18 + u_p3/90 - 43*u_m1/36 + 13*u_m2/180, -13*u_0/8 + 11*u_p1/8 - 7*u_p2/16 + u_p3/16 + 9*u_m1/16 + u_m2/16, -2*u_0/9 - u_p1/9 + 5*u_p2/36 - u_p3/36 + 13*u_m1/36 - 5*u_m2/36, 11*u_0/24 - 3*u_p1/8 + 7*u_p2/48 - u_p3/48 - 13*u_m1/48 + u_m2/16, -u_0/12 + u_p1/12 - u_p2/24 + u_p3/120 + u_m1/24 - u_m2/120)
        end
    elseif idx == -1
        for i in eachindex(uval)
            uval[i] = @evalpoly(ξ[i], 29*u_0/20 - 21*u_p1/20 + 37*u_p2/60 - 13*u_p3/60 + u_p4/30 + u_m1/6, u_0/18 + 53*u_p1/36 - 41*u_p2/36 + 4*u_p3/9 - 13*u_p4/180 - 137*u_m1/180, -41*u_0/16 + 21*u_p1/8 - 11*u_p2/8 + 7*u_p3/16 - u_p4/16 + 15*u_m1/16, 67*u_0/36 - 26*u_p1/9 + 20*u_p2/9 - 31*u_p3/36 + 5*u_p4/36 - 17*u_m1/36, -23*u_0/48 + 7*u_p1/8 - 19*u_p2/24 + 17*u_p3/48 - u_p4/16 + 5*u_m1/48, u_0/24 - u_p1/12 + u_p2/12 - u_p3/24 + u_p4/120 - u_m1/120)
        end
    else # idx == 0
        for i in eachindex(uval)
            uval[i] = @evalpoly(ξ[i], 49*u_0/20 - 71*u_p1/20 + 79*u_p2/20 - 163*u_p3/60 + 31*u_p4/30 - u_p5/6, -203*u_0/45 + 116*u_p1/9 - 589*u_p2/36 + 427*u_p3/36 - 167*u_p4/36 + 137*u_p5/180, 49*u_0/16 - 183*u_p1/16 + 139*u_p2/8 - 109*u_p3/8 + 89*u_p4/16 - 15*u_p5/16, -35*u_0/36 + 151*u_p1/36 - 65*u_p2/9 + 56*u_p3/9 - 97*u_p4/36 + 17*u_p5/36, 7*u_0/48 - 11*u_p1/16 + 31*u_p2/24 - 29*u_p3/24 + 9*u_p4/16 - 5*u_p5/48, -u_0/120 + u_p1/24 - u_p2/12 + u_p3/12 - u_p4/24 + u_p5/120)
        end
    end
    nothing
end

function interpolate!(uval, ξ, u_m6, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, u_p6, balance_law, mesh, eno::ModifiedENO{7})
    idx = eno.choose_stencil(u_m6, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, u_p6)

    if idx == -6
        for i in eachindex(uval)
            uval[i] = @evalpoly(ξ[i], u_0/7 + 223*u_m1/140 - 197*u_m2/140 + 153*u_m3/140 - 241*u_m4/420 + 37*u_m5/210 - u_m6/42, 7*u_0/10 + 14*u_m1/45 - 43*u_m2/18 + 85*u_m3/36 - 49*u_m4/36 + 79*u_m5/180 - 11*u_m6/180, 29*u_0/30 - 219*u_m1/80 + 49*u_m2/16 - 47*u_m3/24 + 7*u_m4/8 - 19*u_m5/80 + 7*u_m6/240, 7*u_0/12 - 91*u_m1/36 + 41*u_m2/9 - 40*u_m3/9 + 91*u_m4/36 - 29*u_m5/36 + u_m6/9, 25*u_0/144 - 43*u_m1/48 + 23*u_m2/12 - 157*u_m3/72 + 67*u_m4/48 - 23*u_m5/48 + 5*u_m6/72, u_0/40 - 17*u_m1/120 + u_m2/3 - 5*u_m3/12 + 7*u_m4/24 - 13*u_m5/120 + u_m6/60, u_0/720 - u_m1/120 + u_m2/48 - u_m3/36 + u_m4/48 - u_m5/120 + u_m6/720)
        end
    elseif idx == -5
        for i in eachindex(uval)
            uval[i] = @evalpoly(ξ[i], 13*u_0/42 - u_p1/42 + 153*u_m1/140 - 241*u_m2/420 + 109*u_m3/420 - 31*u_m4/420 + u_m5/105, 203*u_0/180 - 11*u_p1/180 - 35*u_m1/36 - u_m2/4 + 2*u_m3/9 - 7*u_m4/90 + u_m5/90, 61*u_0/80 + 7*u_p1/240 - 17*u_m1/8 + 49*u_m2/24 - 15*u_m3/16 + 21*u_m4/80 - u_m5/30, -7*u_0/36 + u_p1/9 - 7*u_m1/36 + 2*u_m2/3 - 5*u_m3/9 + 7*u_m4/36 - u_m5/36, -5*u_0/16 + 5*u_p1/72 + 9*u_m1/16 - 37*u_m2/72 + u_m3/4 - u_m4/16 + u_m5/144, -11*u_0/120 + u_p1/60 + 5*u_m1/24 - u_m2/4 + u_m3/6 - 7*u_m4/120 + u_m5/120, -u_0/120 + u_p1/720 + u_m1/48 - u_m2/36 + u_m3/48 - u_m4/120 + u_m5/720)
        end
    elseif idx == -4
        for i in eachindex(uval)
            uval[i] = @evalpoly(ξ[i], 107*u_0/210 - 19*u_p1/210 + u_p2/105 + 319*u_m1/420 - 101*u_m2/420 + 5*u_m3/84 - u_m4/140, 49*u_0/36 - 5*u_p1/36 + u_p2/90 - 49*u_m1/36 + 5*u_m2/36 - u_m3/90, u_0/16 + 21*u_p1/80 - u_p2/30 - 23*u_m1/24 + 7*u_m2/8 - 19*u_m3/80 + 7*u_m4/240, -7*u_0/9 + 11*u_p1/36 - u_p2/36 + 7*u_m1/9 - 11*u_m2/36 + u_m3/36, -u_0/6 + u_p1/48 + u_p2/144 + 23*u_m1/72 - 13*u_m2/48 + 5*u_m3/48 - u_m4/72, u_0/12 - u_p1/24 + u_p2/120 - u_m1/12 + u_m2/24 - u_m3/120, u_0/48 - u_p1/120 + u_p2/720 - u_m1/36 + u_m2/48 - u_m3/120 + u_m4/720)
        end
    elseif idx == -3
        for i in eachindex(uval)
            uval[i] = @evalpoly(ξ[i], 319*u_0/420 - 101*u_p1/420 + 5*u_p2/84 - u_p3/140 + 107*u_m1/210 - 19*u_m2/210 + u_m3/105, 49*u_0/36 - 5*u_p1/36 + u_p2/90 - 49*u_m1/36 + 5*u_m2/36 - u_m3/90, -23*u_0/24 + 7*u_p1/8 - 19*u_p2/80 + 7*u_p3/240 + u_m1/16 + 21*u_m2/80 - u_m3/30, -7*u_0/9 + 11*u_p1/36 - u_p2/36 + 7*u_m1/9 - 11*u_m2/36 + u_m3/36, 23*u_0/72 - 13*u_p1/48 + 5*u_p2/48 - u_p3/72 - u_m1/6 + u_m2/48 + u_m3/144, u_0/12 - u_p1/24 + u_p2/120 - u_m1/12 + u_m2/24 - u_m3/120, -u_0/36 + u_p1/48 - u_p2/120 + u_p3/720 + u_m1/48 - u_m2/120 + u_m3/720)
        end
    elseif idx == -2
        for i in eachindex(uval)
            uval[i] = @evalpoly(ξ[i], 153*u_0/140 - 241*u_p1/420 + 109*u_p2/420 - 31*u_p3/420 + u_p4/105 + 13*u_m1/42 - u_m2/42, 35*u_0/36 + u_p1/4 - 2*u_p2/9 + 7*u_p3/90 - u_p4/90 - 203*u_m1/180 + 11*u_m2/180, -17*u_0/8 + 49*u_p1/24 - 15*u_p2/16 + 21*u_p3/80 - u_p4/30 + 61*u_m1/80 + 7*u_m2/240, 7*u_0/36 - 2*u_p1/3 + 5*u_p2/9 - 7*u_p3/36 + u_p4/36 + 7*u_m1/36 - u_m2/9, 9*u_0/16 - 37*u_p1/72 + u_p2/4 - u_p3/16 + u_p4/144 - 5*u_m1/16 + 5*u_m2/72, -5*u_0/24 + u_p1/4 - u_p2/6 + 7*u_p3/120 - u_p4/120 + 11*u_m1/120 - u_m2/60, u_0/48 - u_p1/36 + u_p2/48 - u_p3/120 + u_p4/720 - u_m1/120 + u_m2/720)
        end
    elseif idx == -1
        for i in eachindex(uval)
            uval[i] = @evalpoly(ξ[i], 223*u_0/140 - 197*u_p1/140 + 153*u_p2/140 - 241*u_p3/420 + 37*u_p4/210 - u_p5/42 + u_m1/7, -14*u_0/45 + 43*u_p1/18 - 85*u_p2/36 + 49*u_p3/36 - 79*u_p4/180 + 11*u_p5/180 - 7*u_m1/10, -219*u_0/80 + 49*u_p1/16 - 47*u_p2/24 + 7*u_p3/8 - 19*u_p4/80 + 7*u_p5/240 + 29*u_m1/30, 91*u_0/36 - 41*u_p1/9 + 40*u_p2/9 - 91*u_p3/36 + 29*u_p4/36 - u_p5/9 - 7*u_m1/12, -43*u_0/48 + 23*u_p1/12 - 157*u_p2/72 + 67*u_p3/48 - 23*u_p4/48 + 5*u_p5/72 + 25*u_m1/144, 17*u_0/120 - u_p1/3 + 5*u_p2/12 - 7*u_p3/24 + 13*u_p4/120 - u_p5/60 - u_m1/40, -u_0/120 + u_p1/48 - u_p2/36 + u_p3/48 - u_p4/120 + u_p5/720 + u_m1/720)
        end
    else # idx == 0
        for i in eachindex(uval)
            uval[i] = @evalpoly(ξ[i], 363*u_0/140 - 617*u_p1/140 + 853*u_p2/140 - 2341*u_p3/420 + 667*u_p4/210 - 43*u_p5/42 + u_p6/7, -469*u_0/90 + 769*u_p1/45 - 967*u_p2/36 + 931*u_p3/36 - 545*u_p4/36 + 893*u_p5/180 - 7*u_p6/10, 967*u_0/240 - 1379*u_p1/80 + 255*u_p2/8 - 791*u_p3/24 + 321*u_p4/16 - 539*u_p5/80 + 29*u_p6/30, -14*u_0/9 + 277*u_p1/36 - 575*u_p2/36 + 161*u_p3/9 - 103*u_p4/9 + 143*u_p5/36 - 7*u_p6/12, 23*u_0/72 - 83*u_p1/48 + 187*u_p2/48 - 337*u_p3/72 + 19*u_p4/6 - 55*u_p5/48 + 25*u_p6/144, -u_0/30 + 23*u_p1/120 - 11*u_p2/24 + 7*u_p3/12 - 5*u_p4/12 + 19*u_p5/120 - u_p6/40, u_0/720 - u_p1/120 + u_p2/48 - u_p3/36 + u_p4/48 - u_p5/120 + u_p6/720)
        end
    end
    nothing
end

function interpolate!(uval, ξ, u_m7, u_m6, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, u_p6, u_p7, balance_law, mesh, eno::ModifiedENO{8})
    idx = eno.choose_stencil(u_m7, u_m6, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, u_p6, u_p7)

    if idx == -7
        for i in eachindex(uval)
            uval[i] = @evalpoly(ξ[i], u_0/8 + 481*u_m1/280 - 499*u_m2/280 + 481*u_m3/280 - 1007*u_m4/840 + 463*u_m5/840 - 25*u_m6/168 + u_m7/56, 363*u_0/560 + 97*u_m1/144 - 2503*u_m2/720 + 601*u_m3/144 - 457*u_m4/144 + 1099*u_m5/720 - 61*u_m6/144 + 29*u_m7/560, 469*u_0/480 - 1349*u_m1/480 + 105*u_m2/32 - 223*u_m3/96 + 119*u_m4/96 - 73*u_m5/160 + 49*u_m6/480 - u_m7/96, 967*u_0/1440 - 4529*u_m1/1440 + 9227*u_m2/1440 - 241*u_m3/32 + 539*u_m4/96 - 3827*u_m5/1440 + 1049*u_m6/1440 - 127*u_m7/1440, 35*u_0/144 - 199*u_m1/144 + 27*u_m2/8 - 83*u_m3/18 + 551*u_m4/144 - 31*u_m5/16 + 5*u_m6/9 - 5*u_m7/72, 23*u_0/480 - 29*u_m1/96 + 391*u_m2/480 - 39*u_m3/32 + 35*u_m4/32 - 283*u_m5/480 + 17*u_m6/96 - 11*u_m7/480, 7*u_0/1440 - 47*u_m1/1440 + 3*u_m2/32 - 43*u_m3/288 + 41*u_m4/288 - 13*u_m5/160 + 37*u_m6/1440 - u_m7/288, u_0/5040 - u_m1/720 + u_m2/240 - u_m3/144 + u_m4/144 - u_m5/240 + u_m6/720 - u_m7/5040)
        end
    elseif idx == -6
        for i in eachindex(uval)
            uval[i] = @evalpoly(ξ[i], 15*u_0/56 - u_p1/56 + 341*u_m1/280 - 219*u_m2/280 + 131*u_m3/280 - 167*u_m4/840 + 43*u_m5/840 - u_m6/168, 17*u_0/16 - 29*u_p1/560 - 559*u_m1/720 - 83*u_m2/144 + 79*u_m3/144 - 197*u_m4/720 + 11*u_m5/144 - 47*u_m6/5040, 143*u_0/160 + u_p1/96 - 403*u_m1/160 + 259*u_m2/96 - 51*u_m3/32 + 21*u_m4/32 - 79*u_m5/480 + 3*u_m6/160, -49*u_0/1440 + 127*u_p1/1440 - 973*u_m1/1440 + 47*u_m2/32 - 391*u_m3/288 + 973*u_m4/1440 - 271*u_m5/1440 + 11*u_m6/480, -5*u_0/16 + 5*u_p1/72 + 9*u_m1/16 - 37*u_m2/72 + u_m3/4 - u_m4/16 + u_m5/144, -13*u_0/96 + 11*u_p1/480 + 163*u_m1/480 - 15*u_m2/32 + 37*u_m3/96 - 91*u_m4/480 + 5*u_m5/96 - u_m6/160, -11*u_0/480 + u_p1/288 + 31*u_m1/480 - 29*u_m2/288 + 3*u_m3/32 - 5*u_m4/96 + 23*u_m5/1440 - u_m6/480, -u_0/720 + u_p1/5040 + u_m1/240 - u_m2/144 + u_m3/144 - u_m4/240 + u_m5/720 - u_m6/5040)
        end
    elseif idx == -5
        for i in eachindex(uval)
            uval[i] = @evalpoly(ξ[i], 73*u_0/168 - 11*u_p1/168 + u_p2/168 + 743*u_m1/840 - 307*u_m2/840 + 113*u_m3/840 - 9*u_m4/280 + u_m5/280, 953*u_0/720 - 91*u_p1/720 + 47*u_p2/5040 - 187*u_m1/144 + 11*u_m2/144 + 19*u_m3/720 - u_m4/80 + u_m5/560, 59*u_0/160 + 77*u_p1/480 - 3*u_p2/160 - 47*u_m1/32 + 133*u_m2/96 - 87*u_m3/160 + 21*u_m4/160 - 7*u_m5/480, -973*u_0/1440 + 391*u_p1/1440 - 11*u_p2/480 + 175*u_m1/288 - 13*u_m2/96 - 107*u_m3/1440 + 49*u_m4/1440 - 7*u_m5/1440, -5*u_0/16 + 5*u_p1/72 + 9*u_m1/16 - 37*u_m2/72 + u_m3/4 - u_m4/16 + u_m5/144, 19*u_0/480 - 13*u_p1/480 + u_p2/160 - u_m1/96 - u_m2/32 + 17*u_m3/480 - 7*u_m4/480 + u_m5/480, 17*u_0/480 - 19*u_p1/1440 + u_p2/480 - 5*u_m1/96 + 13*u_m2/288 - 11*u_m3/480 + u_m4/160 - u_m5/1440, u_0/240 - u_p1/720 + u_p2/5040 - u_m1/144 + u_m2/144 - u_m3/240 + u_m4/720 - u_m5/5040)
        end
    elseif idx == -4
        for i in eachindex(uval)
            uval[i] = @evalpoly(ξ[i], 533*u_0/840 - 139*u_p1/840 + 29*u_p2/840 - u_p3/280 + 533*u_m1/840 - 139*u_m2/840 + 29*u_m3/840 - u_m4/280, 205*u_0/144 - 127*u_p1/720 + 17*u_p2/720 - u_p3/560 - 205*u_m1/144 + 127*u_m2/720 - 17*u_m3/720 + u_m4/560, -43*u_0/96 + 91*u_p1/160 - 13*u_p2/96 + 7*u_p3/480 - 43*u_m1/96 + 91*u_m2/160 - 13*u_m3/96 + 7*u_m4/480, -91*u_0/96 + 587*u_p1/1440 - 89*u_p2/1440 + 7*u_p3/1440 + 91*u_m1/96 - 587*u_m2/1440 + 89*u_m3/1440 - 7*u_m4/1440, 11*u_0/144 - u_p1/8 + u_p2/18 - u_p3/144 + 11*u_m1/144 - u_m2/8 + u_m3/18 - u_m4/144, 5*u_0/32 - 41*u_p1/480 + 11*u_p2/480 - u_p3/480 - 5*u_m1/32 + 41*u_m2/480 - 11*u_m3/480 + u_m4/480, -u_0/288 + u_p1/160 - u_p2/288 + u_p3/1440 - u_m1/288 + u_m2/160 - u_m3/288 + u_m4/1440, -u_0/144 + u_p1/240 - u_p2/720 + u_p3/5040 + u_m1/144 - u_m2/240 + u_m3/720 - u_m4/5040)
        end
    elseif idx == -3
        for i in eachindex(uval)
            uval[i] = @evalpoly(ξ[i], 743*u_0/840 - 307*u_p1/840 + 113*u_p2/840 - 9*u_p3/280 + u_p4/280 + 73*u_m1/168 - 11*u_m2/168 + u_m3/168, 187*u_0/144 - 11*u_p1/144 - 19*u_p2/720 + u_p3/80 - u_p4/560 - 953*u_m1/720 + 91*u_m2/720 - 47*u_m3/5040, -47*u_0/32 + 133*u_p1/96 - 87*u_p2/160 + 21*u_p3/160 - 7*u_p4/480 + 59*u_m1/160 + 77*u_m2/480 - 3*u_m3/160, -175*u_0/288 + 13*u_p1/96 + 107*u_p2/1440 - 49*u_p3/1440 + 7*u_p4/1440 + 973*u_m1/1440 - 391*u_m2/1440 + 11*u_m3/480, 9*u_0/16 - 37*u_p1/72 + u_p2/4 - u_p3/16 + u_p4/144 - 5*u_m1/16 + 5*u_m2/72, u_0/96 + u_p1/32 - 17*u_p2/480 + 7*u_p3/480 - u_p4/480 - 19*u_m1/480 + 13*u_m2/480 - u_m3/160, -5*u_0/96 + 13*u_p1/288 - 11*u_p2/480 + u_p3/160 - u_p4/1440 + 17*u_m1/480 - 19*u_m2/1440 + u_m3/480, u_0/144 - u_p1/144 + u_p2/240 - u_p3/720 + u_p4/5040 - u_m1/240 + u_m2/720 - u_m3/5040)
        end
    elseif idx == -2
        for i in eachindex(uval)
            uval[i] = @evalpoly(ξ[i], 341*u_0/280 - 219*u_p1/280 + 131*u_p2/280 - 167*u_p3/840 + 43*u_p4/840 - u_p5/168 + 15*u_m1/56 - u_m2/56, 559*u_0/720 + 83*u_p1/144 - 79*u_p2/144 + 197*u_p3/720 - 11*u_p4/144 + 47*u_p5/5040 - 17*u_m1/16 + 29*u_m2/560, -403*u_0/160 + 259*u_p1/96 - 51*u_p2/32 + 21*u_p3/32 - 79*u_p4/480 + 3*u_p5/160 + 143*u_m1/160 + u_m2/96, 973*u_0/1440 - 47*u_p1/32 + 391*u_p2/288 - 973*u_p3/1440 + 271*u_p4/1440 - 11*u_p5/480 + 49*u_m1/1440 - 127*u_m2/1440, 9*u_0/16 - 37*u_p1/72 + u_p2/4 - u_p3/16 + u_p4/144 - 5*u_m1/16 + 5*u_m2/72, -163*u_0/480 + 15*u_p1/32 - 37*u_p2/96 + 91*u_p3/480 - 5*u_p4/96 + u_p5/160 + 13*u_m1/96 - 11*u_m2/480, 31*u_0/480 - 29*u_p1/288 + 3*u_p2/32 - 5*u_p3/96 + 23*u_p4/1440 - u_p5/480 - 11*u_m1/480 + u_m2/288, -u_0/240 + u_p1/144 - u_p2/144 + u_p3/240 - u_p4/720 + u_p5/5040 + u_m1/720 - u_m2/5040)
        end
    elseif idx == -1
        for i in eachindex(uval)
            uval[i] = @evalpoly(ξ[i], 481*u_0/280 - 499*u_p1/280 + 481*u_p2/280 - 1007*u_p3/840 + 463*u_p4/840 - 25*u_p5/168 + u_p6/56 + u_m1/8, -97*u_0/144 + 2503*u_p1/720 - 601*u_p2/144 + 457*u_p3/144 - 1099*u_p4/720 + 61*u_p5/144 - 29*u_p6/560 - 363*u_m1/560, -1349*u_0/480 + 105*u_p1/32 - 223*u_p2/96 + 119*u_p3/96 - 73*u_p4/160 + 49*u_p5/480 - u_p6/96 + 469*u_m1/480, 4529*u_0/1440 - 9227*u_p1/1440 + 241*u_p2/32 - 539*u_p3/96 + 3827*u_p4/1440 - 1049*u_p5/1440 + 127*u_p6/1440 - 967*u_m1/1440, -199*u_0/144 + 27*u_p1/8 - 83*u_p2/18 + 551*u_p3/144 - 31*u_p4/16 + 5*u_p5/9 - 5*u_p6/72 + 35*u_m1/144, 29*u_0/96 - 391*u_p1/480 + 39*u_p2/32 - 35*u_p3/32 + 283*u_p4/480 - 17*u_p5/96 + 11*u_p6/480 - 23*u_m1/480, -47*u_0/1440 + 3*u_p1/32 - 43*u_p2/288 + 41*u_p3/288 - 13*u_p4/160 + 37*u_p5/1440 - u_p6/288 + 7*u_m1/1440, u_0/720 - u_p1/240 + u_p2/144 - u_p3/144 + u_p4/240 - u_p5/720 + u_p6/5040 - u_m1/5040)
        end
    else # idx == 0
        for i in eachindex(uval)
            uval[i] = @evalpoly(ξ[i], 761*u_0/280 - 1479*u_p1/280 + 2441*u_p2/280 - 8357*u_p3/840 + 6343*u_p4/840 - 613*u_p5/168 + 57*u_p6/56 - u_p7/8, -29531*u_0/5040 + 15571*u_p1/720 - 29141*u_p2/720 + 6991*u_p3/144 - 5447*u_p4/144 + 13373*u_p5/720 - 419*u_p6/80 + 363*u_p7/560, 801*u_0/160 - 11557*u_p1/480 + 8383*u_p2/160 - 2149*u_p3/32 + 5209*u_p4/96 - 4361*u_p5/160 + 1249*u_p6/160 - 469*u_p7/480, -1069*u_0/480 + 17849*u_p1/1440 - 43307*u_p2/1440 + 11921*u_p3/288 - 3355*u_p4/96 + 26027*u_p5/1440 - 7609*u_p6/1440 + 967*u_p7/1440, 9*u_0/16 - 247*u_p1/72 + 9*u_p2 - 211*u_p3/16 + 1681*u_p4/144 - 25*u_p5/4 + 15*u_p6/8 - 35*u_p7/144, -13*u_0/160 + 253*u_p1/480 - 703*u_p2/480 + 217*u_p3/96 - 67*u_p4/32 + 559*u_p5/480 - 173*u_p6/480 + 23*u_p7/480, u_0/160 - 61*u_p1/1440 + 59*u_p2/480 - 19*u_p3/96 + 55*u_p4/288 - 53*u_p5/480 + 17*u_p6/480 - 7*u_p7/1440, -u_0/5040 + u_p1/720 - u_p2/240 + u_p3/144 - u_p4/144 + u_p5/240 - u_p6/720 + u_p7/5040)
        end
    end
    nothing
end


################################################################################

"""
    ClassicalChoiceENO

Classical choice of the stencil in an ENO reconstruction, minimising the absolute
values of the divided differences lexicographically.
"""
struct ClassicalChoiceENO end

@inline function (::ClassicalChoiceENO)(u_m1, u_0, u_p1)
    if abs(-u_m1 + u_0) < abs(-u_0 + u_p1)
        idx = -1
    else
        idx = 0
    end

    idx
end

@inline function (::ClassicalChoiceENO)(u_m2, u_m1, u_0, u_p1, u_p2)
    if abs(-u_m1 + u_0) < abs(-u_0 + u_p1)
        if abs(-2*u_m1 + u_m2 + u_0) < abs(u_m1 - 2*u_0 + u_p1)
            idx = -2
        else
            idx = -1
        end
    else
        if abs(u_m1 - 2*u_0 + u_p1) < abs(u_0 - 2*u_p1 + u_p2)
            idx = -1
        else
            idx = 0
        end
    end

    idx
end

@inline function (::ClassicalChoiceENO)(u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3)
    if abs(-u_m1 + u_0) < abs(-u_0 + u_p1)
        if abs(-2*u_m1 + u_m2 + u_0) < abs(u_m1 - 2*u_0 + u_p1)
            if abs(-3*u_m1 + 3*u_m2 - u_m3 + u_0) < abs(3*u_m1 - u_m2 - 3*u_0 + u_p1)
                idx = -3
            else
                idx = -2
            end
        else
            if abs(3*u_m1 - u_m2 - 3*u_0 + u_p1) < abs(-u_m1 + 3*u_0 - 3*u_p1 + u_p2)
                idx = -2
            else
                idx = -1
            end
        end
    else
        if abs(u_m1 - 2*u_0 + u_p1) < abs(u_0 - 2*u_p1 + u_p2)
            if abs(3*u_m1 - u_m2 - 3*u_0 + u_p1) < abs(-u_m1 + 3*u_0 - 3*u_p1 + u_p2)
                idx = -2
            else
                idx = -1
            end
        else
            if abs(-u_m1 + 3*u_0 - 3*u_p1 + u_p2) < abs(-u_0 + 3*u_p1 - 3*u_p2 + u_p3)
                idx = -1
            else
                idx = 0
            end
        end
    end

    idx
end

@inline function (::ClassicalChoiceENO)(u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4)
    if abs(-u_m1 + u_0) < abs(-u_0 + u_p1)
        if abs(-2*u_m1 + u_m2 + u_0) < abs(u_m1 - 2*u_0 + u_p1)
            if abs(-3*u_m1 + 3*u_m2 - u_m3 + u_0) < abs(3*u_m1 - u_m2 - 3*u_0 + u_p1)
                if abs(-4*u_m1 + 6*u_m2 - 4*u_m3 + u_m4 + u_0) < abs(6*u_m1 - 4*u_m2 + u_m3 - 4*u_0 + u_p1)
                    idx = -4
                else
                    idx = -3
                end
            else
                if abs(6*u_m1 - 4*u_m2 + u_m3 - 4*u_0 + u_p1) < abs(-4*u_m1 + u_m2 + 6*u_0 - 4*u_p1 + u_p2)
                    idx = -3
                else
                    idx = -2
                end
            end
        else
            if abs(3*u_m1 - u_m2 - 3*u_0 + u_p1) < abs(-u_m1 + 3*u_0 - 3*u_p1 + u_p2)
                if abs(6*u_m1 - 4*u_m2 + u_m3 - 4*u_0 + u_p1) < abs(-4*u_m1 + u_m2 + 6*u_0 - 4*u_p1 + u_p2)
                    idx = -3
                else
                    idx = -2
                end
            else
                if abs(-4*u_m1 + u_m2 + 6*u_0 - 4*u_p1 + u_p2) < abs(u_m1 - 4*u_0 + 6*u_p1 - 4*u_p2 + u_p3)
                    idx = -2
                else
                    idx = -1
                end
            end
        end
    else
        if abs(u_m1 - 2*u_0 + u_p1) < abs(u_0 - 2*u_p1 + u_p2)
            if abs(3*u_m1 - u_m2 - 3*u_0 + u_p1) < abs(-u_m1 + 3*u_0 - 3*u_p1 + u_p2)
                if abs(6*u_m1 - 4*u_m2 + u_m3 - 4*u_0 + u_p1) < abs(-4*u_m1 + u_m2 + 6*u_0 - 4*u_p1 + u_p2)
                    idx = -3
                else
                    idx = -2
                end
            else
                if abs(-4*u_m1 + u_m2 + 6*u_0 - 4*u_p1 + u_p2) < abs(u_m1 - 4*u_0 + 6*u_p1 - 4*u_p2 + u_p3)
                    idx = -2
                else
                    idx = -1
                end
            end
        else
            if abs(-u_m1 + 3*u_0 - 3*u_p1 + u_p2) < abs(-u_0 + 3*u_p1 - 3*u_p2 + u_p3)
                if abs(-4*u_m1 + u_m2 + 6*u_0 - 4*u_p1 + u_p2) < abs(u_m1 - 4*u_0 + 6*u_p1 - 4*u_p2 + u_p3)
                    idx = -2
                else
                    idx = -1
                end
            else
                if abs(u_m1 - 4*u_0 + 6*u_p1 - 4*u_p2 + u_p3) < abs(u_0 - 4*u_p1 + 6*u_p2 - 4*u_p3 + u_p4)
                    idx = -1
                else
                    idx = 0
                end
            end
        end
    end

    idx
end

@inline function (::ClassicalChoiceENO)(u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5)
    if abs(-u_m1 + u_0) < abs(-u_0 + u_p1)
        if abs(-2*u_m1 + u_m2 + u_0) < abs(u_m1 - 2*u_0 + u_p1)
            if abs(-3*u_m1 + 3*u_m2 - u_m3 + u_0) < abs(3*u_m1 - u_m2 - 3*u_0 + u_p1)
                if abs(-4*u_m1 + 6*u_m2 - 4*u_m3 + u_m4 + u_0) < abs(6*u_m1 - 4*u_m2 + u_m3 - 4*u_0 + u_p1)
                    if abs(-5*u_m1 + 10*u_m2 - 10*u_m3 + 5*u_m4 - u_m5 + u_0) < abs(10*u_m1 - 10*u_m2 + 5*u_m3 - u_m4 - 5*u_0 + u_p1)
                        idx = -5
                    else
                        idx = -4
                    end
                else
                    if abs(10*u_m1 - 10*u_m2 + 5*u_m3 - u_m4 - 5*u_0 + u_p1) < abs(-10*u_m1 + 5*u_m2 - u_m3 + 10*u_0 - 5*u_p1 + u_p2)
                        idx = -4
                    else
                        idx = -3
                    end
                end
            else
                if abs(6*u_m1 - 4*u_m2 + u_m3 - 4*u_0 + u_p1) < abs(-4*u_m1 + u_m2 + 6*u_0 - 4*u_p1 + u_p2)
                    if abs(10*u_m1 - 10*u_m2 + 5*u_m3 - u_m4 - 5*u_0 + u_p1) < abs(-10*u_m1 + 5*u_m2 - u_m3 + 10*u_0 - 5*u_p1 + u_p2)
                        idx = -4
                    else
                        idx = -3
                    end
                else
                    if abs(-10*u_m1 + 5*u_m2 - u_m3 + 10*u_0 - 5*u_p1 + u_p2) < abs(5*u_m1 - u_m2 - 10*u_0 + 10*u_p1 - 5*u_p2 + u_p3)
                        idx = -3
                    else
                        idx = -2
                    end
                end
            end
        else
            if abs(3*u_m1 - u_m2 - 3*u_0 + u_p1) < abs(-u_m1 + 3*u_0 - 3*u_p1 + u_p2)
                if abs(6*u_m1 - 4*u_m2 + u_m3 - 4*u_0 + u_p1) < abs(-4*u_m1 + u_m2 + 6*u_0 - 4*u_p1 + u_p2)
                    if abs(10*u_m1 - 10*u_m2 + 5*u_m3 - u_m4 - 5*u_0 + u_p1) < abs(-10*u_m1 + 5*u_m2 - u_m3 + 10*u_0 - 5*u_p1 + u_p2)
                        idx = -4
                    else
                        idx = -3
                    end
                else
                    if abs(-10*u_m1 + 5*u_m2 - u_m3 + 10*u_0 - 5*u_p1 + u_p2) < abs(5*u_m1 - u_m2 - 10*u_0 + 10*u_p1 - 5*u_p2 + u_p3)
                        idx = -3
                    else
                        idx = -2
                    end
                end
            else
                if abs(-4*u_m1 + u_m2 + 6*u_0 - 4*u_p1 + u_p2) < abs(u_m1 - 4*u_0 + 6*u_p1 - 4*u_p2 + u_p3)
                    if abs(-10*u_m1 + 5*u_m2 - u_m3 + 10*u_0 - 5*u_p1 + u_p2) < abs(5*u_m1 - u_m2 - 10*u_0 + 10*u_p1 - 5*u_p2 + u_p3)
                        idx = -3
                    else
                        idx = -2
                    end
                else
                    if abs(5*u_m1 - u_m2 - 10*u_0 + 10*u_p1 - 5*u_p2 + u_p3) < abs(-u_m1 + 5*u_0 - 10*u_p1 + 10*u_p2 - 5*u_p3 + u_p4)
                        idx = -2
                    else
                        idx = -1
                    end
                end
            end
        end
    else
        if abs(u_m1 - 2*u_0 + u_p1) < abs(u_0 - 2*u_p1 + u_p2)
            if abs(3*u_m1 - u_m2 - 3*u_0 + u_p1) < abs(-u_m1 + 3*u_0 - 3*u_p1 + u_p2)
                if abs(6*u_m1 - 4*u_m2 + u_m3 - 4*u_0 + u_p1) < abs(-4*u_m1 + u_m2 + 6*u_0 - 4*u_p1 + u_p2)
                    if abs(10*u_m1 - 10*u_m2 + 5*u_m3 - u_m4 - 5*u_0 + u_p1) < abs(-10*u_m1 + 5*u_m2 - u_m3 + 10*u_0 - 5*u_p1 + u_p2)
                        idx = -4
                    else
                        idx = -3
                    end
                else
                    if abs(-10*u_m1 + 5*u_m2 - u_m3 + 10*u_0 - 5*u_p1 + u_p2) < abs(5*u_m1 - u_m2 - 10*u_0 + 10*u_p1 - 5*u_p2 + u_p3)
                        idx = -3
                    else
                        idx = -2
                    end
                end
            else
                if abs(-4*u_m1 + u_m2 + 6*u_0 - 4*u_p1 + u_p2) < abs(u_m1 - 4*u_0 + 6*u_p1 - 4*u_p2 + u_p3)
                    if abs(-10*u_m1 + 5*u_m2 - u_m3 + 10*u_0 - 5*u_p1 + u_p2) < abs(5*u_m1 - u_m2 - 10*u_0 + 10*u_p1 - 5*u_p2 + u_p3)
                        idx = -3
                    else
                        idx = -2
                    end
                else
                    if abs(5*u_m1 - u_m2 - 10*u_0 + 10*u_p1 - 5*u_p2 + u_p3) < abs(-u_m1 + 5*u_0 - 10*u_p1 + 10*u_p2 - 5*u_p3 + u_p4)
                        idx = -2
                    else
                        idx = -1
                    end
                end
            end
        else
            if abs(-u_m1 + 3*u_0 - 3*u_p1 + u_p2) < abs(-u_0 + 3*u_p1 - 3*u_p2 + u_p3)
                if abs(-4*u_m1 + u_m2 + 6*u_0 - 4*u_p1 + u_p2) < abs(u_m1 - 4*u_0 + 6*u_p1 - 4*u_p2 + u_p3)
                    if abs(-10*u_m1 + 5*u_m2 - u_m3 + 10*u_0 - 5*u_p1 + u_p2) < abs(5*u_m1 - u_m2 - 10*u_0 + 10*u_p1 - 5*u_p2 + u_p3)
                        idx = -3
                    else
                        idx = -2
                    end
                else
                    if abs(5*u_m1 - u_m2 - 10*u_0 + 10*u_p1 - 5*u_p2 + u_p3) < abs(-u_m1 + 5*u_0 - 10*u_p1 + 10*u_p2 - 5*u_p3 + u_p4)
                        idx = -2
                    else
                        idx = -1
                    end
                end
            else
                if abs(u_m1 - 4*u_0 + 6*u_p1 - 4*u_p2 + u_p3) < abs(u_0 - 4*u_p1 + 6*u_p2 - 4*u_p3 + u_p4)
                    if abs(5*u_m1 - u_m2 - 10*u_0 + 10*u_p1 - 5*u_p2 + u_p3) < abs(-u_m1 + 5*u_0 - 10*u_p1 + 10*u_p2 - 5*u_p3 + u_p4)
                        idx = -2
                    else
                        idx = -1
                    end
                else
                    if abs(-u_m1 + 5*u_0 - 10*u_p1 + 10*u_p2 - 5*u_p3 + u_p4) < abs(-u_0 + 5*u_p1 - 10*u_p2 + 10*u_p3 - 5*u_p4 + u_p5)
                        idx = -1
                    else
                        idx = 0
                    end
                end
            end
        end
    end

    idx
end

@inline function (::ClassicalChoiceENO)(u_m6, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, u_p6)
    if abs(-u_m1 + u_0) < abs(-u_0 + u_p1)
        if abs(-2*u_m1 + u_m2 + u_0) < abs(u_m1 - 2*u_0 + u_p1)
            if abs(-3*u_m1 + 3*u_m2 - u_m3 + u_0) < abs(3*u_m1 - u_m2 - 3*u_0 + u_p1)
                if abs(-4*u_m1 + 6*u_m2 - 4*u_m3 + u_m4 + u_0) < abs(6*u_m1 - 4*u_m2 + u_m3 - 4*u_0 + u_p1)
                    if abs(-5*u_m1 + 10*u_m2 - 10*u_m3 + 5*u_m4 - u_m5 + u_0) < abs(10*u_m1 - 10*u_m2 + 5*u_m3 - u_m4 - 5*u_0 + u_p1)
                        if abs(-6*u_m1 + 15*u_m2 - 20*u_m3 + 15*u_m4 - 6*u_m5 + u_m6 + u_0) < abs(15*u_m1 - 20*u_m2 + 15*u_m3 - 6*u_m4 + u_m5 - 6*u_0 + u_p1)
                            idx = -6
                        else
                            idx = -5
                        end
                    else
                        if abs(15*u_m1 - 20*u_m2 + 15*u_m3 - 6*u_m4 + u_m5 - 6*u_0 + u_p1) < abs(-20*u_m1 + 15*u_m2 - 6*u_m3 + u_m4 + 15*u_0 - 6*u_p1 + u_p2)
                            idx = -5
                        else
                            idx = -4
                        end
                    end
                else
                    if abs(10*u_m1 - 10*u_m2 + 5*u_m3 - u_m4 - 5*u_0 + u_p1) < abs(-10*u_m1 + 5*u_m2 - u_m3 + 10*u_0 - 5*u_p1 + u_p2)
                        if abs(15*u_m1 - 20*u_m2 + 15*u_m3 - 6*u_m4 + u_m5 - 6*u_0 + u_p1) < abs(-20*u_m1 + 15*u_m2 - 6*u_m3 + u_m4 + 15*u_0 - 6*u_p1 + u_p2)
                            idx = -5
                        else
                            idx = -4
                        end
                    else
                        if abs(-20*u_m1 + 15*u_m2 - 6*u_m3 + u_m4 + 15*u_0 - 6*u_p1 + u_p2) < abs(15*u_m1 - 6*u_m2 + u_m3 - 20*u_0 + 15*u_p1 - 6*u_p2 + u_p3)
                            idx = -4
                        else
                            idx = -3
                        end
                    end
                end
            else
                if abs(6*u_m1 - 4*u_m2 + u_m3 - 4*u_0 + u_p1) < abs(-4*u_m1 + u_m2 + 6*u_0 - 4*u_p1 + u_p2)
                    if abs(10*u_m1 - 10*u_m2 + 5*u_m3 - u_m4 - 5*u_0 + u_p1) < abs(-10*u_m1 + 5*u_m2 - u_m3 + 10*u_0 - 5*u_p1 + u_p2)
                        if abs(15*u_m1 - 20*u_m2 + 15*u_m3 - 6*u_m4 + u_m5 - 6*u_0 + u_p1) < abs(-20*u_m1 + 15*u_m2 - 6*u_m3 + u_m4 + 15*u_0 - 6*u_p1 + u_p2)
                            idx = -5
                        else
                            idx = -4
                        end
                    else
                        if abs(-20*u_m1 + 15*u_m2 - 6*u_m3 + u_m4 + 15*u_0 - 6*u_p1 + u_p2) < abs(15*u_m1 - 6*u_m2 + u_m3 - 20*u_0 + 15*u_p1 - 6*u_p2 + u_p3)
                            idx = -4
                        else
                            idx = -3
                        end
                    end
                else
                    if abs(-10*u_m1 + 5*u_m2 - u_m3 + 10*u_0 - 5*u_p1 + u_p2) < abs(5*u_m1 - u_m2 - 10*u_0 + 10*u_p1 - 5*u_p2 + u_p3)
                        if abs(-20*u_m1 + 15*u_m2 - 6*u_m3 + u_m4 + 15*u_0 - 6*u_p1 + u_p2) < abs(15*u_m1 - 6*u_m2 + u_m3 - 20*u_0 + 15*u_p1 - 6*u_p2 + u_p3)
                            idx = -4
                        else
                            idx = -3
                        end
                    else
                        if abs(15*u_m1 - 6*u_m2 + u_m3 - 20*u_0 + 15*u_p1 - 6*u_p2 + u_p3) < abs(-6*u_m1 + u_m2 + 15*u_0 - 20*u_p1 + 15*u_p2 - 6*u_p3 + u_p4)
                            idx = -3
                        else
                            idx = -2
                        end
                    end
                end
            end
        else
            if abs(3*u_m1 - u_m2 - 3*u_0 + u_p1) < abs(-u_m1 + 3*u_0 - 3*u_p1 + u_p2)
                if abs(6*u_m1 - 4*u_m2 + u_m3 - 4*u_0 + u_p1) < abs(-4*u_m1 + u_m2 + 6*u_0 - 4*u_p1 + u_p2)
                    if abs(10*u_m1 - 10*u_m2 + 5*u_m3 - u_m4 - 5*u_0 + u_p1) < abs(-10*u_m1 + 5*u_m2 - u_m3 + 10*u_0 - 5*u_p1 + u_p2)
                        if abs(15*u_m1 - 20*u_m2 + 15*u_m3 - 6*u_m4 + u_m5 - 6*u_0 + u_p1) < abs(-20*u_m1 + 15*u_m2 - 6*u_m3 + u_m4 + 15*u_0 - 6*u_p1 + u_p2)
                            idx = -5
                        else
                            idx = -4
                        end
                    else
                        if abs(-20*u_m1 + 15*u_m2 - 6*u_m3 + u_m4 + 15*u_0 - 6*u_p1 + u_p2) < abs(15*u_m1 - 6*u_m2 + u_m3 - 20*u_0 + 15*u_p1 - 6*u_p2 + u_p3)
                            idx = -4
                        else
                            idx = -3
                        end
                    end
                else
                    if abs(-10*u_m1 + 5*u_m2 - u_m3 + 10*u_0 - 5*u_p1 + u_p2) < abs(5*u_m1 - u_m2 - 10*u_0 + 10*u_p1 - 5*u_p2 + u_p3)
                        if abs(-20*u_m1 + 15*u_m2 - 6*u_m3 + u_m4 + 15*u_0 - 6*u_p1 + u_p2) < abs(15*u_m1 - 6*u_m2 + u_m3 - 20*u_0 + 15*u_p1 - 6*u_p2 + u_p3)
                            idx = -4
                        else
                            idx = -3
                        end
                    else
                        if abs(15*u_m1 - 6*u_m2 + u_m3 - 20*u_0 + 15*u_p1 - 6*u_p2 + u_p3) < abs(-6*u_m1 + u_m2 + 15*u_0 - 20*u_p1 + 15*u_p2 - 6*u_p3 + u_p4)
                            idx = -3
                        else
                            idx = -2
                        end
                    end
                end
            else
                if abs(-4*u_m1 + u_m2 + 6*u_0 - 4*u_p1 + u_p2) < abs(u_m1 - 4*u_0 + 6*u_p1 - 4*u_p2 + u_p3)
                    if abs(-10*u_m1 + 5*u_m2 - u_m3 + 10*u_0 - 5*u_p1 + u_p2) < abs(5*u_m1 - u_m2 - 10*u_0 + 10*u_p1 - 5*u_p2 + u_p3)
                        if abs(-20*u_m1 + 15*u_m2 - 6*u_m3 + u_m4 + 15*u_0 - 6*u_p1 + u_p2) < abs(15*u_m1 - 6*u_m2 + u_m3 - 20*u_0 + 15*u_p1 - 6*u_p2 + u_p3)
                            idx = -4
                        else
                            idx = -3
                        end
                    else
                        if abs(15*u_m1 - 6*u_m2 + u_m3 - 20*u_0 + 15*u_p1 - 6*u_p2 + u_p3) < abs(-6*u_m1 + u_m2 + 15*u_0 - 20*u_p1 + 15*u_p2 - 6*u_p3 + u_p4)
                            idx = -3
                        else
                            idx = -2
                        end
                    end
                else
                    if abs(5*u_m1 - u_m2 - 10*u_0 + 10*u_p1 - 5*u_p2 + u_p3) < abs(-u_m1 + 5*u_0 - 10*u_p1 + 10*u_p2 - 5*u_p3 + u_p4)
                        if abs(15*u_m1 - 6*u_m2 + u_m3 - 20*u_0 + 15*u_p1 - 6*u_p2 + u_p3) < abs(-6*u_m1 + u_m2 + 15*u_0 - 20*u_p1 + 15*u_p2 - 6*u_p3 + u_p4)
                            idx = -3
                        else
                            idx = -2
                        end
                    else
                        if abs(-6*u_m1 + u_m2 + 15*u_0 - 20*u_p1 + 15*u_p2 - 6*u_p3 + u_p4) < abs(u_m1 - 6*u_0 + 15*u_p1 - 20*u_p2 + 15*u_p3 - 6*u_p4 + u_p5)
                            idx = -2
                        else
                            idx = -1
                        end
                    end
                end
            end
        end
    else
        if abs(u_m1 - 2*u_0 + u_p1) < abs(u_0 - 2*u_p1 + u_p2)
            if abs(3*u_m1 - u_m2 - 3*u_0 + u_p1) < abs(-u_m1 + 3*u_0 - 3*u_p1 + u_p2)
                if abs(6*u_m1 - 4*u_m2 + u_m3 - 4*u_0 + u_p1) < abs(-4*u_m1 + u_m2 + 6*u_0 - 4*u_p1 + u_p2)
                    if abs(10*u_m1 - 10*u_m2 + 5*u_m3 - u_m4 - 5*u_0 + u_p1) < abs(-10*u_m1 + 5*u_m2 - u_m3 + 10*u_0 - 5*u_p1 + u_p2)
                        if abs(15*u_m1 - 20*u_m2 + 15*u_m3 - 6*u_m4 + u_m5 - 6*u_0 + u_p1) < abs(-20*u_m1 + 15*u_m2 - 6*u_m3 + u_m4 + 15*u_0 - 6*u_p1 + u_p2)
                            idx = -5
                        else
                            idx = -4
                        end
                    else
                        if abs(-20*u_m1 + 15*u_m2 - 6*u_m3 + u_m4 + 15*u_0 - 6*u_p1 + u_p2) < abs(15*u_m1 - 6*u_m2 + u_m3 - 20*u_0 + 15*u_p1 - 6*u_p2 + u_p3)
                            idx = -4
                        else
                            idx = -3
                        end
                    end
                else
                    if abs(-10*u_m1 + 5*u_m2 - u_m3 + 10*u_0 - 5*u_p1 + u_p2) < abs(5*u_m1 - u_m2 - 10*u_0 + 10*u_p1 - 5*u_p2 + u_p3)
                        if abs(-20*u_m1 + 15*u_m2 - 6*u_m3 + u_m4 + 15*u_0 - 6*u_p1 + u_p2) < abs(15*u_m1 - 6*u_m2 + u_m3 - 20*u_0 + 15*u_p1 - 6*u_p2 + u_p3)
                            idx = -4
                        else
                            idx = -3
                        end
                    else
                        if abs(15*u_m1 - 6*u_m2 + u_m3 - 20*u_0 + 15*u_p1 - 6*u_p2 + u_p3) < abs(-6*u_m1 + u_m2 + 15*u_0 - 20*u_p1 + 15*u_p2 - 6*u_p3 + u_p4)
                            idx = -3
                        else
                            idx = -2
                        end
                    end
                end
            else
                if abs(-4*u_m1 + u_m2 + 6*u_0 - 4*u_p1 + u_p2) < abs(u_m1 - 4*u_0 + 6*u_p1 - 4*u_p2 + u_p3)
                    if abs(-10*u_m1 + 5*u_m2 - u_m3 + 10*u_0 - 5*u_p1 + u_p2) < abs(5*u_m1 - u_m2 - 10*u_0 + 10*u_p1 - 5*u_p2 + u_p3)
                        if abs(-20*u_m1 + 15*u_m2 - 6*u_m3 + u_m4 + 15*u_0 - 6*u_p1 + u_p2) < abs(15*u_m1 - 6*u_m2 + u_m3 - 20*u_0 + 15*u_p1 - 6*u_p2 + u_p3)
                            idx = -4
                        else
                            idx = -3
                        end
                    else
                        if abs(15*u_m1 - 6*u_m2 + u_m3 - 20*u_0 + 15*u_p1 - 6*u_p2 + u_p3) < abs(-6*u_m1 + u_m2 + 15*u_0 - 20*u_p1 + 15*u_p2 - 6*u_p3 + u_p4)
                            idx = -3
                        else
                            idx = -2
                        end
                    end
                else
                    if abs(5*u_m1 - u_m2 - 10*u_0 + 10*u_p1 - 5*u_p2 + u_p3) < abs(-u_m1 + 5*u_0 - 10*u_p1 + 10*u_p2 - 5*u_p3 + u_p4)
                        if abs(15*u_m1 - 6*u_m2 + u_m3 - 20*u_0 + 15*u_p1 - 6*u_p2 + u_p3) < abs(-6*u_m1 + u_m2 + 15*u_0 - 20*u_p1 + 15*u_p2 - 6*u_p3 + u_p4)
                            idx = -3
                        else
                            idx = -2
                        end
                    else
                        if abs(-6*u_m1 + u_m2 + 15*u_0 - 20*u_p1 + 15*u_p2 - 6*u_p3 + u_p4) < abs(u_m1 - 6*u_0 + 15*u_p1 - 20*u_p2 + 15*u_p3 - 6*u_p4 + u_p5)
                            idx = -2
                        else
                            idx = -1
                        end
                    end
                end
            end
        else
            if abs(-u_m1 + 3*u_0 - 3*u_p1 + u_p2) < abs(-u_0 + 3*u_p1 - 3*u_p2 + u_p3)
                if abs(-4*u_m1 + u_m2 + 6*u_0 - 4*u_p1 + u_p2) < abs(u_m1 - 4*u_0 + 6*u_p1 - 4*u_p2 + u_p3)
                    if abs(-10*u_m1 + 5*u_m2 - u_m3 + 10*u_0 - 5*u_p1 + u_p2) < abs(5*u_m1 - u_m2 - 10*u_0 + 10*u_p1 - 5*u_p2 + u_p3)
                        if abs(-20*u_m1 + 15*u_m2 - 6*u_m3 + u_m4 + 15*u_0 - 6*u_p1 + u_p2) < abs(15*u_m1 - 6*u_m2 + u_m3 - 20*u_0 + 15*u_p1 - 6*u_p2 + u_p3)
                            idx = -4
                        else
                            idx = -3
                        end
                    else
                        if abs(15*u_m1 - 6*u_m2 + u_m3 - 20*u_0 + 15*u_p1 - 6*u_p2 + u_p3) < abs(-6*u_m1 + u_m2 + 15*u_0 - 20*u_p1 + 15*u_p2 - 6*u_p3 + u_p4)
                            idx = -3
                        else
                            idx = -2
                        end
                    end
                else
                    if abs(5*u_m1 - u_m2 - 10*u_0 + 10*u_p1 - 5*u_p2 + u_p3) < abs(-u_m1 + 5*u_0 - 10*u_p1 + 10*u_p2 - 5*u_p3 + u_p4)
                        if abs(15*u_m1 - 6*u_m2 + u_m3 - 20*u_0 + 15*u_p1 - 6*u_p2 + u_p3) < abs(-6*u_m1 + u_m2 + 15*u_0 - 20*u_p1 + 15*u_p2 - 6*u_p3 + u_p4)
                            idx = -3
                        else
                            idx = -2
                        end
                    else
                        if abs(-6*u_m1 + u_m2 + 15*u_0 - 20*u_p1 + 15*u_p2 - 6*u_p3 + u_p4) < abs(u_m1 - 6*u_0 + 15*u_p1 - 20*u_p2 + 15*u_p3 - 6*u_p4 + u_p5)
                            idx = -2
                        else
                            idx = -1
                        end
                    end
                end
            else
                if abs(u_m1 - 4*u_0 + 6*u_p1 - 4*u_p2 + u_p3) < abs(u_0 - 4*u_p1 + 6*u_p2 - 4*u_p3 + u_p4)
                    if abs(5*u_m1 - u_m2 - 10*u_0 + 10*u_p1 - 5*u_p2 + u_p3) < abs(-u_m1 + 5*u_0 - 10*u_p1 + 10*u_p2 - 5*u_p3 + u_p4)
                        if abs(15*u_m1 - 6*u_m2 + u_m3 - 20*u_0 + 15*u_p1 - 6*u_p2 + u_p3) < abs(-6*u_m1 + u_m2 + 15*u_0 - 20*u_p1 + 15*u_p2 - 6*u_p3 + u_p4)
                            idx = -3
                        else
                            idx = -2
                        end
                    else
                        if abs(-6*u_m1 + u_m2 + 15*u_0 - 20*u_p1 + 15*u_p2 - 6*u_p3 + u_p4) < abs(u_m1 - 6*u_0 + 15*u_p1 - 20*u_p2 + 15*u_p3 - 6*u_p4 + u_p5)
                            idx = -2
                        else
                            idx = -1
                        end
                    end
                else
                    if abs(-u_m1 + 5*u_0 - 10*u_p1 + 10*u_p2 - 5*u_p3 + u_p4) < abs(-u_0 + 5*u_p1 - 10*u_p2 + 10*u_p3 - 5*u_p4 + u_p5)
                        if abs(-6*u_m1 + u_m2 + 15*u_0 - 20*u_p1 + 15*u_p2 - 6*u_p3 + u_p4) < abs(u_m1 - 6*u_0 + 15*u_p1 - 20*u_p2 + 15*u_p3 - 6*u_p4 + u_p5)
                            idx = -2
                        else
                            idx = -1
                        end
                    else
                        if abs(u_m1 - 6*u_0 + 15*u_p1 - 20*u_p2 + 15*u_p3 - 6*u_p4 + u_p5) < abs(u_0 - 6*u_p1 + 15*u_p2 - 20*u_p3 + 15*u_p4 - 6*u_p5 + u_p6)
                            idx = -1
                        else
                            idx = 0
                        end
                    end
                end
            end
        end
    end

    idx
end

@inline function (::ClassicalChoiceENO)(u_m7, u_m6, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, u_p6, u_p7)
    if abs(-u_m1 + u_0) < abs(-u_0 + u_p1)
        if abs(-2*u_m1 + u_m2 + u_0) < abs(u_m1 - 2*u_0 + u_p1)
            if abs(-3*u_m1 + 3*u_m2 - u_m3 + u_0) < abs(3*u_m1 - u_m2 - 3*u_0 + u_p1)
                if abs(-4*u_m1 + 6*u_m2 - 4*u_m3 + u_m4 + u_0) < abs(6*u_m1 - 4*u_m2 + u_m3 - 4*u_0 + u_p1)
                    if abs(-5*u_m1 + 10*u_m2 - 10*u_m3 + 5*u_m4 - u_m5 + u_0) < abs(10*u_m1 - 10*u_m2 + 5*u_m3 - u_m4 - 5*u_0 + u_p1)
                        if abs(-6*u_m1 + 15*u_m2 - 20*u_m3 + 15*u_m4 - 6*u_m5 + u_m6 + u_0) < abs(15*u_m1 - 20*u_m2 + 15*u_m3 - 6*u_m4 + u_m5 - 6*u_0 + u_p1)
                            if abs(-7*u_m1 + 21*u_m2 - 35*u_m3 + 35*u_m4 - 21*u_m5 + 7*u_m6 - u_m7 + u_0) < abs(21*u_m1 - 35*u_m2 + 35*u_m3 - 21*u_m4 + 7*u_m5 - u_m6 - 7*u_0 + u_p1)
                                idx = -7
                            else
                                idx = -6
                            end
                        else
                            if abs(21*u_m1 - 35*u_m2 + 35*u_m3 - 21*u_m4 + 7*u_m5 - u_m6 - 7*u_0 + u_p1) < abs(-35*u_m1 + 35*u_m2 - 21*u_m3 + 7*u_m4 - u_m5 + 21*u_0 - 7*u_p1 + u_p2)
                                idx = -6
                            else
                                idx = -5
                            end
                        end
                    else
                        if abs(15*u_m1 - 20*u_m2 + 15*u_m3 - 6*u_m4 + u_m5 - 6*u_0 + u_p1) < abs(-20*u_m1 + 15*u_m2 - 6*u_m3 + u_m4 + 15*u_0 - 6*u_p1 + u_p2)
                            if abs(21*u_m1 - 35*u_m2 + 35*u_m3 - 21*u_m4 + 7*u_m5 - u_m6 - 7*u_0 + u_p1) < abs(-35*u_m1 + 35*u_m2 - 21*u_m3 + 7*u_m4 - u_m5 + 21*u_0 - 7*u_p1 + u_p2)
                                idx = -6
                            else
                                idx = -5
                            end
                        else
                            if abs(-35*u_m1 + 35*u_m2 - 21*u_m3 + 7*u_m4 - u_m5 + 21*u_0 - 7*u_p1 + u_p2) < abs(35*u_m1 - 21*u_m2 + 7*u_m3 - u_m4 - 35*u_0 + 21*u_p1 - 7*u_p2 + u_p3)
                                idx = -5
                            else
                                idx = -4
                            end
                        end
                    end
                else
                    if abs(10*u_m1 - 10*u_m2 + 5*u_m3 - u_m4 - 5*u_0 + u_p1) < abs(-10*u_m1 + 5*u_m2 - u_m3 + 10*u_0 - 5*u_p1 + u_p2)
                        if abs(15*u_m1 - 20*u_m2 + 15*u_m3 - 6*u_m4 + u_m5 - 6*u_0 + u_p1) < abs(-20*u_m1 + 15*u_m2 - 6*u_m3 + u_m4 + 15*u_0 - 6*u_p1 + u_p2)
                            if abs(21*u_m1 - 35*u_m2 + 35*u_m3 - 21*u_m4 + 7*u_m5 - u_m6 - 7*u_0 + u_p1) < abs(-35*u_m1 + 35*u_m2 - 21*u_m3 + 7*u_m4 - u_m5 + 21*u_0 - 7*u_p1 + u_p2)
                                idx = -6
                            else
                                idx = -5
                            end
                        else
                            if abs(-35*u_m1 + 35*u_m2 - 21*u_m3 + 7*u_m4 - u_m5 + 21*u_0 - 7*u_p1 + u_p2) < abs(35*u_m1 - 21*u_m2 + 7*u_m3 - u_m4 - 35*u_0 + 21*u_p1 - 7*u_p2 + u_p3)
                                idx = -5
                            else
                                idx = -4
                            end
                        end
                    else
                        if abs(-20*u_m1 + 15*u_m2 - 6*u_m3 + u_m4 + 15*u_0 - 6*u_p1 + u_p2) < abs(15*u_m1 - 6*u_m2 + u_m3 - 20*u_0 + 15*u_p1 - 6*u_p2 + u_p3)
                            if abs(-35*u_m1 + 35*u_m2 - 21*u_m3 + 7*u_m4 - u_m5 + 21*u_0 - 7*u_p1 + u_p2) < abs(35*u_m1 - 21*u_m2 + 7*u_m3 - u_m4 - 35*u_0 + 21*u_p1 - 7*u_p2 + u_p3)
                                idx = -5
                            else
                                idx = -4
                            end
                        else
                            if abs(35*u_m1 - 21*u_m2 + 7*u_m3 - u_m4 - 35*u_0 + 21*u_p1 - 7*u_p2 + u_p3) < abs(-21*u_m1 + 7*u_m2 - u_m3 + 35*u_0 - 35*u_p1 + 21*u_p2 - 7*u_p3 + u_p4)
                                idx = -4
                            else
                                idx = -3
                            end
                        end
                    end
                end
            else
                if abs(6*u_m1 - 4*u_m2 + u_m3 - 4*u_0 + u_p1) < abs(-4*u_m1 + u_m2 + 6*u_0 - 4*u_p1 + u_p2)
                    if abs(10*u_m1 - 10*u_m2 + 5*u_m3 - u_m4 - 5*u_0 + u_p1) < abs(-10*u_m1 + 5*u_m2 - u_m3 + 10*u_0 - 5*u_p1 + u_p2)
                        if abs(15*u_m1 - 20*u_m2 + 15*u_m3 - 6*u_m4 + u_m5 - 6*u_0 + u_p1) < abs(-20*u_m1 + 15*u_m2 - 6*u_m3 + u_m4 + 15*u_0 - 6*u_p1 + u_p2)
                            if abs(21*u_m1 - 35*u_m2 + 35*u_m3 - 21*u_m4 + 7*u_m5 - u_m6 - 7*u_0 + u_p1) < abs(-35*u_m1 + 35*u_m2 - 21*u_m3 + 7*u_m4 - u_m5 + 21*u_0 - 7*u_p1 + u_p2)
                                idx = -6
                            else
                                idx = -5
                            end
                        else
                            if abs(-35*u_m1 + 35*u_m2 - 21*u_m3 + 7*u_m4 - u_m5 + 21*u_0 - 7*u_p1 + u_p2) < abs(35*u_m1 - 21*u_m2 + 7*u_m3 - u_m4 - 35*u_0 + 21*u_p1 - 7*u_p2 + u_p3)
                                idx = -5
                            else
                                idx = -4
                            end
                        end
                    else
                        if abs(-20*u_m1 + 15*u_m2 - 6*u_m3 + u_m4 + 15*u_0 - 6*u_p1 + u_p2) < abs(15*u_m1 - 6*u_m2 + u_m3 - 20*u_0 + 15*u_p1 - 6*u_p2 + u_p3)
                            if abs(-35*u_m1 + 35*u_m2 - 21*u_m3 + 7*u_m4 - u_m5 + 21*u_0 - 7*u_p1 + u_p2) < abs(35*u_m1 - 21*u_m2 + 7*u_m3 - u_m4 - 35*u_0 + 21*u_p1 - 7*u_p2 + u_p3)
                                idx = -5
                            else
                                idx = -4
                            end
                        else
                            if abs(35*u_m1 - 21*u_m2 + 7*u_m3 - u_m4 - 35*u_0 + 21*u_p1 - 7*u_p2 + u_p3) < abs(-21*u_m1 + 7*u_m2 - u_m3 + 35*u_0 - 35*u_p1 + 21*u_p2 - 7*u_p3 + u_p4)
                                idx = -4
                            else
                                idx = -3
                            end
                        end
                    end
                else
                    if abs(-10*u_m1 + 5*u_m2 - u_m3 + 10*u_0 - 5*u_p1 + u_p2) < abs(5*u_m1 - u_m2 - 10*u_0 + 10*u_p1 - 5*u_p2 + u_p3)
                        if abs(-20*u_m1 + 15*u_m2 - 6*u_m3 + u_m4 + 15*u_0 - 6*u_p1 + u_p2) < abs(15*u_m1 - 6*u_m2 + u_m3 - 20*u_0 + 15*u_p1 - 6*u_p2 + u_p3)
                            if abs(-35*u_m1 + 35*u_m2 - 21*u_m3 + 7*u_m4 - u_m5 + 21*u_0 - 7*u_p1 + u_p2) < abs(35*u_m1 - 21*u_m2 + 7*u_m3 - u_m4 - 35*u_0 + 21*u_p1 - 7*u_p2 + u_p3)
                                idx = -5
                            else
                                idx = -4
                            end
                        else
                            if abs(35*u_m1 - 21*u_m2 + 7*u_m3 - u_m4 - 35*u_0 + 21*u_p1 - 7*u_p2 + u_p3) < abs(-21*u_m1 + 7*u_m2 - u_m3 + 35*u_0 - 35*u_p1 + 21*u_p2 - 7*u_p3 + u_p4)
                                idx = -4
                            else
                                idx = -3
                            end
                        end
                    else
                        if abs(15*u_m1 - 6*u_m2 + u_m3 - 20*u_0 + 15*u_p1 - 6*u_p2 + u_p3) < abs(-6*u_m1 + u_m2 + 15*u_0 - 20*u_p1 + 15*u_p2 - 6*u_p3 + u_p4)
                            if abs(35*u_m1 - 21*u_m2 + 7*u_m3 - u_m4 - 35*u_0 + 21*u_p1 - 7*u_p2 + u_p3) < abs(-21*u_m1 + 7*u_m2 - u_m3 + 35*u_0 - 35*u_p1 + 21*u_p2 - 7*u_p3 + u_p4)
                                idx = -4
                            else
                                idx = -3
                            end
                        else
                            if abs(-21*u_m1 + 7*u_m2 - u_m3 + 35*u_0 - 35*u_p1 + 21*u_p2 - 7*u_p3 + u_p4) < abs(7*u_m1 - u_m2 - 21*u_0 + 35*u_p1 - 35*u_p2 + 21*u_p3 - 7*u_p4 + u_p5)
                                idx = -3
                            else
                                idx = -2
                            end
                        end
                    end
                end
            end
        else
            if abs(3*u_m1 - u_m2 - 3*u_0 + u_p1) < abs(-u_m1 + 3*u_0 - 3*u_p1 + u_p2)
                if abs(6*u_m1 - 4*u_m2 + u_m3 - 4*u_0 + u_p1) < abs(-4*u_m1 + u_m2 + 6*u_0 - 4*u_p1 + u_p2)
                    if abs(10*u_m1 - 10*u_m2 + 5*u_m3 - u_m4 - 5*u_0 + u_p1) < abs(-10*u_m1 + 5*u_m2 - u_m3 + 10*u_0 - 5*u_p1 + u_p2)
                        if abs(15*u_m1 - 20*u_m2 + 15*u_m3 - 6*u_m4 + u_m5 - 6*u_0 + u_p1) < abs(-20*u_m1 + 15*u_m2 - 6*u_m3 + u_m4 + 15*u_0 - 6*u_p1 + u_p2)
                            if abs(21*u_m1 - 35*u_m2 + 35*u_m3 - 21*u_m4 + 7*u_m5 - u_m6 - 7*u_0 + u_p1) < abs(-35*u_m1 + 35*u_m2 - 21*u_m3 + 7*u_m4 - u_m5 + 21*u_0 - 7*u_p1 + u_p2)
                                idx = -6
                            else
                                idx = -5
                            end
                        else
                            if abs(-35*u_m1 + 35*u_m2 - 21*u_m3 + 7*u_m4 - u_m5 + 21*u_0 - 7*u_p1 + u_p2) < abs(35*u_m1 - 21*u_m2 + 7*u_m3 - u_m4 - 35*u_0 + 21*u_p1 - 7*u_p2 + u_p3)
                                idx = -5
                            else
                                idx = -4
                            end
                        end
                    else
                        if abs(-20*u_m1 + 15*u_m2 - 6*u_m3 + u_m4 + 15*u_0 - 6*u_p1 + u_p2) < abs(15*u_m1 - 6*u_m2 + u_m3 - 20*u_0 + 15*u_p1 - 6*u_p2 + u_p3)
                            if abs(-35*u_m1 + 35*u_m2 - 21*u_m3 + 7*u_m4 - u_m5 + 21*u_0 - 7*u_p1 + u_p2) < abs(35*u_m1 - 21*u_m2 + 7*u_m3 - u_m4 - 35*u_0 + 21*u_p1 - 7*u_p2 + u_p3)
                                idx = -5
                            else
                                idx = -4
                            end
                        else
                            if abs(35*u_m1 - 21*u_m2 + 7*u_m3 - u_m4 - 35*u_0 + 21*u_p1 - 7*u_p2 + u_p3) < abs(-21*u_m1 + 7*u_m2 - u_m3 + 35*u_0 - 35*u_p1 + 21*u_p2 - 7*u_p3 + u_p4)
                                idx = -4
                            else
                                idx = -3
                            end
                        end
                    end
                else
                    if abs(-10*u_m1 + 5*u_m2 - u_m3 + 10*u_0 - 5*u_p1 + u_p2) < abs(5*u_m1 - u_m2 - 10*u_0 + 10*u_p1 - 5*u_p2 + u_p3)
                        if abs(-20*u_m1 + 15*u_m2 - 6*u_m3 + u_m4 + 15*u_0 - 6*u_p1 + u_p2) < abs(15*u_m1 - 6*u_m2 + u_m3 - 20*u_0 + 15*u_p1 - 6*u_p2 + u_p3)
                            if abs(-35*u_m1 + 35*u_m2 - 21*u_m3 + 7*u_m4 - u_m5 + 21*u_0 - 7*u_p1 + u_p2) < abs(35*u_m1 - 21*u_m2 + 7*u_m3 - u_m4 - 35*u_0 + 21*u_p1 - 7*u_p2 + u_p3)
                                idx = -5
                            else
                                idx = -4
                            end
                        else
                            if abs(35*u_m1 - 21*u_m2 + 7*u_m3 - u_m4 - 35*u_0 + 21*u_p1 - 7*u_p2 + u_p3) < abs(-21*u_m1 + 7*u_m2 - u_m3 + 35*u_0 - 35*u_p1 + 21*u_p2 - 7*u_p3 + u_p4)
                                idx = -4
                            else
                                idx = -3
                            end
                        end
                    else
                        if abs(15*u_m1 - 6*u_m2 + u_m3 - 20*u_0 + 15*u_p1 - 6*u_p2 + u_p3) < abs(-6*u_m1 + u_m2 + 15*u_0 - 20*u_p1 + 15*u_p2 - 6*u_p3 + u_p4)
                            if abs(35*u_m1 - 21*u_m2 + 7*u_m3 - u_m4 - 35*u_0 + 21*u_p1 - 7*u_p2 + u_p3) < abs(-21*u_m1 + 7*u_m2 - u_m3 + 35*u_0 - 35*u_p1 + 21*u_p2 - 7*u_p3 + u_p4)
                                idx = -4
                            else
                                idx = -3
                            end
                        else
                            if abs(-21*u_m1 + 7*u_m2 - u_m3 + 35*u_0 - 35*u_p1 + 21*u_p2 - 7*u_p3 + u_p4) < abs(7*u_m1 - u_m2 - 21*u_0 + 35*u_p1 - 35*u_p2 + 21*u_p3 - 7*u_p4 + u_p5)
                                idx = -3
                            else
                                idx = -2
                            end
                        end
                    end
                end
            else
                if abs(-4*u_m1 + u_m2 + 6*u_0 - 4*u_p1 + u_p2) < abs(u_m1 - 4*u_0 + 6*u_p1 - 4*u_p2 + u_p3)
                    if abs(-10*u_m1 + 5*u_m2 - u_m3 + 10*u_0 - 5*u_p1 + u_p2) < abs(5*u_m1 - u_m2 - 10*u_0 + 10*u_p1 - 5*u_p2 + u_p3)
                        if abs(-20*u_m1 + 15*u_m2 - 6*u_m3 + u_m4 + 15*u_0 - 6*u_p1 + u_p2) < abs(15*u_m1 - 6*u_m2 + u_m3 - 20*u_0 + 15*u_p1 - 6*u_p2 + u_p3)
                            if abs(-35*u_m1 + 35*u_m2 - 21*u_m3 + 7*u_m4 - u_m5 + 21*u_0 - 7*u_p1 + u_p2) < abs(35*u_m1 - 21*u_m2 + 7*u_m3 - u_m4 - 35*u_0 + 21*u_p1 - 7*u_p2 + u_p3)
                                idx = -5
                            else
                                idx = -4
                            end
                        else
                            if abs(35*u_m1 - 21*u_m2 + 7*u_m3 - u_m4 - 35*u_0 + 21*u_p1 - 7*u_p2 + u_p3) < abs(-21*u_m1 + 7*u_m2 - u_m3 + 35*u_0 - 35*u_p1 + 21*u_p2 - 7*u_p3 + u_p4)
                                idx = -4
                            else
                                idx = -3
                            end
                        end
                    else
                        if abs(15*u_m1 - 6*u_m2 + u_m3 - 20*u_0 + 15*u_p1 - 6*u_p2 + u_p3) < abs(-6*u_m1 + u_m2 + 15*u_0 - 20*u_p1 + 15*u_p2 - 6*u_p3 + u_p4)
                            if abs(35*u_m1 - 21*u_m2 + 7*u_m3 - u_m4 - 35*u_0 + 21*u_p1 - 7*u_p2 + u_p3) < abs(-21*u_m1 + 7*u_m2 - u_m3 + 35*u_0 - 35*u_p1 + 21*u_p2 - 7*u_p3 + u_p4)
                                idx = -4
                            else
                                idx = -3
                            end
                        else
                            if abs(-21*u_m1 + 7*u_m2 - u_m3 + 35*u_0 - 35*u_p1 + 21*u_p2 - 7*u_p3 + u_p4) < abs(7*u_m1 - u_m2 - 21*u_0 + 35*u_p1 - 35*u_p2 + 21*u_p3 - 7*u_p4 + u_p5)
                                idx = -3
                            else
                                idx = -2
                            end
                        end
                    end
                else
                    if abs(5*u_m1 - u_m2 - 10*u_0 + 10*u_p1 - 5*u_p2 + u_p3) < abs(-u_m1 + 5*u_0 - 10*u_p1 + 10*u_p2 - 5*u_p3 + u_p4)
                        if abs(15*u_m1 - 6*u_m2 + u_m3 - 20*u_0 + 15*u_p1 - 6*u_p2 + u_p3) < abs(-6*u_m1 + u_m2 + 15*u_0 - 20*u_p1 + 15*u_p2 - 6*u_p3 + u_p4)
                            if abs(35*u_m1 - 21*u_m2 + 7*u_m3 - u_m4 - 35*u_0 + 21*u_p1 - 7*u_p2 + u_p3) < abs(-21*u_m1 + 7*u_m2 - u_m3 + 35*u_0 - 35*u_p1 + 21*u_p2 - 7*u_p3 + u_p4)
                                idx = -4
                            else
                                idx = -3
                            end
                        else
                            if abs(-21*u_m1 + 7*u_m2 - u_m3 + 35*u_0 - 35*u_p1 + 21*u_p2 - 7*u_p3 + u_p4) < abs(7*u_m1 - u_m2 - 21*u_0 + 35*u_p1 - 35*u_p2 + 21*u_p3 - 7*u_p4 + u_p5)
                                idx = -3
                            else
                                idx = -2
                            end
                        end
                    else
                        if abs(-6*u_m1 + u_m2 + 15*u_0 - 20*u_p1 + 15*u_p2 - 6*u_p3 + u_p4) < abs(u_m1 - 6*u_0 + 15*u_p1 - 20*u_p2 + 15*u_p3 - 6*u_p4 + u_p5)
                            if abs(-21*u_m1 + 7*u_m2 - u_m3 + 35*u_0 - 35*u_p1 + 21*u_p2 - 7*u_p3 + u_p4) < abs(7*u_m1 - u_m2 - 21*u_0 + 35*u_p1 - 35*u_p2 + 21*u_p3 - 7*u_p4 + u_p5)
                                idx = -3
                            else
                                idx = -2
                            end
                        else
                            if abs(7*u_m1 - u_m2 - 21*u_0 + 35*u_p1 - 35*u_p2 + 21*u_p3 - 7*u_p4 + u_p5) < abs(-u_m1 + 7*u_0 - 21*u_p1 + 35*u_p2 - 35*u_p3 + 21*u_p4 - 7*u_p5 + u_p6)
                                idx = -2
                            else
                                idx = -1
                            end
                        end
                    end
                end
            end
        end
    else
        if abs(u_m1 - 2*u_0 + u_p1) < abs(u_0 - 2*u_p1 + u_p2)
            if abs(3*u_m1 - u_m2 - 3*u_0 + u_p1) < abs(-u_m1 + 3*u_0 - 3*u_p1 + u_p2)
                if abs(6*u_m1 - 4*u_m2 + u_m3 - 4*u_0 + u_p1) < abs(-4*u_m1 + u_m2 + 6*u_0 - 4*u_p1 + u_p2)
                    if abs(10*u_m1 - 10*u_m2 + 5*u_m3 - u_m4 - 5*u_0 + u_p1) < abs(-10*u_m1 + 5*u_m2 - u_m3 + 10*u_0 - 5*u_p1 + u_p2)
                        if abs(15*u_m1 - 20*u_m2 + 15*u_m3 - 6*u_m4 + u_m5 - 6*u_0 + u_p1) < abs(-20*u_m1 + 15*u_m2 - 6*u_m3 + u_m4 + 15*u_0 - 6*u_p1 + u_p2)
                            if abs(21*u_m1 - 35*u_m2 + 35*u_m3 - 21*u_m4 + 7*u_m5 - u_m6 - 7*u_0 + u_p1) < abs(-35*u_m1 + 35*u_m2 - 21*u_m3 + 7*u_m4 - u_m5 + 21*u_0 - 7*u_p1 + u_p2)
                                idx = -6
                            else
                                idx = -5
                            end
                        else
                            if abs(-35*u_m1 + 35*u_m2 - 21*u_m3 + 7*u_m4 - u_m5 + 21*u_0 - 7*u_p1 + u_p2) < abs(35*u_m1 - 21*u_m2 + 7*u_m3 - u_m4 - 35*u_0 + 21*u_p1 - 7*u_p2 + u_p3)
                                idx = -5
                            else
                                idx = -4
                            end
                        end
                    else
                        if abs(-20*u_m1 + 15*u_m2 - 6*u_m3 + u_m4 + 15*u_0 - 6*u_p1 + u_p2) < abs(15*u_m1 - 6*u_m2 + u_m3 - 20*u_0 + 15*u_p1 - 6*u_p2 + u_p3)
                            if abs(-35*u_m1 + 35*u_m2 - 21*u_m3 + 7*u_m4 - u_m5 + 21*u_0 - 7*u_p1 + u_p2) < abs(35*u_m1 - 21*u_m2 + 7*u_m3 - u_m4 - 35*u_0 + 21*u_p1 - 7*u_p2 + u_p3)
                                idx = -5
                            else
                                idx = -4
                            end
                        else
                            if abs(35*u_m1 - 21*u_m2 + 7*u_m3 - u_m4 - 35*u_0 + 21*u_p1 - 7*u_p2 + u_p3) < abs(-21*u_m1 + 7*u_m2 - u_m3 + 35*u_0 - 35*u_p1 + 21*u_p2 - 7*u_p3 + u_p4)
                                idx = -4
                            else
                                idx = -3
                            end
                        end
                    end
                else
                    if abs(-10*u_m1 + 5*u_m2 - u_m3 + 10*u_0 - 5*u_p1 + u_p2) < abs(5*u_m1 - u_m2 - 10*u_0 + 10*u_p1 - 5*u_p2 + u_p3)
                        if abs(-20*u_m1 + 15*u_m2 - 6*u_m3 + u_m4 + 15*u_0 - 6*u_p1 + u_p2) < abs(15*u_m1 - 6*u_m2 + u_m3 - 20*u_0 + 15*u_p1 - 6*u_p2 + u_p3)
                            if abs(-35*u_m1 + 35*u_m2 - 21*u_m3 + 7*u_m4 - u_m5 + 21*u_0 - 7*u_p1 + u_p2) < abs(35*u_m1 - 21*u_m2 + 7*u_m3 - u_m4 - 35*u_0 + 21*u_p1 - 7*u_p2 + u_p3)
                                idx = -5
                            else
                                idx = -4
                            end
                        else
                            if abs(35*u_m1 - 21*u_m2 + 7*u_m3 - u_m4 - 35*u_0 + 21*u_p1 - 7*u_p2 + u_p3) < abs(-21*u_m1 + 7*u_m2 - u_m3 + 35*u_0 - 35*u_p1 + 21*u_p2 - 7*u_p3 + u_p4)
                                idx = -4
                            else
                                idx = -3
                            end
                        end
                    else
                        if abs(15*u_m1 - 6*u_m2 + u_m3 - 20*u_0 + 15*u_p1 - 6*u_p2 + u_p3) < abs(-6*u_m1 + u_m2 + 15*u_0 - 20*u_p1 + 15*u_p2 - 6*u_p3 + u_p4)
                            if abs(35*u_m1 - 21*u_m2 + 7*u_m3 - u_m4 - 35*u_0 + 21*u_p1 - 7*u_p2 + u_p3) < abs(-21*u_m1 + 7*u_m2 - u_m3 + 35*u_0 - 35*u_p1 + 21*u_p2 - 7*u_p3 + u_p4)
                                idx = -4
                            else
                                idx = -3
                            end
                        else
                            if abs(-21*u_m1 + 7*u_m2 - u_m3 + 35*u_0 - 35*u_p1 + 21*u_p2 - 7*u_p3 + u_p4) < abs(7*u_m1 - u_m2 - 21*u_0 + 35*u_p1 - 35*u_p2 + 21*u_p3 - 7*u_p4 + u_p5)
                                idx = -3
                            else
                                idx = -2
                            end
                        end
                    end
                end
            else
                if abs(-4*u_m1 + u_m2 + 6*u_0 - 4*u_p1 + u_p2) < abs(u_m1 - 4*u_0 + 6*u_p1 - 4*u_p2 + u_p3)
                    if abs(-10*u_m1 + 5*u_m2 - u_m3 + 10*u_0 - 5*u_p1 + u_p2) < abs(5*u_m1 - u_m2 - 10*u_0 + 10*u_p1 - 5*u_p2 + u_p3)
                        if abs(-20*u_m1 + 15*u_m2 - 6*u_m3 + u_m4 + 15*u_0 - 6*u_p1 + u_p2) < abs(15*u_m1 - 6*u_m2 + u_m3 - 20*u_0 + 15*u_p1 - 6*u_p2 + u_p3)
                            if abs(-35*u_m1 + 35*u_m2 - 21*u_m3 + 7*u_m4 - u_m5 + 21*u_0 - 7*u_p1 + u_p2) < abs(35*u_m1 - 21*u_m2 + 7*u_m3 - u_m4 - 35*u_0 + 21*u_p1 - 7*u_p2 + u_p3)
                                idx = -5
                            else
                                idx = -4
                            end
                        else
                            if abs(35*u_m1 - 21*u_m2 + 7*u_m3 - u_m4 - 35*u_0 + 21*u_p1 - 7*u_p2 + u_p3) < abs(-21*u_m1 + 7*u_m2 - u_m3 + 35*u_0 - 35*u_p1 + 21*u_p2 - 7*u_p3 + u_p4)
                                idx = -4
                            else
                                idx = -3
                            end
                        end
                    else
                        if abs(15*u_m1 - 6*u_m2 + u_m3 - 20*u_0 + 15*u_p1 - 6*u_p2 + u_p3) < abs(-6*u_m1 + u_m2 + 15*u_0 - 20*u_p1 + 15*u_p2 - 6*u_p3 + u_p4)
                            if abs(35*u_m1 - 21*u_m2 + 7*u_m3 - u_m4 - 35*u_0 + 21*u_p1 - 7*u_p2 + u_p3) < abs(-21*u_m1 + 7*u_m2 - u_m3 + 35*u_0 - 35*u_p1 + 21*u_p2 - 7*u_p3 + u_p4)
                                idx = -4
                            else
                                idx = -3
                            end
                        else
                            if abs(-21*u_m1 + 7*u_m2 - u_m3 + 35*u_0 - 35*u_p1 + 21*u_p2 - 7*u_p3 + u_p4) < abs(7*u_m1 - u_m2 - 21*u_0 + 35*u_p1 - 35*u_p2 + 21*u_p3 - 7*u_p4 + u_p5)
                                idx = -3
                            else
                                idx = -2
                            end
                        end
                    end
                else
                    if abs(5*u_m1 - u_m2 - 10*u_0 + 10*u_p1 - 5*u_p2 + u_p3) < abs(-u_m1 + 5*u_0 - 10*u_p1 + 10*u_p2 - 5*u_p3 + u_p4)
                        if abs(15*u_m1 - 6*u_m2 + u_m3 - 20*u_0 + 15*u_p1 - 6*u_p2 + u_p3) < abs(-6*u_m1 + u_m2 + 15*u_0 - 20*u_p1 + 15*u_p2 - 6*u_p3 + u_p4)
                            if abs(35*u_m1 - 21*u_m2 + 7*u_m3 - u_m4 - 35*u_0 + 21*u_p1 - 7*u_p2 + u_p3) < abs(-21*u_m1 + 7*u_m2 - u_m3 + 35*u_0 - 35*u_p1 + 21*u_p2 - 7*u_p3 + u_p4)
                                idx = -4
                            else
                                idx = -3
                            end
                        else
                            if abs(-21*u_m1 + 7*u_m2 - u_m3 + 35*u_0 - 35*u_p1 + 21*u_p2 - 7*u_p3 + u_p4) < abs(7*u_m1 - u_m2 - 21*u_0 + 35*u_p1 - 35*u_p2 + 21*u_p3 - 7*u_p4 + u_p5)
                                idx = -3
                            else
                                idx = -2
                            end
                        end
                    else
                        if abs(-6*u_m1 + u_m2 + 15*u_0 - 20*u_p1 + 15*u_p2 - 6*u_p3 + u_p4) < abs(u_m1 - 6*u_0 + 15*u_p1 - 20*u_p2 + 15*u_p3 - 6*u_p4 + u_p5)
                            if abs(-21*u_m1 + 7*u_m2 - u_m3 + 35*u_0 - 35*u_p1 + 21*u_p2 - 7*u_p3 + u_p4) < abs(7*u_m1 - u_m2 - 21*u_0 + 35*u_p1 - 35*u_p2 + 21*u_p3 - 7*u_p4 + u_p5)
                                idx = -3
                            else
                                idx = -2
                            end
                        else
                            if abs(7*u_m1 - u_m2 - 21*u_0 + 35*u_p1 - 35*u_p2 + 21*u_p3 - 7*u_p4 + u_p5) < abs(-u_m1 + 7*u_0 - 21*u_p1 + 35*u_p2 - 35*u_p3 + 21*u_p4 - 7*u_p5 + u_p6)
                                idx = -2
                            else
                                idx = -1
                            end
                        end
                    end
                end
            end
        else
            if abs(-u_m1 + 3*u_0 - 3*u_p1 + u_p2) < abs(-u_0 + 3*u_p1 - 3*u_p2 + u_p3)
                if abs(-4*u_m1 + u_m2 + 6*u_0 - 4*u_p1 + u_p2) < abs(u_m1 - 4*u_0 + 6*u_p1 - 4*u_p2 + u_p3)
                    if abs(-10*u_m1 + 5*u_m2 - u_m3 + 10*u_0 - 5*u_p1 + u_p2) < abs(5*u_m1 - u_m2 - 10*u_0 + 10*u_p1 - 5*u_p2 + u_p3)
                        if abs(-20*u_m1 + 15*u_m2 - 6*u_m3 + u_m4 + 15*u_0 - 6*u_p1 + u_p2) < abs(15*u_m1 - 6*u_m2 + u_m3 - 20*u_0 + 15*u_p1 - 6*u_p2 + u_p3)
                            if abs(-35*u_m1 + 35*u_m2 - 21*u_m3 + 7*u_m4 - u_m5 + 21*u_0 - 7*u_p1 + u_p2) < abs(35*u_m1 - 21*u_m2 + 7*u_m3 - u_m4 - 35*u_0 + 21*u_p1 - 7*u_p2 + u_p3)
                                idx = -5
                            else
                                idx = -4
                            end
                        else
                            if abs(35*u_m1 - 21*u_m2 + 7*u_m3 - u_m4 - 35*u_0 + 21*u_p1 - 7*u_p2 + u_p3) < abs(-21*u_m1 + 7*u_m2 - u_m3 + 35*u_0 - 35*u_p1 + 21*u_p2 - 7*u_p3 + u_p4)
                                idx = -4
                            else
                                idx = -3
                            end
                        end
                    else
                        if abs(15*u_m1 - 6*u_m2 + u_m3 - 20*u_0 + 15*u_p1 - 6*u_p2 + u_p3) < abs(-6*u_m1 + u_m2 + 15*u_0 - 20*u_p1 + 15*u_p2 - 6*u_p3 + u_p4)
                            if abs(35*u_m1 - 21*u_m2 + 7*u_m3 - u_m4 - 35*u_0 + 21*u_p1 - 7*u_p2 + u_p3) < abs(-21*u_m1 + 7*u_m2 - u_m3 + 35*u_0 - 35*u_p1 + 21*u_p2 - 7*u_p3 + u_p4)
                                idx = -4
                            else
                                idx = -3
                            end
                        else
                            if abs(-21*u_m1 + 7*u_m2 - u_m3 + 35*u_0 - 35*u_p1 + 21*u_p2 - 7*u_p3 + u_p4) < abs(7*u_m1 - u_m2 - 21*u_0 + 35*u_p1 - 35*u_p2 + 21*u_p3 - 7*u_p4 + u_p5)
                                idx = -3
                            else
                                idx = -2
                            end
                        end
                    end
                else
                    if abs(5*u_m1 - u_m2 - 10*u_0 + 10*u_p1 - 5*u_p2 + u_p3) < abs(-u_m1 + 5*u_0 - 10*u_p1 + 10*u_p2 - 5*u_p3 + u_p4)
                        if abs(15*u_m1 - 6*u_m2 + u_m3 - 20*u_0 + 15*u_p1 - 6*u_p2 + u_p3) < abs(-6*u_m1 + u_m2 + 15*u_0 - 20*u_p1 + 15*u_p2 - 6*u_p3 + u_p4)
                            if abs(35*u_m1 - 21*u_m2 + 7*u_m3 - u_m4 - 35*u_0 + 21*u_p1 - 7*u_p2 + u_p3) < abs(-21*u_m1 + 7*u_m2 - u_m3 + 35*u_0 - 35*u_p1 + 21*u_p2 - 7*u_p3 + u_p4)
                                idx = -4
                            else
                                idx = -3
                            end
                        else
                            if abs(-21*u_m1 + 7*u_m2 - u_m3 + 35*u_0 - 35*u_p1 + 21*u_p2 - 7*u_p3 + u_p4) < abs(7*u_m1 - u_m2 - 21*u_0 + 35*u_p1 - 35*u_p2 + 21*u_p3 - 7*u_p4 + u_p5)
                                idx = -3
                            else
                                idx = -2
                            end
                        end
                    else
                        if abs(-6*u_m1 + u_m2 + 15*u_0 - 20*u_p1 + 15*u_p2 - 6*u_p3 + u_p4) < abs(u_m1 - 6*u_0 + 15*u_p1 - 20*u_p2 + 15*u_p3 - 6*u_p4 + u_p5)
                            if abs(-21*u_m1 + 7*u_m2 - u_m3 + 35*u_0 - 35*u_p1 + 21*u_p2 - 7*u_p3 + u_p4) < abs(7*u_m1 - u_m2 - 21*u_0 + 35*u_p1 - 35*u_p2 + 21*u_p3 - 7*u_p4 + u_p5)
                                idx = -3
                            else
                                idx = -2
                            end
                        else
                            if abs(7*u_m1 - u_m2 - 21*u_0 + 35*u_p1 - 35*u_p2 + 21*u_p3 - 7*u_p4 + u_p5) < abs(-u_m1 + 7*u_0 - 21*u_p1 + 35*u_p2 - 35*u_p3 + 21*u_p4 - 7*u_p5 + u_p6)
                                idx = -2
                            else
                                idx = -1
                            end
                        end
                    end
                end
            else
                if abs(u_m1 - 4*u_0 + 6*u_p1 - 4*u_p2 + u_p3) < abs(u_0 - 4*u_p1 + 6*u_p2 - 4*u_p3 + u_p4)
                    if abs(5*u_m1 - u_m2 - 10*u_0 + 10*u_p1 - 5*u_p2 + u_p3) < abs(-u_m1 + 5*u_0 - 10*u_p1 + 10*u_p2 - 5*u_p3 + u_p4)
                        if abs(15*u_m1 - 6*u_m2 + u_m3 - 20*u_0 + 15*u_p1 - 6*u_p2 + u_p3) < abs(-6*u_m1 + u_m2 + 15*u_0 - 20*u_p1 + 15*u_p2 - 6*u_p3 + u_p4)
                            if abs(35*u_m1 - 21*u_m2 + 7*u_m3 - u_m4 - 35*u_0 + 21*u_p1 - 7*u_p2 + u_p3) < abs(-21*u_m1 + 7*u_m2 - u_m3 + 35*u_0 - 35*u_p1 + 21*u_p2 - 7*u_p3 + u_p4)
                                idx = -4
                            else
                                idx = -3
                            end
                        else
                            if abs(-21*u_m1 + 7*u_m2 - u_m3 + 35*u_0 - 35*u_p1 + 21*u_p2 - 7*u_p3 + u_p4) < abs(7*u_m1 - u_m2 - 21*u_0 + 35*u_p1 - 35*u_p2 + 21*u_p3 - 7*u_p4 + u_p5)
                                idx = -3
                            else
                                idx = -2
                            end
                        end
                    else
                        if abs(-6*u_m1 + u_m2 + 15*u_0 - 20*u_p1 + 15*u_p2 - 6*u_p3 + u_p4) < abs(u_m1 - 6*u_0 + 15*u_p1 - 20*u_p2 + 15*u_p3 - 6*u_p4 + u_p5)
                            if abs(-21*u_m1 + 7*u_m2 - u_m3 + 35*u_0 - 35*u_p1 + 21*u_p2 - 7*u_p3 + u_p4) < abs(7*u_m1 - u_m2 - 21*u_0 + 35*u_p1 - 35*u_p2 + 21*u_p3 - 7*u_p4 + u_p5)
                                idx = -3
                            else
                                idx = -2
                            end
                        else
                            if abs(7*u_m1 - u_m2 - 21*u_0 + 35*u_p1 - 35*u_p2 + 21*u_p3 - 7*u_p4 + u_p5) < abs(-u_m1 + 7*u_0 - 21*u_p1 + 35*u_p2 - 35*u_p3 + 21*u_p4 - 7*u_p5 + u_p6)
                                idx = -2
                            else
                                idx = -1
                            end
                        end
                    end
                else
                    if abs(-u_m1 + 5*u_0 - 10*u_p1 + 10*u_p2 - 5*u_p3 + u_p4) < abs(-u_0 + 5*u_p1 - 10*u_p2 + 10*u_p3 - 5*u_p4 + u_p5)
                        if abs(-6*u_m1 + u_m2 + 15*u_0 - 20*u_p1 + 15*u_p2 - 6*u_p3 + u_p4) < abs(u_m1 - 6*u_0 + 15*u_p1 - 20*u_p2 + 15*u_p3 - 6*u_p4 + u_p5)
                            if abs(-21*u_m1 + 7*u_m2 - u_m3 + 35*u_0 - 35*u_p1 + 21*u_p2 - 7*u_p3 + u_p4) < abs(7*u_m1 - u_m2 - 21*u_0 + 35*u_p1 - 35*u_p2 + 21*u_p3 - 7*u_p4 + u_p5)
                                idx = -3
                            else
                                idx = -2
                            end
                        else
                            if abs(7*u_m1 - u_m2 - 21*u_0 + 35*u_p1 - 35*u_p2 + 21*u_p3 - 7*u_p4 + u_p5) < abs(-u_m1 + 7*u_0 - 21*u_p1 + 35*u_p2 - 35*u_p3 + 21*u_p4 - 7*u_p5 + u_p6)
                                idx = -2
                            else
                                idx = -1
                            end
                        end
                    else
                        if abs(u_m1 - 6*u_0 + 15*u_p1 - 20*u_p2 + 15*u_p3 - 6*u_p4 + u_p5) < abs(u_0 - 6*u_p1 + 15*u_p2 - 20*u_p3 + 15*u_p4 - 6*u_p5 + u_p6)
                            if abs(7*u_m1 - u_m2 - 21*u_0 + 35*u_p1 - 35*u_p2 + 21*u_p3 - 7*u_p4 + u_p5) < abs(-u_m1 + 7*u_0 - 21*u_p1 + 35*u_p2 - 35*u_p3 + 21*u_p4 - 7*u_p5 + u_p6)
                                idx = -2
                            else
                                idx = -1
                            end
                        else
                            if abs(-u_m1 + 7*u_0 - 21*u_p1 + 35*u_p2 - 35*u_p3 + 21*u_p4 - 7*u_p5 + u_p6) < abs(-u_0 + 7*u_p1 - 21*u_p2 + 35*u_p3 - 35*u_p4 + 21*u_p5 - 7*u_p6 + u_p7)
                                idx = -1
                            else
                                idx = 0
                            end
                        end
                    end
                end
            end
        end
    end

    idx
end


################################################################################

"""
    BiasedENOChoice

Biased choice of the stencil in an ENO reconstruction, minimising the absolute
values of the divided differences lexicographically with a bias towards the
central stencil(s).
"""
struct BiasedENOChoice end

@inline function (::BiasedENOChoice)(u_m1, u_0, u_p1)
    if abs(-u_m1 + u_0) < abs(-u_0 + u_p1)
        idx = -1
    else
        idx = 0
    end

    idx
end

@inline function (::BiasedENOChoice)(u_m2, u_m1, u_0, u_p1, u_p2)
    if abs(-u_m1 + u_0) < abs(-u_0 + u_p1)
        if 2*abs(-2*u_m1 + u_m2 + u_0) < abs(u_m1 - 2*u_0 + u_p1)
            idx = -2
        else
            idx = -1
        end
    else
        if abs(u_m1 - 2*u_0 + u_p1) < 2*abs(u_0 - 2*u_p1 + u_p2)
            idx = -1
        else
            idx = 0
        end
    end

    idx
end

@inline function (::BiasedENOChoice)(u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3)
    if abs(-u_m1 + u_0) < abs(-u_0 + u_p1)
        if 2*abs(-2*u_m1 + u_m2 + u_0) < abs(u_m1 - 2*u_0 + u_p1)
            if 2*abs(-3*u_m1 + 3*u_m2 - u_m3 + u_0) < abs(3*u_m1 - u_m2 - 3*u_0 + u_p1)
                idx = -3
            else
                idx = -2
            end
        else
            if abs(3*u_m1 - u_m2 - 3*u_0 + u_p1) < abs(-u_m1 + 3*u_0 - 3*u_p1 + u_p2)
                idx = -2
            else
                idx = -1
            end
        end
    else
        if abs(u_m1 - 2*u_0 + u_p1) < 2*abs(u_0 - 2*u_p1 + u_p2)
            if abs(3*u_m1 - u_m2 - 3*u_0 + u_p1) < abs(-u_m1 + 3*u_0 - 3*u_p1 + u_p2)
                idx = -2
            else
                idx = -1
            end
        else
            if abs(-u_m1 + 3*u_0 - 3*u_p1 + u_p2) < 2*abs(-u_0 + 3*u_p1 - 3*u_p2 + u_p3)
                idx = -1
            else
                idx = 0
            end
        end
    end

    idx
end

@inline function (::BiasedENOChoice)(u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4)
    if abs(-u_m1 + u_0) < abs(-u_0 + u_p1)
        if 2*abs(-2*u_m1 + u_m2 + u_0) < abs(u_m1 - 2*u_0 + u_p1)
            if 2*abs(-3*u_m1 + 3*u_m2 - u_m3 + u_0) < abs(3*u_m1 - u_m2 - 3*u_0 + u_p1)
                if 2*abs(-4*u_m1 + 6*u_m2 - 4*u_m3 + u_m4 + u_0) < abs(6*u_m1 - 4*u_m2 + u_m3 - 4*u_0 + u_p1)
                    idx = -4
                else
                    idx = -3
                end
            else
                if 2*abs(6*u_m1 - 4*u_m2 + u_m3 - 4*u_0 + u_p1) < abs(-4*u_m1 + u_m2 + 6*u_0 - 4*u_p1 + u_p2)
                    idx = -3
                else
                    idx = -2
                end
            end
        else
            if abs(3*u_m1 - u_m2 - 3*u_0 + u_p1) < abs(-u_m1 + 3*u_0 - 3*u_p1 + u_p2)
                if 2*abs(6*u_m1 - 4*u_m2 + u_m3 - 4*u_0 + u_p1) < abs(-4*u_m1 + u_m2 + 6*u_0 - 4*u_p1 + u_p2)
                    idx = -3
                else
                    idx = -2
                end
            else
                if abs(-4*u_m1 + u_m2 + 6*u_0 - 4*u_p1 + u_p2) < 2*abs(u_m1 - 4*u_0 + 6*u_p1 - 4*u_p2 + u_p3)
                    idx = -2
                else
                    idx = -1
                end
            end
        end
    else
        if abs(u_m1 - 2*u_0 + u_p1) < 2*abs(u_0 - 2*u_p1 + u_p2)
            if abs(3*u_m1 - u_m2 - 3*u_0 + u_p1) < abs(-u_m1 + 3*u_0 - 3*u_p1 + u_p2)
                if 2*abs(6*u_m1 - 4*u_m2 + u_m3 - 4*u_0 + u_p1) < abs(-4*u_m1 + u_m2 + 6*u_0 - 4*u_p1 + u_p2)
                    idx = -3
                else
                    idx = -2
                end
            else
                if abs(-4*u_m1 + u_m2 + 6*u_0 - 4*u_p1 + u_p2) < 2*abs(u_m1 - 4*u_0 + 6*u_p1 - 4*u_p2 + u_p3)
                    idx = -2
                else
                    idx = -1
                end
            end
        else
            if abs(-u_m1 + 3*u_0 - 3*u_p1 + u_p2) < 2*abs(-u_0 + 3*u_p1 - 3*u_p2 + u_p3)
                if abs(-4*u_m1 + u_m2 + 6*u_0 - 4*u_p1 + u_p2) < 2*abs(u_m1 - 4*u_0 + 6*u_p1 - 4*u_p2 + u_p3)
                    idx = -2
                else
                    idx = -1
                end
            else
                if abs(u_m1 - 4*u_0 + 6*u_p1 - 4*u_p2 + u_p3) < 2*abs(u_0 - 4*u_p1 + 6*u_p2 - 4*u_p3 + u_p4)
                    idx = -1
                else
                    idx = 0
                end
            end
        end
    end

    idx
end


################################################################################

"""
    MinL2Choice

Choice of the stencil in a modified ENO reconstruction minimising the L² norm of
the reconstructed polynomial.
"""
struct MinL2Choice end

@inline function (::MinL2Choice)(u_m1, u_0, u_p1)
    idx = -1
    val_old = (-u_m1 + u_0)^2 # + 4*u_0^2

    val_new = (-u_0 + u_p1)^2 # + 4*u_0^2
    if val_new < val_old
        idx = 0
        val_old = val_new
    end

    idx
end

@inline function (::MinL2Choice)(u_m2, u_m1, u_0, u_p1, u_p2)
    idx = -2
    val_old = (-2*u_m1 + u_m2 + u_0)^2 + (-12*u_m1 + 3*u_m2 + 9*u_0)^2 # + 144*u_0^2

    val_new = (u_m1 - 2*u_0 + u_p1)^2 + (-3*u_m1 + 3*u_p1)^2 # + 144*u_0^2
    if val_new < val_old
        idx = -1
        val_old = val_new
    end

    val_new = (u_0 - 2*u_p1 + u_p2)^2 + (-9*u_0 + 12*u_p1 - 3*u_p2)^2 # + 144*u_0^2
    if val_new < val_old
        idx = 0
        val_old = val_new
    end

    idx
end

@inline function (::MinL2Choice)(u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3)
    idx = -3
    val_old = (-3*u_m1 + 3*u_m2 - u_m3 + u_0)^2 + (-50*u_m1 + 40*u_m2 - 10*u_m3 + 20*u_0)^2 + (-177*u_m1 + 87*u_m2 - 19*u_m3 + 109*u_0)^2 # + 14400*u_0^2

    val_new = (3*u_m1 - u_m2 - 3*u_0 + u_p1)^2 + (10*u_m1 - 20*u_0 + 10*u_p1)^2 + (-63*u_m1 + 11*u_m2 + 33*u_0 + 19*u_p1)^2 # + 14400*u_0^2
    if val_new < val_old
        idx = -2
        val_old = val_new
    end

    val_new = (-u_m1 + 3*u_0 - 3*u_p1 + u_p2)^2 + (10*u_m1 - 20*u_0 + 10*u_p1)^2 + (-19*u_m1 - 33*u_0 + 63*u_p1 - 11*u_p2)^2 # + 14400*u_0^2
    if val_new < val_old
        idx = -1
        val_old = val_new
    end

    val_new = (-u_0 + 3*u_p1 - 3*u_p2 + u_p3)^2 + (20*u_0 - 50*u_p1 + 40*u_p2 - 10*u_p3)^2 + (-109*u_0 + 177*u_p1 - 87*u_p2 + 19*u_p3)^2 # + 14400*u_0^2
    if val_new < val_old
        idx = 0
        val_old = val_new
    end

    idx
end

@inline function (::MinL2Choice)(u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4)
    idx = -4
    val_old = (-4*u_m1 + 6*u_m2 - 4*u_m3 + u_m4 + u_0)^2 + (-126*u_m1 + 168*u_m2 - 98*u_m3 + 21*u_m4 + 35*u_0)^2 + (-1200*u_m1 + 1310*u_m2 - 640*u_m3 + 125*u_m4 + 405*u_0)^2 + (-3234*u_m1 + 2352*u_m2 - 1022*u_m3 + 189*u_m4 + 1715*u_0)^2 # + 2822400*u_0^2

    val_new = (6*u_m1 - 4*u_m2 + u_m3 - 4*u_0 + u_p1)^2 + (84*u_m1 - 42*u_m2 + 7*u_m3 - 70*u_0 + 21*u_p1)^2 + (50*u_m1 + 60*u_m2 - 15*u_m3 - 220*u_0 + 125*u_p1)^2 + (-1344*u_m1 + 462*u_m2 - 77*u_m3 + 770*u_0 + 189*u_p1)^2 # + 2822400*u_0^2
    if val_new < val_old
        idx = -3
        val_old = val_new
    end

    val_new = (-4*u_m1 + u_m2 + 6*u_0 - 4*u_p1 + u_p2)^2 + (14*u_m1 - 7*u_m2 - 14*u_p1 + 7*u_p2)^2 + (200*u_m1 - 15*u_m2 - 370*u_0 + 200*u_p1 - 15*u_p2)^2 + (-574*u_m1 + 77*u_m2 + 574*u_p1 - 77*u_p2)^2 # + 2822400*u_0^2
    if val_new < val_old
        idx = -2
        val_old = val_new
    end

    val_new = (u_m1 - 4*u_0 + 6*u_p1 - 4*u_p2 + u_p3)^2 + (-21*u_m1 + 70*u_0 - 84*u_p1 + 42*u_p2 - 7*u_p3)^2 + (125*u_m1 - 220*u_0 + 50*u_p1 + 60*u_p2 - 15*u_p3)^2 + (-189*u_m1 - 770*u_0 + 1344*u_p1 - 462*u_p2 + 77*u_p3)^2 # + 2822400*u_0^2
    if val_new < val_old
        idx = -1
        val_old = val_new
    end

    val_new = (u_0 - 4*u_p1 + 6*u_p2 - 4*u_p3 + u_p4)^2 + (-35*u_0 + 126*u_p1 - 168*u_p2 + 98*u_p3 - 21*u_p4)^2 + (405*u_0 - 1200*u_p1 + 1310*u_p2 - 640*u_p3 + 125*u_p4)^2 + (-1715*u_0 + 3234*u_p1 - 2352*u_p2 + 1022*u_p3 - 189*u_p4)^2 # + 2822400*u_0^2
    if val_new < val_old
        idx = 0
        val_old = val_new
    end

    idx
end

@inline function (::MinL2Choice)(u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5)
    idx = -5
    val_old = (-5*u_m1 + 10*u_m2 - 10*u_m3 + 5*u_m4 - u_m5 + u_0)^2 + (-252*u_m1 + 468*u_m2 - 432*u_m3 + 198*u_m4 - 36*u_m5 + 54*u_0)^2 + (-4438*u_m1 + 7364*u_m2 - 6104*u_m3 + 2548*u_m4 - 434*u_m5 + 1064*u_0)^2 + (-31500*u_m1 + 43380*u_m2 - 31320*u_m3 + 12150*u_m4 - 1980*u_m5 + 9270*u_0)^2 + (-71157*u_m1 + 68226*u_m2 - 44286*u_m3 + 16347*u_m4 - 2589*u_m5 + 33459*u_0)^2 # + 914457600*u_0^2

    val_new = (10*u_m1 - 10*u_m2 + 5*u_m3 - u_m4 - 5*u_0 + u_p1)^2 + (288*u_m1 - 252*u_m2 + 108*u_m3 - 18*u_m4 - 162*u_0 + 36*u_p1)^2 + (2072*u_m1 - 1316*u_m2 + 406*u_m3 - 56*u_m4 - 1540*u_0 + 434*u_p1)^2 + (-1800*u_m1 + 3780*u_m2 - 1620*u_m3 + 270*u_m4 - 2610*u_0 + 1980*u_p1)^2 + (-32322*u_m1 + 16446*u_m2 - 5451*u_m3 + 813*u_m4 + 17925*u_0 + 2589*u_p1)^2 # + 914457600*u_0^2
    if val_new < val_old
        idx = -4
        val_old = val_new
    end

    val_new = (-10*u_m1 + 5*u_m2 - u_m3 + 10*u_0 - 5*u_p1 + u_p2)^2 + (-72*u_m1 + 18*u_m2 + 108*u_0 - 72*u_p1 + 18*u_p2)^2 + (952*u_m1 - 476*u_m2 + 70*u_m3 - 700*u_0 + 98*u_p1 + 56*u_p2)^2 + (3600*u_m1 - 270*u_m2 - 6660*u_0 + 3600*u_p1 - 270*u_p2)^2 + (-16062*u_m1 + 4251*u_m2 - 573*u_m3 + 5730*u_0 + 7467*u_p1 - 813*u_p2)^2 # + 914457600*u_0^2
    if val_new < val_old
        idx = -3
        val_old = val_new
    end

    val_new = (5*u_m1 - u_m2 - 10*u_0 + 10*u_p1 - 5*u_p2 + u_p3)^2 + (-72*u_m1 + 18*u_m2 + 108*u_0 - 72*u_p1 + 18*u_p2)^2 + (-98*u_m1 - 56*u_m2 + 700*u_0 - 952*u_p1 + 476*u_p2 - 70*u_p3)^2 + (3600*u_m1 - 270*u_m2 - 6660*u_0 + 3600*u_p1 - 270*u_p2)^2 + (-7467*u_m1 + 813*u_m2 - 5730*u_0 + 16062*u_p1 - 4251*u_p2 + 573*u_p3)^2 # + 914457600*u_0^2
    if val_new < val_old
        idx = -2
        val_old = val_new
    end

    val_new = (-u_m1 + 5*u_0 - 10*u_p1 + 10*u_p2 - 5*u_p3 + u_p4)^2 + (36*u_m1 - 162*u_0 + 288*u_p1 - 252*u_p2 + 108*u_p3 - 18*u_p4)^2 + (-434*u_m1 + 1540*u_0 - 2072*u_p1 + 1316*u_p2 - 406*u_p3 + 56*u_p4)^2 + (1980*u_m1 - 2610*u_0 - 1800*u_p1 + 3780*u_p2 - 1620*u_p3 + 270*u_p4)^2 + (-2589*u_m1 - 17925*u_0 + 32322*u_p1 - 16446*u_p2 + 5451*u_p3 - 813*u_p4)^2 # + 914457600*u_0^2
    if val_new < val_old
        idx = -1
        val_old = val_new
    end

    val_new = (-u_0 + 5*u_p1 - 10*u_p2 + 10*u_p3 - 5*u_p4 + u_p5)^2 + (54*u_0 - 252*u_p1 + 468*u_p2 - 432*u_p3 + 198*u_p4 - 36*u_p5)^2 + (-1064*u_0 + 4438*u_p1 - 7364*u_p2 + 6104*u_p3 - 2548*u_p4 + 434*u_p5)^2 + (9270*u_0 - 31500*u_p1 + 43380*u_p2 - 31320*u_p3 + 12150*u_p4 - 1980*u_p5)^2 + (-33459*u_0 + 71157*u_p1 - 68226*u_p2 + 44286*u_p3 - 16347*u_p4 + 2589*u_p5)^2 # + 914457600*u_0^2
    if val_new < val_old
        idx = 0
        val_old = val_new
    end

    idx
end

@inline function (::MinL2Choice)(u_m6, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, u_p6)
    idx = -6
    val_old = (-6*u_m1 + 15*u_m2 - 20*u_m3 + 15*u_m4 - 6*u_m5 + u_m6 + u_0)^2 + (-440*u_m1 + 1045*u_m2 - 1320*u_m3 + 935*u_m4 - 352*u_m5 + 55*u_m6 + 77*u_0)^2 + (-12204*u_m1 + 26946*u_m2 - 31704*u_m3 + 21006*u_m4 - 7452*u_m5 + 1110*u_m6 + 2298*u_0)^2 + (-157696*u_m1 + 312158*u_m2 - 334488*u_m3 + 206206*u_m4 - 69608*u_m5 + 10010*u_m6 + 33418*u_0)^2 + (-923934*u_m1 + 1531695*u_m2 - 1458820*u_m3 + 844635*u_m4 - 274494*u_m5 + 38489*u_m6 + 242429*u_0)^2 + (-1837704*u_m1 + 2181597*u_m2 - 1881792*u_m3 + 1040259*u_m4 - 329208*u_m5 + 45375*u_m6 + 781473*u_0)^2 # + 442597478400*u_0^2

    val_new = (15*u_m1 - 20*u_m2 + 15*u_m3 - 6*u_m4 + u_m5 - 6*u_0 + u_p1)^2 + (715*u_m1 - 880*u_m2 + 605*u_m3 - 220*u_m4 + 33*u_m5 - 308*u_0 + 55*u_p1)^2 + (11106*u_m1 - 11904*u_m2 + 7146*u_m3 - 2304*u_m4 + 318*u_m5 - 5472*u_0 + 1110*u_p1)^2 + (52514*u_m1 - 38192*u_m2 + 15862*u_m3 - 4004*u_m4 + 462*u_m5 - 36652*u_0 + 10010*u_p1)^2 + (-115665*u_m1 + 184580*u_m2 - 111705*u_m3 + 36366*u_m4 - 5071*u_m5 - 26994*u_0 + 38489*u_p1)^2 + (-884829*u_m1 + 593472*u_m2 - 293667*u_m3 + 87384*u_m4 - 11583*u_m5 + 463848*u_0 + 45375*u_p1)^2 # + 442597478400*u_0^2
    if val_new < val_old
        idx = -5
        val_old = val_new
    end

    val_new = (-20*u_m1 + 15*u_m2 - 6*u_m3 + u_m4 + 15*u_0 - 6*u_p1 + u_p2)^2 + (-440*u_m1 + 275*u_m2 - 88*u_m3 + 11*u_m4 + 385*u_0 - 176*u_p1 + 33*u_p2)^2 + (-24*u_m1 - 774*u_m2 + 468*u_m3 - 78*u_m4 + 1206*u_0 - 1116*u_p1 + 318*u_p2)^2 + (36344*u_m1 - 22022*u_m2 + 6160*u_m3 - 770*u_m4 - 26950*u_0 + 6776*u_p1 + 462*u_p2)^2 + (61820*u_m1 + 7095*u_m2 - 5214*u_m3 + 869*u_m4 - 133485*u_0 + 73986*u_p1 - 5071*u_p2)^2 + (-479424*u_m1 + 188067*u_m2 - 50424*u_m3 + 6303*u_m4 + 220605*u_0 + 126456*u_p1 - 11583*u_p2)^2 # + 442597478400*u_0^2
    if val_new < val_old
        idx = -4
        val_old = val_new
    end

    val_new = (15*u_m1 - 6*u_m2 + u_m3 - 20*u_0 + 15*u_p1 - 6*u_p2 + u_p3)^2 + (-55*u_m1 + 44*u_m2 - 11*u_m3 + 55*u_p1 - 44*u_p2 + 11*u_p3)^2 + (-2754*u_m1 + 864*u_m2 - 78*u_m3 + 3936*u_0 - 2754*u_p1 + 864*u_p2 - 78*u_p3)^2 + (9394*u_m1 - 5852*u_m2 + 770*u_m3 - 9394*u_p1 + 5852*u_p2 - 770*u_p3)^2 + (92235*u_m1 - 11154*u_m2 + 869*u_m3 - 163900*u_0 + 92235*u_p1 - 11154*u_p2 + 869*u_p3)^2 + (-258819*u_m1 + 55704*u_m2 - 6303*u_m3 + 258819*u_p1 - 55704*u_p2 + 6303*u_p3)^2 # + 442597478400*u_0^2
    if val_new < val_old
        idx = -3
        val_old = val_new
    end

    val_new = (-6*u_m1 + u_m2 + 15*u_0 - 20*u_p1 + 15*u_p2 - 6*u_p3 + u_p4)^2 + (176*u_m1 - 33*u_m2 - 385*u_0 + 440*u_p1 - 275*u_p2 + 88*u_p3 - 11*u_p4)^2 + (-1116*u_m1 + 318*u_m2 + 1206*u_0 - 24*u_p1 - 774*u_p2 + 468*u_p3 - 78*u_p4)^2 + (-6776*u_m1 - 462*u_m2 + 26950*u_0 - 36344*u_p1 + 22022*u_p2 - 6160*u_p3 + 770*u_p4)^2 + (73986*u_m1 - 5071*u_m2 - 133485*u_0 + 61820*u_p1 + 7095*u_p2 - 5214*u_p3 + 869*u_p4)^2 + (-126456*u_m1 + 11583*u_m2 - 220605*u_0 + 479424*u_p1 - 188067*u_p2 + 50424*u_p3 - 6303*u_p4)^2 # + 442597478400*u_0^2
    if val_new < val_old
        idx = -2
        val_old = val_new
    end

    val_new = (u_m1 - 6*u_0 + 15*u_p1 - 20*u_p2 + 15*u_p3 - 6*u_p4 + u_p5)^2 + (-55*u_m1 + 308*u_0 - 715*u_p1 + 880*u_p2 - 605*u_p3 + 220*u_p4 - 33*u_p5)^2 + (1110*u_m1 - 5472*u_0 + 11106*u_p1 - 11904*u_p2 + 7146*u_p3 - 2304*u_p4 + 318*u_p5)^2 + (-10010*u_m1 + 36652*u_0 - 52514*u_p1 + 38192*u_p2 - 15862*u_p3 + 4004*u_p4 - 462*u_p5)^2 + (38489*u_m1 - 26994*u_0 - 115665*u_p1 + 184580*u_p2 - 111705*u_p3 + 36366*u_p4 - 5071*u_p5)^2 + (-45375*u_m1 - 463848*u_0 + 884829*u_p1 - 593472*u_p2 + 293667*u_p3 - 87384*u_p4 + 11583*u_p5)^2 # + 442597478400*u_0^2
    if val_new < val_old
        idx = -1
        val_old = val_new
    end

    val_new = (u_0 - 6*u_p1 + 15*u_p2 - 20*u_p3 + 15*u_p4 - 6*u_p5 + u_p6)^2 + (-77*u_0 + 440*u_p1 - 1045*u_p2 + 1320*u_p3 - 935*u_p4 + 352*u_p5 - 55*u_p6)^2 + (2298*u_0 - 12204*u_p1 + 26946*u_p2 - 31704*u_p3 + 21006*u_p4 - 7452*u_p5 + 1110*u_p6)^2 + (-33418*u_0 + 157696*u_p1 - 312158*u_p2 + 334488*u_p3 - 206206*u_p4 + 69608*u_p5 - 10010*u_p6)^2 + (242429*u_0 - 923934*u_p1 + 1531695*u_p2 - 1458820*u_p3 + 844635*u_p4 - 274494*u_p5 + 38489*u_p6)^2 + (-781473*u_0 + 1837704*u_p1 - 2181597*u_p2 + 1881792*u_p3 - 1040259*u_p4 + 329208*u_p5 - 45375*u_p6)^2 # + 442597478400*u_0^2
    if val_new < val_old
        idx = 0
        val_old = val_new
    end

    idx
end

@inline function (::MinL2Choice)(u_m7, u_m6, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, u_p6, u_p7)
    idx = -7
    val_old = (-7*u_m1 + 21*u_m2 - 35*u_m3 + 35*u_m4 - 21*u_m5 + 7*u_m6 - u_m7 + u_0)^2 + (-702*u_m1 + 2028*u_m2 - 3250*u_m3 + 3120*u_m4 - 1794*u_m5 + 572*u_m6 - 78*u_m7 + 104*u_0)^2 + (-27995*u_m1 + 76835*u_m2 - 117095*u_m3 + 107085*u_m4 - 58817*u_m5 + 17985*u_m6 - 2365*u_m7 + 4367*u_0)^2 + (-563004*u_m1 + 1437696*u_m2 - 2052804*u_m3 + 1774656*u_m4 - 930852*u_m5 + 274560*u_m6 - 35100*u_m7 + 94848*u_0)^2 + (-29673917*u_m1/5 + 68100851*u_m2/5 - 17870125*u_m3 + 14534793*u_m4 - 36569351*u_m5/5 + 10474737*u_m6/5 - 1310491*u_m7/5 + 5654831*u_0/5)^2 + (-30262518*u_m1 + 58544772*u_m2 - 69130490*u_m3 + 53161680*u_m4 - 25857546*u_m5 + 7240948*u_m6 - 891462*u_m7 + 7194616*u_0)^2 + (-272888473*u_m1/5 + 385568469*u_m2/5 - 82913545*u_m3 + 61033687*u_m4 - 144757899*u_m5/5 + 39885703*u_m6/5 - 4855279*u_m7/5 + 106446769*u_0/5)^2 # + 299195895398400*u_0^2

    val_new = (21*u_m1 - 35*u_m2 + 35*u_m3 - 21*u_m4 + 7*u_m5 - u_m6 - 7*u_0 + u_p1)^2 + (1482*u_m1 - 2340*u_m2 + 2210*u_m3 - 1248*u_m4 + 390*u_m5 - 52*u_m6 - 520*u_0 + 78*u_p1)^2 + (38225*u_m1 - 55605*u_m2 + 48455*u_m3 - 25355*u_m4 + 7403*u_m5 - 935*u_m6 - 14553*u_0 + 2365*u_p1)^2 + (419796*u_m1 - 527904*u_m2 + 404196*u_m3 - 190944*u_m4 + 51948*u_m5 - 6240*u_m6 - 185952*u_0 + 35100*u_p1)^2 + (7019831*u_m1/5 - 1057329*u_m2 + 476749*u_m3 - 713531*u_m4/5 + 124397*u_m5/5 - 9191*u_m6/5 - 4829097*u_0/5 + 1310491*u_p1/5)^2 + (-5301582*u_m1 + 8622900*u_m2 - 6728150*u_m3 + 3239808*u_m4 - 896610*u_m5 + 109252*u_m6 + 62920*u_0 + 891462*u_p1)^2 + (-136940661*u_m1/5 + 22734569*u_m2 - 14939639*u_m3 + 33272811*u_m4/5 - 8810087*u_m5/5 + 1043471*u_m6/5 + 67604537*u_0/5 + 4855279*u_p1/5)^2 # + 299195895398400*u_0^2
    if val_new < val_old
        idx = -6
        val_old = val_new
    end

    val_new = (-35*u_m1 + 35*u_m2 - 21*u_m3 + 7*u_m4 - u_m5 + 21*u_0 - 7*u_p1 + u_p2)^2 + (-1430*u_m1 + 1300*u_m2 - 702*u_m3 + 208*u_m4 - 26*u_m5 + 936*u_0 - 338*u_p1 + 52*u_p2)^2 + (-14135*u_m1 + 9845*u_m2 - 3905*u_m3 + 825*u_m4 - 77*u_m5 + 11627*u_0 - 5115*u_p1 + 935*u_p2)^2 + (70356*u_m1 - 91104*u_m2 + 54756*u_m3 - 16224*u_m4 + 2028*u_m5 - 11232*u_0 - 14820*u_p1 + 6240*u_p2)^2 + (1301027*u_m1 - 928655*u_m2 + 1869049*u_m3/5 - 456183*u_m4/5 + 50869*u_m5/5 - 4571749*u_0/5 + 1236963*u_p1/5 + 9191*u_p2/5)^2 + (816530*u_m1 + 975260*u_m2 - 610038*u_m3 + 180752*u_m4 - 22594*u_m5 - 2996136*u_0 + 1765478*u_p1 - 109252*u_p2)^2 + (-15701257*u_m1 + 8125975*u_m2 - 16263819*u_m3/5 + 4055623*u_m4/5 - 462319*u_m5/5 + 38387349*u_0/5 + 13203047*u_p1/5 - 1043471*u_p2/5)^2 # + 299195895398400*u_0^2
    if val_new < val_old
        idx = -5
        val_old = val_new
    end

    val_new = (35*u_m1 - 21*u_m2 + 7*u_m3 - u_m4 - 35*u_0 + 21*u_p1 - 7*u_p2 + u_p3)^2 + (390*u_m1 - 156*u_m2 + 26*u_m3 - 520*u_0 + 390*u_p1 - 156*u_p2 + 26*u_p3)^2 + (-8745*u_m1 + 5533*u_m2 - 1749*u_m3 + 209*u_m4 + 7315*u_0 - 2959*u_p1 + 319*u_p2 + 77*u_p3)^2 + (-71604*u_m1 + 22464*u_m2 - 2028*u_m3 + 102336*u_0 - 71604*u_p1 + 22464*u_p2 - 2028*u_p3)^2 + (588861*u_m1 - 1794611*u_m2/5 + 444717*u_m3/5 - 49231*u_m4/5 - 344617*u_0 - 187369*u_p1/5 + 416143*u_p2/5 - 50869*u_p3/5)^2 + (2398110*u_m1 - 290004*u_m2 + 22594*u_m3 - 4261400*u_0 + 2398110*u_p1 - 290004*u_p2 + 22594*u_p3)^2 + (-9228791*u_m1 + 14740011*u_m2/5 - 3318887*u_m3/5 + 357071*u_m4/5 + 2499497*u_0 + 26147979*u_p1/5 - 4742023*u_p2/5 + 462319*u_p3/5)^2 # + 299195895398400*u_0^2
    if val_new < val_old
        idx = -4
        val_old = val_new
    end

    val_new = (-21*u_m1 + 7*u_m2 - u_m3 + 35*u_0 - 35*u_p1 + 21*u_p2 - 7*u_p3 + u_p4)^2 + (390*u_m1 - 156*u_m2 + 26*u_m3 - 520*u_0 + 390*u_p1 - 156*u_p2 + 26*u_p3)^2 + (2959*u_m1 - 319*u_m2 - 77*u_m3 - 7315*u_0 + 8745*u_p1 - 5533*u_p2 + 1749*u_p3 - 209*u_p4)^2 + (-71604*u_m1 + 22464*u_m2 - 2028*u_m3 + 102336*u_0 - 71604*u_p1 + 22464*u_p2 - 2028*u_p3)^2 + (187369*u_m1/5 - 416143*u_m2/5 + 50869*u_m3/5 + 344617*u_0 - 588861*u_p1 + 1794611*u_p2/5 - 444717*u_p3/5 + 49231*u_p4/5)^2 + (2398110*u_m1 - 290004*u_m2 + 22594*u_m3 - 4261400*u_0 + 2398110*u_p1 - 290004*u_p2 + 22594*u_p3)^2 + (-26147979*u_m1/5 + 4742023*u_m2/5 - 462319*u_m3/5 - 2499497*u_0 + 9228791*u_p1 - 14740011*u_p2/5 + 3318887*u_p3/5 - 357071*u_p4/5)^2 # + 299195895398400*u_0^2
    if val_new < val_old
        idx = -3
        val_old = val_new
    end

    val_new = (7*u_m1 - u_m2 - 21*u_0 + 35*u_p1 - 35*u_p2 + 21*u_p3 - 7*u_p4 + u_p5)^2 + (-338*u_m1 + 52*u_m2 + 936*u_0 - 1430*u_p1 + 1300*u_p2 - 702*u_p3 + 208*u_p4 - 26*u_p5)^2 + (5115*u_m1 - 935*u_m2 - 11627*u_0 + 14135*u_p1 - 9845*u_p2 + 3905*u_p3 - 825*u_p4 + 77*u_p5)^2 + (-14820*u_m1 + 6240*u_m2 - 11232*u_0 + 70356*u_p1 - 91104*u_p2 + 54756*u_p3 - 16224*u_p4 + 2028*u_p5)^2 + (-1236963*u_m1/5 - 9191*u_m2/5 + 4571749*u_0/5 - 1301027*u_p1 + 928655*u_p2 - 1869049*u_p3/5 + 456183*u_p4/5 - 50869*u_p5/5)^2 + (1765478*u_m1 - 109252*u_m2 - 2996136*u_0 + 816530*u_p1 + 975260*u_p2 - 610038*u_p3 + 180752*u_p4 - 22594*u_p5)^2 + (-13203047*u_m1/5 + 1043471*u_m2/5 - 38387349*u_0/5 + 15701257*u_p1 - 8125975*u_p2 + 16263819*u_p3/5 - 4055623*u_p4/5 + 462319*u_p5/5)^2 # + 299195895398400*u_0^2
    if val_new < val_old
        idx = -2
        val_old = val_new
    end

    val_new = (-u_m1 + 7*u_0 - 21*u_p1 + 35*u_p2 - 35*u_p3 + 21*u_p4 - 7*u_p5 + u_p6)^2 + (78*u_m1 - 520*u_0 + 1482*u_p1 - 2340*u_p2 + 2210*u_p3 - 1248*u_p4 + 390*u_p5 - 52*u_p6)^2 + (-2365*u_m1 + 14553*u_0 - 38225*u_p1 + 55605*u_p2 - 48455*u_p3 + 25355*u_p4 - 7403*u_p5 + 935*u_p6)^2 + (35100*u_m1 - 185952*u_0 + 419796*u_p1 - 527904*u_p2 + 404196*u_p3 - 190944*u_p4 + 51948*u_p5 - 6240*u_p6)^2 + (-1310491*u_m1/5 + 4829097*u_0/5 - 7019831*u_p1/5 + 1057329*u_p2 - 476749*u_p3 + 713531*u_p4/5 - 124397*u_p5/5 + 9191*u_p6/5)^2 + (891462*u_m1 + 62920*u_0 - 5301582*u_p1 + 8622900*u_p2 - 6728150*u_p3 + 3239808*u_p4 - 896610*u_p5 + 109252*u_p6)^2 + (-4855279*u_m1/5 - 67604537*u_0/5 + 136940661*u_p1/5 - 22734569*u_p2 + 14939639*u_p3 - 33272811*u_p4/5 + 8810087*u_p5/5 - 1043471*u_p6/5)^2 # + 299195895398400*u_0^2
    if val_new < val_old
        idx = -1
        val_old = val_new
    end

    val_new = (-u_0 + 7*u_p1 - 21*u_p2 + 35*u_p3 - 35*u_p4 + 21*u_p5 - 7*u_p6 + u_p7)^2 + (104*u_0 - 702*u_p1 + 2028*u_p2 - 3250*u_p3 + 3120*u_p4 - 1794*u_p5 + 572*u_p6 - 78*u_p7)^2 + (-4367*u_0 + 27995*u_p1 - 76835*u_p2 + 117095*u_p3 - 107085*u_p4 + 58817*u_p5 - 17985*u_p6 + 2365*u_p7)^2 + (94848*u_0 - 563004*u_p1 + 1437696*u_p2 - 2052804*u_p3 + 1774656*u_p4 - 930852*u_p5 + 274560*u_p6 - 35100*u_p7)^2 + (-5654831*u_0/5 + 29673917*u_p1/5 - 68100851*u_p2/5 + 17870125*u_p3 - 14534793*u_p4 + 36569351*u_p5/5 - 10474737*u_p6/5 + 1310491*u_p7/5)^2 + (7194616*u_0 - 30262518*u_p1 + 58544772*u_p2 - 69130490*u_p3 + 53161680*u_p4 - 25857546*u_p5 + 7240948*u_p6 - 891462*u_p7)^2 + (-106446769*u_0/5 + 272888473*u_p1/5 - 385568469*u_p2/5 + 82913545*u_p3 - 61033687*u_p4 + 144757899*u_p5/5 - 39885703*u_p6/5 + 4855279*u_p7/5)^2 # + 299195895398400*u_0^2
    if val_new < val_old
        idx = 0
        val_old = val_new
    end

    idx
end


################################################################################

"""
    BiasedMinL2Choice

Choice of the stencil in a modified ENO reconstruction minimising the L² norm of
the reconstructed polynomial with a bias towards the central reconstruction.
"""
struct BiasedMinL2Choice end

@inline function (::BiasedMinL2Choice)(u_m1, u_0, u_p1)
    idx = -1
    val_old = (-u_m1 + u_0)^2 # + 4*u_0^2

    val_new = (-u_0 + u_p1)^2 # + 4*u_0^2
    if val_new < val_old
        idx = 0
        val_old = val_new
    end

    idx
end

@inline function (::BiasedMinL2Choice)(u_m2, u_m1, u_0, u_p1, u_p2)
    idx = -2
    val_old = (-2*u_m1 + u_m2 + u_0)^2 + (-12*u_m1 + 3*u_m2 + 9*u_0)^2 # + 144*u_0^2
    val_old *= 2

    val_new = (u_m1 - 2*u_0 + u_p1)^2 + (-3*u_m1 + 3*u_p1)^2 # + 144*u_0^2
    if val_new < val_old
        idx = -1
        val_old = val_new
    end

    val_new = (u_0 - 2*u_p1 + u_p2)^2 + (-9*u_0 + 12*u_p1 - 3*u_p2)^2 # + 144*u_0^2
    val_new *= 2
    if val_new < val_old
        idx = 0
        val_old = val_new
    end

    idx
end

@inline function (::BiasedMinL2Choice)(u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3)
    idx = -3
    val_old = (-3*u_m1 + 3*u_m2 - u_m3 + u_0)^2 + (-50*u_m1 + 40*u_m2 - 10*u_m3 + 20*u_0)^2 + (-177*u_m1 + 87*u_m2 - 19*u_m3 + 109*u_0)^2 # + 14400*u_0^2
    val_old *= 2

    val_new = (3*u_m1 - u_m2 - 3*u_0 + u_p1)^2 + (10*u_m1 - 20*u_0 + 10*u_p1)^2 + (-63*u_m1 + 11*u_m2 + 33*u_0 + 19*u_p1)^2 # + 14400*u_0^2
    if val_new < val_old
        idx = -2
        val_old = val_new
    end

    val_new = (-u_m1 + 3*u_0 - 3*u_p1 + u_p2)^2 + (10*u_m1 - 20*u_0 + 10*u_p1)^2 + (-19*u_m1 - 33*u_0 + 63*u_p1 - 11*u_p2)^2 # + 14400*u_0^2
    if val_new < val_old
        idx = -1
        val_old = val_new
    end

    val_new = (-u_0 + 3*u_p1 - 3*u_p2 + u_p3)^2 + (20*u_0 - 50*u_p1 + 40*u_p2 - 10*u_p3)^2 + (-109*u_0 + 177*u_p1 - 87*u_p2 + 19*u_p3)^2 # + 14400*u_0^2
    val_new *= 2
    if val_new < val_old
        idx = 0
        val_old = val_new
    end

    idx
end

@inline function (::BiasedMinL2Choice)(u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4)
    idx = -4
    val_old = (-4*u_m1 + 6*u_m2 - 4*u_m3 + u_m4 + u_0)^2 + (-126*u_m1 + 168*u_m2 - 98*u_m3 + 21*u_m4 + 35*u_0)^2 + (-1200*u_m1 + 1310*u_m2 - 640*u_m3 + 125*u_m4 + 405*u_0)^2 + (-3234*u_m1 + 2352*u_m2 - 1022*u_m3 + 189*u_m4 + 1715*u_0)^2 # + 2822400*u_0^2
    val_old *= 4

    val_new = (6*u_m1 - 4*u_m2 + u_m3 - 4*u_0 + u_p1)^2 + (84*u_m1 - 42*u_m2 + 7*u_m3 - 70*u_0 + 21*u_p1)^2 + (50*u_m1 + 60*u_m2 - 15*u_m3 - 220*u_0 + 125*u_p1)^2 + (-1344*u_m1 + 462*u_m2 - 77*u_m3 + 770*u_0 + 189*u_p1)^2 # + 2822400*u_0^2
    val_new *= 2
    if val_new < val_old
        idx = -3
        val_old = val_new
    end

    val_new = (-4*u_m1 + u_m2 + 6*u_0 - 4*u_p1 + u_p2)^2 + (14*u_m1 - 7*u_m2 - 14*u_p1 + 7*u_p2)^2 + (200*u_m1 - 15*u_m2 - 370*u_0 + 200*u_p1 - 15*u_p2)^2 + (-574*u_m1 + 77*u_m2 + 574*u_p1 - 77*u_p2)^2 # + 2822400*u_0^2
    if val_new < val_old
        idx = -2
        val_old = val_new
    end

    val_new = (u_m1 - 4*u_0 + 6*u_p1 - 4*u_p2 + u_p3)^2 + (-21*u_m1 + 70*u_0 - 84*u_p1 + 42*u_p2 - 7*u_p3)^2 + (125*u_m1 - 220*u_0 + 50*u_p1 + 60*u_p2 - 15*u_p3)^2 + (-189*u_m1 - 770*u_0 + 1344*u_p1 - 462*u_p2 + 77*u_p3)^2 # + 2822400*u_0^2
    val_new *= 2
    if val_new < val_old
        idx = -1
        val_old = val_new
    end

    val_new = (u_0 - 4*u_p1 + 6*u_p2 - 4*u_p3 + u_p4)^2 + (-35*u_0 + 126*u_p1 - 168*u_p2 + 98*u_p3 - 21*u_p4)^2 + (405*u_0 - 1200*u_p1 + 1310*u_p2 - 640*u_p3 + 125*u_p4)^2 + (-1715*u_0 + 3234*u_p1 - 2352*u_p2 + 1022*u_p3 - 189*u_p4)^2 # + 2822400*u_0^2
    val_new *= 4
    if val_new < val_old
        idx = 0
        val_old = val_new
    end

    idx
end


################################################################################

"""
    LexMinLegendreChoice

Choice of the stencil in a modified ENO reconstruction minimising the Legendre
coefficients of the reconstructed polynomial lexicographically.
"""
struct LexMinLegendreChoice end

@inline function (::LexMinLegendreChoice)(u_m1, u_0, u_p1)
    idx = -1
    val_old = ( (-u_m1 + u_0)^2,)
    # first coefficient 4*u_0^2

    val_new = ( (-u_0 + u_p1)^2,)
    # first coefficient 4*u_0^2
    if val_new < val_old
        idx = 0
        val_old = val_new
    end

    idx
end

@inline function (::LexMinLegendreChoice)(u_m2, u_m1, u_0, u_p1, u_p2)
    idx = -2
    val_old = ( (-12*u_m1 + 3*u_m2 + 9*u_0)^2, (-2*u_m1 + u_m2 + u_0)^2,)
    # first coefficient 144*u_0^2

    val_new = ( (-3*u_m1 + 3*u_p1)^2, (u_m1 - 2*u_0 + u_p1)^2,)
    # first coefficient 144*u_0^2
    if val_new < val_old
        idx = -1
        val_old = val_new
    end

    val_new = ( (-9*u_0 + 12*u_p1 - 3*u_p2)^2, (u_0 - 2*u_p1 + u_p2)^2,)
    # first coefficient 144*u_0^2
    if val_new < val_old
        idx = 0
        val_old = val_new
    end

    idx
end

@inline function (::LexMinLegendreChoice)(u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3)
    idx = -3
    val_old = ( (-177*u_m1 + 87*u_m2 - 19*u_m3 + 109*u_0)^2, (-50*u_m1 + 40*u_m2 - 10*u_m3 + 20*u_0)^2, (-3*u_m1 + 3*u_m2 - u_m3 + u_0)^2,)
    # first coefficient 14400*u_0^2

    val_new = ( (-63*u_m1 + 11*u_m2 + 33*u_0 + 19*u_p1)^2, (10*u_m1 - 20*u_0 + 10*u_p1)^2, (3*u_m1 - u_m2 - 3*u_0 + u_p1)^2,)
    # first coefficient 14400*u_0^2
    if val_new < val_old
        idx = -2
        val_old = val_new
    end

    val_new = ( (-19*u_m1 - 33*u_0 + 63*u_p1 - 11*u_p2)^2, (10*u_m1 - 20*u_0 + 10*u_p1)^2, (-u_m1 + 3*u_0 - 3*u_p1 + u_p2)^2,)
    # first coefficient 14400*u_0^2
    if val_new < val_old
        idx = -1
        val_old = val_new
    end

    val_new = ( (-109*u_0 + 177*u_p1 - 87*u_p2 + 19*u_p3)^2, (20*u_0 - 50*u_p1 + 40*u_p2 - 10*u_p3)^2, (-u_0 + 3*u_p1 - 3*u_p2 + u_p3)^2,)
    # first coefficient 14400*u_0^2
    if val_new < val_old
        idx = 0
        val_old = val_new
    end

    idx
end

@inline function (::LexMinLegendreChoice)(u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4)
    idx = -4
    val_old = ( (-3234*u_m1 + 2352*u_m2 - 1022*u_m3 + 189*u_m4 + 1715*u_0)^2, (-1200*u_m1 + 1310*u_m2 - 640*u_m3 + 125*u_m4 + 405*u_0)^2, (-126*u_m1 + 168*u_m2 - 98*u_m3 + 21*u_m4 + 35*u_0)^2, (-4*u_m1 + 6*u_m2 - 4*u_m3 + u_m4 + u_0)^2,)
    # first coefficient 2822400*u_0^2

    val_new = ( (-1344*u_m1 + 462*u_m2 - 77*u_m3 + 770*u_0 + 189*u_p1)^2, (50*u_m1 + 60*u_m2 - 15*u_m3 - 220*u_0 + 125*u_p1)^2, (84*u_m1 - 42*u_m2 + 7*u_m3 - 70*u_0 + 21*u_p1)^2, (6*u_m1 - 4*u_m2 + u_m3 - 4*u_0 + u_p1)^2,)
    # first coefficient 2822400*u_0^2
    if val_new < val_old
        idx = -3
        val_old = val_new
    end

    val_new = ( (-574*u_m1 + 77*u_m2 + 574*u_p1 - 77*u_p2)^2, (200*u_m1 - 15*u_m2 - 370*u_0 + 200*u_p1 - 15*u_p2)^2, (14*u_m1 - 7*u_m2 - 14*u_p1 + 7*u_p2)^2, (-4*u_m1 + u_m2 + 6*u_0 - 4*u_p1 + u_p2)^2,)
    # first coefficient 2822400*u_0^2
    if val_new < val_old
        idx = -2
        val_old = val_new
    end

    val_new = ( (-189*u_m1 - 770*u_0 + 1344*u_p1 - 462*u_p2 + 77*u_p3)^2, (125*u_m1 - 220*u_0 + 50*u_p1 + 60*u_p2 - 15*u_p3)^2, (-21*u_m1 + 70*u_0 - 84*u_p1 + 42*u_p2 - 7*u_p3)^2, (u_m1 - 4*u_0 + 6*u_p1 - 4*u_p2 + u_p3)^2,)
    # first coefficient 2822400*u_0^2
    if val_new < val_old
        idx = -1
        val_old = val_new
    end

    val_new = ( (-1715*u_0 + 3234*u_p1 - 2352*u_p2 + 1022*u_p3 - 189*u_p4)^2, (405*u_0 - 1200*u_p1 + 1310*u_p2 - 640*u_p3 + 125*u_p4)^2, (-35*u_0 + 126*u_p1 - 168*u_p2 + 98*u_p3 - 21*u_p4)^2, (u_0 - 4*u_p1 + 6*u_p2 - 4*u_p3 + u_p4)^2,)
    # first coefficient 2822400*u_0^2
    if val_new < val_old
        idx = 0
        val_old = val_new
    end

    idx
end
