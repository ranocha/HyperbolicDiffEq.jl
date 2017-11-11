
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
        edge_u[1,cell] = u_0/2 + u_m1/2
        edge_u[2,cell] = 3*u_0/2 - u_m1/2
    else # idx == 0
        edge_u[1,cell] = 3*u_0/2 - u_p1/2
        edge_u[2,cell] = u_0/2 + u_p1/2
    end
    nothing
end

function (eno::ModifiedENO{3})(edge_u, cell, u_m2, u_m1, u_0, u_p1, u_p2, balance_law, meshx)
    idx = eno.choose_stencil(u_m2, u_m1, u_0, u_p1, u_p2)

    @inbounds if idx == -2
        edge_u[1,cell] = u_0/3 + 5*u_m1/6 - u_m2/6
        edge_u[2,cell] = 11*u_0/6 - 7*u_m1/6 + u_m2/3
    elseif idx == -1
        edge_u[1,cell] = 5*u_0/6 - u_p1/6 + u_m1/3
        edge_u[2,cell] = 5*u_0/6 + u_p1/3 - u_m1/6
    else # idx == 0
        edge_u[1,cell] = 11*u_0/6 - 7*u_p1/6 + u_p2/3
        edge_u[2,cell] = u_0/3 + 5*u_p1/6 - u_p2/6
    end
    nothing
end

function (eno::ModifiedENO{4})(edge_u, cell, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, balance_law, meshx)
    idx = eno.choose_stencil(u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3)

    @inbounds if idx == -3
        edge_u[1,cell] = u_0/4 + 13*u_m1/12 - 5*u_m2/12 + u_m3/12
        edge_u[2,cell] = 25*u_0/12 - 23*u_m1/12 + 13*u_m2/12 - u_m3/4
    elseif idx == -2
        edge_u[1,cell] = 7*u_0/12 - u_p1/12 + 7*u_m1/12 - u_m2/12
        edge_u[2,cell] = 13*u_0/12 + u_p1/4 - 5*u_m1/12 + u_m2/12
    elseif idx == -1
        edge_u[1,cell] = 13*u_0/12 - 5*u_p1/12 + u_p2/12 + u_m1/4
        edge_u[2,cell] = 7*u_0/12 + 7*u_p1/12 - u_p2/12 - u_m1/12
    else # idx == 0
        edge_u[1,cell] = 25*u_0/12 - 23*u_p1/12 + 13*u_p2/12 - u_p3/4
        edge_u[2,cell] = u_0/4 + 13*u_p1/12 - 5*u_p2/12 + u_p3/12
    end
    nothing
end

function (eno::ModifiedENO{5})(edge_u, cell, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, balance_law, meshx)
    idx = eno.choose_stencil(u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4)

    @inbounds if idx == -4
        edge_u[1,cell] = u_0/5 + 77*u_m1/60 - 43*u_m2/60 + 17*u_m3/60 - u_m4/20
        edge_u[2,cell] = 137*u_0/60 - 163*u_m1/60 + 137*u_m2/60 - 21*u_m3/20 + u_m4/5
    elseif idx == -3
        edge_u[1,cell] = 9*u_0/20 - u_p1/20 + 47*u_m1/60 - 13*u_m2/60 + u_m3/30
        edge_u[2,cell] = 77*u_0/60 + u_p1/5 - 43*u_m1/60 + 17*u_m2/60 - u_m3/20
    elseif idx == -2
        edge_u[1,cell] = 47*u_0/60 - 13*u_p1/60 + u_p2/30 + 9*u_m1/20 - u_m2/20
        edge_u[2,cell] = 47*u_0/60 + 9*u_p1/20 - u_p2/20 - 13*u_m1/60 + u_m2/30
    elseif idx == -1
        edge_u[1,cell] = 77*u_0/60 - 43*u_p1/60 + 17*u_p2/60 - u_p3/20 + u_m1/5
        edge_u[2,cell] = 9*u_0/20 + 47*u_p1/60 - 13*u_p2/60 + u_p3/30 - u_m1/20
    else # idx == 0
        edge_u[1,cell] = 137*u_0/60 - 163*u_p1/60 + 137*u_p2/60 - 21*u_p3/20 + u_p4/5
        edge_u[2,cell] = u_0/5 + 77*u_p1/60 - 43*u_p2/60 + 17*u_p3/60 - u_p4/20
    end
    nothing
end

function (eno::ModifiedENO{6})(edge_u, cell, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, balance_law, meshx)
    idx = eno.choose_stencil(u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5)

    @inbounds if idx == -5
        edge_u[1,cell] = u_0/6 + 29*u_m1/20 - 21*u_m2/20 + 37*u_m3/60 - 13*u_m4/60 + u_m5/30
        edge_u[2,cell] = 49*u_0/20 - 71*u_m1/20 + 79*u_m2/20 - 163*u_m3/60 + 31*u_m4/30 - u_m5/6
    elseif idx == -4
        edge_u[1,cell] = 11*u_0/30 - u_p1/30 + 19*u_m1/20 - 23*u_m2/60 + 7*u_m3/60 - u_m4/60
        edge_u[2,cell] = 29*u_0/20 + u_p1/6 - 21*u_m1/20 + 37*u_m2/60 - 13*u_m3/60 + u_m4/30
    elseif idx == -3
        edge_u[1,cell] = 37*u_0/60 - 2*u_p1/15 + u_p2/60 + 37*u_m1/60 - 2*u_m2/15 + u_m3/60
        edge_u[2,cell] = 19*u_0/20 + 11*u_p1/30 - u_p2/30 - 23*u_m1/60 + 7*u_m2/60 - u_m3/60
    elseif idx == -2
        edge_u[1,cell] = 19*u_0/20 - 23*u_p1/60 + 7*u_p2/60 - u_p3/60 + 11*u_m1/30 - u_m2/30
        edge_u[2,cell] = 37*u_0/60 + 37*u_p1/60 - 2*u_p2/15 + u_p3/60 - 2*u_m1/15 + u_m2/60
    elseif idx == -1
        edge_u[1,cell] = 29*u_0/20 - 21*u_p1/20 + 37*u_p2/60 - 13*u_p3/60 + u_p4/30 + u_m1/6
        edge_u[2,cell] = 11*u_0/30 + 19*u_p1/20 - 23*u_p2/60 + 7*u_p3/60 - u_p4/60 - u_m1/30
    else # idx == 0
        edge_u[1,cell] = 49*u_0/20 - 71*u_p1/20 + 79*u_p2/20 - 163*u_p3/60 + 31*u_p4/30 - u_p5/6
        edge_u[2,cell] = u_0/6 + 29*u_p1/20 - 21*u_p2/20 + 37*u_p3/60 - 13*u_p4/60 + u_p5/30
    end

    nothing
end

function (eno::ModifiedENO{7})(edge_u, cell, u_m6, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, u_p6, balance_law, meshx)
    idx = eno.choose_stencil(u_m6, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, u_p6)

    @inbounds if idx == -6
            edge_u[1,cell] = u_0/7 + 223*u_m1/140 - 197*u_m2/140 + 153*u_m3/140 - 241*u_m4/420 + 37*u_m5/210 - u_m6/42
            edge_u[2,cell] = 363*u_0/140 - 617*u_m1/140 + 853*u_m2/140 - 2341*u_m3/420 + 667*u_m4/210 - 43*u_m5/42 + u_m6/7
        elseif idx == -5
            edge_u[1,cell] = 13*u_0/42 - u_p1/42 + 153*u_m1/140 - 241*u_m2/420 + 109*u_m3/420 - 31*u_m4/420 + u_m5/105
            edge_u[2,cell] = 223*u_0/140 + u_p1/7 - 197*u_m1/140 + 153*u_m2/140 - 241*u_m3/420 + 37*u_m4/210 - u_m5/42
        elseif idx == -4
            edge_u[1,cell] = 107*u_0/210 - 19*u_p1/210 + u_p2/105 + 319*u_m1/420 - 101*u_m2/420 + 5*u_m3/84 - u_m4/140
            edge_u[2,cell] = 153*u_0/140 + 13*u_p1/42 - u_p2/42 - 241*u_m1/420 + 109*u_m2/420 - 31*u_m3/420 + u_m4/105
        elseif idx == -3
            edge_u[1,cell] = 319*u_0/420 - 101*u_p1/420 + 5*u_p2/84 - u_p3/140 + 107*u_m1/210 - 19*u_m2/210 + u_m3/105
            edge_u[2,cell] = 319*u_0/420 + 107*u_p1/210 - 19*u_p2/210 + u_p3/105 - 101*u_m1/420 + 5*u_m2/84 - u_m3/140
        elseif idx == -2
            edge_u[1,cell] = 153*u_0/140 - 241*u_p1/420 + 109*u_p2/420 - 31*u_p3/420 + u_p4/105 + 13*u_m1/42 - u_m2/42
            edge_u[2,cell] = 107*u_0/210 + 319*u_p1/420 - 101*u_p2/420 + 5*u_p3/84 - u_p4/140 - 19*u_m1/210 + u_m2/105
        elseif idx == -1
            edge_u[1,cell] = 223*u_0/140 - 197*u_p1/140 + 153*u_p2/140 - 241*u_p3/420 + 37*u_p4/210 - u_p5/42 + u_m1/7
            edge_u[2,cell] = 13*u_0/42 + 153*u_p1/140 - 241*u_p2/420 + 109*u_p3/420 - 31*u_p4/420 + u_p5/105 - u_m1/42
        else # idx == 0
            edge_u[1,cell] = 363*u_0/140 - 617*u_p1/140 + 853*u_p2/140 - 2341*u_p3/420 + 667*u_p4/210 - 43*u_p5/42 + u_p6/7
            edge_u[2,cell] = u_0/7 + 223*u_p1/140 - 197*u_p2/140 + 153*u_p3/140 - 241*u_p4/420 + 37*u_p5/210 - u_p6/42
        end

    nothing
end

function (eno::ModifiedENO{8})(edge_u, cell, u_m7, u_m6, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, u_p6, u_p7, balance_law, meshx)
    idx = eno.choose_stencil(u_m7, u_m6, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, u_p6, u_p7)

    @inbounds if idx == -7
        edge_u[1,cell] = u_0/8 + 481*u_m1/280 - 499*u_m2/280 + 481*u_m3/280 - 1007*u_m4/840 + 463*u_m5/840 - 25*u_m6/168 + u_m7/56
        edge_u[2,cell] = 761*u_0/280 - 1479*u_m1/280 + 2441*u_m2/280 - 8357*u_m3/840 + 6343*u_m4/840 - 613*u_m5/168 + 57*u_m6/56 - u_m7/8
    elseif idx == -6
        edge_u[1,cell] = 15*u_0/56 - u_p1/56 + 341*u_m1/280 - 219*u_m2/280 + 131*u_m3/280 - 167*u_m4/840 + 43*u_m5/840 - u_m6/168
        edge_u[2,cell] = 481*u_0/280 + u_p1/8 - 499*u_m1/280 + 481*u_m2/280 - 1007*u_m3/840 + 463*u_m4/840 - 25*u_m5/168 + u_m6/56
    elseif idx == -5
        edge_u[1,cell] = 73*u_0/168 - 11*u_p1/168 + u_p2/168 + 743*u_m1/840 - 307*u_m2/840 + 113*u_m3/840 - 9*u_m4/280 + u_m5/280
        edge_u[2,cell] = 341*u_0/280 + 15*u_p1/56 - u_p2/56 - 219*u_m1/280 + 131*u_m2/280 - 167*u_m3/840 + 43*u_m4/840 - u_m5/168
    elseif idx == -4
        edge_u[1,cell] = 533*u_0/840 - 139*u_p1/840 + 29*u_p2/840 - u_p3/280 + 533*u_m1/840 - 139*u_m2/840 + 29*u_m3/840 - u_m4/280
        edge_u[2,cell] = 743*u_0/840 + 73*u_p1/168 - 11*u_p2/168 + u_p3/168 - 307*u_m1/840 + 113*u_m2/840 - 9*u_m3/280 + u_m4/280
    elseif idx == -3
        edge_u[1,cell] = 743*u_0/840 - 307*u_p1/840 + 113*u_p2/840 - 9*u_p3/280 + u_p4/280 + 73*u_m1/168 - 11*u_m2/168 + u_m3/168
        edge_u[2,cell] = 533*u_0/840 + 533*u_p1/840 - 139*u_p2/840 + 29*u_p3/840 - u_p4/280 - 139*u_m1/840 + 29*u_m2/840 - u_m3/280
    elseif idx == -2
        edge_u[1,cell] = 341*u_0/280 - 219*u_p1/280 + 131*u_p2/280 - 167*u_p3/840 + 43*u_p4/840 - u_p5/168 + 15*u_m1/56 - u_m2/56
        edge_u[2,cell] = 73*u_0/168 + 743*u_p1/840 - 307*u_p2/840 + 113*u_p3/840 - 9*u_p4/280 + u_p5/280 - 11*u_m1/168 + u_m2/168
    elseif idx == -1
        edge_u[1,cell] = 481*u_0/280 - 499*u_p1/280 + 481*u_p2/280 - 1007*u_p3/840 + 463*u_p4/840 - 25*u_p5/168 + u_p6/56 + u_m1/8
        edge_u[2,cell] = 15*u_0/56 + 341*u_p1/280 - 219*u_p2/280 + 131*u_p3/280 - 167*u_p4/840 + 43*u_p5/840 - u_p6/168 - u_m1/56
    else # idx == 0
        edge_u[1,cell] = 761*u_0/280 - 1479*u_p1/280 + 2441*u_p2/280 - 8357*u_p3/840 + 6343*u_p4/840 - 613*u_p5/168 + 57*u_p6/56 - u_p7/8
        edge_u[2,cell] = u_0/8 + 481*u_p1/280 - 499*u_p2/280 + 481*u_p3/280 - 1007*u_p4/840 + 463*u_p5/840 - 25*u_p6/168 + u_p7/56
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
