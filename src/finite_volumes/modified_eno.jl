
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
