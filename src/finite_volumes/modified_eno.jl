
"""
    ModifiedENO{K, ChooseStencil}

A reconstruction using the same polynomials as the essentially non-oscillatory
(ENO) reconstruction on `K` cells. The stencil is chosen by `choose_stencil`.
"""
struct ModifiedENO{K, ChooseStencil} <: AbstractReconstruction
    choose_stencil::ChooseStencil
end

function ModifiedENO{K}(::Val{K}=Val{2}(), choose_stencil=MinAbsDiffLex())
    ModifiedENO{K,typeof(choose_stencil)}(choose_stencil)
end

function ModifiedENO(K::Integer, choose_stencil=MinAbsDiffLex())
    ModifiedENO{K,typeof(choose_stencil)}(choose_stencil)
end


# ENO{K} uses K different polynomials of degree K-1
# -> the total stencil width is 2K-1
@inline stencil_width(::ModifiedENO{1}) = Val{1}()
@inline stencil_width(::ModifiedENO{2}) = Val{3}()
@inline stencil_width(::ModifiedENO{3}) = Val{5}()
@inline stencil_width(::ModifiedENO{4}) = Val{7}()
@inline stencil_width(::ModifiedENO{5}) = Val{9}()


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

#=
#TODO
function (eno::ModifiedENO{6})(edge_u, cell, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, balance_law, meshx)
    idx = eno.choose_stencil(u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5)


    nothing
end
=#

################################################################################

"""
Classical choice of the stencil in an ENO reconstruction, minimising the absolute
values of the divided differences lexicographically.
"""
struct MinAbsDiffLex end


@inline function (choose_stencil::MinAbsDiffLex)(u_m1, u_0, u_p1)
  if abs(u_m1 - u_0) < abs(u_0 - u_p1)
    idx = -1
  else
    idx = 0
  end

  idx
end


@inline function (choose_stencil::MinAbsDiffLex)(u_m2, u_m1, u_0, u_p1, u_p2)
  if abs(u_m1 - u_0) < abs(u_0 - u_p1)
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


@inline function (choose_stencil::MinAbsDiffLex)(u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3)
  if abs(u_m1 - u_0) < abs(u_0 - u_p1)
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


@inline function (choose_stencil::MinAbsDiffLex)(u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4)
  if abs(u_m1 - u_0) < abs(u_0 - u_p1)
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


@inline function (choose_stencil::MinAbsDiffLex)(u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5)
  if abs(u_m1 - u_0) < abs(u_0 - u_p1)
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
