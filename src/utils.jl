
"""
    logmean(a, b)

The logarithmic mean value (a - b) / ( log(a) - log(b) ).
"""
@inline function logmean(a, b)
    if a == b
        a
    else
        (a-b) / (log(a)-log(b))
    end
end

Base.@pure function logmean(a::Float64, b::Float64)
  # see Ismail, Roe (2009): Affordable, entropy consistent...
  const ɛ = 0.01
  ζ = a/b
  f = (ζ-1)/(ζ+1)
  u = f*f
  if u < ɛ
    F = 1 + u/3 + u*u/5 + u*u*u/7
  else
    F = log(ζ)/(2f)
  end
  (a+b)/(2F)
end


"""
    integrate(func, u::AbstractArray{U,4}, meshx, meshy, basis::NodalBasis)

Map the function `func` to the coefficients `u` and integrate with respect to
the meshes `meshx, meshy` and the polynomial basis `basis`.
"""
function integrate(func, u::AbstractArray{U,4}, meshx, meshy, basis::NodalBasis) where U
    @unpack weights = basis
    Pp1 = length(weights)
    @boundscheck begin
        @assert Pp1 == size(u,1) == size(u,2)
        @assert numcells(meshx) == size(u,3)
        @assert numcells(meshy) == size(u,4)
    end

    res = zero(func(first(u)))
    @inbounds for iy in cell_indices(meshy), ix in cell_indices(meshx)
        jac = volume(ix, meshx) * volume(iy, meshy) / 2
        for ny in 1:Pp1, nx in 1:Pp1
            res += func(u[nx,ny,ix,iy]) * weights[nx] * weights[ny] * jac
        end
    end
    res
end

"""
    integrate(u::AbstractArray{U,4}, meshx, meshy, basis::NodalBasis)

Integrate the coefficients `u`with respect to the meshes `meshx, meshy` and the
polynomial basis `basis`.
"""
function integrate(u::AbstractArray{U,4}, meshx, meshy, basis) where U
    integrate(identity, u, meshx, meshy, basis)
end
