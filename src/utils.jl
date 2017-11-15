
"""
    logmean(a::Real, b::Real)

The logarithmic mean value (a - b) / ( log(a) - log(b) ).
"""
function logmean(_a::T, _b::T) where T<:Real
    a, b = minmax(_a,_b)
    if a < 0
        throw(DomainError("Cannot compute logarithmic mean of negative valuess."))
    end

    if a == b
        a
    else
        (a-b) / (log(a/b))
    end
end
logmean(a::Real, b::Real) = logmean(promote(a,b)...)

function logmean(_a::Float64, _b::Float64)
    a, b = minmax(_a,_b)
    if a < 0
        throw(DomainError("Cannot compute logarithmic mean of negative valuess."))
    end

    if a == b
        return a
    end

    # see Ismail, Roe (2009): Affordable, entropy consistent...
    ζ = a / b
    f = (ζ-1) / (ζ+1)
    u = f*f
    if u < 1.0e-2
        F = @evalpoly(u, 1, 1/3, 1/5, 1/7)
    else
        F = log(ζ)/(2f)
    end
    (a+b)/(2F)
end

"""
    integrate(func, u::AbstractArray{U,1}, meshx::AbstractMesh1D)

Map the function `func` to the coefficients `u` and integrate with respect to
the mesh `meshx`.
"""
function integrate(func, u::AbstractArray{U,1}, meshx::AbstractMesh1D) where U
    @boundscheck begin
        @assert numcells(meshx) == size(u,1)
    end

    res = zero(func(first(u)))
    @inbounds for ix in cell_indices(meshx)
        res += func(u[ix]) * volume(ix, meshx)
    end
    res
end

"""
    integrate(u::AbstractArray{U,1}, meshx::AbstractMesh1D)

Integrate the coefficients `u` with respect to the mesh `meshx`.
"""
function integrate(u::AbstractArray{U,1}, meshx::AbstractMesh1D) where U
    integrate(identity, u, meshx)
end


"""
    integrate(func, u::AbstractArray{U,2}, meshx::AbstractMesh1D, basis::NodalBasis)

Map the function `func` to the coefficients `u` and integrate with respect to
the mesh `meshx` and the polynomial basis `basis`.
"""
function integrate(func, u::AbstractArray{U,2}, meshx::AbstractMesh1D, basis::NodalBasis) where U
    @unpack weights = basis
    Pp1 = length(weights)
    @boundscheck begin
        @assert Pp1 == size(u,1)
        @assert numcells(meshx) == size(u,2)
    end

    res = zero(func(first(u)))
    @inbounds for ix in cell_indices(meshx)
        jac = volume(ix, meshx) / 2
        for nx in 1:Pp1
            res += func(u[nx,ix]) * weights[nx] * jac
        end
    end
    res
end

"""
    integrate(u::AbstractArray{U,2}, meshx::AbstractMesh1D, basis::NodalBasis)

Integrate the coefficients `u` with respect to the mesh `meshx` and the
polynomial basis `basis`.
"""
function integrate(u::AbstractArray{U,2}, meshx::AbstractMesh1D, basis::NodalBasis) where U
    integrate(identity, u, meshx, basis)
end


"""
    integrate(func, u::AbstractArray{U,4}, meshx::AbstractMesh1D, meshy::AbstractMesh1D, basis::NodalBasis)

Map the function `func` to the coefficients `u` and integrate with respect to
the meshes `meshx, meshy` and the polynomial basis `basis`.
"""
function integrate(func, u::AbstractArray{U,4}, meshx::AbstractMesh1D, meshy::AbstractMesh1D, basis::NodalBasis) where U
    @unpack weights = basis
    Pp1 = length(weights)
    @boundscheck begin
        @assert Pp1 == size(u,1) == size(u,2)
        @assert numcells(meshx) == size(u,3)
        @assert numcells(meshy) == size(u,4)
    end

    res = zero(func(first(u)))
    @inbounds for iy in cell_indices(meshy), ix in cell_indices(meshx)
        jac = volume(ix, meshx) * volume(iy, meshy) / 4
        for ny in 1:Pp1, nx in 1:Pp1
            res += func(u[nx,ny,ix,iy]) * weights[nx] * weights[ny] * jac
        end
    end
    res
end

"""
    integrate(u::AbstractArray{U,4}, meshx::AbstractMesh1D, meshy::AbstractMesh1D, basis::NodalBasis)

Integrate the coefficients `u` with respect to the meshes `meshx, meshy` and the
polynomial basis `basis`.
"""
function integrate(u::AbstractArray{U,4}, meshx::AbstractMesh1D, meshy::AbstractMesh1D, basis::NodalBasis) where U
    integrate(identity, u, meshx, meshy, basis)
end


"""
    integrate(func, u::AbstractArray{U,6}, meshx, meshy, meshz, basis::NodalBasis)

Map the function `func` to the coefficients `u` and integrate with respect to
the meshes `meshx, meshy, meshz` and the polynomial basis `basis`.
"""
function integrate(func, u::AbstractArray{U,6}, meshx::AbstractMesh1D, meshy::AbstractMesh1D, meshz::AbstractMesh1D, basis::NodalBasis) where U
    @unpack weights = basis
    Pp1 = length(weights)
    @boundscheck begin
        @assert Pp1 == size(u,1) == size(u,2) == size(u,3)
        @assert numcells(meshx) == size(u,4)
        @assert numcells(meshy) == size(u,5)
        @assert numcells(meshz) == size(u,6)
    end

    res = zero(func(first(u)))
    @inbounds for iz in cell_indices(meshz), iy in cell_indices(meshy), ix in cell_indices(meshx)
        jac = volume(ix, meshx) * volume(iy, meshy) * volume(iz, meshz)  / 8
        for nz in 1:Pp1, ny in 1:Pp1, nx in 1:Pp1
            res += func(u[nx,ny,nz,ix,iy,iz]) * weights[nx] * weights[ny] * weights[nz] * jac
        end
    end
    res
end

"""
    integrate(u::AbstractArray{U,6}, meshx, meshy, meshz, basis::NodalBasis)

Integrate the coefficients `u` with respect to the meshes `meshx, meshy, meshz`
and the polynomial basis `basis`.
"""
function integrate(u::AbstractArray{U,6}, meshx::AbstractMesh1D, meshy::AbstractMesh1D, meshz::AbstractMesh1D, basis::NodalBasis) where U
    integrate(identity, u, meshx, meshy, meshz, basis)
end
