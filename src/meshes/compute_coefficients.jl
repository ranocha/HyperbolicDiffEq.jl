
"""
    compute_coefficients(u, mesh::AbstractMesh1D)

Computse the coefficients of the function `u` in a first order finite volume
method on the `mesh`.
"""
function compute_coefficients(u, mesh::AbstractMesh1D, order::Integer=1)
    xmin, xmax = bounds(mesh)
    uval = zeros(typeof(u((xmin+xmax)/2)), numcells(mesh))
    compute_coefficients!(uval, u, mesh, order)
    uval
end

"""
    compute_coefficients!(uval, u, mesh::AbstractMesh1D)

Computse the coefficients of the function `u` in a first order finite volume
method on the `mesh` and stores them in `uval`.
"""
function compute_coefficients!(uval, u, mesh::AbstractMesh1D, order::Integer=1)
    if eltype(uval) <: AbstractArray
        basis = GaussLegendre(order-1, eltype(eltype(uval)))
    else
        basis = GaussLegendre(order-1, eltype(uval))
    end
    x = similar(basis.nodes)

    for cell in cell_indices(mesh)
        xmin, xmax = bounds(cell, mesh)
        map_from_canonical!(x, basis.nodes, xmin, xmax, basis)
        uval[cell] = integrate(u, x, basis.weights) / 2
    end
    nothing
end


"""
    evaluate_coefficients(u, mesh::AbstractMesh1D)

Evaluates the coefficients `u` in a first order finite volume method on the `mesh`.
Returns `xplot, uplot` as vectors that can be used for plots.
"""
function evaluate_coefficients(u, mesh::AbstractMesh1D)
    xplot = zeros(2*numcells(mesh))
    uplot = zeros(eltype(u), 2*numcells(mesh))

    evaluate_coefficients!(xplot, uplot, u, mesh)
end

"""
    evaluate_coefficients!(xplot, uplot, u, mesh::AbstractMesh1D)

Evaluates the coefficients `u` in a first order finite volume method on the `mesh`.
Returns `xplot, uplot` as vectors that can be used for plots.
"""
function evaluate_coefficients!(xplot, uplot, u, mesh::AbstractMesh1D)
  @assert length(xplot) == 2*numcells(mesh)
  @assert length(uplot) == 2*numcells(mesh)

  idx = 1
  @inbounds for cell in cell_indices(mesh)
    xmin, xmax = bounds(cell, mesh)

    xplot[idx  ] = xmin
    xplot[idx+1] = xmax
    uplot[idx] = uplot[idx+1] = u[cell]

    idx = idx+2
  end

  xplot, uplot
end


"""
    evaluate_coefficients(u, balance_law, mesh::AbstractMesh1D,
                          reconstruction::AbstractReconstruction,
                          npoints=2*order(reconstruction))

Evaluates the coefficients `u` in a finite volume method on `mesh` using
`reconstruction`.
Returns `xplot, uplot` as vectors that can be used for plots.
"""
function evaluate_coefficients(u, balance_law, meshx::AbstractMesh1D,
                               reconstruction::AbstractReconstruction,
                               npoints=2*order(reconstruction))
    xplot = zeros(npoints*numcells(meshx))
    uplot = zeros(eltype(u), npoints*numcells(meshx))

    evaluate_coefficients!(xplot, uplot, u, balance_law, meshx, reconstruction)
end

"""
    evaluate_coefficients!(xplot, uplot, u, balance_law, mesh::AbstractMesh1D,
                           reconstruction::AbstractReconstruction)

Evaluates the coefficients `u` in a finite volume method on `mesh` using
`reconstruction`.
Returns `xplot, uplot` as vectors that can be used for plots.
"""
function evaluate_coefficients!(xplot, uplot, u, balance_law, meshx::AbstractMesh1D,
                                reconstruction::AbstractReconstruction)
    npoints = length(xplot) ÷ numcells(meshx)
    @assert length(uplot) == npoints*numcells(meshx)

    ξ = linspace(0+eps(), 1-eps(), npoints) |> collect
    x = zeros(ξ)

    uval = zeros(eltype(u), npoints)
    @inbounds for cell in cell_indices(meshx)
        xmin, xmax = bounds(cell, meshx)
        interpolate!(uval, ξ, cell, u, balance_law, meshx, reconstruction, stencil_width_val(reconstruction))
        for nx in 1:npoints
            xplot[(cell-1)*npoints+(nx-1)+1] = xmin + ξ[nx]*(xmax-xmin)
            uplot[(cell-1)*npoints+(nx-1)+1] = uval[nx]
        end
    end

    xplot, uplot
end


################################################################################

function compute_coefficients(u, mesh::AbstractMesh1D, basis::NodalBasis, parallel=Val{:serial}())
    xmin, xmax = bounds(mesh)
    uval = Array{typeof(u((xmin+xmax)/2))}(length(basis.nodes), numcells(mesh))
    compute_coefficients!(uval, u, mesh, basis)
    uval
end

function compute_coefficients!(uval, u, mesh::AbstractMesh1D, basis::NodalBasis)
    @unpack nodes = basis
    Pp1 = length(nodes)

    for cell in cell_indices(mesh)
        xmin, xmax = bounds(cell, mesh)
        for n in 1:Pp1
            x = map_from_canonical(nodes[n], xmin, xmax, basis)
            uval[n,cell] = u(x)
        end
    end

    nothing
end

function evaluate_coefficients(u, meshx::AbstractMesh1D, basis::NodalBasis,
                                npoints=2*length(basis.nodes))
    xplot = zeros(npoints*numcells(meshx))
    uplot = zeros(eltype(u), npoints*numcells(meshx))

    evaluate_coefficients!(xplot, uplot, u, meshx, basis)
end

function evaluate_coefficients!(xplot, uplot, u, meshx::AbstractMesh1D, basis::NodalBasis)
    npoints = length(xplot) ÷ numcells(meshx)

    @assert length(xplot) == npoints*numcells(meshx)
    @assert length(xplot) == length(uplot)

    ξ = linspace(-1+eps(), 1-eps(), npoints) |> collect
    x = zeros(ξ)

    Nx  = numcells(meshx)
    Pp1 = length(basis.nodes)
    uval = zeros(eltype(u), npoints)
    for ix in 1:Nx
        xmin, xmax = bounds(ix, meshx)
        map_from_canonical!(x, ξ, xmin, xmax, basis)
        interpolate!(uval, ξ, u[:,ix], basis)

        for jx in 1:npoints
            xplot[(ix-1)*npoints+(jx-1)+1] = x[jx]
            uplot[(ix-1)*npoints+(jx-1)+1] = uval[jx]
        end
    end

    xplot, uplot
end


function compute_coefficients(u, meshx::AbstractMesh1D, meshy::AbstractMesh1D, basis::NodalBasis)
    xmin, xmax = bounds(meshx)
    ymin, ymax = bounds(meshy)
    Pp1 = length(basis.nodes)
    uval = zeros(typeof(u((xmin+xmax)/2,(ymin+ymax)/2)), Pp1, Pp1, numcells(meshx), numcells(meshy))
    compute_coefficients!(uval, u, meshx, meshy, basis)
    uval
end

function compute_coefficients!(uval, u, meshx::AbstractMesh1D, meshy::AbstractMesh1D, basis::NodalBasis)
    @unpack nodes = basis
    Pp1 = length(nodes)

    for iy in cell_indices(meshy), ix in cell_indices(meshx)
        xmin, xmax = bounds(ix, meshx)
        ymin, ymax = bounds(iy, meshy)
        for ny in 1:Pp1, nx in 1:Pp1
            x = map_from_canonical(nodes[nx], xmin, xmax, basis)
            y = map_from_canonical(nodes[ny], ymin, ymax, basis)
            uval[nx,ny,ix,iy] = u(x,y)
        end
    end

    nothing
end

function evaluate_coefficients(u, meshx::AbstractMesh1D, meshy::AbstractMesh1D,
                                basis::NodalBasis, npoints=2*length(basis.nodes))
    xplot = zeros(npoints*numcells(meshx))
    yplot = zeros(npoints*numcells(meshy))
    uplot = zeros(eltype(u), npoints*numcells(meshy), npoints*numcells(meshx))

    evaluate_coefficients!(xplot, yplot, uplot, u, meshx, meshy, basis)
end

function evaluate_coefficients!(xplot, yplot, uplot, u, meshx::AbstractMesh1D,
                                meshy::AbstractMesh1D, basis::NodalBasis)
    npoints = length(xplot) ÷ numcells(meshx)

    @assert length(xplot) == npoints*numcells(meshx)
    @assert length(yplot) == npoints*numcells(meshy)
    @assert length(xplot) == size(uplot,2)
    @assert length(yplot) == size(uplot,1)

    ξ = linspace(-1+eps(),1-eps(),npoints) |> collect
    η = linspace(-1+eps(),1-eps(),npoints) |> collect
    intX = Jacobi.interp_mat(ξ, basis.nodes)
    intY = Jacobi.interp_mat(η, basis.nodes)

    Nx = numcells(meshx)
    Ny = numcells(meshy)
    Pp1 = length(basis.nodes)
    tmp = zeros(eltype(u), Pp1, npoints)
    uval = zeros(eltype(u), npoints, npoints)
    for iy in 1:Ny, ix in 1:Nx
        xmin, xmax = bounds(ix, meshx)
        ymin, ymax = bounds(iy, meshy)
        x = map_from_canonical(ξ, xmin, xmax, basis)
        y = map_from_canonical(η, ymin, ymax, basis)

        for nx in 1:Pp1
            A_mul_B!(view(tmp,nx,:), intY, view(u,nx,:,ix,iy))
        end
        A_mul_B!(uval, intX, tmp)

        for jy in 1:npoints, jx in 1:npoints
            xplot[(ix-1)*npoints+(jx-1)+1] = x[jx]
            yplot[(iy-1)*npoints+(jy-1)+1] = y[jy]
            uplot[(iy-1)*npoints+(jy-1)+1, (ix-1)*npoints+(jx-1)+1] = uval[jx,jy]
        end
    end

    xplot, yplot, uplot
end


function compute_coefficients(u, meshx::AbstractMesh1D, meshy::AbstractMesh1D,
                                 meshz::AbstractMesh1D, basis::NodalBasis)
    xmin, xmax = bounds(meshx)
    ymin, ymax = bounds(meshy)
    zmin, zmax = bounds(meshz)
    Pp1 = length(basis.nodes)
    uval = zeros(typeof(u((xmin+xmax)/2,(ymin+ymax)/2,(zmin+zmax)/2)), Pp1, Pp1, Pp1,
                    numcells(meshx), numcells(meshy), numcells(meshz))
    compute_coefficients!(uval, u, meshx, meshy, meshz, basis)
    uval
end

function compute_coefficients!(uval, u, meshx::AbstractMesh1D, meshy::AbstractMesh1D,
                                        meshz::AbstractMesh1D, basis::NodalBasis)
    @unpack nodes = basis
    Pp1 = length(nodes)

    for iz in cell_indices(meshz), iy in cell_indices(meshy), ix in cell_indices(meshx)
        xmin, xmax = bounds(ix, meshx)
        ymin, ymax = bounds(iy, meshy)
        zmin, zmax = bounds(iz, meshz)
        for nz in 1:Pp1, ny in 1:Pp1, nx in 1:Pp1
            x = map_from_canonical(nodes[nx], xmin, xmax, basis)
            y = map_from_canonical(nodes[ny], ymin, ymax, basis)
            z = map_from_canonical(nodes[nz], zmin, zmax, basis)
            uval[nx,ny,nz,ix,iy,iz] = u(x,y,z)
        end
    end

    nothing
end
