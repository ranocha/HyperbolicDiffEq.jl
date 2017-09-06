
"""
    compute_coefficients(u, mesh::AbstractMesh1D)

Computse the coefficients of the function `u` in a first order finite volume
method on the `mesh`.
"""
function compute_coefficients(u, mesh::AbstractMesh1D)
    xmin, xmax = bounds(mesh)
    uval = zeros(typeof(u((xmin+xmax)/2)), numcells(mesh))
    compute_coefficients!(uval, u, mesh)
    uval
end

"""
    compute_coefficients!(uval, u, mesh::AbstractMesh1D)

Computse the coefficients of the function `u` in a first order finite volume
method on the `mesh` and stores them in `uval`.
"""
function compute_coefficients!(uval, u, mesh::AbstractMesh1D)
    for cell in cell_indices(mesh)
        xmin, xmax = bounds(cell, mesh)
        uval[cell] = u((xmin+xmax)/2)
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


################################################################################

function compute_coefficients(u, mesh::AbstractMesh1D, basis::NodalBasis)
    xmin, xmax = bounds(mesh)
    uval = zeros(typeof(u((xmin+xmax)/2)), length(basis.nodes), numcells(mesh))
    compute_coefficients!(uval, u, mesh)
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
