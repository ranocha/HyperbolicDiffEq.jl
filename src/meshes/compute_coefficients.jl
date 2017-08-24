
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
