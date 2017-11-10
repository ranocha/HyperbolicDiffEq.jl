"""
An abstract mesh in one space dimension.
"""
abstract type AbstractMesh1D <: AbstractMesh end

"""
An abstract periodic mesh in one space dimension.
"""
abstract type AbstractPeriodicMesh1D <: AbstractMesh1D end


isperiodic(mesh) = false
isperiodic(mesh::AbstractPeriodicMesh1D) = true

@inline numedges(mesh::AbstractMesh1D) = numcells(mesh)+1
@inline edge_indices(mesh::AbstractMesh1D) = 1:numedges(mesh)


"""
    left_edge(cell::Int, mesh::AbstractMesh1D)

The index of the edge to the left of `cell` in `mesh`.
"""
@inline function left_edge(cell::Int, mesh::AbstractMesh1D)
    @boundscheck begin
        assert(1 <= cell <= numcells(mesh))
    end
    cell
end

"""
    right_edge(cell::Int, mesh::AbstractMesh1D)

The index of the edge to the right of `cell` in `mesh`.
"""
@inline function right_edge(cell::Int, mesh::AbstractMesh1D)
    @boundscheck begin
        assert(1 <= cell <= numcells(mesh))
    end
    cell+1
end


"""
    left_cell(edge::Int, mesh::AbstractMesh1D)

The index of the cell to the left of `edge` in `mesh`.
"""
@inline function left_cell(edge::Int, mesh::AbstractMesh1D)
    @boundscheck begin
        assert(1 <= edge <= numedges(mesh))
    end
    if isperiodic(mesh)
        edge == 1 ? numcells(mesh) : edge-1
    else
        error("To be implemented.")
    end
end

"""
    right_cell(edge::Int, mesh::AbstractMesh1D)

The index of the cell to the right of `edge` in `mesh`.
"""
@inline function right_cell(edge::Int, mesh::AbstractMesh1D)
    @boundscheck begin
        assert(1 <= edge <= numedges(mesh))
    end
    if isperiodic(mesh)
        edge == numedges(mesh) ? 1 : edge
    else
        error("To be implemented.")
    end
end



function Base.show(io::IO, mesh::AbstractMesh1D)
  print(io,
    typeof(mesh), " with ", numcells(mesh), " cells in ", bounds(mesh))
end


@recipe function f(mesh::AbstractMesh1D; add_marker = false)
  markershape --> (add_marker ? :vline : :none)
  delete!(d, :add_marker)

  legend := false
  ylim := (-1,1)
  yticks --> false
  grid --> false
  border --> false
  xlabel --> L"x"
  size --> (600,50)
  color --> [:blue :orange]

  xlims = zeros(2, numcells(mesh))
  for cell in 1:numcells(mesh)
    xlims[1,cell], xlims[2,cell] = bounds(cell,mesh)
  end

  xlims, zeros(xlims)
end



################################################################################

"""
A uniform periodic mesh in one space dimension of `Nx` cells between
`xmin` and `xmax`.
"""
immutable UniformPeriodicMesh1D{T<:Real} <: AbstractPeriodicMesh1D
  xmin::T
  xmax::T
  Nx::Int

  Δx::T

  function UniformPeriodicMesh1D{T}(xmin::T, xmax::T, Nx::Int) where T
    Nx > 0 || error("The number of elements `Nx` must be positive [`Nx == $Nx`].")
    xmin < xmax || error("`xmin` must be smaller than `xmax` [`xmin == $xmin, xmax == $xmax`].")

    new(xmin, xmax, Nx, (xmax-xmin)/Nx)
  end
end

function UniformPeriodicMesh1D(xmin::Real, xmax::Real, Nx::Integer)
  xmin, xmax = promote(xmin, xmax)
  UniformPeriodicMesh1D{typeof(xmin)}(xmin, xmax, Int(Nx))
end

@inline numcells(mesh::UniformPeriodicMesh1D) = mesh.Nx
@inline cell_indices(mesh::UniformPeriodicMesh1D) = 1:numcells(mesh)
@inline volume(cell::Int, mesh::UniformPeriodicMesh1D) = mesh.Δx

function bounds(cell::Int, mesh::UniformPeriodicMesh1D)
  @unpack Δx = mesh

  xmin = mesh.xmin + (cell-1)*Δx
  xmin = xmin + eps(xmin)

  xmax = mesh.xmin + cell*Δx
  xmax = xmax - eps(xmax)

  xmin, xmax
end

@inline bounds(mesh::UniformPeriodicMesh1D) = mesh.xmin, mesh.xmax



################################################################################

"""
A uniform mesh in one space dimension of `Nx` cells between `xmin` and `xmax`.
"""
immutable UniformMesh1D{T<:Real} <: AbstractMesh1D
  xmin::T
  xmax::T
  Nx::Int

  Δx::T

  function UniformMesh1D{T}(xmin::T, xmax::T, Nx::Int) where T
    Nx > 0 || error("The number of elements `Nx` must be positive [`Nx == $Nx`].")
    xmin < xmax || error("`xmin` must be smaller than `xmax` [`xmin == $xmin, xmax == $xmax`].")

    new(xmin, xmax, Nx, (xmax-xmin)/Nx)
  end
end

function UniformMesh1D(xmin::Real, xmax::Real, Nx::Integer)
  xmin, xmax = promote(xmin, xmax)
  UniformMesh1D{typeof(xmin)}(xmin, xmax, Int(Nx))
end

@inline numcells(mesh::UniformMesh1D) = mesh.Nx
@inline cell_indices(mesh::UniformMesh1D) = 1:numcells(mesh)
@inline volume(cell::Int, mesh::UniformMesh1D) = mesh.Δx

function bounds(cell::Int, mesh::UniformMesh1D)
  @unpack Δx = mesh

  xmin = mesh.xmin + (cell-1)*Δx
  xmin = xmin + eps(xmin)

  xmax = mesh.xmin + cell*Δx
  xmax = xmax - eps(xmax)

  xmin, xmax
end

@inline bounds(mesh::UniformMesh1D) = mesh.xmin, mesh.xmax
