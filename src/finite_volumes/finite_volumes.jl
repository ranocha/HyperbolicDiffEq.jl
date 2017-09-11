
abstract type FiniteVolumeSemidiscretisation <: AbstractSemidiscretisation end

struct FirstOrderFV{BalanceLaw, Mesh, NumFlux,
                    Fluxes, UseThreads} <: FiniteVolumeSemidiscretisation
    balance_law::BalanceLaw
    mesh::Mesh
    numflux::NumFlux

    fluxes::Fluxes
end

function FirstOrderFV(balance_law::AbstractBalanceLaw{1}, mesh::AbstractMesh1D,
                        numflux, usethreads=false)
    if !isperiodic(mesh)
        error("The mesh must be periodic in order to infer boundary conditions.")
    end

    U = variables(balance_law)
    u = zero(U)

    F = typeof(numflux(u,u,balance_law))
    fluxes = zeros(F, numedges(mesh))

    FirstOrderFV{typeof(balance_law), typeof(mesh), typeof(numflux),
                    typeof(fluxes), Val{usethreads}}(
        balance_law, mesh, numflux, fluxes
    )
end


@inline function usethreads{BalanceLaw, Mesh, NumFLux, Fluxes}(
    fv::FirstOrderFV{BalanceLaw, Mesh, NumFLux, Fluxes, Val{false}})
  false
end

@inline function usethreads{BalanceLaw, Mesh, NumFLux, Fluxes}(
    fv::FirstOrderFV{BalanceLaw, Mesh, NumFLux, Fluxes, Val{true}})
  true
end

function show(io::IO, fv::FirstOrderFV)
  print(io, "First order finite volume method ",
        "\n  Balance law:   ", fv.balance_law,
        "\n  Mesh:          ", fv.mesh,
        "\n  Numerical flux: ", fv.numflux)
end



"""
    (fv::FirstOrderFV)(t, u, du)

Apply a first order finite volume semidiscretisation.
"""
function (fv::FirstOrderFV)(t, u, du)
  @boundscheck begin
    if length(u) != length(du)
      error("length(u) = $(length(u)) != $(length(du)) = length(du)")
    end
    length(u) != numcells(fv.mesh) && error("length(u) != numcells(fv.mesh)")
    if eltype(u) != variables(fv.balance_law)
      error("eltype(u) == $(eltype(u)) != $(variables(fv.balance_law)) == variables(fv.balance_law)")
    end
  end

  @unpack balance_law, mesh, numflux, fluxes = fv

  compute_fluxes!(fluxes, numflux, u, mesh, balance_law, Val{usethreads(fv)})
  compute_du!(du, fluxes, mesh, Val{usethreads(fv)})

  nothing
end


function compute_fluxes!(fluxes, numflux, u, mesh::AbstractMesh1D, balance_law, ::Type{Val{true}})
    Threads.@threads for edge in edge_indices(mesh)
        @inbounds left = left_cell(edge, mesh)
        @inbounds right = right_cell(edge, mesh)
        @inbounds fluxes[edge] = numflux(u[left], u[right], balance_law)
    end
end

function compute_fluxes!(fluxes, numflux, u, mesh::AbstractMesh1D, balance_law, ::Type{Val{false}})
    for edge in edge_indices(mesh)
        @inbounds left = left_cell(edge, mesh)
        @inbounds right = right_cell(edge, mesh)
        @inbounds fluxes[edge] = numflux(u[left], u[right], balance_law)
    end
end


function compute_du!(du, fluxes, mesh::AbstractMesh1D, ::Type{Val{true}})
    Threads.@threads for cell in cell_indices(mesh)
        @inbounds left = left_edge(cell, mesh)
        @inbounds right = right_edge(cell, mesh)
        @inbounds du[cell] = -( fluxes[right] - fluxes[left] ) / volume(cell, mesh)
    end
end

function compute_du!(du, fluxes, mesh::AbstractMesh1D, ::Type{Val{false}})
    for cell in cell_indices(mesh)
        @inbounds left = left_edge(cell, mesh)
        @inbounds right = right_edge(cell, mesh)
        @inbounds du[cell] = -( fluxes[right] - fluxes[left] ) / volume(cell, mesh)
    end
end

################################################################################

function semidiscretise(fv::AbstractSemidiscretisation, u₀func, tspan)
  u₀ = compute_coefficients(u₀func, fv.mesh)

  ODEProblem(fv, u₀, tspan)
end

"""
    max_dt(t, u, fv::FirstOrderFV, cfl=0.5)

Compute the maximal time step `dt` satisfying the CFL condition `dt <= cfl * dx / 2λ`,
where `dx` is the length of a cell and `λ` the greatest absolute value of the
speed in this cell.
"""
function max_dt(t, u, fv::FirstOrderFV)
    @unpack balance_law, mesh = fv
    dt = mapreduce(cell->volume(cell,mesh) / max_abs_speed(u[cell],balance_law),
                    min, typemax(t), cell_indices(mesh))
    if dt == Inf
        dt = mapreduce(cell->volume(cell,mesh), min, dt, cell_indices(mesh))
    end
    cfl * dt
end
