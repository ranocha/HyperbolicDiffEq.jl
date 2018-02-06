
struct UniformPeriodicReconstructedFV1D{BalanceLaw, T, Fnum, Reconstruction,
                                        EdgeU, Fluxes, Parallel} <: FiniteVolumeSemidiscretisation
    balance_law::BalanceLaw
    meshx::UniformPeriodicMesh1D{T}
    fnum::Fnum
    reconstruction::Reconstruction

    edge_u::EdgeU
    fluxes::Fluxes
    parallel::Parallel

    function UniformPeriodicReconstructedFV1D(
            balance_law::BalanceLaw, meshx::UniformPeriodicMesh1D{T},
            fnum::Fnum, reconstruction::Reconstruction,
            edge_u::EdgeU, fluxes::Fluxes, parallel::Parallel) where {BalanceLaw, T, Fnum, Reconstruction, EdgeU, Fluxes, Parallel}
        @assert typeof(fluxes) <: AbstractArray{variables(balance_law), 1}
        @assert typeof(edge_u) <: AbstractArray{variables(balance_law), 2}
        @assert size(fluxes, 1) == size(edge_u, 2) == numcells(meshx)
        @assert size(edge_u, 1) == 2
        new{BalanceLaw, T, Fnum, Reconstruction, EdgeU, Fluxes, Parallel}(
            balance_law, meshx, fnum, reconstruction, edge_u, fluxes, parallel)
    end
end

function UniformPeriodicReconstructedFV1D(balance_law, meshx, fnum, reconstruction,
                                          parallel=Val{:serial}())
    fluxes = Array{variables(balance_law)}(numcells(meshx))
    edge_u = Array{variables(balance_law)}(2, numcells(meshx))
    UniformPeriodicReconstructedFV1D(balance_law, meshx, fnum, reconstruction,
                                     edge_u, fluxes, parallel)
end

function Base.show(io::IO, fv::UniformPeriodicReconstructedFV1D)
    print(io, "Finite volume method ",
            "\n  Balance law:    ", fv.balance_law,
            "\n  Mesh:           ", fv.meshx,
            "\n  Numerical flux: ", fv.fnum,
            "\n  Reconstruction: ", fv.reconstruction)
end

function evaluate_coefficients(u, fv::UniformPeriodicReconstructedFV1D,
                               npoints=2*order(fv.reconstruction))
    evaluate_coefficients(u, fv.balance_law, fv.meshx, fv.reconstruction, npoints)
end


@noinline function (semidisc::UniformPeriodicReconstructedFV1D)(du, u, p, t)
    @boundscheck begin
        if size(u) != size(du)
            error("size(u) = $(size(u)) != $(size(du)) = size(du)")
        end
        @assert length(u) == numcells(semidisc.meshx)
        @assert length(u) >= stencil_width(semidisc.reconstruction) - 1
        if eltype(u) != variables(semidisc.balance_law)
            error("eltype(u) == $(eltype(u)) != $(variables(semidisc.balance_law)) == variables(semidisc.balance_law)")
        end
    end

    fill!(du, zero(eltype(du)))
    add_numerical_fluxes!(du, u, semidisc, t)

    nothing
end

function add_numerical_fluxes!(du, u, fv::UniformPeriodicReconstructedFV1D, t)
    @unpack balance_law, meshx, fnum, reconstruction, edge_u, fluxes, parallel = fv

    reconstruct!(edge_u, u, balance_law, meshx, reconstruction, stencil_width_val(reconstruction), parallel)
    compute_numerical_fluxes!(fluxes, edge_u, balance_law, meshx, fnum, parallel)
    compute_du!(du, fluxes, meshx, parallel)
end


function compute_numerical_fluxes!(fluxes, edge_u, balance_law, meshx::UniformPeriodicMesh1D,
                                   fnum, parallel)
    fluxes[1] = fnum(edge_u[2,end], edge_u[1,1], balance_law)
    for edge in 2:length(fluxes)
        fluxes[edge] = fnum(edge_u[2,edge-1], edge_u[1,edge], balance_law)
    end
    nothing
end

function compute_numerical_fluxes!(fluxes, edge_u, balance_law, meshx::UniformPeriodicMesh1D,
                                   fnum, parallel::Val{:threads})
    @inbounds fluxes[1] = fnum(edge_u[2,end], edge_u[1,1], balance_law)
    @inbounds Threads.@threads for edge in 2:length(fluxes)
        fluxes[edge] = fnum(edge_u[2,edge-1], edge_u[1,edge], balance_law)
    end
    nothing
end


function compute_du!(du, fluxes, meshx::UniformPeriodicMesh1D, parallel::Val{:serial})
    @inbounds for cell in 1:numcells(meshx)-1
        du[cell] -= ( fluxes[cell+1] - fluxes[cell] ) / volume(cell, meshx)
    end
    @inbounds du[end] -= ( fluxes[1] - fluxes[end] ) / volume(numcells(meshx), meshx)
    nothing
end

function compute_du!(du, fluxes, meshx::UniformPeriodicMesh1D, parallel::Val{:threads})
    @inbounds Threads.@threads for cell in 1:numcells(meshx)-1
        du[cell] -= ( fluxes[cell+1] - fluxes[cell] ) / volume(cell, meshx)
    end
    @inbounds du[end] -= ( fluxes[1] - fluxes[end] ) / volume(numcells(meshx), meshx)
    nothing
end

################################################################################

function semidiscretise(fv::UniformPeriodicReconstructedFV1D, u₀func, tspan)
    u₀ = compute_coefficients(u₀func, fv.meshx, order(fv.reconstruction))

    ODEProblem(fv, u₀, tspan)
end

"""
    max_dt(t, u, fv::UniformPeriodicReconstructedFV1D, cfl=0.5)

Compute the maximal time step `dt` satisfying the CFL condition `dt <= cfl * dx / 2λ`,
where `dx` is the length of a cell and `λ` the greatest absolute value of the
speed in this cell.
"""
function max_dt(t, u, fv::UniformPeriodicReconstructedFV1D, cfl=0.5)
    @unpack balance_law, meshx = fv
    dt = mapreduce(cell->volume(cell,meshx) / max_abs_speed(u[cell],balance_law),
                    min, typemax(t), cell_indices(meshx))
    if dt == Inf
        dt = mapreduce(cell->volume(cell,meshx), min, dt, cell_indices(meshx))
    end
    cfl * dt
end
