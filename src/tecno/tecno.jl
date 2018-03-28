
struct ScalarUniformPeriodicTecno1D{BalanceLaw, T, Fnum, Diss, Reconstruction,
                                    ComputeEntropyVar, EntropyVar, EdgeEntropyVar, 
                                    Fluxes, Parallel} <: FiniteVolumeSemidiscretisation
    balance_law::BalanceLaw
    meshx::UniformPeriodicMesh1D{T}
    fnum::Fnum
    diss::Diss
    reconstruction::Reconstruction

    compute_entropy_var::ComputeEntropyVar
    entropy_var::EntropyVar
    edge_entropy_var::EdgeEntropyVar
    fluxes::Fluxes
    parallel::Parallel

    function ScalarUniformPeriodicTecno1D(
            balance_law::BalanceLaw, meshx::UniformPeriodicMesh1D{T},
            fnum::Fnum, diss::Diss, reconstruction::Reconstruction,
            compute_entropy_var::ComputeEntropyVar,
            entropy_var::EntropyVar, edge_entropy_var::EdgeEntropyVar,
            fluxes::Fluxes, parallel::Parallel) where {BalanceLaw, T, Fnum, Diss, Reconstruction, ComputeEntropyVar, EntropyVar, EdgeEntropyVar, Fluxes, Parallel}
        @assert typeof(fluxes) <: AbstractArray{variables(balance_law), 1}
        @assert typeof(entropy_var) <: AbstractArray{variables(balance_law), 1}
        @assert typeof(edge_entropy_var) <: AbstractArray{variables(balance_law), 2}
        @assert size(fluxes, 1) == size(entropy_var, 1) == size(edge_entropy_var, 2) == numcells(meshx)
        @assert size(edge_entropy_var, 1) == 2
        new{BalanceLaw, T, Fnum, Diss, Reconstruction, ComputeEntropyVar, EntropyVar, EdgeEntropyVar, Fluxes, Parallel}(
            balance_law, meshx, fnum, diss, reconstruction, compute_entropy_var, entropy_var, edge_entropy_var, fluxes, parallel)
    end
end

function ScalarUniformPeriodicTecno1D(balance_law, meshx, fnum, diss, reconstruction,
                                      compute_entropy_var, parallel=Val{:serial}())
    fluxes = Array{variables(balance_law)}(numcells(meshx))
    entropy_var = Array{variables(balance_law)}(numcells(meshx))
    edge_entropy_var = Array{variables(balance_law)}(2, numcells(meshx))
    ScalarUniformPeriodicTecno1D(balance_law, meshx, fnum, diss, reconstruction,
                                 compute_entropy_var, entropy_var, edge_entropy_var, fluxes, parallel)
end

function Base.show(io::IO, fv::ScalarUniformPeriodicTecno1D)
    print(io, "TeCNO finite volume method ",
            "\n  Balance law:    ", fv.balance_law,
            "\n  Mesh:           ", fv.meshx,
            "\n  Numerical flux: ", fv.fnum,
            "\n  Dissipation:    ", fv.diss,
            "\n  Reconstruction: ", fv.reconstruction)
end

function evaluate_coefficients(u, fv::ScalarUniformPeriodicTecno1D,
                               npoints=2*order(fv.reconstruction))
    evaluate_coefficients(u, fv.balance_law, fv.meshx, fv.reconstruction, npoints)
end


@noinline function (semidisc::ScalarUniformPeriodicTecno1D)(du, u, p, t)
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

function add_numerical_fluxes!(du, u, fv::ScalarUniformPeriodicTecno1D, t)
    @unpack balance_law, meshx, fnum, diss, reconstruction, compute_entropy_var, entropy_var, edge_entropy_var, fluxes, parallel = fv

    # L² entropy
    entropy_var .= compute_entropy_var.(u, balance_law)
    reconstruct!(edge_entropy_var, entropy_var, balance_law, meshx, reconstruction, stencil_width_val(reconstruction), parallel)
    compute_EC_fluxes!(fluxes, u, balance_law, meshx, fnum, stencil_width_val(reconstruction), parallel)
    compute_ES_fluxes!(fluxes, edge_entropy_var, balance_law, meshx, diss, parallel)
    compute_du!(du, fluxes, meshx, parallel)
end


# Entropy conservative fluxes
function compute_EC_fluxes!(fluxes, u, balance_law, meshx::UniformPeriodicMesh1D,
                            fnum, ::Val{1}, parallel)
    @inbounds fluxes[1] = fnum(u[end], u[1], balance_law)
    @inbounds for edge in 2:length(fluxes)
        fluxes[edge] = fnum(u[edge-1], u[edge], balance_law)
    end
    nothing
end

function compute_EC_fluxes!(fluxes, u, balance_law, meshx::UniformPeriodicMesh1D,
                            fnum, ::Val{3}, parallel)
    T = eltype(u)
    c1 = T(4) / 3
    c2 = T(-1) / 6

    @inbounds begin
        fluxes[1] = c1 * fnum(u[end], u[1], balance_law) +
                        c2 * fnum(u[end-1], u[1], balance_law) + 
                        c2 * fnum(u[end], u[2], balance_law)
        fluxes[2] = c1 * fnum(u[1], u[2], balance_law) +
                        c2 * fnum(u[end], u[2], balance_law) + 
                        c2 * fnum(u[1], u[3], balance_law)
        fluxes[end] = c1 * fnum(u[end-1], u[end], balance_law) +
                        c2 * fnum(u[end-2], u[end], balance_law) + 
                        c2 * fnum(u[end-1], u[1], balance_law)
    end
    @inbounds for edge in 3:length(fluxes)-1
        fluxes[edge] = c1 * fnum(u[edge-1], u[edge], balance_law) +
                        c2 * fnum(u[edge-2], u[edge], balance_law) + 
                        c2 * fnum(u[edge-1], u[edge+1], balance_law)
    end
    nothing
end

function compute_EC_fluxes!(fluxes, u, balance_law, meshx::UniformPeriodicMesh1D,
                            fnum, ::Val{5}, parallel)
    T = eltype(u)
    c1 = T(3) / 2
    c2 = T(-3) / 10
    c3 = T(1) / 30

    @inbounds begin
        fluxes[1] = c1 * fnum(u[end], u[1], balance_law) +
                        c2 * fnum(u[end-1], u[1], balance_law) + 
                        c2 * fnum(u[end], u[2], balance_law) + 
                        c3 * fnum(u[end-2], u[1], balance_law) + 
                        c3 * fnum(u[end-1], u[2], balance_law) + 
                        c3 * fnum(u[end], u[3], balance_law)
        fluxes[2] = c1 * fnum(u[1], u[2], balance_law) +
                        c2 * fnum(u[end], u[2], balance_law) + 
                        c2 * fnum(u[1], u[3], balance_law) + 
                        c3 * fnum(u[end-1], u[2], balance_law) + 
                        c3 * fnum(u[end], u[3], balance_law) + 
                        c3 * fnum(u[1], u[4], balance_law)
        fluxes[3] = c1 * fnum(u[2], u[3], balance_law) +
                        c2 * fnum(u[1], u[3], balance_law) + 
                        c2 * fnum(u[2], u[4], balance_law) + 
                        c3 * fnum(u[end], u[3], balance_law) + 
                        c3 * fnum(u[1], u[4], balance_law) + 
                        c3 * fnum(u[2], u[5], balance_law)
        fluxes[end-1] = c1 * fnum(u[end-2], u[end-1], balance_law) +
                        c2 * fnum(u[end-3], u[end-1], balance_law) + 
                        c2 * fnum(u[end-2], u[end], balance_law) + 
                        c3 * fnum(u[end-4], u[end-1], balance_law) + 
                        c3 * fnum(u[end-3], u[end], balance_law) + 
                        c3 * fnum(u[end-2], u[1], balance_law)
        fluxes[end] = c1 * fnum(u[end-1], u[end], balance_law) +
                        c2 * fnum(u[end-2], u[end], balance_law) + 
                        c2 * fnum(u[end-1], u[1], balance_law) + 
                        c3 * fnum(u[end-3], u[end], balance_law) + 
                        c3 * fnum(u[end-2], u[1], balance_law) + 
                        c3 * fnum(u[end-1], u[2], balance_law)
    end
    @inbounds for edge in 4:length(fluxes)-2
        fluxes[edge] = c1 * fnum(u[edge-1], u[edge], balance_law) +
                        c2 * fnum(u[edge-2], u[edge], balance_law) + 
                        c2 * fnum(u[edge-1], u[edge+1], balance_law) + 
                        c3 * fnum(u[edge-3], u[edge], balance_law) + 
                        c3 * fnum(u[edge-2], u[edge+1], balance_law) + 
                        c3 * fnum(u[edge-1], u[edge+2], balance_law)
    end
    nothing
end

# Add dissipation operators to get entropy stable fluxes
function compute_ES_fluxes!(fluxes, edge_entropy_var, balance_law, meshx::UniformPeriodicMesh1D,
                            diss, parallel)
    @inbounds fluxes[1] += diss(edge_entropy_var[2,end], edge_entropy_var[1,1], balance_law)
    @inbounds for edge in 2:length(fluxes)
        fluxes[edge] += diss(edge_entropy_var[2,edge-1], edge_entropy_var[1,edge], balance_law)
    end
    nothing
end

function compute_ES_fluxes!(fluxes, edge_entropy_var, balance_law, meshx::UniformPeriodicMesh1D,
                            diss, parallel::Val{:threads})
    @inbounds fluxes[1] += diss(edge_entropy_var[2,end], edge_entropy_var[1,1], balance_law)
    @inbounds Threads.@threads for edge in 2:length(fluxes)
        fluxes[edge] += diss(edge_entropy_var[2,edge-1], edge_entropy_var[1,edge], balance_law)
    end
    nothing
end

################################################################################

function semidiscretise(fv::ScalarUniformPeriodicTecno1D, u₀func, tspan)
    u₀ = compute_coefficients(u₀func, fv.meshx, order(fv.reconstruction))

    ODEProblem(fv, u₀, tspan)
end

"""
    max_dt(t, u, fv::ScalarUniformPeriodicTecno1D, cfl=0.5)

Compute the maximal time step `dt` satisfying the CFL condition `dt <= cfl * dx / 2λ`,
where `dx` is the length of a cell and `λ` the greatest absolute value of the
speed in this cell.
"""
function max_dt(t, u, fv::ScalarUniformPeriodicTecno1D, cfl=0.5)
    @unpack balance_law, meshx = fv
    dt = mapreduce(cell->volume(cell,meshx) / max_abs_speed(u[cell],balance_law),
                    min, typemax(t), cell_indices(meshx))
    if dt == Inf
        dt = mapreduce(cell->volume(cell,meshx), min, dt, cell_indices(meshx))
    end
    cfl * dt
end

function integrate(func, u, fv::ScalarUniformPeriodicTecno1D)
    integrate(func, u, fv.meshx)
end