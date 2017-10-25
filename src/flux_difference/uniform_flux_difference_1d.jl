
struct UniformFluxDiffDisc1D{BalanceLaw, T, Basis, Fvol, Fnumint, Fnumext, LeftBC, RightBC,
                                Fluxes, Parallel} <: AbstractSemidiscretisation
    balance_law::BalanceLaw
    meshx::UniformMesh1D{T}
    basis::Basis
    fvol::Fvol
    fnumint::Fnumint
    fnumext::Fnumext
    left_bc::LeftBC
    right_bc::RightBC

    fluxes::Fluxes
    parallel::Parallel

    function UniformFluxDiffDisc1D(balance_law::BalanceLaw, meshx::UniformMesh1D{T},
                                    basis::Basis, fvol::Fvol, fnumint::Fnumint,
                                    fnumext::Fnumext, left_bc::LeftBC, right_bc::RightBC,
                                    fluxes::Fluxes, parallel::Parallel) where
                                        {BalanceLaw, T, Basis, Fvol, Fnumint, Fnumext,
                                        LeftBC, RightBC, Fluxes, Parallel}
        @assert typeof(fluxes) <: AbstractArray{variables(balance_law), 1}
        @assert size(fluxes, 1) == numedges(meshx)
        new{BalanceLaw, T, Basis, Fvol, Fnumint, Fnumext, LeftBC, RightBC, Fluxes, Parallel}(
            balance_law, meshx, basis, fvol, fnumint, fnumext, left_bc, right_bc,
            fluxes, parallel)
    end
end

function UniformFluxDiffDisc1D(balance_law, meshx, basis, fvol, fnumint, fnumext,
                                left_bc, right_bc, parallel=Val{:serial}())
    fluxes = Array{variables(balance_law)}(numedges(meshx))
    UniformFluxDiffDisc1D(balance_law, meshx, basis, fvol, fnumint, fnumext,
                            left_bc, right_bc, fluxes, parallel)
end


@noinline function (semidisc::UniformFluxDiffDisc1D)(t, u, du)
  @boundscheck begin
    if size(u) != size(du)
      error("size(u) = $(size(u)) != $(size(du)) = size(du)")
    end
    @assert size(u,1) == length(semidisc.basis.nodes)
    @assert size(u,2) == numcells(semidisc.meshx)
    if eltype(u) != variables(semidisc.balance_law)
      error("eltype(u) == $(eltype(u)) != $(variables(semidisc.balance_law)) == variables(semidisc.balance_law)")
    end
  end

  fill!(du, zero(eltype(du)))
  add_flux_differences!(du, u, semidisc)
  add_numerical_fluxes!(du, u, semidisc, t)

  nothing
end


function add_flux_differences!(du, u, semidisc::UniformFluxDiffDisc1D)
    Nx  = size(u, 2)
    Pp1 = size(u, 1)

    @unpack balance_law, fvol, basis, parallel = semidisc
    @unpack D = basis
    @boundscheck @assert Pp1 == size(D,1) == size(D,2)
    jacx = 2 / semidisc.meshx.Δx

    if Pp1 <= 1
        return nothing
    end

    add_flux_differences_inner_loop!(du, u, balance_law, fvol, Nx, Pp1, D, jacx, parallel)
end

@inline function add_flux_differences_inner_loop!(du, u, balance_law, fvol, Nx, Pp1,
                                                    D, jacx, parallel)
    dims = (Pp1, Nx)

    @inbounds for ix in Base.OneTo(Nx)
        for nx in Base.OneTo(Pp1)
            idx = sub2ind(dims, nx, ix)
            u_idx = u[idx]
            # compute x derivative
            # at first for different indices
            for k in (nx+1):Pp1
                idxk = sub2ind(dims, k, ix)
                f = fvol(u_idx, u[idxk], balance_law)
                du[idx]  -= 2*jacx*D[nx,k] * f
                du[idxk] -= 2*jacx*D[k,nx] * f
            end
            # then for the diagonal element, using the consistency of the flux
            du[idx] -= 2*jacx*D[nx,nx] * flux(u_idx, balance_law)
        end
    end

    nothing
end
@inline function add_flux_differences_inner_loop!(du, u, balance_law, fvol, Nx, Pp1,
                                                    D, jacx, parallel::Val{:threads})
    dims = (Pp1, Nx)

     @inbounds Threads.@threads for ix in Base.OneTo(Nx)
        for nx in Base.OneTo(Pp1)
            idx = sub2ind(dims, nx, ix)
            u_idx = u[idx]
            # compute x derivative
            # at first for different indices
            for k in (nx+1):Pp1
                idxk = sub2ind(dims, k, ix)
                f = fvol(u_idx, u[idxk], balance_law)
                du[idx]  -= 2*jacx*D[nx,k] * f
                du[idxk] -= 2*jacx*D[k,nx] * f
            end
            # then for the diagonal element, using the consistency of the flux
            du[idx] -= 2*jacx*D[nx,nx] * flux(u_idx, balance_law)
        end
    end

    nothing
end


function add_numerical_fluxes!(du, u, semidisc::UniformFluxDiffDisc1D, t)
    Nx  = size(u, 2)
    Pp1 = size(u, 1)

    @unpack balance_law, basis, fnumint, fnumext, fluxes, parallel, left_bc, right_bc = semidisc
    @boundscheck @assert Pp1 == length(basis.weights)

    jacx = 2 / semidisc.meshx.Δx

    add_numerical_fluxes_inner_loop1!(du, fluxes, u, t, balance_law, fnumint,
                                        fnumext, left_bc, right_bc, Nx, basis, jacx, parallel)
    add_numerical_fluxes_inner_loop2!(du, fluxes, u, balance_law, Nx, basis, jacx, parallel)
end

@inline function add_numerical_fluxes_inner_loop1!(du, fluxes, u, t, balance_law,
                                                    fnumint, fnumext, left_bc, right_bc,
                                                    Nx, basis::LobattoLegendre,
                                                    jacx, parallel)
    # calculate external numerical fluxes
    @inbounds fluxes[1] = fnumext(left_bc(t), u[1,1], balance_law)
    @inbounds fluxes[end] = fnumext(u[end,end], right_bc(t), balance_law)

    # calculate internal numerical fluxes
    @inbounds for ix in 2:Nx
        # flux x - left
        fluxes[ix] = fnumint(u[end,ix-1], u[1,ix], balance_law)
    end

    nothing
end
@inline function add_numerical_fluxes_inner_loop1!(du, fluxes, u, t, balance_law,
                                                    fnumint, fnumext, left_bc, right_bc,
                                                    Nx, basis::LobattoLegendre,
                                                    jacx, parallel::Val{:threads})
    # calculate external numerical fluxes
    @inbounds fluxes[1] = fnumext(left_bc(t), u[1,1], balance_law)
    @inbounds fluxes[end] = fnumext(u[end,end], right_bc(t), balance_law)

    # calculate internal numerical fluxes
    @inbounds Threads.@threads for ix in 2:Nx
        # flux x - left
        fluxes[ix] = fnumint(u[end,ix-1], u[1,ix], balance_law)
    end

    nothing
end

@inline function add_numerical_fluxes_inner_loop2!(du, fluxes, u, balance_law,
                                                    Nx, basis::LobattoLegendre,
                                                    jacx, parallel)
    @inbounds i_ω1 = 1 / basis.weights[1]
    @inbounds i_ωend = 1 / basis.weights[end]
    Pp1 = length(basis.nodes)
    dims = (Pp1, Nx)

    # add numerical fluxes
    @inbounds for ix in Base.OneTo(Nx)
        # flux x - left
        idx = sub2ind(dims, 1, ix)
        du[idx] += ( fluxes[ix] - flux(u[idx], balance_law) ) * jacx * i_ω1

        # flux x - right
        idx = sub2ind(dims, Pp1, ix)
        du[idx] -= ( fluxes[ix+1] - flux(u[idx], balance_law) ) * jacx * i_ωend
    end

    nothing
end
@inline function add_numerical_fluxes_inner_loop2!(du, fluxes, u, balance_law,
                                                    Nx, basis::LobattoLegendre,
                                                    jacx, parallel::Val{:threads})
    @inbounds i_ω1 = 1 / basis.weights[1]
    @inbounds i_ωend = 1 / basis.weights[end]
    Pp1 = length(basis.nodes)
    dims = (Pp1, Nx)

    # add numerical fluxes
    @inbounds Threads.@threads for ix in Base.OneTo(Nx)
        # flux x - left
        idx = sub2ind(dims, 1, ix)
        du[idx] += ( fluxes[ix] - flux(u[idx], balance_law) ) * jacx * i_ω1

        # flux x - right
        idx = sub2ind(dims, Pp1, ix)
        du[idx] -= ( fluxes[ix+1] - flux(u[idx], balance_law) ) * jacx * i_ωend
    end

    nothing
end


function semidiscretise(semidisc::UniformFluxDiffDisc1D, u₀func, tspan)
    @unpack meshx, basis, parallel = semidisc
    u₀ = compute_coefficients(u₀func, meshx, basis, parallel)

    ODEProblem(semidisc, u₀, tspan)
end

"""
    max_dt(t, u, semidisc::UniformFluxDiffDisc1D, cfl)

Compute the maximal time step `dt` satisfying the CFL condition `dt <= cfl * dx / ((2p+1)*λ)`,
where `dx` is the length of a cell, `λ` the greatest absolute value of the
speed in this cell, and `p` the polynomial degree.
"""
function max_dt(t, u, semidisc::UniformFluxDiffDisc1D, cfl=0.5)
    @unpack balance_law, meshx, basis = semidisc
    @unpack Δx = meshx
    p = length(basis.nodes)-1
    factor = Δx / (2p+1)

    dt = mapreduce(u->factor/max_abs_speed(u, balance_law), min, typemax(t), u)
    if dt == Inf
        dt = factor
    end
    cfl * dt
end
