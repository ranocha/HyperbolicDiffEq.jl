
struct UniformPeriodicFluxDiffDisc2D{BalanceLaw, T, Fvol, Fnum, Fluxes, UseThreads<:Union{Val{true},Val{false}}} <: AbstractSemidiscretisation
    balance_law::BalanceLaw
    meshx::UniformPeriodicMesh1D{T}
    meshy::UniformPeriodicMesh1D{T}
    basis::LobattoLegendre{T}
    fvol::Fvol
    fnum::Fnum

    fluxes::Fluxes
    usethreads::UseThreads

    function UniformPeriodicFluxDiffDisc2D(balance_law::BalanceLaw, meshx::UniformPeriodicMesh1D{T}, meshy::UniformPeriodicMesh1D{T}, basis::LobattoLegendre{T}, fvol::Fvol, fnum::Fnum, fluxes::Fluxes, usethreads::UseThreads) where {BalanceLaw, T, Fvol, Fnum, Fluxes, UseThreads}
        @assert typeof(fluxes) <: AbstractArray{variables(balance_law), 4}
        @assert size(fluxes, 1) == length(basis.nodes)
        @assert size(fluxes, 2) == 2
        @assert size(fluxes, 3) == numedges(meshx)
        @assert size(fluxes, 4) == numedges(meshy)
        new{BalanceLaw, T, Fvol, Fnum, Fluxes, UseThreads}(balance_law, meshx, meshy, basis, fvol, fnum, fluxes, usethreads)
    end
end

function UniformPeriodicFluxDiffDisc2D(balance_law, meshx, meshy, basis, fvol, fnum, usethreads::Bool=false)
    fluxes = Array{variables(balance_law)}(length(basis.nodes), 2, numedges(meshx), numedges(meshy))
    UniformPeriodicFluxDiffDisc2D(balance_law, meshx, meshy, basis, fvol, fnum, fluxes, Val{usethreads}())
end


@noinline function (semidisc::UniformPeriodicFluxDiffDisc2D)(t, u, du)
  @boundscheck begin
    if size(u) != size(du)
      error("size(u) = $(size(u)) != $(size(du)) = size(du)")
    end
    @assert size(u,1) == size(u,2) == length(semidisc.basis.nodes)
    @assert size(u,3) == numcells(semidisc.meshx)
    @assert size(u,4) == numcells(semidisc.meshy)
    if eltype(u) != variables(semidisc.balance_law)
      error("eltype(u) == $(eltype(u)) != $(variables(semidisc.balance_law)) == variables(semidisc.balance_law)")
    end
  end

  fill!(du, zero(eltype(du)))
  add_flux_differences!(du, u, semidisc)
  add_numerical_fluxes!(du, u, semidisc)

  nothing
end


function add_flux_differences!(du, u, semidisc::UniformPeriodicFluxDiffDisc2D)
    Nx = size(u, 3)
    Ny = size(u, 4)
    Pp1 = size(u, 1)
    @boundscheck @assert Pp1 == size(u, 2)

    @unpack balance_law, fvol, basis, usethreads = semidisc
    @unpack D = basis
    @boundscheck @assert Pp1 == size(D,1) == size(D,2)
    jacx = 2 / semidisc.meshx.Δx
    jacy = 2 / semidisc.meshy.Δx

    if Pp1 <= 1
        return nothing
    end

    add_flux_differences_inner_loop!(du, u, balance_law, fvol, Nx, Ny, Pp1, D, jacx, jacy, usethreads)
end

@inline function add_flux_differences_inner_loop!(du, u, balance_law, fvol, Nx, Ny, Pp1,
                                                    D, jacx, jacy, ::Val{true})
    dirx = Val{:x}()
    diry = Val{:y}()
    NxNy = Nx*Ny
    dims = (Pp1,Pp1,Nx*Ny)
    @inbounds Threads.@threads for ixy in Base.OneTo(NxNy)
        for ny in Base.OneTo(Pp1), nx in Base.OneTo(Pp1)
            idx = sub2ind(dims, nx, ny, ixy)
            u_idx = u[idx]
            # compute x derivative
            # at first for different indices
            for k in (nx+1):Pp1
                idxk = sub2ind(dims, k, ny, ixy)
                f = fvol(u_idx, u[idxk], balance_law, dirx)
                du[idx]  -= 2*jacx*D[nx,k] * f
                du[idxk] -= 2*jacx*D[k,nx] * f
            end
            # then for the diagonal element, using the consistency of the flux
            du[idx] -= 2*jacx*D[nx,nx] * flux(u_idx, balance_law, dirx)

            # compute y derivative
            # at first for different indices
            for k in (ny+1):Pp1
                idxk = sub2ind(dims, nx, k, ixy)
                f = fvol(u_idx, u[idxk], balance_law, diry)
                du[idx]  -= 2*jacy*D[ny,k] * f
                du[idxk] -= 2*jacy*D[k,ny] * f
            end
            # then for the diagonal element, using the consistency of the flux
            du[idx] -= 2*jacy*D[ny,ny] * flux(u_idx, balance_law, diry)
        end
    end
    nothing
end

@inline function add_flux_differences_inner_loop!(du, u, balance_law, fvol, Nx, Ny, Pp1,
                                                    D, jacx, jacy, ::Val{false})
    dirx = Val{:x}()
    diry = Val{:y}()
    NxNy = Nx*Ny
    dims = (Pp1,Pp1,Nx*Ny)
    @inbounds for ixy in Base.OneTo(NxNy)
        for ny in Base.OneTo(Pp1), nx in Base.OneTo(Pp1)
            idx = sub2ind(dims, nx, ny, ixy)
            u_idx = u[idx]
            # compute x derivative
            # at first for different indices
            for k in (nx+1):Pp1
                idxk = sub2ind(dims, k, ny, ixy)
                f = fvol(u_idx, u[idxk], balance_law, dirx)
                du[idx]  -= 2*jacx*D[nx,k] * f
                du[idxk] -= 2*jacx*D[k,nx] * f
            end
            # then for the diagonal element, using the consistency of the flux
            du[idx] -= 2*jacx*D[nx,nx] * flux(u_idx, balance_law, dirx)

            # compute y derivative
            # at first for different indices
            for k in (ny+1):Pp1
                idxk = sub2ind(dims, nx, k, ixy)
                f = fvol(u_idx, u[idxk], balance_law, diry)
                du[idx]  -= 2*jacy*D[ny,k] * f
                du[idxk] -= 2*jacy*D[k,ny] * f
            end
            # then for the diagonal element, using the consistency of the flux
            du[idx] -= 2*jacy*D[ny,ny] * flux(u_idx, balance_law, diry)
        end
    end
    nothing
end


function add_numerical_fluxes!(du, u, semidisc::UniformPeriodicFluxDiffDisc2D)
    Nx = size(u, 3)
    Ny = size(u, 4)
    Pp1 = size(u, 1)
    @boundscheck @assert Pp1 == size(u, 2)

    @unpack balance_law, meshx, meshy, fnum, fluxes, usethreads = semidisc
    ω = semidisc.basis.weights
    @boundscheck @assert Pp1 == length(ω)

    jacx = 2 / semidisc.meshx.Δx
    jacy = 2 / semidisc.meshy.Δx

    add_numerical_fluxes_inner_loop1!(du, fluxes, u, balance_law, fnum, Nx, Ny, Pp1, ω, jacx, jacy, usethreads)
    add_numerical_fluxes_inner_loop2!(du, fluxes, u, balance_law, fnum, Nx, Ny, Pp1, ω, jacx, jacy, usethreads)
end

@inline function add_numerical_fluxes_inner_loop1!(du, fluxes, u, balance_law, fnum,
                                                    Nx, Ny, Pp1, ω, jacx, jacy, ::Val{true})
    dirx = Val{:x}()
    diry = Val{:y}()
    dims = (Pp1, Pp1, Nx, Ny)
    # calculate numerical fluxes
    @inbounds Threads.@threads for ixy in Base.OneTo(Nx*Ny)
        ix, iy = ind2sub((Nx,Ny), ixy)
        ixm1 = ix ==  1 ? Nx : ix-1
        iym1 = iy ==  1 ? Ny : iy-1

        # flux x
        for ny in 1:Pp1
            # flux x - left
            fluxes[ny,1,ix,iy] = fnum(u[end,ny,ixm1,iy], u[1,ny,ix,iy], balance_law, dirx)
        end

        # flux y
        for nx in 1:Pp1
            # flux y - bottom
            fluxes[nx,2,ix,iy] = fnum(u[nx,end,ix,iym1], u[nx,1,ix,iy], balance_law, diry)
        end
    end
    nothing
end
@inline function add_numerical_fluxes_inner_loop2!(du, fluxes, u, balance_law, fnum,
                                                    Nx, Ny, Pp1, ω, jacx, jacy, ::Val{true})
    dirx = Val{:x}()
    diry = Val{:y}()
    @inbounds i_ω1 = 1 / ω[1]
    @inbounds i_ωend = 1 / ω[end]
    dims = (Pp1, Pp1, Nx, Ny)
    # add numerical fluxes
    @inbounds Threads.@threads for ixy in Base.OneTo(Nx*Ny)
        ix, iy = ind2sub((Nx,Ny), ixy)
        ixp1 = ix == Nx ?  1 : ix+1
        iyp1 = iy == Ny ?  1 : iy+1

        # flux x
        for ny in 1:Pp1
            # flux x - left
            idx = sub2ind(dims, 1, ny, ix, iy)
            du[idx] += ( fluxes[ny,1,ix,iy] - flux(u[idx], balance_law, dirx) ) * jacx * i_ω1

            # flux x - right
            idx = sub2ind(dims, Pp1, ny, ix, iy)
            du[idx] -= ( fluxes[ny,1,ixp1,iy] - flux(u[idx], balance_law, dirx) ) * jacx * i_ωend
        end

        # flux y
        for nx in 1:Pp1
            # flux y - bottom
            idx = sub2ind(dims, nx, 1, ix, iy)
            du[idx] += ( fluxes[nx,2,ix,iy] - flux(u[idx], balance_law, diry) ) * jacy * i_ω1

            # flux y - top
            idx = sub2ind(dims, nx, Pp1, ix, iy)
            du[idx] -= ( fluxes[nx,2,ix,iyp1] - flux(u[idx], balance_law, diry) ) * jacy * i_ωend
        end
    end
    nothing
end

@inline function add_numerical_fluxes_inner_loop1!(du, fluxes, u, balance_law, fnum,
                                                    Nx, Ny, Pp1, ω, jacx, jacy, ::Val{false})
    dirx = Val{:x}()
    diry = Val{:y}()
    dims = (Pp1, Pp1, Nx, Ny)
    # calculate numerical fluxes
    @inbounds for ixy in Base.OneTo(Nx*Ny)
        ix, iy = ind2sub((Nx,Ny), ixy)
        ixm1 = ix ==  1 ? Nx : ix-1
        iym1 = iy ==  1 ? Ny : iy-1

        # flux x
        for ny in 1:Pp1
            # flux x - left
            fluxes[ny,1,ix,iy] = fnum(u[end,ny,ixm1,iy], u[1,ny,ix,iy], balance_law, dirx)
        end

        # flux y
        for nx in 1:Pp1
            # flux y - bottom
            fluxes[nx,2,ix,iy] = fnum(u[nx,end,ix,iym1], u[nx,1,ix,iy], balance_law, diry)
        end
    end
    nothing
end
@inline function add_numerical_fluxes_inner_loop2!(du, fluxes, u, balance_law, fnum,
                                                    Nx, Ny, Pp1, ω, jacx, jacy, ::Val{false})
    dirx = Val{:x}()
    diry = Val{:y}()
    @inbounds i_ω1 = 1 / ω[1]
    @inbounds i_ωend = 1 / ω[end]
    dims = (Pp1, Pp1, Nx, Ny)
    # add numerical fluxes
    @inbounds for ixy in Base.OneTo(Nx*Ny)
        ix, iy = ind2sub((Nx,Ny), ixy)
        ixp1 = ix == Nx ?  1 : ix+1
        iyp1 = iy == Ny ?  1 : iy+1

        # flux x
        for ny in 1:Pp1
            # flux x - left
            idx = sub2ind(dims, 1, ny, ix, iy)
            du[idx] += ( fluxes[ny,1,ix,iy] - flux(u[idx], balance_law, dirx) ) * jacx * i_ω1

            # flux x - right
            idx = sub2ind(dims, Pp1, ny, ix, iy)
            du[idx] -= ( fluxes[ny,1,ixp1,iy] - flux(u[idx], balance_law, dirx) ) * jacx * i_ωend
        end

        # flux y
        for nx in 1:Pp1
            # flux y - bottom
            idx = sub2ind(dims, nx, 1, ix, iy)
            du[idx] += ( fluxes[nx,2,ix,iy] - flux(u[idx], balance_law, diry) ) * jacy * i_ω1

            # flux y - top
            idx = sub2ind(dims, nx, Pp1, ix, iy)
            du[idx] -= ( fluxes[nx,2,ix,iyp1] - flux(u[idx], balance_law, diry) ) * jacy * i_ωend
        end
    end
    nothing
end


function semidiscretise(semidisc::UniformPeriodicFluxDiffDisc2D, u₀func, tspan)
    @unpack meshx, meshy, basis = semidisc
    u₀ = compute_coefficients(u₀func, meshx, meshy, basis)

    ODEProblem(semidisc, u₀, tspan)
end



################################################################################


struct UniformPeriodicFluxDiffDisc3D{BalanceLaw, T, Fvol, Fnum, UseThreads<:Union{Val{true},Val{false}}} <: AbstractSemidiscretisation
    balance_law::BalanceLaw
    meshx::UniformPeriodicMesh1D{T}
    meshy::UniformPeriodicMesh1D{T}
    meshz::UniformPeriodicMesh1D{T}
    basis::LobattoLegendre{T}
    fvol::Fvol
    fnum::Fnum

    usethreads::UseThreads
end

function UniformPeriodicFluxDiffDisc3D(balance_law, meshx, meshy, meshz, basis, fvol, fnum, usethreads::Bool=false)
    UniformPeriodicFluxDiffDisc3D(balance_law, meshx, meshy, meshz, basis, fvol, fnum, Val{usethreads}())
end


@noinline function (semidisc::UniformPeriodicFluxDiffDisc3D)(t, u, du)
  @boundscheck begin
    if size(u) != size(du)
      error("size(u) = $(size(u)) != $(size(du)) = size(du)")
    end
    @assert size(u,1) == size(u,2) == size(u,3) == length(semidisc.basis.nodes)
    @assert size(u,4) == numcells(semidisc.meshx)
    @assert size(u,5) == numcells(semidisc.meshy)
    @assert size(u,6) == numcells(semidisc.meshz)
    if eltype(u) != variables(semidisc.balance_law)
      error("eltype(u) == $(eltype(u)) != $(variables(semidisc.balance_law)) == variables(semidisc.balance_law)")
    end
  end

  fill!(du, zero(eltype(du)))
  add_flux_differences!(du, u, semidisc)
  add_numerical_fluxes!(du, u, semidisc)

  nothing
end


function add_flux_differences!(du, u, semidisc::UniformPeriodicFluxDiffDisc3D)
    Nx = size(u, 4)
    Ny = size(u, 5)
    Nz = size(u, 6)
    Pp1 = size(u, 1)
    @boundscheck @assert Pp1 == size(u, 2) == size(u, 3)

    @unpack balance_law, fvol, basis, usethreads = semidisc
    @unpack D = basis
    @boundscheck @assert Pp1 == size(D,1) == size(D,2)
    jacx = 2 / semidisc.meshx.Δx
    jacy = 2 / semidisc.meshy.Δx
    jacz = 2 / semidisc.meshz.Δx

    if Pp1 <= 1
        return nothing
    end

    add_flux_differences_inner_loop!(du, u, balance_law, fvol, Nx, Ny, Nz, Pp1, D, jacx, jacy, jacz, usethreads)
end

@inline function add_flux_differences_inner_loop!(du, u, balance_law, fvol, Nx, Ny, Nz,
                                                    Pp1, D, jacx, jacy, jacz, ::Val{true})
    dirx = Val{:x}()
    diry = Val{:y}()
    dirz = Val{:z}()
    @inbounds Threads.@threads for ixyz in Base.OneTo(Nx*Ny*Nz)
        ix, iy, iz = ind2sub((Nx,Ny,Nz), ixyz)
        for nz in Base.OneTo(Pp1), ny in Base.OneTo(Pp1), nx in Base.OneTo(Pp1)
            # compute x derivative
            for k in Base.OneTo(Pp1)
                du[nx,ny,nz,ix,iy,iz] -= 2*jacx*D[nx,k] * fvol( u[nx,ny,nz,ix,iy,iz],
                                                                u[k ,ny,nz,ix,iy,iz],
                                                                balance_law,
                                                                dirx)
            end

            # compute y derivative
            for k in Base.OneTo(Pp1)
                du[nx,ny,nz,ix,iy,iz] -= 2*jacx*D[ny,k] * fvol( u[nx,ny,nz,ix,iy,iz],
                                                                u[nx,k ,nz,ix,iy,iz],
                                                                balance_law,
                                                                diry)
            end

            # compute z derivative
            for k in Base.OneTo(Pp1)
                du[nx,ny,nz,ix,iy,iz] -= 2*jacx*D[nz,k] * fvol( u[nx,ny,nz,ix,iy,iz],
                                                                u[nx,ny,k ,ix,iy,iz],
                                                                balance_law,
                                                                dirz)
            end
        end
    end
    nothing
end

@inline function add_flux_differences_inner_loop!(du, u, balance_law, fvol, Nx, Ny, Nz,
                                                    Pp1, D, jacx, jacy, jacz, ::Val{false})
    dirx = Val{:x}()
    diry = Val{:y}()
    dirz = Val{:z}()
    @inbounds for ixyz in Base.OneTo(Nx*Ny*Nz)
        ix, iy, iz = ind2sub((Nx,Ny,Nz), ixyz)
        for nz in Base.OneTo(Pp1), ny in Base.OneTo(Pp1), nx in Base.OneTo(Pp1)
            # compute x derivative
            for k in Base.OneTo(Pp1)
                du[nx,ny,nz,ix,iy,iz] -= 2*jacx*D[nx,k] * fvol( u[nx,ny,nz,ix,iy,iz],
                                                                u[k ,ny,nz,ix,iy,iz],
                                                                balance_law,
                                                                dirx)
            end

            # compute y derivative
            for k in Base.OneTo(Pp1)
                du[nx,ny,nz,ix,iy,iz] -= 2*jacx*D[ny,k] * fvol( u[nx,ny,nz,ix,iy,iz],
                                                                u[nx,k ,nz,ix,iy,iz],
                                                                balance_law,
                                                                diry)
            end

            # compute z derivative
            for k in Base.OneTo(Pp1)
                du[nx,ny,nz,ix,iy,iz] -= 2*jacx*D[nz,k] * fvol( u[nx,ny,nz,ix,iy,iz],
                                                                u[nx,ny,k ,ix,iy,iz],
                                                                balance_law,
                                                                dirz)
            end
        end
    end
    nothing
end


function add_numerical_fluxes!(du, u, semidisc::UniformPeriodicFluxDiffDisc3D)
    Nx = size(u, 4)
    Ny = size(u, 5)
    Nz = size(u, 6)
    Pp1 = size(u, 1)
    @boundscheck @assert Pp1 == size(u, 2) == size(u, 3)

    @unpack balance_law, meshx, meshy, fnum, usethreads = semidisc
    ω = semidisc.basis.weights
    @boundscheck @assert Pp1 == length(ω)

    jacx = 2 / semidisc.meshx.Δx
    jacy = 2 / semidisc.meshy.Δx
    jacz = 2 / semidisc.meshz.Δx

    add_numerical_fluxes_inner_loop!(du, u, balance_law, fnum, Nx, Ny, Nz, Pp1, ω, jacx, jacy, jacz, usethreads)
end

@inline function add_numerical_fluxes_inner_loop!(du, u, balance_law, fnum, Nx, Ny, Nz,
                                                    Pp1, ω, jacx, jacy, jacz, ::Val{true})
    dirx = Val{:x}()
    diry = Val{:y}()
    dirz = Val{:y}()
    @inbounds Threads.@threads for ixyz in Base.OneTo(Nx*Ny*Nz)
        ix, iy, iz = ind2sub((Nx,Ny,Nz), ixyz)
        ixm1 = ix ==  1 ? Nx : ix-1
        ixp1 = ix == Nx ?  1 : ix+1
        iym1 = iy ==  1 ? Ny : iy-1
        iyp1 = iy == Ny ?  1 : iy+1
        izm1 = iz ==  1 ? Nz : iz-1
        izp1 = iz == Nz ?  1 : iz+1

        # flux x
        for nz in Base.OneTo(Pp1), ny in Base.OneTo(Pp1)
            # flux x - left
            du[1,ny,nz,ix,iy,iz] += (
                fnum(u[end,ny,nz,ixm1,iy,iz], u[1,ny,nz,ix,iy,iz], balance_law, dirx)
                - flux(u[1,ny,nz,ix,iy,iz], balance_law, dirx) ) * jacx / ω[1]

            # flux x - right
            du[end,ny,nz,ix,iy,iz] -= (
                fnum(u[end,ny,nz,ix,iy,iz], u[1,ny,nz,ixp1,iy,iz], balance_law, dirx)
                - flux(u[end,ny,nz,ix,iy,iz], balance_law, dirx) ) * jacx / ω[end]
        end

        # flux y
        for nz in Base.OneTo(Pp1), nx in Base.OneTo(Pp1)
            # flux y - bottom
            du[nx,1,nz,ix,iy,iz] += (
                fnum(u[nx,end,nz,ix,iym1,iz], u[nx,1,nz,ix,iy,iz], balance_law, diry)
                - flux(u[nx,1,nz,ix,iy,iz], balance_law, diry) ) * jacy / ω[1]

            # flux y - top
            du[nx,end,nz,ix,iy,iz] -= (
                fnum(u[nx,end,nz,ix,iy,iz], u[nx,1,nz,ix,iyp1,iz], balance_law, diry)
                - flux(u[nx,end,nz,ix,iy,iz], balance_law, diry) ) * jacy / ω[end]
        end

        # flux z
        for ny in Base.OneTo(Pp1), nx in Base.OneTo(Pp1)
            # flux z - back
            du[nx,ny,1,ix,iy,iz] += (
                fnum(u[nx,ny,end,ix,iy,izm1], u[nx,ny,1,ix,iy,iz], balance_law, dirz)
                - flux(u[nx,ny,1,ix,iy,iz], balance_law, dirz) ) * jacz / ω[1]

            # flux z - front
            du[nx,ny,end,ix,iy,iz] -= (
                fnum(u[nx,ny,end,ix,iy,iz], u[nx,ny,1,ix,iy,izp1], balance_law, dirz)
                - flux(u[nx,ny,end,ix,iy,iz], balance_law, dirz) ) * jacz / ω[end]
        end
    end
    nothing
end

@inline function add_numerical_fluxes_inner_loop!(du, u, balance_law, fnum, Nx, Ny, Nz,
                                                    Pp1, ω, jacx, jacy, jacz, ::Val{false})
    dirx = Val{:x}()
    diry = Val{:y}()
    dirz = Val{:y}()
    @inbounds for ixyz in Base.OneTo(Nx*Ny*Nz)
        ix, iy, iz = ind2sub((Nx,Ny,Nz), ixyz)
        ixm1 = ix ==  1 ? Nx : ix-1
        ixp1 = ix == Nx ?  1 : ix+1
        iym1 = iy ==  1 ? Ny : iy-1
        iyp1 = iy == Ny ?  1 : iy+1
        izm1 = iz ==  1 ? Nz : iz-1
        izp1 = iz == Nz ?  1 : iz+1

        # flux x
        for nz in Base.OneTo(Pp1), ny in Base.OneTo(Pp1)
            # flux x - left
            du[1,ny,nz,ix,iy,iz] += (
                fnum(u[end,ny,nz,ixm1,iy,iz], u[1,ny,nz,ix,iy,iz], balance_law, dirx)
                - flux(u[1,ny,nz,ix,iy,iz], balance_law, dirx) ) * jacx / ω[1]

            # flux x - right
            du[end,ny,nz,ix,iy,iz] -= (
                fnum(u[end,ny,nz,ix,iy,iz], u[1,ny,nz,ixp1,iy,iz], balance_law, dirx)
                - flux(u[end,ny,nz,ix,iy,iz], balance_law, dirx) ) * jacx / ω[end]
        end

        # flux y
        for nz in Base.OneTo(Pp1), nx in Base.OneTo(Pp1)
            # flux y - bottom
            du[nx,1,nz,ix,iy,iz] += (
                fnum(u[nx,end,nz,ix,iym1,iz], u[nx,1,nz,ix,iy,iz], balance_law, diry)
                - flux(u[nx,1,nz,ix,iy,iz], balance_law, diry) ) * jacy / ω[1]

            # flux y - top
            du[nx,end,nz,ix,iy,iz] -= (
                fnum(u[nx,end,nz,ix,iy,iz], u[nx,1,nz,ix,iyp1,iz], balance_law, diry)
                - flux(u[nx,end,nz,ix,iy,iz], balance_law, diry) ) * jacy / ω[end]
        end

        # flux z
        for ny in Base.OneTo(Pp1), nx in Base.OneTo(Pp1)
            # flux z - back
            du[nx,ny,1,ix,iy,iz] += (
                fnum(u[nx,ny,end,ix,iy,izm1], u[nx,ny,1,ix,iy,iz], balance_law, dirz)
                - flux(u[nx,ny,1,ix,iy,iz], balance_law, dirz) ) * jacz / ω[1]

            # flux z - front
            du[nx,ny,end,ix,iy,iz] -= (
                fnum(u[nx,ny,end,ix,iy,iz], u[nx,ny,1,ix,iy,izp1], balance_law, dirz)
                - flux(u[nx,ny,end,ix,iy,iz], balance_law, dirz) ) * jacz / ω[end]
        end
    end
    nothing
end


function semidiscretise(semidisc::UniformPeriodicFluxDiffDisc3D, u₀func, tspan)
    @unpack meshx, meshy, meshz, basis = semidisc
    u₀ = compute_coefficients(u₀func, meshx, meshy, meshz, basis)

    ODEProblem(semidisc, u₀, tspan)
end
