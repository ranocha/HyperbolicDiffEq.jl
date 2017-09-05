
struct UniformPeriodicFluxDiffDisc2D{BalanceLaw, T, Fvol, Fnum, UseThreads<:Union{Val{true},Val{false}}} <: AbstractSemidiscretisation
    balance_law::BalanceLaw
    meshx::UniformPeriodicMesh1D{T}
    meshy::UniformPeriodicMesh1D{T}
    basis::LobattoLegendre{T}
    fvol::Fvol
    fnum::Fnum

    usethreads::UseThreads
end

function UniformPeriodicFluxDiffDisc2D(balance_law, meshx, meshy, basis, fvol, fnum, usethreads::Bool=false)
    UniformPeriodicFluxDiffDisc2D(balance_law, meshx, meshy, basis, fvol, fnum, Val{usethreads}())
end


function (semidisc::UniformPeriodicFluxDiffDisc2D)(t, u, du)
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

@inline function add_flux_differences_inner_loop!(du, u, balance_law, fvol, Nx, Ny, Pp1, D, jacx, jacy, ::Val{true})
    dirx = Val{:x}()
    diry = Val{:y}()
    @inbounds Threads.@threads for ixy in Base.OneTo(Nx*Ny)
        iy, ix = divrem(ixy-1, Nx) .+ 1
        for ny in Base.OneTo(Pp1), nx in Base.OneTo(Pp1)
            # compute x derivative
            for k in Base.OneTo(Pp1)
                du[nx,ny,ix,iy] -= 2*jacx*D[nx,k] * fvol(u[nx,ny,ix,iy],
                                                         u[k ,ny,ix,iy],
                                                         balance_law,
                                                         dirx)
            end

            # compute y derivative
            for k in Base.OneTo(Pp1)
                du[nx,ny,ix,iy] -= 2*jacy*D[ny,k] * fvol(u[nx,ny,ix,iy],
                                                         u[nx,k ,ix,iy],
                                                         balance_law,
                                                         diry)
            end
        end
    end
    nothing
end

@inline function add_flux_differences_inner_loop!(du, u, balance_law, fvol, Nx, Ny, Pp1,
                                            D, jacx, jacy, ::Val{false})
    dirx = Val{:x}()
    diry = Val{:y}()
    @inbounds for iy in Base.OneTo(Ny), ix in Base.OneTo(Nx)
        for ny in Base.OneTo(Pp1), nx in Base.OneTo(Pp1)
            # compute x derivative
            for k in Base.OneTo(Pp1)
                du[nx,ny,ix,iy] -= 2*jacx*D[nx,k] * fvol(u[nx,ny,ix,iy],
                                                         u[k ,ny,ix,iy],
                                                         balance_law,
                                                         dirx)
            end

            # compute y derivative
            for k in Base.OneTo(Pp1)
                du[nx,ny,ix,iy] -= 2*jacy*D[ny,k] * fvol(u[nx,ny,ix,iy],
                                                         u[nx,k ,ix,iy],
                                                         balance_law,
                                                         diry)
            end
        end
    end
    nothing
end


function add_numerical_fluxes!(du, u, semidisc::UniformPeriodicFluxDiffDisc2D)
    Nx = size(u, 3)
    Ny = size(u, 4)
    Pp1 = size(u, 1)
    @boundscheck @assert Pp1 == size(u, 2)

    @unpack balance_law, meshx, meshy, fnum, usethreads = semidisc
    ω = semidisc.basis.weights
    @boundscheck @assert Pp1 == length(ω)

    jacx = 2 / semidisc.meshx.Δx
    jacy = 2 / semidisc.meshy.Δx

    add_numerical_fluxes_inner_loop!(du, u, balance_law, fnum, Nx, Ny, Pp1, ω, jacx, jacy, usethreads)
end

@inline function add_numerical_fluxes_inner_loop!(du, u, balance_law, fnum, Nx, Ny, Pp1,
                                            ω, jacx, jacy, ::Val{true})
    dirx = Val{:x}()
    diry = Val{:y}()
    @inbounds Threads.@threads for ixy in Base.OneTo(Nx*Ny)
        iy, ix = divrem(ixy-1, Nx) .+ 1
        ixm1 = ix ==  1 ? Nx : ix-1
        ixp1 = ix == Nx ?  1 : ix+1
        iym1 = iy ==  1 ? Ny : iy-1
        iyp1 = iy == Ny ?  1 : iy+1

        # flux x
        for ny in 1:Pp1
            # flux x - left
            du[1,ny,ix,iy] += (
                fnum(u[end,ny,ixm1,iy], u[1,ny,ix,iy], balance_law, dirx)
                - flux(u[1,ny,ix,iy], balance_law, dirx) ) * jacx / ω[1]

            # flux x - right
            du[end,ny,ix,iy] -= (
                fnum(u[end,ny,ix,iy], u[1,ny,ixp1,iy], balance_law, dirx)
                - flux(u[end,ny,ix,iy], balance_law, dirx) ) * jacx / ω[end]
        end

        # flux y
        for nx in 1:Pp1
            # flux y - bottom
            du[nx,1,ix,iy] += (
                fnum(u[nx,end,ix,iym1], u[nx,1,ix,iy], balance_law, diry)
                - flux(u[nx,1,ix,iy], balance_law, diry) ) * jacy / ω[1]

            # flux y - top
            du[nx,end,ix,iy] -= (
                fnum(u[nx,end,ix,iy], u[nx,1,ix,iyp1], balance_law, diry)
                - flux(u[nx,end,ix,iy], balance_law, diry) ) * jacy / ω[end]
        end
    end
    nothing
end

@inline function add_numerical_fluxes_inner_loop!(du, u, balance_law, fnum, Nx, Ny, Pp1,
                                            ω, jacx, jacy, ::Val{false})
    dirx = Val{:x}()
    diry = Val{:y}()
    @inbounds for iy in Base.OneTo(Ny), ix in Base.OneTo(Nx)
        ixm1 = ix ==  1 ? Nx : ix-1
        ixp1 = ix == Nx ?  1 : ix+1
        iym1 = iy ==  1 ? Ny : iy-1
        iyp1 = iy == Ny ?  1 : iy+1

        # flux x
        for ny in 1:Pp1
            # flux x - left
            du[1,ny,ix,iy] += (
                fnum(u[end,ny,ixm1,iy], u[1,ny,ix,iy], balance_law, dirx)
                - flux(u[1,ny,ix,iy], balance_law, dirx) ) * jacx / ω[1]

            # flux x - right
            du[end,ny,ix,iy] -= (
                fnum(u[end,ny,ix,iy], u[1,ny,ixp1,iy], balance_law, dirx)
                - flux(u[end,ny,ix,iy], balance_law, dirx) ) * jacx / ω[end]
        end

        # flux y
        for nx in 1:Pp1
            # flux y - bottom
            du[nx,1,ix,iy] += (
                fnum(u[nx,end,ix,iym1], u[nx,1,ix,iy], balance_law, diry)
                - flux(u[nx,1,ix,iy], balance_law, diry) ) * jacy / ω[1]

            # flux y - top
            du[nx,end,ix,iy] -= (
                fnum(u[nx,end,ix,iy], u[nx,1,ix,iyp1], balance_law, diry)
                - flux(u[nx,end,ix,iy], balance_law, diry) ) * jacy / ω[end]
        end
    end
    nothing
end


function semidiscretise(semidisc::UniformPeriodicFluxDiffDisc2D, u₀func, tspan)
    @unpack meshx, meshy, basis = semidisc
    u₀ = compute_coefficients(u₀func, meshx, meshy, basis)

    ODEProblem(semidisc, u₀, tspan)
end
