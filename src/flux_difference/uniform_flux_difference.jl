
struct UniformPeriodicFluxDiffDisc2D{BalanceLaw, T, Fvol, Fnum, Fluxes, UseThreads} <: AbstractSemidiscretisation
    balance_law::BalanceLaw
    meshx::UniformPeriodicMesh1D{T}
    meshy::UniformPeriodicMesh1D{T}
    basis::LobattoLegendre{T}
    fvol::Fvol
    fnum::Fnum

    usethreads::UseThreads
end

function UniformPeriodicFluxDiffDisc2D(balance_law::AbstractBalanceLaw{2}, meshx, meshy,
                                        basis, fvol, fnum, usethreads::Bool=false)
    UniformPeriodicFluxDiffDisc2D(balance_law, meshx, meshy, basis, fvol, fnum, fluxes,
                                    Val{usethreads}())
end


function (semidisc::UniformPeriodicFluxDiffDisc2D)(t, u, du)
  @boundscheck begin
    if size(u) != size(du)
      error("size(u) = $(size(u)) != $(size(du)) = size(du)")
    end
    size(u,1) != numcells(semidisc.meshx) && error("size(u,1) != numcells(semidisc.meshx)")
    size(u,2) != numcells(semidisc.meshy) && error("size(u,2) != numcells(semidisc.meshy)")
    if eltype(u) != variables(semidisc.balance_law)
      error("eltype(u) == $(eltype(u)) != $(variables(semidisc.balance_law)) == variables(semidisc.balance_law)")
    end
  end

  du .= 0
  add_flux_differences!(du, u, semidisc, semidisc.usethreads)
  add_numerical_fluxes!(du, u, semidisc, semidisc.usethreads)

  nothing
end


function add_flux_differences!(du, u, semidisc::UniformPeriodicFluxDiffDisc2D, usethreads)
    Nx = size(u, 3)
    Ny = size(u, 4)
    Pp1 = size(u, 1)
    assert(Pp1 == size(u, 2))

    @unpack balance_law, fvol, basis = semidisc
    @unpack D = basis
    jacx = 2 / semidisc.meshx.Δx
    jacy = 2 / semidisc.meshy.Δx

    for iy in 1:Ny
        for ix in 1:Nx
            for ny in 1:Pp1
                for nx in 1:Pp1
                    # compute x derivative
                    if Pp1 > 1
                        for k in 1:Pp1
                            du[nx,ny,ix,iy] -= 2*jacx*D[nx,k] * fvol(u[nx,ny,ix,iy],
                                                                     u[k ,ny,ix,iy],
                                                                     balance_law,
                                                                     Val{:x}())
                        end
                    end
                    # compute y derivative
                    if Pp1 > 1
                        for k in 1:Pp1
                            du[nx,ny,ix,iy] -= 2*jacy*D[ny,k] * fvol(u[nx,ny,ix,iy],
                                                                     u[nx,k ,ix,iy],
                                                                     balance_law,
                                                                     Val{:y}())
                        end
                    end
                end
            end
        end
    end

    nothing
end


function add_numerical_fluxes!(du, u, semidisc::UniformPeriodicFluxDiffDisc2D, usethreads)
    Nx = size(u, 3)
    Ny = size(u, 4)
    Pp1 = size(u, 1)
    assert(Pp1 == size(u, 2))

    @unpack balance_law, meshx, meshy, fnum = semidisc
    ω = semidisc.basis.weights

    jacx = 2 / semidisc.meshx.Δx
    jacy = 2 / semidisc.meshy.Δx

    for iy in 1:Ny
        for ix in 1:Nx
            ixm1 = ix ==  1 ? Nx : ix-1
            ixp1 = ix == Nx ?  1 : ix+1
            iym1 = iy ==  1 ? Ny : iy-1
            iyp1 = iy == Ny ?  1 : iy+1

            # flux x
            direction = Val{:x}()
            for ny in 1:Pp1
                # flux x - left
                du[1,ny,ix,iy] += (
                    fnum(u[end,ny,ixm1,iy], u[1,ny,ix,iy], balance_law, direction)
                    - flux(u[1,ny,ix,iy], balance_law, direction) ) * jacx / ω[1]

                # flux x - right
                du[end,ny,ix,iy] += (
                    fnum(u[end,ny,ix,iy], u[1,ny,ix+1,iy], balance_law, direction)
                    - flux(u[end,ny,ix,iy], balance_law, direction) ) * jacx / ω[end]
            end

            # flux y
            for ny in 1:Pp1
                # flux y - bottom
                du[nx,1,ix,iy] += (
                    fnum(u[nx,end,ix,iym1], u[nx,1,ix,iy], balance_law, direction)
                    - flux(u[1,ny,ix,iy], balance_law, direction) ) * jacy / ω[1]

                # flux y - top
                du[nx,end,ix,iy] += (
                    fnum(u[nx,end,ix,iy], u[nx,1,ix,iyp1], balance_law, direction)
                    - flux(u[end,ny,ix,iy], balance_law, direction) ) * jacy / ω[end]
            end
        end
    end

    nothing
end


function semidiscretise(semidisc::UniformPeriodicFluxDiffDisc2D, u₀func, tspan)
    @unpack meshx, meshy, basis = semidisc
    u₀ = compute_coefficients(u₀func, meshx, meshy, basis)

    ODEProblem(semidisc, u₀, tspan)
end
