
"""
    Euler{T<:Real,Dim}

The compᵣessible Euler equations in `Dim` dimensions using `T` as scalar type.
"""
struct Euler{T<:Real,Dim} <: AbstractBalanceLaw{Dim}
    γ::T
end

function Euler(γ=1., Dim=1)
  Euler{typeof(γ), Dim}(γ)
end

function show{T,Dim}(io::IO, model::Euler{T,Dim})
  pᵣint(io, "Compᵣessible Euler equations with γ=", model.γ, " {T=", T, ", Dim=", Dim, "}")
end
"""
Conserved variables (ϱ, ϱvx, ϱvy, ϱe) of the Euler equations in two space dimensions.
"""
struct EulerVar2D{T} <: FieldVector{4,T}
  ϱ  ::T
  ϱvx::T
  ϱvy::T
  ϱe ::T
end

function (::Type{EulerVar2D{T}}){T}(val::Real)
  EulerVar2D{T}(val, val, val, val)
end

function similar_type{T}(::EulerVar2D{T})
  EulerVar2D{T}
end

@inline variables{T}(model::Euler{T,2}) = EulerVar2D{T}

@inline function primitive_variables(u::EulerVar2D, model::Euler)
    @unpack ϱ, ϱvx, ϱvy, ϱe = u
    @unpack γ = model
    if ϱ ≈ 0
        vx = zero(ϱ)
        vy = zero(ϱ)
    else
        vx = ϱvx / ϱ
        vy = ϱvy / ϱ
    end
    p = (γ-1) * (ϱe - ϱvx*vx/2 - ϱvy*vy/2)

    ϱ, vx, vy, p
end

@inline function conserved_variables(ϱ, vx, vy, p, model::Euler)
    @unpack γ = model

    variables(model)(ϱ, ϱ*vx, ϱ*vy, (ϱ*vx^2+ϱ*vy^2)/2+(γ-1)*p)
end

@inline function kinetic_energy(u::EulerVar2D, model::Euler)
    ϱ, vx, vy, p = primitive_variables(u, model)
    !
    ϱ * (vx^2 + vy^2) / 2
end

@inline function flux{T}(u::EulerVar2D{T}, model::Euler, dir::Val{:x})
    @unpack ϱvx, ϱvy, ϱe = u
    ϱ, vx, vy, p = primitive_variables(u, model)

    SVector{4,T}(ϱvx, ϱvx*vx + p, ϱvx*vy, (ϱe+p)*vx)
end

@inline function flux{T}(u::EulerVar2D{T}, model::Euler, dir::Val{:y})
    @unpack ϱvx, ϱvy, ϱe = u
    ϱ, vx, vy, p = primitive_variables(u, model)

    SVector{4,T}(ϱvy, ϱvy*vx, ϱvy*vy + p, (ϱe+p)*vy)
end

@inline function max_abs_speed(u::EulerVar2D, model::Euler)
    @unpack γ = model
    ϱ, vx, vy, p = primitive_variables(u, model)

    max(abs(vx),abs(vy)) + sqrt(γ * p / ϱ)
end

@inline function max_abs_speed(u::EulerVar2D, model::Euler, dir::Val{:x})
    @unpack γ = model
    ϱ, vx, vy, p = primitive_variables(u, model)

    abs(vx) + sqrt(γ * p / ϱ)
end

@inline function max_abs_speed(u::EulerVar2D, model::Euler, dir::Val{:y})
    @unpack γ = model
    ϱ, vx, vy, p = primitive_variables(u, model)

    abs(vy) + sqrt(γ * p / ϱ)
end


@recipe function f{X,Y,T}(xyumodel::Tuple{X,Y,AbstractArray{EulerVar2D{T}},Euler})
    x, y, u, model = xyumodel

    ϱ  = mappedarray(u->primitive_variables(u,model)[1], u)
    vx = mappedarray(u->primitive_variables(u,model)[2], u)
    vy = mappedarray(u->primitive_variables(u,model)[3], u)
    p  = mappedarray(u->primitive_variables(u,model)[4], u)

    layout --> (2,2)
    legend --> false

    @series begin
        subplot := 1
        title  --> L"\varrho"
        x, y, ϱ
    end

    @series begin
        subplot := 2
        title  --> L"p"
        x, y, p
    end

    @series begin
        subplot := 3
        title  --> L"v_x"
        x, y, vx
    end

    @series begin
        subplot := 4
        title  --> L"v_y"
        x, y, vy
    end
end


################################################################################


Base.@pure function (fnum::SuliciuFlux)(uₗ::EulerVar2D{T}, uᵣ::EulerVar2D{T}, model::Euler{T,2}, dir::Val{:x}) where T
    @unpack γ = model
    ϱeₗ = uₗ.ϱe
    ϱeᵣ = uᵣ.ϱe
    ϱₗ, vxₗ, vyₗ, pₗ = primitive_variables(uₗ, model)
    cₗ = sqrt(γ * pₗ / ϱₗ)
    ɛₗ = ϱeₗ/ϱₗ - vxₗ*vxₗ/2 - vyₗ*vyₗ/2
    ϱᵣ, vxᵣ, vyᵣ, pᵣ = primitive_variables(uᵣ, model)
    cᵣ = sqrt(γ * pᵣ / ϱᵣ)
    ɛᵣ = ϱeᵣ/ϱᵣ - vxᵣ*vxᵣ/2 - vyᵣ*vyᵣ/2
    α = (γ+1)/2

    # compute speeds
    cₗ_ϱₗ = zero(T)
    cᵣ_ϱᵣ = zero(T)
    if pₗ <= pᵣ && 0 < pᵣ
      cₗ_ϱₗ = cₗ + α*max( zero(T), (pᵣ-pₗ)/(ϱᵣ*cᵣ) + vxₗ - vxᵣ )
      cᵣ_ϱᵣ = cᵣ + α*max( zero(T), (pₗ-pᵣ)/(cₗ_ϱₗ*ϱₗ) + vxₗ - vxᵣ )
    elseif pᵣ <= pₗ && 0 < pₗ
      cᵣ_ϱᵣ = cᵣ + α*max( zero(T), (pₗ-pᵣ)/(ϱₗ*cₗ) + vxₗ - vxᵣ )
      cₗ_ϱₗ = cₗ + α*max( zero(T), (pᵣ-pₗ)/(cᵣ_ϱᵣ*ϱᵣ) + vxₗ - vxᵣ )
    end

    # compute intermediate values
    cₗ = cₗ_ϱₗ*ϱₗ
    cᵣ = cᵣ_ϱᵣ*ϱᵣ
    vxs = ( cₗ*vxₗ + cᵣ*vxᵣ + pₗ - pᵣ ) / ( cₗ+cᵣ )
    vxs = ifelse(isnan(vxs), zero(T), vxs)
    pₗs = pᵣs = ( cᵣ*pₗ + cₗ*pᵣ - cₗ*cᵣ*(vxᵣ-vxₗ) ) / (cₗ+cᵣ)
    pₗs = pᵣs = ifelse(isnan(pₗs), zero(T), pₗs)
    ϱₗs = 1 / ( 1/ϱₗ + (cᵣ*(vxᵣ-vxₗ) + pₗ - pᵣ ) / (cₗ*(cₗ+cᵣ)) )
    ϱₗs = ifelse(isnan(ϱₗs), zero(T), ϱₗs)
    ϱᵣs = 1 / ( 1/ϱᵣ + (cₗ*(vxᵣ-vxₗ) + pᵣ - pₗ ) / (cᵣ*(cₗ+cᵣ)) )
    ϱᵣs = ifelse(isnan(ϱᵣs), zero(T), ϱᵣs)
    ɛₗs = ɛₗ + (pₗs*pₗs - pₗ*pₗ)/(2*cₗ*cₗ)
    ɛᵣs = ɛᵣ + (pᵣs*pᵣs - pᵣ*pᵣ)/(2*cᵣ*cᵣ)

    # compute fluxes
    fϱ   = zero(T)
    fϱvx = zero(T)
    fϱvy = zero(T)
    fϱe  = zero(T)
    if 0 <= vxₗ-cₗ_ϱₗ
        fϱ   = ϱₗ*vxₗ
        fϱvx = ϱₗ*vxₗ*vxₗ + pₗ
        fϱvy = ϱₗ*vxₗ*vyₗ
        fϱe  = (ϱₗ*vxₗ*vxₗ/2 + ϱₗ*vyₗ*vyₗ/2 + ϱₗ*ɛₗ + pₗ)*vxₗ
    elseif 0 <= vxs
        fϱ   = ϱₗs*vxs
        fϱvx = ϱₗs*vxs*vxs + pₗs
        fϱvy = ϱₗ*vxs*vyₗ
        fϱe  = (ϱₗs*vxs*vxs/2 + ϱₗs*vyₗ*vyₗ/2 + ϱₗs*ɛₗs + pₗs)*vxs
    elseif 0 <= vxᵣ+cᵣ_ϱᵣ
        fϱ   = ϱᵣs*vxs
        fϱvx = ϱᵣs*vxs*vxs + pᵣs
        fϱvy = ϱₗ*vxs*vyᵣ
        fϱe  = (ϱᵣs*vxs*vxs/2 + ϱᵣs*vyᵣ*vyᵣ/2 + ϱᵣs*ɛᵣs + pᵣs)*vxs
    else
        fϱ   = ϱᵣ*vxᵣ
        fϱvx = ϱᵣ*vxᵣ*vxᵣ + pᵣ
        fϱvy = ϱₗ*vxᵣ*vyᵣ
        fϱe  = (ϱᵣ*vxᵣ*vxᵣ/2 + ϱᵣ*vyᵣ*vyᵣ/2 + ϱᵣ*ɛᵣ + pᵣ)*vxᵣ
    end

    SVector(fϱ, fϱvx, fϱvy, fϱe)
end

Base.@pure function (fnum::SuliciuFlux)(uₗ::EulerVar2D{T}, uᵣ::EulerVar2D{T}, model::Euler{T,2}, dir::Val{:y}) where T
    @unpack γ = model
    ϱeₗ = uₗ.ϱe
    ϱeᵣ = uᵣ.ϱe
    ϱₗ, vxₗ, vyₗ, pₗ = primitive_variables(uₗ, model)
    cₗ = sqrt(γ * pₗ / ϱₗ)
    ɛₗ = ϱeₗ/ϱₗ - vxₗ*vxₗ/2 - vyₗ*vyₗ/2
    ϱᵣ, vxᵣ, vyᵣ, pᵣ = primitive_variables(uᵣ, model)
    cᵣ = sqrt(γ * pᵣ / ϱᵣ)
    ɛᵣ = ϱeᵣ/ϱᵣ - vxᵣ*vxᵣ/2 - vyᵣ*vyᵣ/2
    α = (γ+1)/2

    # compute speeds
    cₗ_ϱₗ = zero(T)
    cᵣ_ϱᵣ = zero(T)
    if pₗ <= pᵣ && 0 < pᵣ
      cₗ_ϱₗ = cₗ + α*max( zero(T), (pᵣ-pₗ)/(ϱᵣ*cᵣ) + vyₗ - vyᵣ )
      cᵣ_ϱᵣ = cᵣ + α*max( zero(T), (pₗ-pᵣ)/(cₗ_ϱₗ*ϱₗ) + vyₗ - vyᵣ )
    elseif pᵣ <= pₗ && 0 < pₗ
      cᵣ_ϱᵣ = cᵣ + α*max( zero(T), (pₗ-pᵣ)/(ϱₗ*cₗ) + vyₗ - vyᵣ )
      cₗ_ϱₗ = cₗ + α*max( zero(T), (pᵣ-pₗ)/(cᵣ_ϱᵣ*ϱᵣ) + vyₗ - vyᵣ )
    end

    # compute intermediate values
    cₗ = cₗ_ϱₗ*ϱₗ
    cᵣ = cᵣ_ϱᵣ*ϱᵣ
    vys = ( cₗ*vyₗ + cᵣ*vyᵣ + pₗ - pᵣ ) / ( cₗ+cᵣ )
    vys = ifelse(isnan(vys), zero(T), vys)
    pₗs = pᵣs = ( cᵣ*pₗ + cₗ*pᵣ - cₗ*cᵣ*(vyᵣ-vyₗ) ) / (cₗ+cᵣ)
    pₗs = pᵣs = ifelse(isnan(pₗs), zero(T), pₗs)
    ϱₗs = 1 / ( 1/ϱₗ + (cᵣ*(vyᵣ-vyₗ) + pₗ - pᵣ ) / (cₗ*(cₗ+cᵣ)) )
    ϱₗs = ifelse(isnan(ϱₗs), zero(T), ϱₗs)
    ϱᵣs = 1 / ( 1/ϱᵣ + (cₗ*(vyᵣ-vyₗ) + pᵣ - pₗ ) / (cᵣ*(cₗ+cᵣ)) )
    ϱᵣs = ifelse(isnan(ϱᵣs), zero(T), ϱᵣs)
    ɛₗs = ɛₗ + (pₗs*pₗs - pₗ*pₗ)/(2*cₗ*cₗ)
    ɛᵣs = ɛᵣ + (pᵣs*pᵣs - pᵣ*pᵣ)/(2*cᵣ*cᵣ)

    # compute fluxes
    fϱ   = zero(T)
    fϱvx = zero(T)
    fϱvy = zero(T)
    fϱe  = zero(T)
    if 0 <= vyₗ-cₗ_ϱₗ
        fϱ   = ϱₗ*vyₗ
        fϱvx = ϱₗ*vxₗ*vyₗ
        fϱvy = ϱₗ*vyₗ*vyₗ + pₗ
        fϱe  = (ϱₗ*vxₗ*vxₗ/2 + ϱₗ*vyₗ*vyₗ/2 + ϱₗ*ɛₗ + pₗ)*vyₗ
    elseif 0 <= vys
        fϱ   = ϱₗs*vys
        fϱvx = ϱₗs*vxₗ*vys
        fϱvy = ϱₗ*vys*vys + pₗs
        fϱe  = (ϱₗs*vxₗ*vxₗ/2 + ϱₗs*vys*vys/2 + ϱₗs*ɛₗs + pₗs)*vys
    elseif 0 <= vyᵣ+cᵣ_ϱᵣ
        fϱ   = ϱᵣs*vys
        fϱvx = ϱᵣs*vxᵣ*vys
        fϱvy = ϱₗ*vys*vys + pᵣs
        fϱe  = (ϱᵣs*vxᵣ*vxᵣ/2 + ϱᵣs*vys*vys/2 + ϱᵣs*ɛᵣs + pᵣs)*vys
    else
        fϱ   = ϱᵣ*vyᵣ
        fϱvx = ϱᵣ*vxᵣ*vyᵣ
        fϱvy = ϱₗ*vyᵣ*vyᵣ + pᵣ
        fϱe  = (ϱᵣ*vxᵣ*vxᵣ/2 + ϱᵣ*vyᵣ*vyᵣ/2 + ϱᵣ*ɛᵣ + pᵣ)*vyᵣ
    end

    SVector(fϱ, fϱvx, fϱvy, fϱe)
end


Base.@pure function (fnum::ChandrashekarFluxEC)(uₗ::EulerVar2D{T}, uᵣ::EulerVar2D{T},
                                                model::Euler{T,2}, dir::Val{:x}) where T
    @unpack γ = model
    ϱₗ, vxₗ, vyₗ, pₗ = primitive_variables(uₗ, model)
    βₗ = ϱₗ / 2pₗ
    ϱᵣ, vxᵣ, vyᵣ, pᵣ = primitive_variables(uᵣ, model)
    βᵣ = ϱᵣ / 2pᵣ

    ϱ    = (ϱₗ + ϱᵣ) / 2
    ϱlog = logmean(ϱₗ, ϱᵣ)
    vx   = (vxₗ + vxᵣ) / 2
    vy   = (vyₗ + vyᵣ) / 2
    v2   = (vxₗ^2+vyₗ^2 + vxᵣ^2+vyᵣ^2) / 2
    β    = (βₗ + βᵣ) / 2
    βlog = logmean(βₗ, βᵣ)

    fϱ   = ϱlog*vx
    fϱvx = vx*fϱ + ϱ/2β
    fϱvy = vy*fϱ
    fϱe  = 1/(2γ-2)*fϱ/βlog - v2*fϱ/2 + vx*fϱvx + vy*fϱvy

    SVector(fϱ, fϱvx, fϱvy, fϱe)
end

Base.@pure function (fnum::ChandrashekarFluxEC)(uₗ::EulerVar2D{T}, uᵣ::EulerVar2D{T},
                                                model::Euler{T,2}, dir::Val{:y}) where T
    @unpack γ = model
    ϱₗ, vxₗ, vyₗ, pₗ = primitive_variables(uₗ, model)
    βₗ = ϱₗ / 2pₗ
    ϱᵣ, vxᵣ, vyᵣ, pᵣ = primitive_variables(uᵣ, model)
    βᵣ = ϱᵣ / 2pᵣ

    ϱ    = (ϱₗ + ϱᵣ) / 2
    ϱlog = logmean(ϱₗ, ϱᵣ)
    vx   = (vxₗ + vxᵣ) / 2
    vy   = (vyₗ + vyᵣ) / 2
    v2   = (vxₗ^2+vyₗ^2 + vxᵣ^2+vyᵣ^2) / 2
    β    = (βₗ + βᵣ) / 2
    βlog = logmean(βₗ, βᵣ)

    fϱ   = ϱlog*vy
    fϱvx = vx*fϱ
    fϱvy = vy*fϱ + ϱ/2β
    fϱe  = 1/(2γ-2)*fϱ/βlog - v2*fϱ/2 + vx*fϱvx + vy*fϱvy

    SVector(fϱ, fϱvx, fϱvy, fϱe)
end



Base.@pure function roe_variables(u::EulerVar2D, model::Euler)
    ϱ, vx, vy, p = primitive_variables(u, model)
    z1 = sqrt(ϱ/p)
    z2 = z1*vx
    z3 = z1*vy
    z5 = p*z1

    z1, z2, z3, z5
end

Base.@pure function (fnum::IsmailRoeFluxEC)(uₗ::EulerVar2D{T}, uᵣ::EulerVar2D{T},
                                            model::Euler{T,2}, dir::Val{:x}) where T
    @unpack γ = model
    z1ₗ, z2ₗ, z3ₗ, z5ₗ = roe_variables(uₗ, model)
    z1ᵣ, z2ᵣ, z3ᵣ, z5ᵣ = roe_variables(uᵣ, model)

    z1    = (z1ₗ + z1ᵣ) / 2
    z1log = logmean(z1ₗ, z1ᵣ)
    z2    = (z2ₗ + z2ᵣ) / 2
    z3    = (z3ₗ + z3ᵣ) / 2
    z5    = (z5ₗ + z5ᵣ) / 2
    z5log = logmean(z5ₗ, z5ᵣ)

    ϱ  = z1*z5log
    vx = z2/z1
    vy = z3/z1
    p1 = z5/z1
    p2 = ( (γ+1)*z5log/z1log + (γ-1)*z5/z1 ) / 2γ
    h  = γ/(γ-1) * p2/ϱ + (vx*vx + vy*vy)/2

    fϱ   = ϱ*vx
    fϱvx = ϱ*vx^2 + p1
    fϱvy = ϱ*vx*vy
    fϱe  = ϱ*vx*h

    SVector(fϱ, fϱvx, fϱvy, fϱe)
end

Base.@pure function (fnum::IsmailRoeFluxEC)(uₗ::EulerVar2D{T}, uᵣ::EulerVar2D{T},
                                            model::Euler{T,2}, dir::Val{:y}) where T
    @unpack γ = model
    z1ₗ, z2ₗ, z3ₗ, z5ₗ = roe_variables(uₗ, model)
    z1ᵣ, z2ᵣ, z3ᵣ, z5ᵣ = roe_variables(uᵣ, model)

    z1    = (z1ₗ + z1ᵣ) / 2
    z1log = logmean(z1ₗ, z1ᵣ)
    z2    = (z2ₗ + z2ᵣ) / 2
    z3    = (z3ₗ + z3ᵣ) / 2
    z5    = (z5ₗ + z5ᵣ) / 2
    z5log = logmean(z5ₗ, z5ᵣ)

    ϱ  = z1*z5log
    vx = z2/z1
    vy = z3/z1
    p1 = z5/z1
    p2 = ( (γ+1)*z5log/z1log + (γ-1)*z5/z1 ) / 2γ
    h  = γ/(γ-1) * p2/ϱ + (vx*vx + vy*vy)/2

    fϱ   = ϱ*vy
    fϱvx = ϱ*vx*vy
    fϱvy = ϱ*vy^2 + p1
    fϱe  = ϱ*vy*h

    SVector(fϱ, fϱvx, fϱvy, fϱe)
end


Base.@pure function (fnum::RanochaFluxECandKEP)(uₗ::EulerVar2D{T}, uᵣ::EulerVar2D{T},
                                                model::Euler{T,2}, dir::Val{:x}) where T
    @unpack γ = model
    ϱₗ, vxₗ, vyₗ, pₗ = primitive_variables(uₗ, model)
    ϱᵣ, vxᵣ, vyᵣ, pᵣ = primitive_variables(uᵣ, model)

    ϱ    = (ϱₗ + ϱᵣ) / 2
    ϱlog = logmean(ϱₗ, ϱᵣ)
    vx   = (vxₗ + vxᵣ) / 2
    vy   = (vyₗ + vyᵣ) / 2
    v2   = (vxₗ^2+vyₗ^2 + vxᵣ^2+vyᵣ^2) / 2
    p    = (pₗ + pᵣ) / 2
    ϱ_p_log = logmean(ϱₗ/pₗ, ϱᵣ/pᵣ)

    #TODO: Test
    fϱ   = ϱlog*vx
    fϱvx = vx*fϱ + p
    fϱvy = vy*fϱ
    fϱe  = vx*fϱvx - ϱlog*v2*vx/2 + ϱlog*vx/((γ-1)*ϱ_p_log) - (pᵣ-pₗ)*(vxᵣ-vxₗ)/4

    SVector(fϱ, fϱvx, fϱvy, fϱe)
end

Base.@pure function (fnum::RanochaFluxECandKEP)(uₗ::EulerVar2D{T}, uᵣ::EulerVar2D{T},
                                                model::Euler{T,2}, dir::Val{:y}) where T
    @unpack γ = model
    ϱₗ, vxₗ, vyₗ, pₗ = primitive_variables(uₗ, model)
    ϱᵣ, vxᵣ, vyᵣ, pᵣ = primitive_variables(uᵣ, model)

    ϱ    = (ϱₗ + ϱᵣ) / 2
    ϱlog = logmean(ϱₗ, ϱᵣ)
    vx   = (vxₗ + vxᵣ) / 2
    vy   = (vyₗ + vyᵣ) / 2
    v2   = (vxₗ^2+vyₗ^2 + vxᵣ^2+vyᵣ^2) / 2
    p    = (pₗ + pᵣ) / 2
    ϱ_p_log = logmean(ϱₗ/pₗ, ϱᵣ/pᵣ)

    # TODO: Test
    fϱ   = ϱlog*vy
    fϱvx = vx*fϱ
    fϱvy = vy*fϱ + p
    fϱe  = vy*fϱvy - ϱlog*v2*vy/2 + ϱlog*vy/((γ-1)*ϱ_p_log) - (pᵣ-pₗ)*(vyᵣ-vyₗ)/4

    SVector(fϱ, fϱvx, fϱvy, fϱe)
end
