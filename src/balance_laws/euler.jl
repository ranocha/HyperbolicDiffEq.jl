
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


function (fnum::SuliciuFlux)(uₗ::EulerVar2D{T}, uᵣ::EulerVar2D{T}, model::Euler{T,2}, dir::Val{:x}) where T
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
    fρ   = zero(T)
    fρvx = zero(T)
    fρvy = zero(T)
    fρe  = zero(T)
    if 0 <= vxₗ-cₗ_ϱₗ
        fρ   = ϱₗ*vxₗ
        fρvx = ϱₗ*vxₗ*vxₗ + pₗ
        fρvy = ϱₗ*vxₗ*vyₗ
        fρe  = (ϱₗ*vxₗ*vxₗ/2 + ϱₗ*vyₗ*vyₗ/2 + ϱₗ*ɛₗ + pₗ)*vxₗ
    elseif 0 <= vxs
        fρ   = ϱₗs*vxs
        fρvx = ϱₗs*vxs*vxs + pₗs
        fρvy = ϱₗ*vxs*vyₗ
        fρe  = (ϱₗs*vxs*vxs/2 + ϱₗs*vyₗ*vyₗ/2 + ϱₗs*ɛₗs + pₗs)*vxs
    elseif 0 <= vxᵣ+cᵣ_ϱᵣ
        fρ   = ϱᵣs*vxs
        fρvx = ϱᵣs*vxs*vxs + pᵣs
        fρvy = ϱₗ*vxs*vyᵣ
        fρe  = (ϱᵣs*vxs*vxs/2 + ϱᵣs*vyᵣ*vyᵣ/2 + ϱᵣs*ɛᵣs + pᵣs)*vxs
    else
        fρ   = ϱᵣ*vxᵣ
        fρvx = ϱᵣ*vxᵣ*vxᵣ + pᵣ
        fρvy = ϱₗ*vxᵣ*vyᵣ
        fρe  = (ϱᵣ*vxᵣ*vxᵣ/2 + ϱᵣ*vyᵣ*vyᵣ/2 + ϱᵣ*ɛᵣ + pᵣ)*vxᵣ
    end

    SVector(fρ, fρvx, fρvy, fρe)
end

function (fnum::SuliciuFlux)(uₗ::EulerVar2D{T}, uᵣ::EulerVar2D{T}, model::Euler{T,2}, dir::Val{:y}) where T
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
    fρ   = zero(T)
    fρvx = zero(T)
    fρvy = zero(T)
    fρe  = zero(T)
    if 0 <= vyₗ-cₗ_ϱₗ
        fρ   = ϱₗ*vyₗ
        fρvx = ϱₗ*vxₗ*vyₗ
        fρvy = ϱₗ*vyₗ*vyₗ + pₗ
        fρe  = (ϱₗ*vxₗ*vxₗ/2 + ϱₗ*vyₗ*vyₗ/2 + ϱₗ*ɛₗ + pₗ)*vyₗ
    elseif 0 <= vys
        fρ   = ϱₗs*vys
        fρvx = ϱₗs*vxₗ*vys
        fρvy = ϱₗ*vys*vys + pₗs
        fρe  = (ϱₗs*vxₗ*vxₗ/2 + ϱₗs*vys*vys/2 + ϱₗs*ɛₗs + pₗs)*vys
    elseif 0 <= vyᵣ+cᵣ_ϱᵣ
        fρ   = ϱᵣs*vys
        fρvx = ϱᵣs*vxᵣ*vys
        fρvy = ϱₗ*vys*vys + pᵣs
        fρe  = (ϱᵣs*vxᵣ*vxᵣ/2 + ϱᵣs*vys*vys/2 + ϱᵣs*ɛᵣs + pᵣs)*vys
    else
        fρ   = ϱᵣ*vyᵣ
        fρvx = ϱᵣ*vxᵣ*vyᵣ
        fρvy = ϱₗ*vyᵣ*vyᵣ + pᵣ
        fρe  = (ϱᵣ*vxᵣ*vxᵣ/2 + ϱᵣ*vyᵣ*vyᵣ/2 + ϱᵣ*ɛᵣ + pᵣ)*vyᵣ
    end

    SVector(fρ, fρvx, fρvy, fρe)
end
