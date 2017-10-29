
"""
    Euler{T<:Real,Dim}

The compressible Euler equations in `Dim` dimensions using `T` as scalar type.
"""
struct Euler{T<:Real,Dim} <: AbstractBalanceLaw{Dim}
    γ::T
end

function Euler(γ=1.4, Dim=1)
  Euler{typeof(γ), Dim}(γ)
end

function show{T,Dim}(io::IO, model::Euler{T,Dim})
  print(io, "Compressible Euler equations with γ=", model.γ, " {T=", T, ", Dim=", Dim, "}")
end


"""
    EulerVar1D{T<:Real}

Conserved variables (ϱ, ϱv, ϱe) of the Euler equations in one space dimension.
"""
struct EulerVar1D{T<:Real} <: FieldVector{3,T}
  ϱ ::T
  ϱv::T
  ϱe::T
end

function (::Type{EulerVar1D{T}}){T}(val::Real)
  EulerVar2D{T}(val, val, val)
end

function similar_type{T}(::EulerVar1D{T})
  EulerVar1D{T}
end

@inline variables{T}(model::Euler{T,1}) = EulerVar1D{T}


"""
    EulerVar2D{T<:Real}

Conserved variables (ϱ, ϱvx, ϱvy, ϱe) of the Euler equations in two space dimensions.
"""
struct EulerVar2D{T<:Real} <: FieldVector{4,T}
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


"""
    EulerVar3D{T<:Real}

Conserved variables (ϱ, ϱvx, ϱvy, ϱvz, ϱe) of the Euler equations in three space dimensions.
"""
struct EulerVar3D{T<:Real} <: FieldVector{5,T}
  ϱ  ::T
  ϱvx::T
  ϱvy::T
  ϱvz::T
  ϱe ::T
end

function (::Type{EulerVar3D{T}}){T}(val::Real)
  EulerVar3D{T}(val, val, val, val, val)
end

function similar_type{T}(::EulerVar3D{T})
  EulerVar3D{T}
end

@inline variables{T}(model::Euler{T,3}) = EulerVar3D{T}


"""
    IntegralQuantitiesEuler{T<:Real}

Some integrated quantities of interest for the Euler equations. Can be used in
callbacks of DifferentialEquations.jl.
"""
struct IntegralQuantitiesEuler{T<:Real} <: FieldVector{2,T}
    kinetic_energy::T
    entropy::T
end

function IntegralQuantitiesEuler(u, model::Euler)
    IntegralQuantitiesEuler(
        kinetic_energy(u, model),
        entropy(u, model)
    )
end


@inline function primitive_variables(u::EulerVar1D, model::Euler)
    @unpack ϱ, ϱv, ϱe = u
    @unpack γ = model
    if ϱ ≈ 0
        v = zero(ϱ)
    else
        v = ϱv / ϱ
    end
    p = (γ-1) * (ϱe - ϱv*v/2)

    ϱ, v, p
end

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

@inline function primitive_variables(u::EulerVar3D, model::Euler)
    @unpack ϱ, ϱvx, ϱvy, ϱvz, ϱe = u
    @unpack γ = model
    if ϱ ≈ 0
        vx = zero(ϱ)
        vy = zero(ϱ)
        vz = zero(ϱ)
    else
        vx = ϱvx / ϱ
        vy = ϱvy / ϱ
        vz = ϱvz / ϱ
    end
    p = (γ-1) * (ϱe - ϱvx*vx/2 - ϱvy*vy/2 - ϱvz*vz/2)

    ϱ, vx, vy, vz, p
end


@inline function conserved_variables(ϱ, v, p, model::Euler)
    @unpack γ = model

    variables(model)(ϱ, ϱ*v, ϱ*v^2/2+p/(γ-1))
end

@inline function conserved_variables(ϱ, vx, vy, p, model::Euler)
    @unpack γ = model

    variables(model)(ϱ, ϱ*vx, ϱ*vy, (ϱ*vx^2+ϱ*vy^2)/2+p/(γ-1))
end

@inline function conserved_variables(ϱ, vx, vy, vz, p, model::Euler)
    @unpack γ = model

    variables(model)(ϱ, ϱ*vx, ϱ*vy, ϱ*vz, (ϱ*vx^2+ϱ*vy^2+ϱ*vz^2)/2+p/(γ-1))
end


@inline function satisfies_physical_constraints(u::EulerVar2D, model::Euler)
    ϱ, vx, vy, p = primitive_variables(u, model)
    ϱ >= 0 && p >= 0
end

@inline function satisfies_physical_constraints(u::EulerVar3D, model::Euler)
    ϱ, vx, vy, vz, p = primitive_variables(u, model)
    ϱ >= 0 && p >= 0
end


Base.@pure function kinetic_energy(u::EulerVar2D, model::Euler)
    ϱ, vx, vy, p = primitive_variables(u, model)

    ϱ * (vx^2 + vy^2) / 2
end

Base.@pure function kinetic_energy(u::EulerVar3D, model::Euler)
    ϱ, vx, vy, vz, p = primitive_variables(u, model)

    ϱ * (vx^2 + vy^2 + vz^2) / 2
end


Base.@pure function entropy(u::EulerVar2D, model::Euler)
    @unpack γ = model
    ϱ, vx, vy, p = primitive_variables(u, model)

    -ϱ*log(p/ϱ^γ) / (γ-1)
end

Base.@pure function entropy(u::EulerVar3D, model::Euler)
    @unpack γ = model
    ϱ, vx, vy, vz, p = primitive_variables(u, model)

    -ϱ*log(p/ϱ^γ) / (γ-1)
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

@inline function flux{T}(u::EulerVar3D{T}, model::Euler, dir::Val{:x})
    @unpack ϱvx, ϱvy, ϱvz, ϱe = u
    ϱ, vx, vy, vz, p = primitive_variables(u, model)

    SVector{5,T}(ϱvx, ϱvx*vx + p, ϱvy*vx, ϱvz*vx, (ϱe+p)*vx)
end

@inline function flux{T}(u::EulerVar3D{T}, model::Euler, dir::Val{:y})
    @unpack ϱvx, ϱvy, ϱvz, ϱe = u
    ϱ, vx, vy, vz, p = primitive_variables(u, model)

    SVector{5,T}(ϱvy, ϱvx*vy, ϱvy*vy + p, ϱvz*vy, (ϱe+p)*vy)
end

@inline function flux{T}(u::EulerVar3D{T}, model::Euler, dir::Val{:z})
    @unpack ϱvx, ϱvy, ϱvz, ϱe = u
    ϱ, vx, vy, vz, p = primitive_variables(u, model)

    SVector{5,T}(ϱvz, ϱvx*vz, ϱvy*vz, ϱvz*vz + p, (ϱe+p)*vz)
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

@inline function max_abs_speed(u::EulerVar3D, model::Euler)
    @unpack γ = model
    ϱ, vx, vy, vz, p = primitive_variables(u, model)

    max(abs(vx),abs(vy),abs(vz)) + sqrt(γ * p / ϱ)
end

@inline function max_abs_speed(u::EulerVar3D, model::Euler, dir::Val{:x})
    @unpack γ = model
    ϱ, vx, vy, vz, p = primitive_variables(u, model)

    abs(vx) + sqrt(γ * p / ϱ)
end

@inline function max_abs_speed(u::EulerVar3D, model::Euler, dir::Val{:y})
    @unpack γ = model
    ϱ, vx, vy, vz, p = primitive_variables(u, model)

    abs(vy) + sqrt(γ * p / ϱ)
end

@inline function max_abs_speed(u::EulerVar3D, model::Euler, dir::Val{:z})
    @unpack γ = model
    ϱ, vx, vy, vz, p = primitive_variables(u, model)

    abs(vz) + sqrt(γ * p / ϱ)
end


@recipe function f{Ξ,T}(ξumodel::Tuple{Ξ,Vector{EulerVar1D{T}},Euler{T,1}})
  ξ, u, model = ξumodel

  ϱ  = mappedarray(u->u.ϱ, u)
  v  = mappedarray(u->primitive_variables(u,model)[2], u)
  p  = mappedarray(u->primitive_variables(u,model)[3], u)

  size   --> (1000,400)
  layout --> (1,3)
  legend --> false

  @series begin
      subplot := 1
      yguide --> L"\varrho"
      ξ, ϱ
  end

  @series begin
      subplot := 2
      yguide --> L"v"
      ξ, v
  end

  @series begin
      subplot := 3
      yguide --> L"p"
      ξ, p
  end
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


@inline function (fnum::SuliciuFlux)(uₗ::EulerVar2D{T}, uᵣ::EulerVar2D{T},
                                        model::Euler{T,2}, dir::Val{:x}) where T
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

@inline function (fnum::SuliciuFlux)(uₗ::EulerVar2D{T}, uᵣ::EulerVar2D{T},
                                        model::Euler{T,2}, dir::Val{:y}) where T
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

# TODO: (fnum::SuliciuFlux)(uₗ::EulerVar3D{T}, uᵣ::EulerVar3D{T}, model::Euler{T,3}, dir::Val{:x}) where T
# TODO: (fnum::SuliciuFlux)(uₗ::EulerVar3D{T}, uᵣ::EulerVar3D{T}, model::Euler{T,3}, dir::Val{:y}) where T
# TODO: (fnum::SuliciuFlux)(uₗ::EulerVar3D{T}, uᵣ::EulerVar3D{T}, model::Euler{T,3}, dir::Val{:z}) where T


@inline function (fnum::ChandrashekarFluxEC)(uₗ::EulerVar2D{T}, uᵣ::EulerVar2D{T},
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
@inline function (fnum::ChandrashekarFluxEC)(uₗ::EulerVar2D{T}, uᵣ::EulerVar2D{T},
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

@inline function (fnum::ChandrashekarFluxEC)(uₗ::EulerVar3D{T}, uᵣ::EulerVar3D{T},
                                                model::Euler{T,3}, dir::Val{:x}) where T
    @unpack γ = model
    ϱₗ, vxₗ, vyₗ, vzₗ, pₗ = primitive_variables(uₗ, model)
    βₗ = ϱₗ / 2pₗ
    ϱᵣ, vxᵣ, vyᵣ, vzᵣ, pᵣ = primitive_variables(uᵣ, model)
    βᵣ = ϱᵣ / 2pᵣ

    ϱ    = (ϱₗ + ϱᵣ) / 2
    ϱlog = logmean(ϱₗ, ϱᵣ)
    vx   = (vxₗ + vxᵣ) / 2
    vy   = (vyₗ + vyᵣ) / 2
    vz   = (vzₗ + vzᵣ) / 2
    v2   = (vxₗ^2+vyₗ^2+vzₗ^2 + vxᵣ^2+vyᵣ^2+vzᵣ^2) / 2
    β    = (βₗ + βᵣ) / 2
    βlog = logmean(βₗ, βᵣ)

    fϱ   = ϱlog*vx
    fϱvx = vx*fϱ + ϱ/2β
    fϱvy = vy*fϱ
    fϱvz = vz*fϱ
    fϱe  = 1/(2γ-2)*fϱ/βlog - v2*fϱ/2 + vx*fϱvx + vy*fϱvy + vz*fϱvz

    SVector(fϱ, fϱvx, fϱvy, fϱvz, fϱe)
end
@inline function (fnum::ChandrashekarFluxEC)(uₗ::EulerVar3D{T}, uᵣ::EulerVar3D{T},
                                                model::Euler{T,3}, dir::Val{:y}) where T
    @unpack γ = model
    ϱₗ, vxₗ, vyₗ, vzₗ, pₗ = primitive_variables(uₗ, model)
    βₗ = ϱₗ / 2pₗ
    ϱᵣ, vxᵣ, vyᵣ, vzᵣ, pᵣ = primitive_variables(uᵣ, model)
    βᵣ = ϱᵣ / 2pᵣ

    ϱ    = (ϱₗ + ϱᵣ) / 2
    ϱlog = logmean(ϱₗ, ϱᵣ)
    vx   = (vxₗ + vxᵣ) / 2
    vy   = (vyₗ + vyᵣ) / 2
    vz   = (vzₗ + vzᵣ) / 2
    v2   = (vxₗ^2+vyₗ^2+vzₗ^2 + vxᵣ^2+vyᵣ^2+vzᵣ^2) / 2
    β    = (βₗ + βᵣ) / 2
    βlog = logmean(βₗ, βᵣ)

    fϱ   = ϱlog*vy
    fϱvx = vx*fϱ
    fϱvy = vy*fϱ + ϱ/2β
    fϱvz = vz*fϱ
    fϱe  = 1/(2γ-2)*fϱ/βlog - v2*fϱ/2 + vx*fϱvx + vy*fϱvy + vz*fϱvz

    SVector(fϱ, fϱvx, fϱvy, fϱvz, fϱe)
end
@inline function (fnum::ChandrashekarFluxEC)(uₗ::EulerVar3D{T}, uᵣ::EulerVar3D{T},
                                                model::Euler{T,3}, dir::Val{:z}) where T
    @unpack γ = model
    ϱₗ, vxₗ, vyₗ, vzₗ, pₗ = primitive_variables(uₗ, model)
    βₗ = ϱₗ / 2pₗ
    ϱᵣ, vxᵣ, vyᵣ, vzᵣ, pᵣ = primitive_variables(uᵣ, model)
    βᵣ = ϱᵣ / 2pᵣ

    ϱ    = (ϱₗ + ϱᵣ) / 2
    ϱlog = logmean(ϱₗ, ϱᵣ)
    vx   = (vxₗ + vxᵣ) / 2
    vy   = (vyₗ + vyᵣ) / 2
    vz   = (vzₗ + vzᵣ) / 2
    v2   = (vxₗ^2+vyₗ^2+vzₗ^2 + vxᵣ^2+vyᵣ^2+vzᵣ^2) / 2
    β    = (βₗ + βᵣ) / 2
    βlog = logmean(βₗ, βᵣ)

    fϱ   = ϱlog*vz
    fϱvx = vx*fϱ
    fϱvy = vy*fϱ
    fϱvz = vz*fϱ + ϱ/2β
    fϱe  = 1/(2γ-2)*fϱ/βlog - v2*fϱ/2 + vx*fϱvx + vy*fϱvy + vz*fϱvz

    SVector(fϱ, fϱvx, fϱvy, fϱvz, fϱe)
end


Base.@pure function roe_variables(u::EulerVar2D, model::Euler)
    ϱ, vx, vy, p = primitive_variables(u, model)
    z1 = sqrt(ϱ/p)
    z2 = z1*vx
    z3 = z1*vy
    z5 = p*z1

    z1, z2, z3, z5
end

Base.@pure function roe_variables(u::EulerVar3D, model::Euler)
    ϱ, vx, vy, vz, p = primitive_variables(u, model)
    z1 = sqrt(ϱ/p)
    z2 = z1*vx
    z3 = z1*vy
    z4 = z1*vz
    z5 = p*z1

    z1, z2, z3, z4, z5
end

@inline function (fnum::IsmailRoeFluxEC)(uₗ::EulerVar2D{T}, uᵣ::EulerVar2D{T},
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
@inline function (fnum::IsmailRoeFluxEC)(uₗ::EulerVar2D{T}, uᵣ::EulerVar2D{T},
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

@inline function (fnum::IsmailRoeFluxEC)(uₗ::EulerVar3D{T}, uᵣ::EulerVar3D{T},
                                            model::Euler{T,3}, dir::Val{:x}) where T
    @unpack γ = model
    z1ₗ, z2ₗ, z3ₗ, z4ₗ, z5ₗ = roe_variables(uₗ, model)
    z1ᵣ, z2ᵣ, z3ᵣ, z4ᵣ, z5ᵣ = roe_variables(uᵣ, model)

    z1    = (z1ₗ + z1ᵣ) / 2
    z1log = logmean(z1ₗ, z1ᵣ)
    z2    = (z2ₗ + z2ᵣ) / 2
    z3    = (z3ₗ + z3ᵣ) / 2
    z4    = (z4ₗ + z4ᵣ) / 2
    z5    = (z5ₗ + z5ᵣ) / 2
    z5log = logmean(z5ₗ, z5ᵣ)

    ϱ  = z1*z5log
    vx = z2/z1
    vy = z3/z1
    vz = z4/z1
    p1 = z5/z1
    p2 = ( (γ+1)*z5log/z1log + (γ-1)*z5/z1 ) / 2γ
    h  = γ/(γ-1) * p2/ϱ + (vx^2+vy^2+vz^2)/2

    fϱ   = ϱ*vx
    fϱvx = ϱ*vx^2 + p1
    fϱvy = ϱ*vx*vy
    fϱvz = ϱ*vx*vz
    fϱe  = ϱ*vx*h

    SVector(fϱ, fϱvx, fϱvy, fϱvz, fϱe)
end
@inline function (fnum::IsmailRoeFluxEC)(uₗ::EulerVar3D{T}, uᵣ::EulerVar3D{T},
                                            model::Euler{T,3}, dir::Val{:y}) where T
    @unpack γ = model
    z1ₗ, z2ₗ, z3ₗ, z4ₗ, z5ₗ = roe_variables(uₗ, model)
    z1ᵣ, z2ᵣ, z3ᵣ, z4ᵣ, z5ᵣ = roe_variables(uᵣ, model)

    z1    = (z1ₗ + z1ᵣ) / 2
    z1log = logmean(z1ₗ, z1ᵣ)
    z2    = (z2ₗ + z2ᵣ) / 2
    z3    = (z3ₗ + z3ᵣ) / 2
    z4    = (z4ₗ + z4ᵣ) / 2
    z5    = (z5ₗ + z5ᵣ) / 2
    z5log = logmean(z5ₗ, z5ᵣ)

    ϱ  = z1*z5log
    vx = z2/z1
    vy = z3/z1
    vz = z4/z1
    p1 = z5/z1
    p2 = ( (γ+1)*z5log/z1log + (γ-1)*z5/z1 ) / 2γ
    h  = γ/(γ-1) * p2/ϱ + (vx^2+vy^2+vz^2)/2

    fϱ   = ϱ*vy
    fϱvx = ϱ*vx*vy
    fϱvy = ϱ*vy^2 + p1
    fϱvz = ϱ*vy*vz
    fϱe  = ϱ*vy*h

    SVector(fϱ, fϱvx, fϱvy, fϱvz, fϱe)
end
@inline function (fnum::IsmailRoeFluxEC)(uₗ::EulerVar3D{T}, uᵣ::EulerVar3D{T},
                                            model::Euler{T,3}, dir::Val{:z}) where T
    @unpack γ = model
    z1ₗ, z2ₗ, z3ₗ, z4ₗ, z5ₗ = roe_variables(uₗ, model)
    z1ᵣ, z2ᵣ, z3ᵣ, z4ᵣ, z5ᵣ = roe_variables(uᵣ, model)

    z1    = (z1ₗ + z1ᵣ) / 2
    z1log = logmean(z1ₗ, z1ᵣ)
    z2    = (z2ₗ + z2ᵣ) / 2
    z3    = (z3ₗ + z3ᵣ) / 2
    z4    = (z4ₗ + z4ᵣ) / 2
    z5    = (z5ₗ + z5ᵣ) / 2
    z5log = logmean(z5ₗ, z5ᵣ)

    ϱ  = z1*z5log
    vx = z2/z1
    vy = z3/z1
    vz = z4/z1
    p1 = z5/z1
    p2 = ( (γ+1)*z5log/z1log + (γ-1)*z5/z1 ) / 2γ
    h  = γ/(γ-1) * p2/ϱ + (vx^2+vy^2+vz^2)/2

    fϱ   = ϱ*vz
    fϱvx = ϱ*vx*vz
    fϱvy = ϱ*vy*vz
    fϱvz = ϱ*vz^2 + p1
    fϱe  = ϱ*vz*h

    SVector(fϱ, fϱvx, fϱvy, fϱvz, fϱe)
end


@inline function (fnum::RanochaFluxECandKEP)(uₗ::EulerVar2D{T}, uᵣ::EulerVar2D{T},
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

    fϱ   = ϱlog*vx
    fϱvx = vx*fϱ + p
    fϱvy = vy*fϱ
    fϱe  = (ϱlog/((γ-1)*ϱ_p_log) + p + ϱlog*(vx^2+vy^2-v2/2)) * vx - (pᵣ-pₗ)*(vxᵣ-vxₗ)/4

    SVector(fϱ, fϱvx, fϱvy, fϱe)
end
@inline function (fnum::RanochaFluxECandKEP)(uₗ::EulerVar2D{T}, uᵣ::EulerVar2D{T},
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

    fϱ   = ϱlog*vy
    fϱvx = vx*fϱ
    fϱvy = vy*fϱ + p
    fϱe  = (ϱlog/((γ-1)*ϱ_p_log) + p + ϱlog*(vx^2+vy^2-v2/2)) * vy - (pᵣ-pₗ)*(vyᵣ-vyₗ)/4

    SVector(fϱ, fϱvx, fϱvy, fϱe)
end

@inline function (fnum::RanochaFluxECandKEP)(uₗ::EulerVar3D{T}, uᵣ::EulerVar3D{T},
                                                model::Euler{T,3}, dir::Val{:x}) where T
    @unpack γ = model
    ϱₗ, vxₗ, vyₗ, vzₗ, pₗ = primitive_variables(uₗ, model)
    ϱᵣ, vxᵣ, vyᵣ, vzᵣ, pᵣ = primitive_variables(uᵣ, model)

    ϱ    = (ϱₗ + ϱᵣ) / 2
    ϱlog = logmean(ϱₗ, ϱᵣ)
    vx   = (vxₗ + vxᵣ) / 2
    vy   = (vyₗ + vyᵣ) / 2
    vz   = (vzₗ + vzᵣ) / 2
    v2   = (vxₗ^2+vyₗ^2+vzₗ^2 + vxᵣ^2+vyᵣ^2+vzᵣ^2) / 2
    p    = (pₗ + pᵣ) / 2
    ϱ_p_log = logmean(ϱₗ/pₗ, ϱᵣ/pᵣ)

    fϱ   = ϱlog*vx
    fϱvx = vx*fϱ + p
    fϱvy = vy*fϱ
    fϱvz = vz*fϱ
    fϱe  = (ϱlog/((γ-1)*ϱ_p_log) + p + ϱlog*(vx^2+vy^2+vz^2-v2/2)) * vx - (pᵣ-pₗ)*(vxᵣ-vxₗ)/4

    SVector(fϱ, fϱvx, fϱvy, fϱvz, fϱe)
end
@inline function (fnum::RanochaFluxECandKEP)(uₗ::EulerVar3D{T}, uᵣ::EulerVar3D{T},
                                                model::Euler{T,3}, dir::Val{:y}) where T
    @unpack γ = model
    ϱₗ, vxₗ, vyₗ, vzₗ, pₗ = primitive_variables(uₗ, model)
    ϱᵣ, vxᵣ, vyᵣ, vzᵣ, pᵣ = primitive_variables(uᵣ, model)

    ϱ    = (ϱₗ + ϱᵣ) / 2
    ϱlog = logmean(ϱₗ, ϱᵣ)
    vx   = (vxₗ + vxᵣ) / 2
    vy   = (vyₗ + vyᵣ) / 2
    vz   = (vzₗ + vzᵣ) / 2
    v2   = (vxₗ^2+vyₗ^2+vzₗ^2 + vxᵣ^2+vyᵣ^2+vzᵣ^2) / 2
    p    = (pₗ + pᵣ) / 2
    ϱ_p_log = logmean(ϱₗ/pₗ, ϱᵣ/pᵣ)

    fϱ   = ϱlog*vy
    fϱvx = vx*fϱ
    fϱvy = vy*fϱ + p
    fϱvz = vz*fϱ
    fϱe  = (ϱlog/((γ-1)*ϱ_p_log) + p + ϱlog*(vx^2+vy^2+vz^2-v2/2)) * vy - (pᵣ-pₗ)*(vyᵣ-vyₗ)/4

    SVector(fϱ, fϱvx, fϱvy, fϱvz, fϱe)
end
@inline function (fnum::RanochaFluxECandKEP)(uₗ::EulerVar3D{T}, uᵣ::EulerVar3D{T},
                                                model::Euler{T,3}, dir::Val{:z}) where T
    @unpack γ = model
    ϱₗ, vxₗ, vyₗ, vzₗ, pₗ = primitive_variables(uₗ, model)
    ϱᵣ, vxᵣ, vyᵣ, vzᵣ, pᵣ = primitive_variables(uᵣ, model)

    ϱ    = (ϱₗ + ϱᵣ) / 2
    ϱlog = logmean(ϱₗ, ϱᵣ)
    vx   = (vxₗ + vxᵣ) / 2
    vy   = (vyₗ + vyᵣ) / 2
    vz   = (vzₗ + vzᵣ) / 2
    v2   = (vxₗ^2+vyₗ^2+vzₗ^2 + vxᵣ^2+vyᵣ^2+vzᵣ^2) / 2
    p    = (pₗ + pᵣ) / 2
    ϱ_p_log = logmean(ϱₗ/pₗ, ϱᵣ/pᵣ)

    fϱ   = ϱlog*vz
    fϱvx = vx*fϱ
    fϱvy = vy*fϱ
    fϱvz = vz*fϱ + p
    fϱe  = (ϱlog/((γ-1)*ϱ_p_log) + p + ϱlog*(vx^2+vy^2+vz^2-v2/2)) * vz - (pᵣ-pₗ)*(vzᵣ-vzₗ)/4

    SVector(fϱ, fϱvx, fϱvy, fϱvz, fϱe)
end


@inline function (fnum::MorinishiFlux)(uₗ::EulerVar2D{T}, uᵣ::EulerVar2D{T},
                                        model::Euler{T,2}, dir::Val{:x}) where T
    @unpack γ = model
    ϱₗ, vxₗ, vyₗ, pₗ = primitive_variables(uₗ, model)
    ϱvxₗ = uₗ.ϱvx
    ϱᵣ, vxᵣ, vyᵣ, pᵣ = primitive_variables(uᵣ, model)
    ϱvxᵣ = uᵣ.ϱvx

    ϱvx  = (ϱvxₗ + ϱvxᵣ) / 2
    vx   = (vxₗ + vxᵣ) / 2
    vy   = (vyₗ + vyᵣ) / 2
    p    = (pₗ + pᵣ) / 2
    pvx  = (pₗ*vxₗ + pᵣ*vxᵣ) / 2
    ϱvxx = (ϱvxₗ*vxₗ + ϱvxᵣ*vxᵣ) / 2
    ϱvxy = (ϱvxₗ*vyₗ + ϱvxᵣ*vyᵣ) / 2
    ϱvxv2= (ϱvxₗ*(vxₗ^2+vyₗ^2) + ϱvxᵣ*(vxᵣ^2+vyᵣ^2)) / 2

    fϱ   = ϱvx
    fϱvx = ϱvx*vx + p
    fϱvy = ϱvx*vy
    fϱe  = γ/(γ-1)*pvx + ϱvxx*vx + ϱvxy*vy - ϱvxv2/2

    SVector(fϱ, fϱvx, fϱvy, fϱe)
end
@inline function (fnum::MorinishiFlux)(uₗ::EulerVar2D{T}, uᵣ::EulerVar2D{T},
                                        model::Euler{T,2}, dir::Val{:y}) where T
    @unpack γ = model
    ϱₗ, vxₗ, vyₗ, pₗ = primitive_variables(uₗ, model)
    ϱvyₗ = uₗ.ϱvy
    ϱᵣ, vxᵣ, vyᵣ, pᵣ = primitive_variables(uᵣ, model)
    ϱvyᵣ = uᵣ.ϱvy

    ϱvy  = (ϱvyₗ + ϱvyᵣ) / 2
    vx   = (vxₗ + vxᵣ) / 2
    vy   = (vyₗ + vyᵣ) / 2
    p    = (pₗ + pᵣ) / 2
    pvy  = (pₗ*vyₗ + pᵣ*vyᵣ) / 2
    ϱvyy = (ϱvyₗ*vyₗ + ϱvyᵣ*vyᵣ) / 2
    ϱvyx = (ϱvyₗ*vxₗ + ϱvyᵣ*vxᵣ) / 2
    ϱvyv2= (ϱvyₗ*(vxₗ^2+vyₗ^2) + ϱvyᵣ*(vxᵣ^2+vyᵣ^2)) / 2

    fϱ   = ϱvy
    fϱvx = ϱvy*vx
    fϱvy = ϱvy*vy + p
    fϱe  = γ/(γ-1)*pvy + ϱvyx*vx + ϱvyy*vy - ϱvyv2/2

    SVector(fϱ, fϱvx, fϱvy, fϱe)
end

@inline function (fnum::MorinishiFlux)(uₗ::EulerVar3D{T}, uᵣ::EulerVar3D{T},
                                        model::Euler{T,3}, dir::Val{:x}) where T
    @unpack γ = model
    ϱₗ, vxₗ, vyₗ, vzₗ, pₗ = primitive_variables(uₗ, model)
    ϱvxₗ = uₗ.ϱvx
    ϱᵣ, vxᵣ, vyᵣ, vzᵣ, pᵣ = primitive_variables(uᵣ, model)
    ϱvxᵣ = uᵣ.ϱvx

    ϱvx  = (ϱvxₗ + ϱvxᵣ) / 2
    vx   = (vxₗ + vxᵣ) / 2
    vy   = (vyₗ + vyᵣ) / 2
    vz   = (vzₗ + vzᵣ) / 2
    p    = (pₗ + pᵣ) / 2
    pvx  = (pₗ*vxₗ + pᵣ*vxᵣ) / 2
    ϱvxx = (ϱvxₗ*vxₗ + ϱvxᵣ*vxᵣ) / 2
    ϱvxy = (ϱvxₗ*vyₗ + ϱvxᵣ*vyᵣ) / 2
    ϱvxz = (ϱvxₗ*vzₗ + ϱvxᵣ*vzᵣ) / 2
    ϱvxv2= (ϱvxₗ*(vxₗ^2+vyₗ^2+vzₗ^2) + ϱvxᵣ*(vxᵣ^2+vyᵣ^2+vzᵣ^2)) / 2

    fϱ   = ϱvx
    fϱvx = ϱvx*vx + p
    fϱvy = ϱvx*vy
    fϱvz = ϱvx*vz
    fϱe  = γ/(γ-1)*pvx + ϱvxx*vx + ϱvxy*vy + ϱvxz*vz - ϱvxv2/2

    SVector(fϱ, fϱvx, fϱvy, fϱvz, fϱe)
end
@inline function (fnum::MorinishiFlux)(uₗ::EulerVar3D{T}, uᵣ::EulerVar3D{T},
                                        model::Euler{T,3}, dir::Val{:y}) where T
    @unpack γ = model
    ϱₗ, vxₗ, vyₗ, vzₗ, pₗ = primitive_variables(uₗ, model)
    ϱvyₗ = uₗ.ϱvy
    ϱᵣ, vxᵣ, vyᵣ, vzᵣ, pᵣ = primitive_variables(uᵣ, model)
    ϱvyᵣ = uᵣ.ϱvy

    ϱvy  = (ϱvyₗ + ϱvyᵣ) / 2
    vx   = (vxₗ + vxᵣ) / 2
    vy   = (vyₗ + vyᵣ) / 2
    vz   = (vzₗ + vzᵣ) / 2
    p    = (pₗ + pᵣ) / 2
    pvy  = (pₗ*vyₗ + pᵣ*vyᵣ) / 2
    ϱvyx = (ϱvyₗ*vxₗ + ϱvyᵣ*vxᵣ) / 2
    ϱvyy = (ϱvyₗ*vyₗ + ϱvyᵣ*vyᵣ) / 2
    ϱvyz = (ϱvyₗ*vzₗ + ϱvyᵣ*vzᵣ) / 2
    ϱvyv2= (ϱvyₗ*(vxₗ^2+vyₗ^2+vzₗ^2) + ϱvyᵣ*(vxᵣ^2+vyᵣ^2+vzᵣ^2)) / 2

    fϱ   = ϱvy
    fϱvx = ϱvy*vx
    fϱvy = ϱvy*vy + p
    fϱvz = ϱvy*vz
    fϱe  = γ/(γ-1)*pvy + ϱvyx*vx + ϱvyy*vy + ϱvyz*vz - ϱvyv2/2

    SVector(fϱ, fϱvx, fϱvy, fϱvz, fϱe)
end
@inline function (fnum::MorinishiFlux)(uₗ::EulerVar3D{T}, uᵣ::EulerVar3D{T},
                                        model::Euler{T,3}, dir::Val{:z}) where T
    @unpack γ = model
    ϱₗ, vxₗ, vyₗ, vzₗ, pₗ = primitive_variables(uₗ, model)
    ϱvzₗ = uₗ.ϱvz
    ϱᵣ, vxᵣ, vyᵣ, vzᵣ, pᵣ = primitive_variables(uᵣ, model)
    ϱvzᵣ = uᵣ.ϱvz

    ϱvz  = (ϱvzₗ + ϱvzᵣ) / 2
    vx   = (vxₗ + vxᵣ) / 2
    vy   = (vyₗ + vyᵣ) / 2
    vz   = (vzₗ + vzᵣ) / 2
    p    = (pₗ + pᵣ) / 2
    pvz  = (pₗ*vzₗ + pᵣ*vzᵣ) / 2
    ϱvzx = (ϱvzₗ*vxₗ + ϱvzᵣ*vxᵣ) / 2
    ϱvzy = (ϱvzₗ*vyₗ + ϱvzᵣ*vyᵣ) / 2
    ϱvzz = (ϱvzₗ*vzₗ + ϱvzᵣ*vzᵣ) / 2
    ϱvzv2= (ϱvzₗ*(vxₗ^2+vyₗ^2+vzₗ^2) + ϱvzᵣ*(vxᵣ^2+vyᵣ^2+vzᵣ^2)) / 2

    fϱ   = ϱvz
    fϱvx = ϱvz*vx
    fϱvy = ϱvz*vy
    fϱvz = ϱvz*vz + p
    fϱe  = γ/(γ-1)*pvz + ϱvzx*vx + ϱvzy*vy + ϱvzz*vz - ϱvzv2/2

    SVector(fϱ, fϱvx, fϱvy, fϱvz, fϱe)
end


@inline function (fnum::DucrosEtAlFlux)(uₗ::EulerVar2D{T}, uᵣ::EulerVar2D{T},
                                        model::Euler{T,2}, dir::Val{:x}) where T
    @unpack γ = model
    ϱₗ, vxₗ, vyₗ, pₗ = primitive_variables(uₗ, model)
    ϱvxₗ, ϱvyₗ, ϱeₗ = uₗ.ϱvx, uₗ.ϱvy, uₗ.ϱe
    ϱᵣ, vxᵣ, vyᵣ, pᵣ = primitive_variables(uᵣ, model)
    ϱvxᵣ, ϱvyᵣ, ϱeᵣ = uᵣ.ϱvx, uᵣ.ϱvy, uᵣ.ϱe

    ϱ   = (ϱₗ + ϱᵣ) / 2
    ϱvx = (ϱvxₗ + ϱvxᵣ) / 2
    ϱvy = (ϱvyₗ + ϱvyᵣ) / 2
    ϱe  = (ϱeₗ + ϱeᵣ) / 2
    vx  = (vxₗ + vxᵣ) / 2
    vy  = (vyₗ + vyᵣ) / 2
    p   = (pₗ + pᵣ) / 2

    fϱ   = ϱ*vx
    fϱvx = ϱvx*vx + p
    fϱvy = ϱvy*vx
    fϱe  = (ϱe+p)*vx

    SVector(fϱ, fϱvx, fϱvy, fϱe)
end
@inline function (fnum::DucrosEtAlFlux)(uₗ::EulerVar2D{T}, uᵣ::EulerVar2D{T},
                                        model::Euler{T,2}, dir::Val{:y}) where T
    @unpack γ = model
    ϱₗ, vxₗ, vyₗ, pₗ = primitive_variables(uₗ, model)
    ϱvxₗ, ϱvyₗ, ϱeₗ = uₗ.ϱvx, uₗ.ϱvy, uₗ.ϱe
    ϱᵣ, vxᵣ, vyᵣ, pᵣ = primitive_variables(uᵣ, model)
    ϱvxᵣ, ϱvyᵣ, ϱeᵣ = uᵣ.ϱvx, uᵣ.ϱvy, uᵣ.ϱe

    ϱ   = (ϱₗ + ϱᵣ) / 2
    ϱvx = (ϱvxₗ + ϱvxᵣ) / 2
    ϱvy = (ϱvyₗ + ϱvyᵣ) / 2
    ϱe  = (ϱeₗ + ϱeᵣ) / 2
    vx  = (vxₗ + vxᵣ) / 2
    vy  = (vyₗ + vyᵣ) / 2
    p   = (pₗ + pᵣ) / 2

    fϱ   = ϱ*vy
    fϱvx = ϱvx*vy
    fϱvy = ϱvy*vy + p
    fϱe  = (ϱe+p)*vy

    SVector(fϱ, fϱvx, fϱvy, fϱe)
end

@inline function (fnum::DucrosEtAlFlux)(uₗ::EulerVar3D{T}, uᵣ::EulerVar3D{T},
                                        model::Euler{T,3}, dir::Val{:x}) where T
    @unpack γ = model
    ϱₗ, vxₗ, vyₗ, vzₗ, pₗ = primitive_variables(uₗ, model)
    ϱvxₗ, ϱvyₗ, ϱvzₗ, ϱeₗ = uₗ.ϱvx, uₗ.ϱvy, uₗ. ϱvz, uₗ.ϱe
    ϱᵣ, vxᵣ, vyᵣ, vzᵣ, pᵣ = primitive_variables(uᵣ, model)
    ϱvxᵣ, ϱvyᵣ, ϱvzᵣ, ϱeᵣ = uᵣ.ϱvx, uᵣ.ϱvy, uᵣ.ϱvz, uᵣ.ϱe

    ϱ   = (ϱₗ + ϱᵣ) / 2
    ϱvx = (ϱvxₗ + ϱvxᵣ) / 2
    ϱvy = (ϱvyₗ + ϱvyᵣ) / 2
    ϱvz = (ϱvzₗ + ϱvzᵣ) / 2
    ϱe  = (ϱeₗ + ϱeᵣ) / 2
    vx  = (vxₗ + vxᵣ) / 2
    vy  = (vyₗ + vyᵣ) / 2
    vz  = (vzₗ + vzᵣ) / 2
    p   = (pₗ + pᵣ) / 2

    fϱ   = ϱ*vx
    fϱvx = ϱvx*vx + p
    fϱvy = ϱvy*vx
    fϱvz = ϱvz*vx
    fϱe  = (ϱe+p)*vx

    SVector(fϱ, fϱvx, fϱvy, fϱvz, fϱe)
end
@inline function (fnum::DucrosEtAlFlux)(uₗ::EulerVar3D{T}, uᵣ::EulerVar3D{T},
                                        model::Euler{T,3}, dir::Val{:y}) where T
    @unpack γ = model
    ϱₗ, vxₗ, vyₗ, vzₗ, pₗ = primitive_variables(uₗ, model)
    ϱvxₗ, ϱvyₗ, ϱvzₗ, ϱeₗ = uₗ.ϱvx, uₗ.ϱvy, uₗ. ϱvz, uₗ.ϱe
    ϱᵣ, vxᵣ, vyᵣ, vzᵣ, pᵣ = primitive_variables(uᵣ, model)
    ϱvxᵣ, ϱvyᵣ, ϱvzᵣ, ϱeᵣ = uᵣ.ϱvx, uᵣ.ϱvy, uᵣ.ϱvz, uᵣ.ϱe

    ϱ   = (ϱₗ + ϱᵣ) / 2
    ϱvx = (ϱvxₗ + ϱvxᵣ) / 2
    ϱvy = (ϱvyₗ + ϱvyᵣ) / 2
    ϱvz = (ϱvzₗ + ϱvzᵣ) / 2
    ϱe  = (ϱeₗ + ϱeᵣ) / 2
    vx  = (vxₗ + vxᵣ) / 2
    vy  = (vyₗ + vyᵣ) / 2
    vz  = (vzₗ + vzᵣ) / 2
    p   = (pₗ + pᵣ) / 2

    fϱ   = ϱ*vy
    fϱvx = ϱvx*vy
    fϱvy = ϱvy*vy + p
    fϱvz = ϱvz*vy
    fϱe  = (ϱe+p)*vy

    SVector(fϱ, fϱvx, fϱvy, fϱvz, fϱe)
end
@inline function (fnum::DucrosEtAlFlux)(uₗ::EulerVar3D{T}, uᵣ::EulerVar3D{T},
                                        model::Euler{T,3}, dir::Val{:z}) where T
    @unpack γ = model
    ϱₗ, vxₗ, vyₗ, vzₗ, pₗ = primitive_variables(uₗ, model)
    ϱvxₗ, ϱvyₗ, ϱvzₗ, ϱeₗ = uₗ.ϱvx, uₗ.ϱvy, uₗ. ϱvz, uₗ.ϱe
    ϱᵣ, vxᵣ, vyᵣ, vzᵣ, pᵣ = primitive_variables(uᵣ, model)
    ϱvxᵣ, ϱvyᵣ, ϱvzᵣ, ϱeᵣ = uᵣ.ϱvx, uᵣ.ϱvy, uᵣ.ϱvz, uᵣ.ϱe

    ϱ   = (ϱₗ + ϱᵣ) / 2
    ϱvx = (ϱvxₗ + ϱvxᵣ) / 2
    ϱvy = (ϱvyₗ + ϱvyᵣ) / 2
    ϱvz = (ϱvzₗ + ϱvzᵣ) / 2
    ϱe  = (ϱeₗ + ϱeᵣ) / 2
    vx  = (vxₗ + vxᵣ) / 2
    vy  = (vyₗ + vyᵣ) / 2
    vz  = (vzₗ + vzᵣ) / 2
    p   = (pₗ + pᵣ) / 2

    fϱ   = ϱ*vz
    fϱvx = ϱvx*vz
    fϱvy = ϱvy*vz
    fϱvz = ϱvz*vz + p
    fϱe  = (ϱe+p)*vz

    SVector(fϱ, fϱvx, fϱvy, fϱvz, fϱe)
end


@inline function (fnum::KennedyGruberFlux)(uₗ::EulerVar2D{T}, uᵣ::EulerVar2D{T},
                                            model::Euler{T,2}, dir::Val{:x}) where T
    @unpack γ = model
    ϱₗ, vxₗ, vyₗ, pₗ = primitive_variables(uₗ, model)
    eₗ = uₗ.ϱe / ϱₗ
    ϱᵣ, vxᵣ, vyᵣ, pᵣ = primitive_variables(uᵣ, model)
    eᵣ = uᵣ.ϱe / ϱᵣ

    ϱ  = (ϱₗ + ϱᵣ) / 2
    vx = (vxₗ + vxᵣ) / 2
    vy = (vyₗ + vyᵣ) / 2
    p  = (pₗ + pᵣ) / 2
    e  = (eₗ + eᵣ) / 2

    fϱ   = ϱ*vx
    fϱvx = ϱ*vx^2 + p
    fϱvy = ϱ*vx*vy
    fϱe  = (ϱ*e+p)*vx

    SVector(fϱ, fϱvx, fϱvy, fϱe)
end
@inline function (fnum::KennedyGruberFlux)(uₗ::EulerVar2D{T}, uᵣ::EulerVar2D{T},
                                            model::Euler{T,2}, dir::Val{:y}) where T
    @unpack γ = model
    ϱₗ, vxₗ, vyₗ, pₗ = primitive_variables(uₗ, model)
    eₗ = uₗ.ϱe / ϱₗ
    ϱᵣ, vxᵣ, vyᵣ, pᵣ = primitive_variables(uᵣ, model)
    eᵣ = uᵣ.ϱe / ϱᵣ

    ϱ  = (ϱₗ + ϱᵣ) / 2
    vx = (vxₗ + vxᵣ) / 2
    vy = (vyₗ + vyᵣ) / 2
    p  = (pₗ + pᵣ) / 2
    e  = (eₗ + eᵣ) / 2

    fϱ   = ϱ*vy
    fϱvx = ϱ*vx*vy
    fϱvy = ϱ*vy^2 + p
    fϱe  = (ϱ*e+p)*vy

    SVector(fϱ, fϱvx, fϱvy, fϱe)
end

@inline function (fnum::KennedyGruberFlux)(uₗ::EulerVar3D{T}, uᵣ::EulerVar3D{T},
                                            model::Euler{T,3}, dir::Val{:x}) where T
    @unpack γ = model
    ϱₗ, vxₗ, vyₗ, vzₗ, pₗ = primitive_variables(uₗ, model)
    eₗ = uₗ.ϱe / ϱₗ
    ϱᵣ, vxᵣ, vyᵣ, vzᵣ, pᵣ = primitive_variables(uᵣ, model)
    eᵣ = uᵣ.ϱe / ϱᵣ

    ϱ  = (ϱₗ + ϱᵣ) / 2
    vx = (vxₗ + vxᵣ) / 2
    vy = (vyₗ + vyᵣ) / 2
    vz = (vzₗ + vzᵣ) / 2
    p  = (pₗ + pᵣ) / 2
    e  = (eₗ + eᵣ) / 2

    fϱ   = ϱ*vx
    fϱvx = ϱ*vx^2 + p
    fϱvy = ϱ*vx*vy
    fϱvz = ϱ*vx*vz
    fϱe  = (ϱ*e+p)*vx

    SVector(fϱ, fϱvx, fϱvy, fϱvz, fϱe)
end
@inline function (fnum::KennedyGruberFlux)(uₗ::EulerVar3D{T}, uᵣ::EulerVar3D{T},
                                            model::Euler{T,3}, dir::Val{:y}) where T
    @unpack γ = model
    ϱₗ, vxₗ, vyₗ, vzₗ, pₗ = primitive_variables(uₗ, model)
    eₗ = uₗ.ϱe / ϱₗ
    ϱᵣ, vxᵣ, vyᵣ, vzᵣ, pᵣ = primitive_variables(uᵣ, model)
    eᵣ = uᵣ.ϱe / ϱᵣ

    ϱ  = (ϱₗ + ϱᵣ) / 2
    vx = (vxₗ + vxᵣ) / 2
    vy = (vyₗ + vyᵣ) / 2
    vz = (vzₗ + vzᵣ) / 2
    p  = (pₗ + pᵣ) / 2
    e  = (eₗ + eᵣ) / 2

    fϱ   = ϱ*vy
    fϱvx = ϱ*vx*vy
    fϱvy = ϱ*vy^2 + p
    fϱvz = ϱ*vy*vz
    fϱe  = (ϱ*e+p)*vy

    SVector(fϱ, fϱvx, fϱvy, fϱvz, fϱe)
end
@inline function (fnum::KennedyGruberFlux)(uₗ::EulerVar3D{T}, uᵣ::EulerVar3D{T},
                                            model::Euler{T,3}, dir::Val{:z}) where T
    @unpack γ = model
    ϱₗ, vxₗ, vyₗ, vzₗ, pₗ = primitive_variables(uₗ, model)
    eₗ = uₗ.ϱe / ϱₗ
    ϱᵣ, vxᵣ, vyᵣ, vzᵣ, pᵣ = primitive_variables(uᵣ, model)
    eᵣ = uᵣ.ϱe / ϱᵣ

    ϱ  = (ϱₗ + ϱᵣ) / 2
    vx = (vxₗ + vxᵣ) / 2
    vy = (vyₗ + vyᵣ) / 2
    vz = (vzₗ + vzᵣ) / 2
    p  = (pₗ + pᵣ) / 2
    e  = (eₗ + eᵣ) / 2

    fϱ   = ϱ*vz
    fϱvx = ϱ*vx*vz
    fϱvy = ϱ*vy*vz
    fϱvz = ϱ*vz^2 + p
    fϱe  = (ϱ*e+p)*vz

    SVector(fϱ, fϱvx, fϱvy, fϱvz, fϱe)
end


@inline function (fnum::PirozzoliFlux)(uₗ::EulerVar2D{T}, uᵣ::EulerVar2D{T},
                                        model::Euler{T,2}, dir::Val{:x}) where T
    @unpack γ = model
    ϱₗ, vxₗ, vyₗ, pₗ = primitive_variables(uₗ, model)
    ϱeₗ = uₗ.ϱe
    ϱᵣ, vxᵣ, vyᵣ, pᵣ = primitive_variables(uᵣ, model)
    ϱeᵣ = uᵣ.ϱe

    ϱ  = (ϱₗ + ϱᵣ) / 2
    vx = (vxₗ + vxᵣ) / 2
    vy = (vyₗ + vyᵣ) / 2
    p  = (pₗ + pᵣ) / 2
    h  = ((ϱeₗ+pₗ)/ϱₗ + (ϱeᵣ+pᵣ)/ϱᵣ) / 2

    fϱ   = ϱ*vx
    fϱvx = ϱ*vx^2 + p
    fϱvy = ϱ*vx*vy
    fϱe  = ϱ*h*vx

    SVector(fϱ, fϱvx, fϱvy, fϱe)
end
@inline function (fnum::PirozzoliFlux)(uₗ::EulerVar2D{T}, uᵣ::EulerVar2D{T},
                                        model::Euler{T,2}, dir::Val{:y}) where T
    @unpack γ = model
    ϱₗ, vxₗ, vyₗ, pₗ = primitive_variables(uₗ, model)
    ϱeₗ = uₗ.ϱe
    ϱᵣ, vxᵣ, vyᵣ, pᵣ = primitive_variables(uᵣ, model)
    ϱeᵣ = uᵣ.ϱe

    ϱ  = (ϱₗ + ϱᵣ) / 2
    vx = (vxₗ + vxᵣ) / 2
    vy = (vyₗ + vyᵣ) / 2
    p  = (pₗ + pᵣ) / 2
    h  = ((ϱeₗ+pₗ)/ϱₗ + (ϱeᵣ+pᵣ)/ϱᵣ) / 2

    fϱ   = ϱ*vy
    fϱvx = ϱ*vx*vy
    fϱvy = ϱ*vy^2 + p
    fϱe  = ϱ*h*vy

    SVector(fϱ, fϱvx, fϱvy, fϱe)
end

@inline function (fnum::PirozzoliFlux)(uₗ::EulerVar3D{T}, uᵣ::EulerVar3D{T},
                                        model::Euler{T,3}, dir::Val{:x}) where T
    @unpack γ = model
    ϱₗ, vxₗ, vyₗ, vzₗ, pₗ = primitive_variables(uₗ, model)
    ϱeₗ = uₗ.ϱe
    ϱᵣ, vxᵣ, vyᵣ, vzᵣ, pᵣ = primitive_variables(uᵣ, model)
    ϱeᵣ = uᵣ.ϱe

    ϱ  = (ϱₗ + ϱᵣ) / 2
    vx = (vxₗ + vxᵣ) / 2
    vy = (vyₗ + vyᵣ) / 2
    vz = (vzₗ + vzᵣ) / 2
    p  = (pₗ + pᵣ) / 2
    h  = ((ϱeₗ+pₗ)/ϱₗ + (ϱeᵣ+pᵣ)/ϱᵣ) / 2

    fϱ   = ϱ*vx
    fϱvx = ϱ*vx^2 + p
    fϱvy = ϱ*vx*vy
    fϱvz = ϱ*vx*vz
    fϱe  = ϱ*h*vx

    SVector(fϱ, fϱvx, fϱvy, fϱvz, fϱe)
end
@inline function (fnum::PirozzoliFlux)(uₗ::EulerVar3D{T}, uᵣ::EulerVar3D{T},
                                        model::Euler{T,3}, dir::Val{:y}) where T
    @unpack γ = model
    ϱₗ, vxₗ, vyₗ, vzₗ, pₗ = primitive_variables(uₗ, model)
    ϱeₗ = uₗ.ϱe
    ϱᵣ, vxᵣ, vyᵣ, vzᵣ, pᵣ = primitive_variables(uᵣ, model)
    ϱeᵣ = uᵣ.ϱe

    ϱ  = (ϱₗ + ϱᵣ) / 2
    vx = (vxₗ + vxᵣ) / 2
    vy = (vyₗ + vyᵣ) / 2
    vz = (vzₗ + vzᵣ) / 2
    p  = (pₗ + pᵣ) / 2
    h  = ((ϱeₗ+pₗ)/ϱₗ + (ϱeᵣ+pᵣ)/ϱᵣ) / 2

    fϱ   = ϱ*vy
    fϱvx = ϱ*vx*vy
    fϱvy = ϱ*vy^2 + p
    fϱvz = ϱ*vy*vz
    fϱe  = ϱ*h*vy

    SVector(fϱ, fϱvx, fϱvy, fϱvz, fϱe)
end
@inline function (fnum::PirozzoliFlux)(uₗ::EulerVar3D{T}, uᵣ::EulerVar3D{T},
                                        model::Euler{T,3}, dir::Val{:z}) where T
    @unpack γ = model
    ϱₗ, vxₗ, vyₗ, vzₗ, pₗ = primitive_variables(uₗ, model)
    ϱeₗ = uₗ.ϱe
    ϱᵣ, vxᵣ, vyᵣ, vzᵣ, pᵣ = primitive_variables(uᵣ, model)
    ϱeᵣ = uᵣ.ϱe

    ϱ  = (ϱₗ + ϱᵣ) / 2
    vx = (vxₗ + vxᵣ) / 2
    vy = (vyₗ + vyᵣ) / 2
    vz = (vzₗ + vzᵣ) / 2
    p  = (pₗ + pᵣ) / 2
    h  = ((ϱeₗ+pₗ)/ϱₗ + (ϱeᵣ+pᵣ)/ϱᵣ) / 2

    fϱ   = ϱ*vy
    fϱvx = ϱ*vx*vz
    fϱvy = ϱ*vy*vz
    fϱvz = ϱ*vz^2 + p
    fϱe  = ϱ*h*vz

    SVector(fϱ, fϱvx, fϱvy, fϱvz, fϱe)
end


################################################################################

struct EulerRiemannSolution{T,T1} <: AbstractRiemannSolution
    prob::RiemannProblem{Euler{T,1},EulerVar1D{T},T1}
    uₘ⁻::EulerVar1D{T} # left middle state
    uₘ⁺::EulerVar1D{T} # left middle state
    σₗ⁻::T # slow speed of left wave
    σₗ⁺::T # fast speed of left wave
    σₘ::T  # speed of the middle (rarefaction) wave
    σᵣ⁻::T # slow speed of second family
    σᵣ⁺::T # fast speed of second family

    function EulerRiemannSolution{T,T1}(prob::RiemannProblem{Euler{T,1},EulerVar1D{T},T1}) where {T,T1}
        uₘ⁻, uₘ⁺, σₗ⁻, σₗ⁺, σₘ, σᵣ⁻, σᵣ⁺ = compute_state_and_speeds(prob.uₗ, prob.uᵣ, prob.model)
        new(prob, uₘ⁻, uₘ⁺, σₗ⁻, σₗ⁺, σₘ, σᵣ⁻, σᵣ⁺)
    end
end

function solve{T,T1}(prob::RiemannProblem{Euler{T,1},EulerVar1D{T},T1})
    EulerRiemannSolution{T,T1}(prob)
end

function minmax_speeds(sol::EulerRiemannSolution)
    sol.σₗ⁻, sol.σᵣ⁺
end

function (sol::EulerRiemannSolution)(ξ::Real)
    @unpack uₘ⁻, uₘ⁺, σₗ⁻, σₗ⁺, σₘ, σᵣ⁻, σᵣ⁺ = sol
    @unpack uₗ, uᵣ, model = sol.prob

    if ξ <= σₗ⁻
        uₗ
    elseif ξ <= σₗ⁺
        rarefaction_wave_1(ξ, uₗ, model)
    elseif ξ <= σₘ
        uₘ⁻
    elseif ξ <= σᵣ⁻
        uₘ⁺
    elseif ξ <= σᵣ⁺
        rarefaction_wave_2(ξ, uᵣ, model)
    else # σᵣ⁺ < ξ
        uᵣ
    end
end

function rarefaction_wave_1{T}(ξ, uₗ::EulerVar1D, model::Euler{T,1})
    @unpack γ = model
    ϱₗ, vₗ, pₗ = primitive_variables(uₗ, model)
    aₗ = sqrt(γ * pₗ / ϱₗ)
    # Toro (2009), Riemann Solvers and Numerical Methods for Fluid Dynamics,
    # equation (4.56)
    ϱ = ϱₗ * ( 2/(γ+1) + (γ-1)/(aₗ*(γ+1))*(vₗ-ξ) )^(2/(γ-1))
    v = (aₗ + (γ-1)/2*vₗ + ξ) * 2/(γ+1)
    p = pₗ * ( 2/(γ+1) + (γ-1)/(aₗ*(γ+1))*(vₗ-ξ) )^(2γ/(γ-1))
    conserved_variables(ϱ, v, p, model)
end

function rarefaction_wave_2{T}(ξ, uᵣ::EulerVar1D, model::Euler{T,1})
    @unpack γ = model
    ϱᵣ, vᵣ, pᵣ = primitive_variables(uᵣ, model)
    aᵣ = sqrt(γ * pᵣ / ϱᵣ)
    # Toro (2009), Riemann Solvers and Numerical Methods for Fluid Dynamics,
    # equation (4.63)
    ϱ = ϱᵣ * ( 2/(γ+1) - (γ-1)/(aᵣ*(γ+1))*(vᵣ-ξ) )^(2/(γ-1))
    v = (-aᵣ + (γ-1)/2*vᵣ + ξ) * 2/(γ+1)
    p = pᵣ * ( 2/(γ+1) - (γ-1)/(aᵣ*(γ+1))*(vᵣ-ξ) )^(2γ/(γ-1))
    conserved_variables(ϱ, v, p, model)
end


function (sol::EulerRiemannSolution)(t::Real, x::Real)
  @unpack x₀, t₀ = sol.prob

  sol((x-x₀)/(t-t₀))
end


function compute_state_and_speeds{T}(uₗ::EulerVar1D{T}, uᵣ::EulerVar1D{T}, model::Euler)
    @unpack γ = model
    ϱₗ, vₗ, pₗ = primitive_variables(uₗ, model)
    aₗ = sqrt(γ * pₗ / ϱₗ)
    ϱᵣ, vᵣ, pᵣ = primitive_variables(uᵣ, model)
    aᵣ = sqrt(γ * pᵣ / ϱᵣ)

    pₘ, vₘ = compute_pₘ_vₘ(ϱₗ, vₗ, pₗ, aₗ, ϱᵣ, vᵣ, pᵣ, aᵣ, γ)

    if pₘ > pₗ
        # left shock wave
        ϱₘ⁻ = ϱₗ * (pₘ/pₗ + (γ-1)/(γ+1)) / (1 + (γ-1)/(γ+1)*pₘ/pₗ)
        σₗ⁻ = σₗ⁺ = vₗ - aₗ*sqrt((γ+1)/2γ * pₘ/pₗ + (γ-1)/2γ)
    else
        # left rarefaction wave
        ϱₘ⁻ = ϱₗ * (pₘ/pₗ)^(1/γ)
        aₘ⁻ = aₗ * (pₘ/pₗ)^((γ-1)/2γ)
        σₗ⁻ = vₗ - aₗ
        σₗ⁺ = vₘ - aₘ⁻
    end

    if pₘ > pᵣ
        # right shock wave
        ϱₘ⁺ = ϱᵣ * (pₘ/pᵣ + (γ-1)/(γ+1)) / (1 + (γ-1)/(γ+1)*pₘ/pᵣ)
        σᵣ⁻ = σᵣ⁺ = vᵣ + aᵣ*sqrt((γ+1)/2γ * pₘ/pᵣ + (γ-1)/2γ)
    else
        # right rarefaction wave
        ϱₘ⁺ = ϱᵣ * (pₘ/pᵣ)^(1/γ)
        aₘ⁺ = aᵣ * (pₘ/pᵣ)^((γ-1)/2γ)
        σᵣ⁻ = vₘ + aₘ⁺
        σᵣ⁺ = vₗ + aₗ
    end

    uₘ⁻ = conserved_variables(ϱₘ⁻, vₘ, pₘ, model)
    uₘ⁺ = conserved_variables(ϱₘ⁺, vₘ, pₘ, model)
    σₘ = vₘ

    uₘ⁻, uₘ⁺, σₗ⁻, σₗ⁺, σₘ, σᵣ⁻, σᵣ⁺
end


function compute_pₘ_vₘ(ϱₗ, vₗ, pₗ, aₗ, ϱᵣ, vᵣ, pᵣ, aᵣ, γ)
    Aₗ = 2/((γ+1)*ϱₗ)
    Bₗ = (γ-1)/(γ+1)*pₗ
    Aᵣ =  2/((γ+1)*ϱᵣ)
    Bᵣ = (γ-1)/(γ+1)*pᵣ

    fₗ = p -> p > pₗ ? (p-pₗ)*sqrt(Aₗ/(p+Bₗ)) : 2aₗ/(γ-1) * ((p/pₗ)^((γ-1)/2γ) - 1)
    fᵣ = p -> p > pᵣ ? (p-pᵣ)*sqrt(Aᵣ/(p+Bᵣ)) : 2aᵣ/(γ-1) * ((p/pᵣ)^((γ-1)/2γ) - 1)

    # two rarefaction approximation
    # Toro (2009), Riemann Solvers and Numerical Methods for Fluid Dynamics,
    # equation (4.46)
    p₀ = ( (aₗ + aᵣ - (γ-1)*(vᵣ-vₗ)/2) / (aₗ/pₗ^((γ-1)/2γ) + aᵣ/pᵣ^((γ-1)/2γ) ) )^(2γ/(γ-1))
    # TODO: Improve; Newton with explicit derivative?
    pₘ = Roots.fzero(p -> fₗ(p)+fᵣ(p)+vᵣ-vₗ, p₀)

    pₘ, (vₗ+vᵣ)/2 + (fᵣ(pₘ) - fₗ(pₘ))/2
end
