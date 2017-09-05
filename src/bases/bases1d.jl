
function map_to_canonical(x, xmin, xmax, ::AbstractBasis)
    ( 2x - (xmax + xmin) ) / (xmax - xmin)
end

function map_from_canonical(ξ, xmin, xmax, ::AbstractBasis)
    ( (xmax + xmin) + ξ * (xmax - xmin) ) / 2
end


struct LobattoLegendre{T<:Real} <: NodalBasis
    nodes::Vector{T}
    weights::Vector{T}
    D::Array{T,2}

    function LobattoLegendre{T}(p::Int) where T
        assert(p >= 0)
        if p == 0
            q = Jacobi.Quadrature(Jacobi.GJ, p+1, 0, 0, T)
            nodes = Jacobi.qzeros(q)
            weights = Jacobi.qweights(q)
            D = Jacobi.qdiff(q)
        else
            q = Jacobi.Quadrature(Jacobi.GLJ, p+1, 0, 0, T)
            nodes = Jacobi.qzeros(q)
            weights = Jacobi.qweights(q)
            D = Jacobi.qdiff(q)
        end

        new(nodes, weights, D)
    end
end

LobattoLegendre(p::Int, T=Float64) = LobattoLegendre{T}(p)
