
function reconstruct!(edge_u, u, balance_law, meshx::UniformPeriodicMesh1D,
                      reconstruction!, stencil_width::Val{1}, parallel)
    # no stencils with overlap; set mean values as reconstructions
    @inbounds for cell in 1:length(u)
        edge_u[1,cell] = edge_u[2,cell] = u[cell]
    end
    nothing
end

function reconstruct!(edge_u, u, balance_law, meshx::UniformPeriodicMesh1D,
                      reconstruction!, stencil_width::Val{3}, parallel)
    # stencils with overlap
    @inbounds begin
        cell = 1
        u_m1 = u[end]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        reconstruction!(edge_u, cell, u_m1, u_0, u_p1, balance_law, meshx)

        cell = length(u)
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[1]
        reconstruction!(edge_u, cell, u_m1, u_0, u_p1, balance_law, meshx)
    end

    @inbounds for cell in 2:length(u)-1
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        reconstruction!(edge_u, cell, u_m1, u_0, u_p1, balance_law, meshx)
    end
    nothing
end

function reconstruct!(edge_u, u, balance_law, meshx::UniformPeriodicMesh1D,
                      reconstruction!, stencil_width::Val{5}, parallel)
    # stencils with overlap
    @inbounds begin
        cell = 1
        u_m2 = u[end-1]
        u_m1 = u[end]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[cell+2]
        reconstruction!(edge_u, cell, u_m2, u_m1, u_0, u_p1, u_p2, balance_law, meshx)

        cell = 2
        u_m2 = u[end]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[cell+2]
        reconstruction!(edge_u, cell, u_m2, u_m1, u_0, u_p1, u_p2, balance_law, meshx)

        cell = length(u)-1
        u_m2 = u[cell-2]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[1]
        reconstruction!(edge_u, cell, u_m2, u_m1, u_0, u_p1, u_p2, balance_law, meshx)

        cell = length(u)
        u_m2 = u[cell-2]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[1]
        u_p2 = u[2]
        reconstruction!(edge_u, cell, u_m2, u_m1, u_0, u_p1, u_p2, balance_law, meshx)
    end

    @inbounds for cell in 2:length(u)-1
        u_m2 = u[cell-2]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[cell+2]
        reconstruction!(edge_u, cell, u_m2, u_m1, u_0, u_p1, u_p2, balance_law, meshx)
    end
    nothing
end

function reconstruct!(edge_u, u, balance_law, meshx::UniformPeriodicMesh1D,
                      reconstruction!, stencil_width::Val{7}, parallel)
    # stencils with overlap
    @inbounds begin
        cell = 1
        u_m3 = u[end-2]
        u_m2 = u[end-1]
        u_m1 = u[end]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[cell+2]
        u_p3 = u[cell+3]
        reconstruction!(edge_u, cell, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, balance_law, meshx)

        cell = 2
        u_m3 = u[end-1]
        u_m2 = u[end]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[cell+2]
        u_p3 = u[cell+3]
        reconstruction!(edge_u, cell, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, balance_law, meshx)

        cell = 3
        u_m3 = u[end]
        u_m2 = u[cell-2]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[cell+2]
        u_p3 = u[cell+3]
        reconstruction!(edge_u, cell, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, balance_law, meshx)

        cell = length(u)-2
        u_m3 = u[cell-3]
        u_m2 = u[cell-2]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[cell+2]
        u_p3 = u[1]
        reconstruction!(edge_u, cell, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, balance_law, meshx)

        cell = length(u)-1
        u_m3 = u[cell-3]
        u_m2 = u[cell-2]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[1]
        u_p3 = u[2]
        reconstruction!(edge_u, cell, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, balance_law, meshx)

        cell = length(u)
        u_m3 = u[cell-3]
        u_m2 = u[cell-2]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[1]
        u_p2 = u[2]
        u_p3 = u[3]
        reconstruction!(edge_u, cell, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, balance_law, meshx)
    end

    @inbounds for cell in 2:length(u)-1
        u_m3 = u[cell-3]
        u_m2 = u[cell-2]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[cell+2]
        u_p3 = u[cell+3]
        reconstruction!(edge_u, cell, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, balance_law, meshx)
    end
    nothing
end
