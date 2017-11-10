
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

    reconstruct_inner_loop!(edge_u, u, balance_law, meshx, reconstruction!, stencil_width, parallel)
end

@inline function reconstruct_inner_loop!(edge_u, u, balance_law, meshx::UniformPeriodicMesh1D,
                                         reconstruction!, stencil_width::Val{3}, parallel)
    # serial
    @inbounds for cell in 2:length(u)-1
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        reconstruction!(edge_u, cell, u_m1, u_0, u_p1, balance_law, meshx)
    end
    nothing
end

@inline function reconstruct_inner_loop!(edge_u, u, balance_law, meshx::UniformPeriodicMesh1D,
                                         reconstruction!, stencil_width::Val{3}, ::Val{:threads})
    # threaded
    @inbounds Threads.@threads for cell in 2:length(u)-1
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

    reconstruct_inner_loop!(edge_u, u, balance_law, meshx, reconstruction!, stencil_width, parallel)
end

@inline function reconstruct_inner_loop!(edge_u, u, balance_law, meshx::UniformPeriodicMesh1D,
                                         reconstruction!, stencil_width::Val{5}, parallel)
    # serial
    @inbounds for cell in 3:length(u)-2
        u_m2 = u[cell-2]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[cell+2]
        reconstruction!(edge_u, cell, u_m2, u_m1, u_0, u_p1, u_p2, balance_law, meshx)
    end
    nothing
end

@inline function reconstruct_inner_loop!(edge_u, u, balance_law, meshx::UniformPeriodicMesh1D,
                                         reconstruction!, stencil_width::Val{5}, ::Val{:threads})
    # threaded
    @inbounds Threads.@threads for cell in 3:length(u)-2
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

    reconstruct_inner_loop!(edge_u, u, balance_law, meshx, reconstruction!, stencil_width, parallel)
end

@inline function reconstruct_inner_loop!(edge_u, u, balance_law, meshx::UniformPeriodicMesh1D,
                                         reconstruction!, stencil_width::Val{7}, parallel)
    # serial
    @inbounds for cell in 4:length(u)-3
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

@inline function reconstruct_inner_loop!(edge_u, u, balance_law, meshx::UniformPeriodicMesh1D,
                                         reconstruction!, stencil_width::Val{7}, ::Val{:threads})
    # threaded
    @inbounds Threads.@threads for cell in 4:length(u)-3
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


function reconstruct!(edge_u, u, balance_law, meshx::UniformPeriodicMesh1D,
                      reconstruction!, stencil_width::Val{9}, parallel)
    # stencils with overlap
    @inbounds begin
        cell = 1
        u_m4 = u[end-3]
        u_m3 = u[end-2]
        u_m2 = u[end-1]
        u_m1 = u[end]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[cell+2]
        u_p3 = u[cell+3]
        u_p4 = u[cell+4]
        reconstruction!(edge_u, cell, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, balance_law, meshx)

        cell = 2
        u_m4 = u[end-2]
        u_m3 = u[end-1]
        u_m2 = u[end]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[cell+2]
        u_p3 = u[cell+3]
        u_p4 = u[cell+4]
        reconstruction!(edge_u, cell, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, balance_law, meshx)

        cell = 3
        u_m4 = u[end-1]
        u_m3 = u[end]
        u_m2 = u[cell-2]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[cell+2]
        u_p3 = u[cell+3]
        u_p4 = u[cell+4]
        reconstruction!(edge_u, cell, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, balance_law, meshx)

        cell = 4
        u_m4 = u[end]
        u_m3 = u[cell-3]
        u_m2 = u[cell-2]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[cell+2]
        u_p3 = u[cell+3]
        u_p4 = u[cell+4]
        reconstruction!(edge_u, cell, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, balance_law, meshx)

        cell = length(u)-3
        u_m4 = u[cell-4]
        u_m3 = u[cell-3]
        u_m2 = u[cell-2]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[cell+2]
        u_p3 = u[cell+3]
        u_p4 = u[1]
        reconstruction!(edge_u, cell, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, balance_law, meshx)

        cell = length(u)-2
        u_m4 = u[cell-4]
        u_m3 = u[cell-3]
        u_m2 = u[cell-2]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[cell+2]
        u_p3 = u[1]
        u_p4 = u[2]
        reconstruction!(edge_u, cell, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, balance_law, meshx)

        cell = length(u)-1
        u_m4 = u[cell-4]
        u_m3 = u[cell-3]
        u_m2 = u[cell-2]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[1]
        u_p3 = u[2]
        u_p4 = u[3]
        reconstruction!(edge_u, cell, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, balance_law, meshx)

        cell = length(u)
        u_m4 = u[cell-4]
        u_m3 = u[cell-3]
        u_m2 = u[cell-2]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[1]
        u_p2 = u[2]
        u_p3 = u[3]
        u_p4 = u[4]
        reconstruction!(edge_u, cell, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, balance_law, meshx)
    end

    reconstruct_inner_loop!(edge_u, u, balance_law, meshx, reconstruction!, stencil_width, parallel)
end

@inline function reconstruct_inner_loop!(edge_u, u, balance_law, meshx::UniformPeriodicMesh1D,
                                         reconstruction!, stencil_width::Val{9}, parallel)
    # serial
    @inbounds for cell in 5:length(u)-4
        u_m4 = u[cell-4]
        u_m3 = u[cell-3]
        u_m2 = u[cell-2]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[cell+2]
        u_p3 = u[cell+3]
        u_p4 = u[cell+4]
        reconstruction!(edge_u, cell, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, balance_law, meshx)
    end
    nothing
end

@inline function reconstruct_inner_loop!(edge_u, u, balance_law, meshx::UniformPeriodicMesh1D,
                                         reconstruction!, stencil_width::Val{9}, ::Val{:threads})
    # threaded
    @inbounds Threads.@threads for cell in 5:length(u)-4
        u_m4 = u[cell-4]
        u_m3 = u[cell-3]
        u_m2 = u[cell-2]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[cell+2]
        u_p3 = u[cell+3]
        u_p4 = u[cell+4]
        reconstruction!(edge_u, cell, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, balance_law, meshx)
    end
    nothing
end


function reconstruct!(edge_u, u, balance_law, meshx::UniformPeriodicMesh1D,
                      reconstruction!, stencil_width::Val{11}, parallel)
    # stencils with overlap
    @inbounds begin
        cell = 1
        u_m5 = u[end-4]
        u_m4 = u[end-3]
        u_m3 = u[end-2]
        u_m2 = u[end-1]
        u_m1 = u[end]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[cell+2]
        u_p3 = u[cell+3]
        u_p4 = u[cell+4]
        u_p5 = u[cell+5]
        reconstruction!(edge_u, cell, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, balance_law, meshx)

        cell = 2
        u_m5 = u[end-3]
        u_m4 = u[end-2]
        u_m3 = u[end-1]
        u_m2 = u[end]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[cell+2]
        u_p3 = u[cell+3]
        u_p4 = u[cell+4]
        u_p5 = u[cell+5]
        reconstruction!(edge_u, cell, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, balance_law, meshx)

        cell = 3
        u_m5 = u[end-2]
        u_m4 = u[end-1]
        u_m3 = u[end]
        u_m2 = u[cell-2]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[cell+2]
        u_p3 = u[cell+3]
        u_p4 = u[cell+4]
        u_p5 = u[cell+5]
        reconstruction!(edge_u, cell, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, balance_law, meshx)

        cell = 4
        u_m5 = u[end-1]
        u_m4 = u[end]
        u_m3 = u[cell-3]
        u_m2 = u[cell-2]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[cell+2]
        u_p3 = u[cell+3]
        u_p4 = u[cell+4]
        u_p5 = u[cell+5]
        reconstruction!(edge_u, cell, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, balance_law, meshx)

        cell = 5
        u_m5 = u[end]
        u_m4 = u[cell-4]
        u_m3 = u[cell-3]
        u_m2 = u[cell-2]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[cell+2]
        u_p3 = u[cell+3]
        u_p4 = u[cell+4]
        u_p5 = u[cell+5]
        reconstruction!(edge_u, cell, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, balance_law, meshx)

        cell = length(u)-4
        u_m5 = u[cell-5]
        u_m4 = u[cell-4]
        u_m3 = u[cell-3]
        u_m2 = u[cell-2]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[cell+2]
        u_p3 = u[cell+3]
        u_p4 = u[cell+4]
        u_p5 = u[1]
        reconstruction!(edge_u, cell, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, balance_law, meshx)

        cell = length(u)-3
        u_m5 = u[cell-5]
        u_m4 = u[cell-4]
        u_m3 = u[cell-3]
        u_m2 = u[cell-2]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[cell+2]
        u_p3 = u[cell+3]
        u_p4 = u[1]
        u_p5 = u[2]
        reconstruction!(edge_u, cell, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, balance_law, meshx)

        cell = length(u)-2
        u_m5 = u[cell-5]
        u_m4 = u[cell-4]
        u_m3 = u[cell-3]
        u_m2 = u[cell-2]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[cell+2]
        u_p3 = u[1]
        u_p4 = u[2]
        u_p5 = u[3]
        reconstruction!(edge_u, cell, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, balance_law, meshx)

        cell = length(u)-1
        u_m5 = u[cell-5]
        u_m4 = u[cell-4]
        u_m3 = u[cell-3]
        u_m2 = u[cell-2]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[1]
        u_p3 = u[2]
        u_p4 = u[3]
        u_p5 = u[4]
        reconstruction!(edge_u, cell, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, balance_law, meshx)

        cell = length(u)
        u_m5 = u[cell-5]
        u_m4 = u[cell-4]
        u_m3 = u[cell-3]
        u_m2 = u[cell-2]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[1]
        u_p2 = u[2]
        u_p3 = u[3]
        u_p4 = u[4]
        u_p5 = u[5]
        reconstruction!(edge_u, cell, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, balance_law, meshx)
    end

    reconstruct_inner_loop!(edge_u, u, balance_law, meshx, reconstruction!, stencil_width, parallel)
end

@inline function reconstruct_inner_loop!(edge_u, u, balance_law, meshx::UniformPeriodicMesh1D,
                                         reconstruction!, stencil_width::Val{11}, parallel)
    # serial
    @inbounds for cell in 6:length(u)-5
        u_m5 = u[cell-5]
        u_m4 = u[cell-4]
        u_m3 = u[cell-3]
        u_m2 = u[cell-2]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[cell+2]
        u_p3 = u[cell+3]
        u_p4 = u[cell+4]
        u_p5 = u[cell+5]
        reconstruction!(edge_u, cell, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, balance_law, meshx)
    end
    nothing
end

@inline function reconstruct_inner_loop!(edge_u, u, balance_law, meshx::UniformPeriodicMesh1D,
                                         reconstruction!, stencil_width::Val{11}, ::Val{:threads})
    # threaded
    @inbounds Threads.@threads for cell in 6:length(u)-5
        u_m5 = u[cell-5]
        u_m4 = u[cell-4]
        u_m3 = u[cell-3]
        u_m2 = u[cell-2]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[cell+2]
        u_p3 = u[cell+3]
        u_p4 = u[cell+4]
        u_p5 = u[cell+5]
        reconstruction!(edge_u, cell, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, balance_law, meshx)
    end
    nothing
end


function reconstruct!(edge_u, u, balance_law, meshx::UniformPeriodicMesh1D,
                      reconstruction!, stencil_width::Val{13}, parallel)
    # stencils with overlap
    @inbounds begin
        cell = 1
        u_m6 = u[end-5]
        u_m5 = u[end-4]
        u_m4 = u[end-3]
        u_m3 = u[end-2]
        u_m2 = u[end-1]
        u_m1 = u[end]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[cell+2]
        u_p3 = u[cell+3]
        u_p4 = u[cell+4]
        u_p5 = u[cell+5]
        u_p6 = u[cell+6]
        reconstruction!(edge_u, cell, u_m6, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, u_p6, balance_law, meshx)

        cell = 2
        u_m6 = u[end-4]
        u_m5 = u[end-3]
        u_m4 = u[end-2]
        u_m3 = u[end-1]
        u_m2 = u[end]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[cell+2]
        u_p3 = u[cell+3]
        u_p4 = u[cell+4]
        u_p5 = u[cell+5]
        u_p6 = u[cell+6]
        reconstruction!(edge_u, cell, u_m6, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, u_p6, balance_law, meshx)

        cell = 3
        u_m6 = u[end-3]
        u_m5 = u[end-2]
        u_m4 = u[end-1]
        u_m3 = u[end]
        u_m2 = u[cell-2]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[cell+2]
        u_p3 = u[cell+3]
        u_p4 = u[cell+4]
        u_p5 = u[cell+5]
        u_p6 = u[cell+6]
        reconstruction!(edge_u, cell, u_m6, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, u_p6, balance_law, meshx)

        cell = 4
        u_m6 = u[end-2]
        u_m5 = u[end-1]
        u_m4 = u[end]
        u_m3 = u[cell-3]
        u_m2 = u[cell-2]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[cell+2]
        u_p3 = u[cell+3]
        u_p4 = u[cell+4]
        u_p5 = u[cell+5]
        u_p6 = u[cell+6]
        reconstruction!(edge_u, cell, u_m6, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, u_p6, balance_law, meshx)

        cell = 5
        u_m6 = u[end-1]
        u_m5 = u[end]
        u_m4 = u[cell-4]
        u_m3 = u[cell-3]
        u_m2 = u[cell-2]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[cell+2]
        u_p3 = u[cell+3]
        u_p4 = u[cell+4]
        u_p5 = u[cell+5]
        u_p6 = u[cell+6]
        reconstruction!(edge_u, cell, u_m6, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, u_p6, balance_law, meshx)

        cell = 6
        u_m6 = u[end]
        u_m5 = u[cell-5]
        u_m4 = u[cell-4]
        u_m3 = u[cell-3]
        u_m2 = u[cell-2]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[cell+2]
        u_p3 = u[cell+3]
        u_p4 = u[cell+4]
        u_p5 = u[cell+5]
        u_p6 = u[cell+6]
        reconstruction!(edge_u, cell, u_m6, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, u_p6, balance_law, meshx)

        cell = length(u)-5
        u_m6 = u[cell-6]
        u_m5 = u[cell-5]
        u_m4 = u[cell-4]
        u_m3 = u[cell-3]
        u_m2 = u[cell-2]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[cell+2]
        u_p3 = u[cell+3]
        u_p4 = u[cell+4]
        u_p5 = u[cell+5]
        u_p6 = u[1]
        reconstruction!(edge_u, cell, u_m6, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, u_p6, balance_law, meshx)

        cell = length(u)-4
        u_m6 = u[cell-6]
        u_m5 = u[cell-5]
        u_m4 = u[cell-4]
        u_m3 = u[cell-3]
        u_m2 = u[cell-2]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[cell+2]
        u_p3 = u[cell+3]
        u_p4 = u[cell+4]
        u_p5 = u[1]
        u_p6 = u[2]
        reconstruction!(edge_u, cell, u_m6, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, u_p6, balance_law, meshx)

        cell = length(u)-3
        u_m6 = u[cell-6]
        u_m5 = u[cell-5]
        u_m4 = u[cell-4]
        u_m3 = u[cell-3]
        u_m2 = u[cell-2]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[cell+2]
        u_p3 = u[cell+3]
        u_p4 = u[1]
        u_p5 = u[2]
        u_p6 = u[3]
        reconstruction!(edge_u, cell, u_m6, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, u_p6, balance_law, meshx)

        cell = length(u)-2
        u_m6 = u[cell-6]
        u_m5 = u[cell-5]
        u_m4 = u[cell-4]
        u_m3 = u[cell-3]
        u_m2 = u[cell-2]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[cell+2]
        u_p3 = u[1]
        u_p4 = u[2]
        u_p5 = u[3]
        u_p6 = u[4]
        reconstruction!(edge_u, cell, u_m6, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, u_p6, balance_law, meshx)

        cell = length(u)-1
        u_m6 = u[cell-6]
        u_m5 = u[cell-5]
        u_m4 = u[cell-4]
        u_m3 = u[cell-3]
        u_m2 = u[cell-2]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[1]
        u_p3 = u[2]
        u_p4 = u[3]
        u_p5 = u[4]
        u_p6 = u[5]
        reconstruction!(edge_u, cell, u_m6, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, u_p6, balance_law, meshx)

        cell = length(u)
        u_m6 = u[cell-6]
        u_m5 = u[cell-5]
        u_m4 = u[cell-4]
        u_m3 = u[cell-3]
        u_m2 = u[cell-2]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[1]
        u_p2 = u[2]
        u_p3 = u[3]
        u_p4 = u[4]
        u_p5 = u[5]
        u_p6 = u[6]
        reconstruction!(edge_u, cell, u_m6, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, u_p6, balance_law, meshx)
    end

    reconstruct_inner_loop!(edge_u, u, balance_law, meshx, reconstruction!, stencil_width, parallel)
end

@inline function reconstruct_inner_loop!(edge_u, u, balance_law, meshx::UniformPeriodicMesh1D,
                                         reconstruction!, stencil_width::Val{13}, parallel)
    # serial
    @inbounds for cell in 7:length(u)-6
        u_m6 = u[cell-6]
        u_m5 = u[cell-5]
        u_m4 = u[cell-4]
        u_m3 = u[cell-3]
        u_m2 = u[cell-2]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[cell+2]
        u_p3 = u[cell+3]
        u_p4 = u[cell+4]
        u_p5 = u[cell+5]
        u_p6 = u[cell+6]
        reconstruction!(edge_u, cell, u_m6, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, u_p6, balance_law, meshx)
    end
    nothing
end

@inline function reconstruct_inner_loop!(edge_u, u, balance_law, meshx::UniformPeriodicMesh1D,
                                         reconstruction!, stencil_width::Val{13}, ::Val{:threads})
    # threaded
    @inbounds Threads.@threads for cell in 7:length(u)-6
        u_m6 = u[cell-6]
        u_m5 = u[cell-5]
        u_m4 = u[cell-4]
        u_m3 = u[cell-3]
        u_m2 = u[cell-2]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[cell+2]
        u_p3 = u[cell+3]
        u_p4 = u[cell+4]
        u_p5 = u[cell+5]
        u_p6 = u[cell+6]
        reconstruction!(edge_u, cell, u_m6, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, u_p6, balance_law, meshx)
    end
    nothing
end


function reconstruct!(edge_u, u, balance_law, meshx::UniformPeriodicMesh1D,
                      reconstruction!, stencil_width::Val{15}, parallel)
    # stencils with overlap
    @inbounds begin
        cell = 1
        u_m7 = u[end-6]
        u_m6 = u[end-5]
        u_m5 = u[end-4]
        u_m4 = u[end-3]
        u_m3 = u[end-2]
        u_m2 = u[end-1]
        u_m1 = u[end]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[cell+2]
        u_p3 = u[cell+3]
        u_p4 = u[cell+4]
        u_p5 = u[cell+5]
        u_p6 = u[cell+6]
        u_p7 = u[cell+7]
        reconstruction!(edge_u, cell, u_m7, u_m6, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, u_p6, u_p7, balance_law, meshx)

        cell = 2
        u_m7 = u[end-5]
        u_m6 = u[end-4]
        u_m5 = u[end-3]
        u_m4 = u[end-2]
        u_m3 = u[end-1]
        u_m2 = u[end]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[cell+2]
        u_p3 = u[cell+3]
        u_p4 = u[cell+4]
        u_p5 = u[cell+5]
        u_p6 = u[cell+6]
        u_p7 = u[cell+7]
        reconstruction!(edge_u, cell, u_m7, u_m6, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, u_p6, u_p7, balance_law, meshx)

        cell = 3
        u_m7 = u[end-4]
        u_m6 = u[end-3]
        u_m5 = u[end-2]
        u_m4 = u[end-1]
        u_m3 = u[end]
        u_m2 = u[cell-2]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[cell+2]
        u_p3 = u[cell+3]
        u_p4 = u[cell+4]
        u_p5 = u[cell+5]
        u_p6 = u[cell+6]
        u_p7 = u[cell+7]
        reconstruction!(edge_u, cell, u_m7, u_m6, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, u_p6, u_p7, balance_law, meshx)

        cell = 4
        u_m7 = u[end-3]
        u_m6 = u[end-2]
        u_m5 = u[end-1]
        u_m4 = u[end]
        u_m3 = u[cell-3]
        u_m2 = u[cell-2]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[cell+2]
        u_p3 = u[cell+3]
        u_p4 = u[cell+4]
        u_p5 = u[cell+5]
        u_p6 = u[cell+6]
        u_p7 = u[cell+7]
        reconstruction!(edge_u, cell, u_m7, u_m6, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, u_p6, u_p7, balance_law, meshx)

        cell = 5
        u_m7 = u[end-2]
        u_m6 = u[end-1]
        u_m5 = u[end]
        u_m4 = u[cell-4]
        u_m3 = u[cell-3]
        u_m2 = u[cell-2]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[cell+2]
        u_p3 = u[cell+3]
        u_p4 = u[cell+4]
        u_p5 = u[cell+5]
        u_p6 = u[cell+6]
        u_p7 = u[cell+7]
        reconstruction!(edge_u, cell, u_m7, u_m6, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, u_p6, u_p7, balance_law, meshx)

        cell = 6
        u_m7 = u[end-1]
        u_m6 = u[end]
        u_m5 = u[cell-5]
        u_m4 = u[cell-4]
        u_m3 = u[cell-3]
        u_m2 = u[cell-2]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[cell+2]
        u_p3 = u[cell+3]
        u_p4 = u[cell+4]
        u_p5 = u[cell+5]
        u_p6 = u[cell+6]
        u_p7 = u[cell+7]
        reconstruction!(edge_u, cell, u_m7, u_m6, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, u_p6, u_p7, balance_law, meshx)

        cell = 7
        u_m7 = u[end]
        u_m6 = u[cell-6]
        u_m5 = u[cell-5]
        u_m4 = u[cell-4]
        u_m3 = u[cell-3]
        u_m2 = u[cell-2]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[cell+2]
        u_p3 = u[cell+3]
        u_p4 = u[cell+4]
        u_p5 = u[cell+5]
        u_p6 = u[cell+6]
        u_p7 = u[cell+7]
        reconstruction!(edge_u, cell, u_m7, u_m6, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, u_p6, u_p7, balance_law, meshx)

        cell = length(u)-6
        u_m7 = u[cell-7]
        u_m6 = u[cell-6]
        u_m5 = u[cell-5]
        u_m4 = u[cell-4]
        u_m3 = u[cell-3]
        u_m2 = u[cell-2]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[cell+2]
        u_p3 = u[cell+3]
        u_p4 = u[cell+4]
        u_p5 = u[cell+5]
        u_p6 = u[cell+6]
        u_p7 = u[1]
        reconstruction!(edge_u, cell, u_m7, u_m6, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, u_p6, u_p7, balance_law, meshx)

        cell = length(u)-5
        u_m7 = u[cell-7]
        u_m6 = u[cell-6]
        u_m5 = u[cell-5]
        u_m4 = u[cell-4]
        u_m3 = u[cell-3]
        u_m2 = u[cell-2]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[cell+2]
        u_p3 = u[cell+3]
        u_p4 = u[cell+4]
        u_p5 = u[cell+5]
        u_p6 = u[1]
        u_p7 = u[2]
        reconstruction!(edge_u, cell, u_m7, u_m6, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, u_p6, u_p7, balance_law, meshx)

        cell = length(u)-4
        u_m7 = u[cell-7]
        u_m6 = u[cell-6]
        u_m5 = u[cell-5]
        u_m4 = u[cell-4]
        u_m3 = u[cell-3]
        u_m2 = u[cell-2]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[cell+2]
        u_p3 = u[cell+3]
        u_p4 = u[cell+4]
        u_p5 = u[1]
        u_p6 = u[2]
        u_p7 = u[3]
        reconstruction!(edge_u, cell, u_m7, u_m6, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, u_p6, u_p7, balance_law, meshx)

        cell = length(u)-3
        u_m7 = u[cell-7]
        u_m6 = u[cell-6]
        u_m5 = u[cell-5]
        u_m4 = u[cell-4]
        u_m3 = u[cell-3]
        u_m2 = u[cell-2]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[cell+2]
        u_p3 = u[cell+3]
        u_p4 = u[1]
        u_p5 = u[2]
        u_p6 = u[3]
        u_p7 = u[4]
        reconstruction!(edge_u, cell, u_m7, u_m6, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, u_p6, u_p7, balance_law, meshx)

        cell = length(u)-2
        u_m7 = u[cell-7]
        u_m6 = u[cell-6]
        u_m5 = u[cell-5]
        u_m4 = u[cell-4]
        u_m3 = u[cell-3]
        u_m2 = u[cell-2]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[cell+2]
        u_p3 = u[1]
        u_p4 = u[2]
        u_p5 = u[3]
        u_p6 = u[4]
        u_p7 = u[5]
        reconstruction!(edge_u, cell, u_m7, u_m6, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, u_p6, u_p7, balance_law, meshx)

        cell = length(u)-1
        u_m7 = u[cell-7]
        u_m6 = u[cell-6]
        u_m5 = u[cell-5]
        u_m4 = u[cell-4]
        u_m3 = u[cell-3]
        u_m2 = u[cell-2]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[1]
        u_p3 = u[2]
        u_p4 = u[3]
        u_p5 = u[4]
        u_p6 = u[5]
        u_p7 = u[6]
        reconstruction!(edge_u, cell, u_m7, u_m6, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, u_p6, u_p7, balance_law, meshx)

        cell = length(u)
        u_m7 = u[cell-7]
        u_m6 = u[cell-6]
        u_m5 = u[cell-5]
        u_m4 = u[cell-4]
        u_m3 = u[cell-3]
        u_m2 = u[cell-2]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[1]
        u_p2 = u[2]
        u_p3 = u[3]
        u_p4 = u[4]
        u_p5 = u[5]
        u_p6 = u[6]
        u_p7 = u[7]
        reconstruction!(edge_u, cell, u_m7, u_m6, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, u_p6, u_p7, balance_law, meshx)
    end

    reconstruct_inner_loop!(edge_u, u, balance_law, meshx, reconstruction!, stencil_width, parallel)
end

@inline function reconstruct_inner_loop!(edge_u, u, balance_law, meshx::UniformPeriodicMesh1D,
                                         reconstruction!, stencil_width::Val{15}, parallel)
    # serial
    @inbounds for cell in 8:length(u)-7
        u_m7 = u[cell-7]
        u_m6 = u[cell-6]
        u_m5 = u[cell-5]
        u_m4 = u[cell-4]
        u_m3 = u[cell-3]
        u_m2 = u[cell-2]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[cell+2]
        u_p3 = u[cell+3]
        u_p4 = u[cell+4]
        u_p5 = u[cell+5]
        u_p6 = u[cell+6]
        u_p7 = u[cell+7]
        reconstruction!(edge_u, cell, u_m7, u_m6, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, u_p6, u_p7, balance_law, meshx)
    end
    nothing
end

@inline function reconstruct_inner_loop!(edge_u, u, balance_law, meshx::UniformPeriodicMesh1D,
                                         reconstruction!, stencil_width::Val{15}, ::Val{:threads})
    # threaded
    @inbounds Threads.@threads for cell in 8:length(u)-7
        u_m7 = u[cell-7]
        u_m6 = u[cell-6]
        u_m5 = u[cell-5]
        u_m4 = u[cell-4]
        u_m3 = u[cell-3]
        u_m2 = u[cell-2]
        u_m1 = u[cell-1]
        u_0  = u[cell]
        u_p1 = u[cell+1]
        u_p2 = u[cell+2]
        u_p3 = u[cell+3]
        u_p4 = u[cell+4]
        u_p5 = u[cell+5]
        u_p6 = u[cell+6]
        u_p7 = u[cell+7]
        reconstruction!(edge_u, cell, u_m7, u_m6, u_m5, u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4, u_p5, u_p6, u_p7, balance_law, meshx)
    end
    nothing
end
