using PM3
using Test
using StructArrays

# Helper function to generate a fresh state for any test
function create_test_state(N_cells=32)
    N_parts = N_cells^3
    
    sim = Simulation(1.0, N_cells) 
    
    C = ΛCDM() 

    particles = StructArray{Particle}(undef, N_parts)
    for i in 1:N_parts
        particles.m[i] = 1.0
        particles.x[i] = rand() * N_cells
        particles.y[i] = rand() * N_cells
        particles.z[i] = rand() * N_cells
        particles.vx[i] = 0.0
        particles.vy[i] = 0.0
        particles.vz[i] = 0.0
        particles.fx[i] = 0.0 
        particles.fy[i] = 0.0
        particles.fz[i] = 0.0
    end
    
    return sim, C, particles
end

@testset "PM3.jl Full Suite" begin
    @testset "Memory Allocation" begin
        include("memory_allocation.jl")
    end
    @testset "Type Stability" begin
        include("type_stability.jl")
    end
    @testset "Physics" begin
        include("test_mesh.jl")
        include("test_poisson.jl")
        include("test_force.jl")
        include("test_leapfrog.jl")
    end
    
    
    # Future test files go here:
    # @testset "Cosmology" begin include("test_cosmology.jl") end
    # @testset "Poisson Solver" begin include("test_poisson.jl") end
end