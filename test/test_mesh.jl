@testset "Mass Assignment (CIC)" begin
    let
        # ---------------------------------------------------------
        # 1. Total Mass Conservation
        # ---------------------------------------------------------
        sim, C, particles = create_test_state(32)
        
        # Give every particle a known mass
        for i in eachindex(particles)
            particles.m[i] = 2.0
        end
        
        assign_density!(particles, sim)
        
        # The sum of the density grid MUST equal the total mass of all particles
        total_particle_mass = sum(particles.m)
        total_grid_mass = sum(sim.density_grid)
        
        @test total_grid_mass ≈ total_particle_mass

        # ---------------------------------------------------------
        # 2. Perfect Symmetry (Fractional Assignment)
        # ---------------------------------------------------------
        sim2, C2, particles2 = create_test_state(32)
        fill!(particles2.m, 0.0) # Clear all random masses
        
        # Place a single 1.0 mass particle precisely at the center of the first cell
        particles2.m[1] = 1.0
        particles2.x[1] = 1.5
        particles2.y[1] = 1.5
        particles2.z[1] = 1.5
        
        assign_density!(particles2, sim2)
        
        # Mass should be perfectly split (0.125) among the 8 surrounding vertices
        @test sim2.density_grid[1, 1, 1] ≈ 0.125
        @test sim2.density_grid[2, 1, 1] ≈ 0.125
        @test sim2.density_grid[1, 2, 1] ≈ 0.125
        @test sim2.density_grid[2, 2, 1] ≈ 0.125
        @test sim2.density_grid[1, 1, 2] ≈ 0.125
        @test sim2.density_grid[2, 1, 2] ≈ 0.125
        @test sim2.density_grid[1, 2, 2] ≈ 0.125
        @test sim2.density_grid[2, 2, 2] ≈ 0.125
        
        # Ensure no ghost mass was created anywhere else on the grid
        @test sum(sim2.density_grid) ≈ 1.0

        # ---------------------------------------------------------
        # 3. Periodic Boundary Wrapping
        # ---------------------------------------------------------
        sim3, C3, particles3 = create_test_state(32)
        fill!(particles3.m, 0.0)
        
        # Place a single particle pushing against the right x-boundary
        particles3.m[1] = 1.0
        particles3.x[1] = 32.75
        particles3.y[1] = 1.0
        particles3.z[1] = 1.0
        
        assign_density!(particles3, sim3)
        
        # Because y and z are exactly 1.0, dr_y and dr_z are 0.0, so the mass stays on the y=1, z=1 plane.
        # However, at x=32.75, 25% of the mass stays at index 32, and 75% wraps around to index 1.
        @test sim3.density_grid[32, 1, 1] ≈ 0.25
        @test sim3.density_grid[1, 1, 1] ≈ 0.75
        
        # Ensure total mass is still exactly 1.0
        @test sum(sim3.density_grid) ≈ 1.0
    end
end