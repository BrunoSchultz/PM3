@testset "Poisson Solver" begin
    let
        # ---------------------------------------------------------
        # 1. DC Mode and Zero-Mean Potential
        # ---------------------------------------------------------
        sim, C, particles = create_test_state(32)
        a = 1.0
        
        # Give particles some uniform mass
        fill!(particles.m, 1.0)
        assign_density!(particles, sim)
        
        solve_poisson!(sim, C, a)
        
    
        # Because the DC mode is zero, the sum of the potential grid must be ~0.0
        # (Using atol because summing a 32^3 grid of floats can leave a microscopic residual)
        total_potential = sum(sim.density_grid)
        @test total_potential ≈ 0.0 atol=1e-10

        # ---------------------------------------------------------
        # 2. Attractive Gravity (Potential Well)
        # ---------------------------------------------------------
        sim2, C2, particles2 = create_test_state(32)
        fill!(particles2.m, 0.0) 
        
        # Place a single massive particle in the middle of the grid
        center_idx = 16
        particles2.m[1] = 100.0
        particles2.x[1] = center_idx + 0.5
        particles2.y[1] = center_idx + 0.5
        particles2.z[1] = center_idx + 0.5
        
        assign_density!(particles2, sim2)
        solve_poisson!(sim2, C2, a)
        
        # The potential at the particle's location must be strongly negative
        core_potential = sim2.density_grid[center_idx, center_idx, center_idx]
        @test core_potential < 0.0
        
        # It should also be the global minimum of the entire grid
        @test core_potential == minimum(sim2.density_grid)

        # ---------------------------------------------------------
        # 3. Linearity of the Poisson Equation
        # ---------------------------------------------------------
        sim3, C3, particles3 = create_test_state(32)
        fill!(particles3.m, 0.0)
        
        # Same particle, but exactly double the mass
        particles3.m[1] = 200.0
        particles3.x[1] = center_idx + 0.5
        particles3.y[1] = center_idx + 0.5
        particles3.z[1] = center_idx + 0.5
        
        assign_density!(particles3, sim3)
        solve_poisson!(sim3, C3, a)
        
        # The potential well should be exactly twice as deep
        doubled_core_potential = sim3.density_grid[center_idx, center_idx, center_idx]
        @test doubled_core_potential ≈ (2.0 * core_potential)
    end
end