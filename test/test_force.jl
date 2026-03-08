@testset "Physics: Force Interpolation (CIC)" begin
    let
        # ---------------------------------------------------------
        # 1. Zero Self-Force
        # ---------------------------------------------------------
        # A particle sitting alone in a void should not exert a net force on itself.
        sim, C, particles = create_test_state(32)
        fill!(particles.m, 0.0)
        
        # Place one massive particle precisely in the center of a cell
        particles.m[1] = 100.0
        particles.x[1] = 16.5
        particles.y[1] = 16.5
        particles.z[1] = 16.5
        
        assign_density!(particles, sim)
        solve_poisson!(sim, C, 1.0)
        calculate_forces!(particles, sim)
        
        # Due to the symmetric CIC weights and central-difference derivatives, 
        # the interpolated force pulling on the particle should perfectly cancel out.
        @test particles.fx[1] ≈ 0.0 atol=1e-10
        @test particles.fy[1] ≈ 0.0 atol=1e-10
        @test particles.fz[1] ≈ 0.0 atol=1e-10

        # ---------------------------------------------------------
        # 2. Directional Gravity (Newtonian Attraction)
        # ---------------------------------------------------------
        sim2, C2, particles2 = create_test_state(32)
        fill!(particles2.m, 0.0)
        
        # Particle 1: The massive central body
        particles2.m[1] = 1000.0
        particles2.x[1] = 16.5
        particles2.y[1] = 16.5
        particles2.z[1] = 16.5
        
        # Particle 2: A light "test" particle sitting directly to the right (+x axis)
        particles2.m[2] = 1.0
        particles2.x[2] = 20.5
        particles2.y[2] = 16.5
        particles2.z[2] = 16.5
        
        assign_density!(particles2, sim2)
        solve_poisson!(sim2, C2, 1.0)
        calculate_forces!(particles2, sim2)
        
        # The force on Particle 2 must pull it LEFT (negative X direction) toward the massive body
        @test particles2.fx[2] < 0.0
        
        # Because it is perfectly aligned on the Y and Z axes, it should feel zero lateral force
        @test particles2.fy[2] ≈ 0.0 atol=1e-10
        @test particles2.fz[2] ≈ 0.0 atol=1e-10
        
        # Action-Reaction: Particle 1 should be pulled RIGHT (+x direction) by Particle 2
        @test particles2.fx[1] > 0.0
    end
end