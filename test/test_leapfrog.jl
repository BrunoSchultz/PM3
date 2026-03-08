@testset "Leapfrog Integrator" begin
    let
        sim, C, particles = create_test_state(32)
        a = 1.0
        da = 0.01
        
        # Manually calculate the expected cosmological factors based on your code
        a_half = a + da / 2.0
        expected_kick = da * f(a, C)
        expected_drift = da / (a_half^2) * f(a_half, C)

        # ---------------------------------------------------------
        # 1. Coasting (Zero Force)
        # ---------------------------------------------------------
        particles.x[1] = 16.0
        particles.vx[1] = 5.0
        particles.fx[1] = 0.0 # No force
        
        step_leapfrog!(particles, sim, C, a, da)
        
        # Velocity must be completely unchanged
        @test particles.vx[1] == 5.0
        # Position must drift exactly by v * drift_factor
        @test particles.x[1] ≈ 16.0 + (5.0 * expected_drift)

        # ---------------------------------------------------------
        # 2. Constant Acceleration
        # ---------------------------------------------------------
        particles.x[2] = 16.0
        particles.vx[2] = 0.0 # Starts stationary
        particles.fx[2] = 10.0 # Apply a strong constant force
        
        step_leapfrog!(particles, sim, C, a, da)
        
        # Velocity must increase by exactly F * kick_factor
        expected_v2 = 10.0 * expected_kick
        @test particles.vx[2] ≈ expected_v2
        
        # Because the code Kicks THEN Drifts, the position uses the NEW velocity
        @test particles.x[2] ≈ 16.0 + (expected_v2 * expected_drift)

        # ---------------------------------------------------------
        # 3. Periodic Boundary Wrapping
        # ---------------------------------------------------------
        # Test A: Pushing past the upper boundary (x > 32)
        particles.x[3] = 31.5
        # Set a velocity that guarantees it moves exactly +2.0 units in position
        particles.vx[3] = 2.0 / expected_drift 
        particles.fx[3] = 0.0
        
        # Test B: Pushing past the lower boundary (x < 0)
        particles.x[4] = 1.5
        # Set a velocity that guarantees it moves exactly -2.0 units in position
        particles.vx[4] = -2.0 / expected_drift 
        particles.fx[4] = 0.0
        
        step_leapfrog!(particles, sim, C, a, da)
        
        # Particle 3 started at 31.5 and moved +2.0 to 33.5. 
        # mod1(33.5, 32) should wrap to exactly 1.5.
        @test particles.x[3] ≈ 1.5
        
        # Particle 4 started at 1.5 and moved -2.0 to -0.5. 
        # mod1(-0.5, 32) should wrap to exactly 31.5.
        @test particles.x[4] ≈ 31.5
    end
end