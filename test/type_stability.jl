# type_stability

@testset "Type Stability & Math Checks" begin
    let
        sim, C, particles = create_test_state(32)
        a = 1.0
        
        # Checks stability AND verifies the physical result
        @test (@inferred Hubble(a, C)) ≈ 67.4
        @test (@inferred f(a, C)) ≈ 1.0
        
        # Dt(C) evaluates to roughly 0.787 for the default ΛCDM parameters
        @test (@inferred Dt(C)) ≈ 0.7870086233060115 
    end
end

@testset "Type Stability: Particles and Core Functions" begin
    let
        sim, C, particles = create_test_state(32)
        a = 1.0
        da = 0.01

        # Test the Particle Constructor
        pos = [1.5, 2.5, 3.5]
        vel = [0.1, -0.2, 0.3]
        mass = 1.0
    
        p = @inferred Particle(pos, vel, mass) 
        
        # Verify the types locked in correctly
        @test p isa Particle
        @test p.x === 1.5
        @test p.vx === 0.1
        @test p.m === 1.0

        @test (@inferred assign_density!(particles, sim)) === nothing
        @test (@inferred solve_poisson!(sim, C, a)) === nothing
        @test (@inferred calculate_forces!(particles, sim)) === nothing
        @test (@inferred step_leapfrog!(particles, sim, C, a, da)) === nothing
    end
end