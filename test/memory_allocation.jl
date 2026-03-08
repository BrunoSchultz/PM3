# memory_allocation.jl
# Note: Test and BenchmarkTools are already loaded by runtests.jl, but it's okay to keep them if you run this file directly.

let
    sim, C, particles = create_test_state(32)
    a = 1.0
    da = 0.01

    # Warm-up to force JIT compilation
    assign_density!(particles, sim)
    solve_poisson!(sim, C, a) 
    calculate_forces!(particles, sim)
    step_leapfrog!(particles, sim, C, a, da)

    @test (@allocated assign_density!(particles, sim)) == 0
    @test (@allocated solve_poisson!(sim, C, a)) <= 16
    @test (@allocated calculate_forces!(particles, sim)) == 0
    @test (@allocated step_leapfrog!(particles, sim, C, a, da)) == 0
end