using PM3
using BenchmarkTools
using StructArrays

# A helper to run exactly one complete time-step
function run_full_step!(particles, sim, C, a, da)
    assign_density!(particles, sim)
    solve_poisson!(sim, C, a)
    calculate_forces!(particles, sim)
    step_leapfrog!(particles, sim, C, a, da)
    return nothing
end

# Setup a physically realistic grid size for benchmarking
const N_cells = 64
const N_parts = N_cells^3
const a = 1.0
const da = 0.01

println("Setting up benchmark state for N_cells = $N_cells, N_parts = $N_parts...")
sim = Simulation(100.0, N_cells)
C = ΛCDM()
particles = StructArray{Particle}(undef, N_parts)

# Initialize particles (simplified for benchmarking)
for i in 1:N_parts
    particles.m[i] = 1.0
    particles.x[i] = rand() * N_cells
    particles.y[i] = rand() * N_cells
    particles.z[i] = rand() * N_cells
    particles.vx[i] = 0.0; particles.vy[i] = 0.0; particles.vz[i] = 0.0
    particles.fx[i] = 0.0; particles.fy[i] = 0.0; particles.fz[i] = 0.0
end

println("Warming up functions...")
run_full_step!(particles, sim, C, a, da)

# ---------------------------------------------------------
# The Benchmarks
# ---------------------------------------------------------

# ---------------------------------------------------------
# The Benchmarks
# ---------------------------------------------------------

println("\n--- Benchmarking Individual Components ---")

print("1. Mass Assignment: ")
@btime assign_density!($particles, $sim)

# Poisson solver modifies the grid in-place, so we give it a fresh grid copy each time
print("2. Poisson Solver: ")
@btime solve_poisson!(s_copy, $C, $a) setup=(s_copy = deepcopy($sim)) evals=1

print("3. Force Interpolation: ")
@btime calculate_forces!($particles, $sim)

# Leapfrog mutates particles, so we give it a fresh particle copy each time
print("4. Leapfrog Step: ")
@btime step_leapfrog!(p_copy, $sim, $C, $a, $da) setup=(p_copy = deepcopy($particles)) evals=1

println("\n--- Benchmarking Full Time-Step ---")
# The full step mutates both, so we copy both!
@btime run_full_step!(p_copy, s_copy, $C, $a, $da) setup=(p_copy = deepcopy($particles); s_copy = deepcopy($sim)) evals=1