using ArgParse
using PM3
using StructArrays
using Serialization

function parse_commandline()
    s = ArgParseSettings(description = "PM3: Particle-Mesh Cosmology Simulation CLI")

    @add_arg_table! s begin
        "--N_cells"
            help = "Number of grid cells in one dimension"
            arg_type = Int
            default = 256
        "--N_parts"
            help = "Total number of particles"
            arg_type = Int
            default = 10_000_000
        "--a_start"
            help = "Starting scale factor"
            arg_type = Float64
            default = 0.01
        "--a_end"
            help = "Ending scale factor"
            arg_type = Float64
            default = 1.0
        "--da"
            help = "Integration step size for scale factor"
            arg_type = Float64
            default = 0.01
        "--out"
            help = "Output file name for the final particle state"
            arg_type = String
            default = "final_state.jls"
    end

    return parse_args(s)
end

function main()
    args = parse_commandline()

    N_cells = args["N_cells"]
    N_parts = args["N_parts"]
    a_start = args["a_start"]
    a_end   = args["a_end"]
    da      = args["da"]
    out     = args["out"]

    println("Initializing PM3 Simulation...")
    sim = Simulation(10000.0, N_cells)
    C = ΛCDM()

    particles = StructArray{Particle}(undef, N_parts)
    @inbounds for i in 1:N_parts
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

    println("Running integration from a = $a_start to $a_end...")
    a = a_start
    while a < a_end
        assign_density!(particles, sim)
        solve_poisson!(sim, C, a)
        calculate_forces!(particles, sim)
        step_leapfrog!(particles, sim, C, a, da)
        
        a += da
    end

    println("Simulation complete. Saving final state to $out...")
    serialize(out, (particles, N_cells))
    println("Done!")
end

main()