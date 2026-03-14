using PM3 
using StructArrays
using Plots

using Plots

function plot_slice(particles, N_cells; slice_center = N_cells/2, thickness = 5.0)
    
    println("Extracting 2D slice at Z = $slice_center ± $(thickness/2)...")
    
    in_slice = @. abs(particles.z - slice_center) < (thickness / 2.0)
    
    # Extract the X and Y coordinates
    x_slice = particles.x[in_slice]
    y_slice = particles.y[in_slice]
    
    println("Found $(length(x_slice)) particles in the slice. Plotting...")

    # Render the Plot
    p = scatter(
        x_slice, y_slice,
        markersize = 1.5,          
        markerstrokewidth = 0,     
        color = :black, 
        alpha = 0.1,               
        legend = false,            
        aspect_ratio = 1.0,        
        xlims = (0, N_cells),      
        ylims = (0, N_cells),
        framestyle = :box,
        grid = false,              
        title = "Particle Slice (Z ≈ $slice_center)",
        xlabel = "X (Code Units)",
        ylabel = "Y (Code Units)",
        dpi = 300                  
    )
    
    return p
end


N_cells = 256
N_parts = 100_000_00

my_box = Simulation(10000.,N_cells)
my_cosmology = ΛCDM()

println("Generating $N_parts random particles...")
particles = StructArray{Particle}(undef, N_parts)
for i in 1:N_parts
    particles.m[i] = 1.0
    
    # Generate random scalars directly
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

a_start = 0.002
a_end = 1.0
steps = 500
da = (a_end - a_start) / steps

# The Engine (Silenced for benchmarking)
function run_simulation!(particles, sim, C, a_start, a_end, da)
    a = a_start
    while a < a_end
        assign_density!(particles, sim)
        solve_poisson!(sim, C, a) 
        calculate_forces!(particles, sim)
        step_leapfrog!(particles, sim, C, a, da)
        a += da
    end 
end

println("Running simulation")
@time run_simulation!(particles, my_box, my_cosmology, a_start, a_end, da)

# Generate the plot
my_plot = plot_slice(particles, N_cells, slice_center = 32.0, thickness = 6.0)

savefig(my_plot, "cosmic_web_slice.png")
println("Saved plot to cosmic_web_slice.png!")