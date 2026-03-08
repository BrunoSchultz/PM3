using StructArrays

function assign_density!(particles::StructVector{Particle}, sim::Simulation)
    N = sim.N_cells
    grid = sim.density_grid

    fill!(grid, 0.0)
    
    @inbounds for i in eachindex(particles)
        px = particles.x[i]
        py = particles.y[i]
        pz = particles.z[i]
        mass = particles.m[i]

        # Get left indices 
        lx = floor(Int, px)
        ly = floor(Int, py)
        lz = floor(Int, pz)

        # Calculate fractions
        dr_x = px - lx;  dl_x = 1.0 - dr_x
        dr_y = py - ly;  dl_y = 1.0 - dr_y
        dr_z = pz - lz;  dl_z = 1.0 - dr_z

        # Enforce periodic boundaries 
        lx_mod = mod1(lx, N); rx_mod = mod1(lx + 1, N)
        ly_mod = mod1(ly, N); ry_mod = mod1(ly + 1, N)
        lz_mod = mod1(lz, N); rz_mod = mod1(lz + 1, N)
        
        w000 = dl_x * dl_y * dl_z
        w100 = dr_x * dl_y * dl_z
        w010 = dl_x * dr_y * dl_z
        w001 = dl_x * dl_y * dr_z
        w110 = dr_x * dr_y * dl_z
        w101 = dr_x * dl_y * dr_z
        w011 = dl_x * dr_y * dr_z
        w111 = dr_x * dr_y * dr_z

        # Assign mass to the grid
        grid[lx_mod, ly_mod, lz_mod] += mass * w000
        grid[rx_mod, ly_mod, lz_mod] += mass * w100
        grid[lx_mod, ry_mod, lz_mod] += mass * w010
        grid[lx_mod, ly_mod, rz_mod] += mass * w001
        grid[rx_mod, ry_mod, lz_mod] += mass * w110
        grid[rx_mod, ly_mod, rz_mod] += mass * w101
        grid[lx_mod, ry_mod, rz_mod] += mass * w011
        grid[rx_mod, ry_mod, rz_mod] += mass * w111
    end 
end
