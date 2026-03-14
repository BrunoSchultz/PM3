# Helper functions to calculate the force
@inline function force_x(phi, i, j, k, N)
    return (phi[mod1(i-1, N), j, k] - phi[mod1(i+1, N), j, k]) * 0.5
end

@inline function force_y(phi, i, j, k, N)
    return (phi[i, mod1(j-1, N), k] - phi[i, mod1(j+1, N), k]) * 0.5
end

@inline function force_z(phi, i, j, k, N)
    return (phi[i, j, mod1(k-1, N)] - phi[i, j, mod1(k+1, N)]) * 0.5
end

function calculate_forces!(particles, sim::Simulation)
    N = sim.N_cells 
    phi = sim.density_grid

    fill!(particles.fx, 0.)
    fill!(particles.fy, 0.)
    fill!(particles.fz, 0.)

    @inbounds Threads.@threads for i in eachindex(particles)
        px = particles.x[i]
        py = particles.y[i]
        pz = particles.z[i]

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
        w110 = dr_x * dr_y * dl_z
        w001 = dl_x * dl_y * dr_z
        w101 = dr_x * dl_y * dr_z
        w011 = dl_x * dr_y * dr_z
        w111 = dr_x * dr_y * dr_z

        # X - Forces
        particles.fx[i] += force_x(phi, lx_mod, ly_mod, lz_mod, N) * w000
        particles.fx[i] += force_x(phi, rx_mod, ly_mod, lz_mod, N) * w100
        particles.fx[i] += force_x(phi, lx_mod, ry_mod, lz_mod, N) * w010
        particles.fx[i] += force_x(phi, rx_mod, ry_mod, lz_mod, N) * w110
        particles.fx[i] += force_x(phi, lx_mod, ly_mod, rz_mod, N) * w001
        particles.fx[i] += force_x(phi, rx_mod, ly_mod, rz_mod, N) * w101
        particles.fx[i] += force_x(phi, lx_mod, ry_mod, rz_mod, N) * w011
        particles.fx[i] += force_x(phi, rx_mod, ry_mod, rz_mod, N) * w111

        # Y - Forces
        particles.fy[i] += force_y(phi, lx_mod, ly_mod, lz_mod, N) * w000
        particles.fy[i] += force_y(phi, rx_mod, ly_mod, lz_mod, N) * w100
        particles.fy[i] += force_y(phi, lx_mod, ry_mod, lz_mod, N) * w010
        particles.fy[i] += force_y(phi, rx_mod, ry_mod, lz_mod, N) * w110
        particles.fy[i] += force_y(phi, lx_mod, ly_mod, rz_mod, N) * w001
        particles.fy[i] += force_y(phi, rx_mod, ly_mod, rz_mod, N) * w101
        particles.fy[i] += force_y(phi, lx_mod, ry_mod, rz_mod, N) * w011
        particles.fy[i] += force_y(phi, rx_mod, ry_mod, rz_mod, N) * w111

        # Z - Forces
        particles.fz[i] += force_z(phi, lx_mod, ly_mod, lz_mod, N) * w000
        particles.fz[i] += force_z(phi, rx_mod, ly_mod, lz_mod, N) * w100
        particles.fz[i] += force_z(phi, lx_mod, ry_mod, lz_mod, N) * w010
        particles.fz[i] += force_z(phi, rx_mod, ry_mod, lz_mod, N) * w110
        particles.fz[i] += force_z(phi, lx_mod, ly_mod, rz_mod, N) * w001
        particles.fz[i] += force_z(phi, rx_mod, ly_mod, rz_mod, N) * w101
        particles.fz[i] += force_z(phi, lx_mod, ry_mod, rz_mod, N) * w011
        particles.fz[i] += force_z(phi, rx_mod, ry_mod, rz_mod, N) * w111
    end 
end