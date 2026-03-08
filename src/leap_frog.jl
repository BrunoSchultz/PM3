function step_leapfrog!(particles, sim, C, a, da)
    N = sim.N_cells
    
    a_half = a + da / 2.0
    
    # Kick is evaluated at current integer step (a_n)
    kick_factor = da * f(a, C) 
    
    # Drift is evaluated at the half-step (a_n+1/2)
    drift_factor = da / (a_half^2) * f(a_half, C) 

    @inbounds for i in eachindex(particles)
        # KICK
        particles.vx[i] += particles.fx[i] * kick_factor
        particles.vy[i] += particles.fy[i] * kick_factor
        particles.vz[i] += particles.fz[i] * kick_factor

        # DRIFT
        particles.x[i] = mod1(particles.x[i] + particles.vx[i] * drift_factor, N) 
        particles.y[i] = mod1(particles.y[i] + particles.vy[i] * drift_factor, N) 
        particles.z[i] = mod1(particles.z[i] + particles.vz[i] * drift_factor, N)  
    end
end