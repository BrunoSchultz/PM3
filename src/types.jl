using StaticArrays
using FFTW

struct Particle
    m::Float64
    x::Float64; y::Float64; z::Float64
    vx::Float64; vy::Float64; vz::Float64
    fx::Float64; fy::Float64; fz::Float64
end

function Particle(pos::AbstractVector, vel::AbstractVector, m::Float64)
    return Particle(
        m, 
        pos[1], pos[2], pos[3], 
        vel[1], vel[2], vel[3], 
        0.0, 0.0, 0.0 # Forces initialize to zero
    )
end

struct Simulation{T_real, T_complex, P_fwd, P_bwd, T_k1D}
    box_size_physical::Float64   # in Mpc or Gpc? Might drop later and handle this when creating particles
    
    N_cells::Int64               # In code units, L_box = N_cells and dx = 1.0

    density_grid::T_real
    density_k_grid::T_complex

    Kx_sq::T_k1D
    Ky_sq::T_k1D
    Kz_sq::T_k1D

    plan_forward::P_fwd
    plan_backward::P_bwd
end 

function Simulation(box_size_physical::Float64, N_cells::Int64)
    N_k = div(N_cells, 2) + 1

    density_grid =  zeros(Float64, N_cells, N_cells, N_cells)
    density_k_grid = zeros(ComplexF64, N_k, N_cells, N_cells)

    plan_forward = plan_rfft(density_grid)
    plan_backward = plan_irfft(density_k_grid, N_cells)

    kx_arr = rfftfreq(N_cells) .* (2 * pi)
    ky_arr = fftfreq(N_cells)  .* (2 * pi)
    kz_arr = fftfreq(N_cells)  .* (2 * pi)

    Kx_sq = @. sin(kx_arr / 2.0)^2
    Ky_sq = @. sin(ky_arr / 2.0)^2
    Kz_sq = @. sin(kz_arr / 2.0)^2
    return Simulation(
        box_size_physical, 
        N_cells, 
        density_grid, 
        density_k_grid,
        Kx_sq,
        Ky_sq,
        Kz_sq,
        plan_forward,
        plan_backward
    )
end




