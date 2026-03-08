using FFTW
using LinearAlgebra

function solve_poisson!(sim::Simulation, C::Cosmology, a::AbstractFloat)
    density_grid = sim.density_grid
    density_grid_k = sim.density_k_grid

    plan_forward = sim.plan_forward
    plan_backward = sim.plan_backward

    mul!(density_grid_k, plan_forward, density_grid)
    apply_greens_function!(sim, C, a) # density_k is now the (fourier transform) of the potential
    mul!(density_grid, plan_backward, density_grid_k)

    return nothing
end

function apply_greens_function!(sim::Simulation, C::Cosmology, a::AbstractFloat)
    N = sim.N_cells
    Kx_sq = sim.Kx_sq
    Ky_sq = sim.Ky_sq
    Kz_sq = sim.Kz_sq
    k_grid = sim.density_k_grid

    prefactor = -(3.0 * C.Ω_m0) / (8.0 * a)

    N_k = div(N, 2) + 1
    @inbounds for k in 1:N
        for j in 1:N
            for i in 1:N_k
                K = Kx_sq[i] + Ky_sq[j] + Kz_sq[k]

                # The singularity at l=m=n=0 handled by setting potential to zero
                if K == 0.0
                    k_grid[i, j, k] = 0.0 + 0.0im
                else
                    # G(k) * delta(k)
                    k_grid[i, j, k] *= prefactor / K
                end
            end
        end
    end
    return nothing
end