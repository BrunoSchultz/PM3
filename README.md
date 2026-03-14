# PM3: A Toy Particle-Mesh Cosmology Simulation

**Note:** This is a hobby side-project built purely for fun and learning! It is *not* heavily optimized for efficiency and lacks several key features found in professional cosmology codes:

* **Proper Initial Conditions (Zel'dovich Approximation):** Currently, `PM3` initializes particles completely at random with zero initial velocity. A real simulation requires generating a primordial density field from a specific power spectrum (e.g., using the Zel'dovich approximation or 2nd-order Lagrangian Perturbation Theory)
* **Adaptive Time-Stepping:** `PM3` uses a fixed step size for the scale factor `da`.
* **Distributed Computing (MPI):** The code is restricted to a single machine's shared memory (using Julia's `Threads.@threads`).

Currently, PM3 operates entirely in code units rather than physical units.

---

## The Particle-Mesh Method

In cosmology, N-body simulations track the gravitational evolution of dark matter over cosmic time. Calculating gravity directly between every pair of particles scales as $O(N^2)$ and becomes computationally unfeasible for millions of particles. 

The **Particle-Mesh (PM)** method solves this by approximating the gravitational field on a discrete 3D grid:
1.  **Mass Assignment:** Particles "smear" their mass onto nearby grid cells.
2.  **Poisson Solver:** Gravity is computed on this grid. By moving to Fourier space (using Fast Fourier Transforms), the differential equations of gravity become simple algebraic multiplications, reducing the computational cost.
3.  **Force Interpolation:** The gravitational forces are calculated on the grid and mapped back to the individual particles.
4.  **Integration:** Particles are nudged forward in time based on these forces.

---

## Code Architecture & Execution Order

My codebase revolves around a few core structures and a straightforward execution loop. 

### Structs
* **`Cosmology`**: Stores cosmological parameters like $\Omega_{m0}$, $\Omega_{\Lambda 0}$, and $H_0$. The default helper function `ΛCDM()` sets up a standard cosmology.
* **`Simulation`**: Holds the pre-allocated density grids, complex k-space grids, squared wave-number arrays (`Kx_sq`, etc.), and the FFT plans (`plan_forward` and `plan_backward`) required for the Poisson solver.
* **`Particle`**: A struct holding mass, 3D positions, 3D velocities, and 3D forces.

### The Simulation Loop
During a simulation run, the engine loop repeatedly calls four primary functions:

1.  **`assign_density!(particles, sim)`**
    * Clears the `density_grid` and maps every particle's mass onto it.
    * Uses Cloud-in-Cell (CIC) interpolation to distribute a particle's mass across the 8 nearest grid cells based on fractional distances (`dr_x`, `dl_x`, etc.), enforcing periodic boundaries.
2.  **`solve_poisson!(sim, C, a)`**
    * Converts the real-space density grid into a gravitational potential grid.
    * It executes a forward FFT (`plan_forward`) to move the density into k-space. Then, `apply_greens_function!` multiplies the k-grid the Green's function. Finally, an inverse FFT (`plan_backward`) transforms the grid back to real space.
3.  **`calculate_forces!(particles, sim)`**
    * Computes the $X$, $Y$, and $Z$ forces acting on each particle.
    * It calculates the gradient of the potential grid using finite differences (`force_x`, `force_y`, `force_z`) and uses the exact same CIC fractional weights (`w000`, `w100`, etc.) to interpolate those grid forces back onto the individual particles.
4.  **`step_leapfrog!(particles, sim, C, a, da)`**
    * Advances the particles forward in time using a symplectic Kick-Drift integrator.
    * It calculates cosmological `kick_factor` and `drift_factor` variables (dependent on the scale factor $a$ and step size `da`). It then updates velocities with the forces (Kick) and positions with the new velocities (Drift), keeping particles within the simulation box boundaries.

---

## Installation Guide

If you've never used Julia before, follow these steps to get the environment set up.

**Step 1: Install Julia**
* Download and install the latest stable release of Julia from the official website: [https://julialang.org/downloads/](https://julialang.org/downloads/).

**Step 2: Clone the Repository**
* Open your terminal and download the codebase:
    ```bash
    git clone [https://github.com/yourusername/PM3.git](https://github.com/yourusername/PM3.git)
    cd PM3
    ```

**Step 3: Setup the Julia Environment**
* Start the Julia REPL by typing `julia` in your terminal from the `PM3` folder.
* Press the `]` key on your keyboard. Your prompt will change to `pkg>`, which means you are in the package manager.
* Run the following commands to activate the project and install all required dependencies:
    ```julia
    pkg> activate .
    pkg> instantiate
    ```
* Once it finishes downloading, press the `Backspace` key to return to the standard `julia>` prompt, then type `exit()` to close it.

---

## Running the Simulation

The workflow is split into two parts: running the computation (`run.jl`) and generating visual plots (`plot.jl`). 

**Important:** You must use the `--project=.` flag when running these scripts from the terminal so that Julia knows to use the local environment you just installed!

### 1. The Compute Script (`run.jl`)
This script initializes the universe, runs the time integration, and saves the final raw data to a `.jls` binary file. The default simulation runs with 256 cells and 10 million particles, which consumes roughly 16 GB of RAM.

**Command Line Options:**
* `--N_cells`: Number of grid cells along one edge of the 3D box (Default: 256).
* `--N_parts`: Total number of particles to simulate (Default: 10000000).
* `--a_start`: The starting cosmic scale factor (Default: 0.1).
* `--a_end`: The ending cosmic scale factor (Default: 1.0).
* `--da`: The integration time-step size (Default: 0.01).
* `--out`: Output file name for the final particle state (Default: `final_state.jls`).

**Example:**
```bash
julia --project=. run.jl --N_cells 128 --N_parts 1000000 --out test_run.jls
```

Once run.jl finishes and saves your .jls file, use plot.jl to visualize the results. This script extracts a 2D slice of the 3D simulation box and renders a scatter plot. Separating this from the main engine allows you to tweak visualization settings without re-running the computation.

Command Line Options:

* `--input`: The .jls data file generated by run.jl (Default: final_state.jls).

* '`--slice_center`': The Z-coordinate center of the visual slice in code units (Default: 32.0).

* `--thickness`: The depth thickness of the slice (Default: 6.0).

* `--out`: Output image filename (Default: cosmic_web_slice.png).

Plot the test run we just generated:

```bash
julia --project=. plot.jl --input test_run.jls --slice_center 64.0 --thickness 10.0 --out my_test_plot.png
```

🤖 AI Usage Disclaimer

The core physics engine and math were written by hand for the sake of learning, the test files (and some of the documentation formatting) were heavily generated using Gemini 3.1 Pro. 