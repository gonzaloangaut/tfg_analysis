import numpy as np
import pandas as pd
import os

def process_and_save_msd(num_cells, max_step, dens, step, rng_seed, max_aspect_ratio):
    num_cells_folder = f"N={num_cells:_}"
    densfolder = f"{dens:.3f}".replace(".", "_")
    out_dir = f"data/{num_cells_folder}/{densfolder}/dat_msd"
    os.makedirs(out_dir, exist_ok=True)
    
    # Calculate the size of the box of the simulation
    L = np.sqrt(num_cells / dens)
    half_L = 0.5 * L

    for seed in rng_seed:
        msd_time_series = []
        
        # Read the first step
        file_t0 = (
            f"data/{num_cells_folder}/{densfolder}/dat/culture_initial_number_of_cells={num_cells}_density={dens}_"
            f"force=Anisotropic_Grosmann_k=3.33_gamma=3_With_Noise_eta=0.033_With_Shrinking_"
            f"rng_seed={seed}_step=00000.dat"
        )
        
        if not os.path.exists(file_t0):
            continue 
            
        df_t0 = pd.read_csv(file_t0)
        
        # Save the position of the "previous" step
        r_prev_wrapped = df_t0[["position_x", "position_y"]].values 
        
        # We copy it
        r_unwrapped = np.copy(r_prev_wrapped)
        
        # Save the first value
        r0_unwrapped = np.copy(r_unwrapped)

        for tic in range(0, max_step + 1, step):
            dat_actual = (
                f"data/{num_cells_folder}/{densfolder}/dat/culture_initial_number_of_cells={num_cells}_density={dens}_"
                f"force=Anisotropic_Grosmann_k=3.33_gamma=3_With_Noise_eta=0.033_With_Shrinking_"
                f"rng_seed={seed}_step={tic:05}.dat"
            )
            
            if os.path.exists(dat_actual):
                df_tic = pd.read_csv(dat_actual)
                r_curr_wrapped = df_tic[["position_x", "position_y"]].values
                
                # Calculate the displacement for the previous frame
                delta_r = r_curr_wrapped - r_prev_wrapped
                
                # Unwrap the periodic boundary conditions
                mask_x = np.abs(delta_r[:, 0]) > half_L
                delta_r[mask_x, 0] -= np.sign(delta_r[mask_x, 0]) * L
                
                mask_y = np.abs(delta_r[:, 1]) > half_L
                delta_r[mask_y, 1] -= np.sign(delta_r[mask_y, 1]) * L
                
                # Sum the real displacement
                if tic > 0:
                    r_unwrapped += delta_r
                
                # Calculate the msd with respect to the origin (GLOBAL)
                sq_disp = np.sum((r_unwrapped - r0_unwrapped)**2, axis=1)
                msd_global = np.mean(sq_disp)
                
                # We still count phenotypes to know the state of the system
                is_round = np.isclose(df_tic["aspect_ratio"], 1.0)
                n_round = is_round.sum()
                n_elongated = num_cells - n_round
                
                msd_time_series.append({
                    "step": tic,
                    "msd_global": msd_global,
                    "n_round": n_round,
                    "n_elongated": n_elongated
                })
                
                # Update memory for next iteration
                r_prev_wrapped = np.copy(r_curr_wrapped)
        
        if msd_time_series:
            df_out = pd.DataFrame(msd_time_series)
            output_file = f"{out_dir}/msd_series_density={dens}_rng_seed={seed}.dat"
            df_out.to_csv(output_file, index=False)


# Main Script
density_list = [0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.71, 0.72, 0.73,
                0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.8, 0.81, 0.82, 0.83, 0.84, 
                0.85, 0.86, 0.87, 0.88, 0.89, 0.9]
nc = 500
max_step = 80_000
step = 100
number_of_realizations = 64
max_aspect_ratio = 5.0
seed_1 = 0x87351080E25CB0FAD77A44A3BE03B491

rng_1 = np.random.default_rng(seed_1)
rng_seed = rng_1.integers(low=2**20, high=2**50, size=number_of_realizations)

for dens in density_list:
    process_and_save_msd(nc, max_step, dens, step, rng_seed, max_aspect_ratio)
    print(f"Completado: densidad {dens}")