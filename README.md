# tfg_analysis

Software for analyzing data from oncostream simulations using the code of the tfg. The repository contains raw simulation data, processing scripts, and visualization notebooks to study clustering, order parameters, and phenotypic transitions.

## Repository structure

- **`data/`**: Contains all the simulation data.  
  Inside, there are folders named after the number of cells. Each contains subdirectories for different densities, following this structure:  
  - **`dat/`**: `.dat` files with the position, orientation, and aspect ratio of each cell for every step and seed.  
  - **`dat_order_parameters/`**: Information about polar and nematic order parameters and the fraction of elongated cells for each step and seed.  
  - Example files are kept only for reference due to the large amount of data.  
  - **`steady_state_images/`**: A folder with images of the steady state for each density is also included.

- **`graphs/`**: Stores all the figures generated during the analysis. Contains figures grouped by the number of cells, as well as comparison plots across different system sizes..

- **`other_scripts/`**: Contains ideas and draft scripts for potential future use.

## Analysis scripts

There are two types of scripts used to analyze the data:

### Python scripts (`.py`)
Generate lighter, intermediate data to simplify the analysis:
- **`dat_msd.py`**:  
  Generates msd data. 
- **`dat_ph_aspect_ratio.py`**:  
  Calculates and saves the quantity of each phenotype for every step and seed.

### Jupyter notebooks (`.ipynb`)
Used to visualize and analyze the processed data:
- **`analysis_order_parameters.ipynb`**:  
  Studies polar and nematic order parameters together with phenotype metrics. Plots their evolution in time and the steady-state values versus density.
- **`analysis_msd.ipynb`**:  
  Studies the msd to identify the distinct macroscopic dynamic regimes.  
