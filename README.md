# Code for "A nonparametric framework for treatment effect modifier discovery in high dimensions"

This repository is organized as follows:

 - `simulations/`: materials for reproducing the simulation study.
 - `application/`: materials for reproducing the FinHER clinical trial analysis.
 
The simulation results can be reproduced by running
`simulations/meals_obs-continuous.R`, `simulations/R/meals_obs-binary.R` and
`simulations/R/meals_rct-survival.R`. We recommend using high performance
computing (HPC) environment. Mock bash scripts are made available in
`simulations/slurm/` for use with an HPC environment using SLURM workload
manager.

The FinHER clinical trial data analysis can be performed by running
`application/R/analysis.Rmd`.
