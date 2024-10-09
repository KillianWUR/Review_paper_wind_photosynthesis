# Leaf photosynthesis modelling with coupled photosynthesis, stomatal conductance and energy balance model, incl. all conductances according to Fick's law.

## Overview

This repository contains scripts which were used as part of the review paper "**[Your Paper Title]**".

- The R script provides an estimation of limitation within a crop canopy using a steady-state model. The vertical profile data used i nthe model are taken from https://doi.org/10.5194/bg-17-4375-2020. Typical parameters for wheat were used (cited in the paper). The script uses a nested iterative procedure to find the steady state photosynthesis rate.
- The Julia script uses an ODE solver to calculate photosynthesis over a time period, including a simple induction and relaxation model. It uses typical values for a tomato crop to show limitations by the boundary layer.


## Repository Contents

- `profile_limitation.R`: The R script that reads data, performs data cleaning, and builds steady-state models to analyze the limitation of gb on leaf photosynthesis within a canopy.
- `daily_limitation.jl`:

## Data Source

The data used in this script is sourced from:

Vilà-Guerau de Arellano, J., Ney, P., Hartogensis, O., De Boer, H., Van Diepen, K., Emin, D., ... & Graf, A. (2020). CloudRoots: integration of advanced instrumental techniques and process modelling of sub-hourly and sub-kilometre land–atmosphere interactions. Biogeosciences, 17(17), 4375-4404.

If you use this script or data, please cite both the original paper and this repository.
