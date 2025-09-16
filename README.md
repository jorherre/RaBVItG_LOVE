RaBVIt-G (adapted) — R code for the paper “Emotional Trade-offs in Successful Romantic Relationships: A Differential Game Approach”

This repository contains a simple adaptation of the RaBVIt-G algorithm tailored to the paper. It keeps the core method and adds a systematic sweep over (alpha, beta) to produce results and figures in a reproducible way.

Requirements

R (version 4.3 or newer).

The packages used by the project (if there is a renv.lock, use renv::restore()).

Quick setup

Open R and run: install.packages("renv"); then run renv::restore()
(If you do not use renv, install the required packages manually when the scripts ask for them.)

How to run

Default run:
Rscript experiments/run.R

With a config file:
Rscript experiments/run.R --config configs/params.yaml

What it does

Explores a grid of (alpha, beta) values.

Computes the quantities of interest and saves tables and figures.

Produces heatmaps with a simple convention: blue = lower values, red = higher values, white ≈ 0. Some cells include numeric labels as references.

Outputs

Numerical results in the “results/” folder.

Figures (heatmaps, etc.) in the “figures/” folder.


Citation

Please cite the main paper and the earlier work where RaBVIt-G was introduced. Add the final DOI/URL to the repository when available.
