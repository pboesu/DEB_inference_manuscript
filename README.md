# Two case studies detailing Bayesian parameter inference for DEB models

This repository contains code that forms the supplementary materials for the manuscript 

Boersch-Supan PH & Johnson LR 2019 ["Two case studies detailing Bayesian parameter inference for dynamic energy budget models"](https://doi.org/10.1016/j.seares.2018.07.014) Journal of Sea Research **143**:57-63

An open access preprint of this manuscript is available on [bioRxiv](https://doi.org/10.1101/259705).

The folder `DEBKiss` contains the DEBKiss case study. The main R script is `DEBKiss_snail.R`. The remaining files contain function definitions and the underlying data.

The folder `StandardDEB` contains the Standard DEB model case study. The main R script is  `earthworm_DEB_inference_fixedEG.R` The subdirectory `data` contains the original earthworm data that the simulations are based on. The subdirectory `src` contains R and C implementations of the scaled standard DEB model.
