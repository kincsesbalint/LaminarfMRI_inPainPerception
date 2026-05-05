<br />
<div align="center">
  <h3 align="center">Laminar fMRI reveals layer‑specific dissociation of sensory and cognitive influences on pain processing </h3>

  <img src="/logo/grouplvlmaineffects_distraction.png" alt="Group-level main effects" width="900">


</div>


## Project content 

This repository contains a full processing pipeline for the Layer-fMRI pain project, from raw image preparation to first-level modeling, layer-specific sampling, and group-level statistics. The raw data has been minimally preprocessed, and first level parameter estimation was performed in the native space.

### Recommended order to reproduce the analysis

1. `preprocessing/`: convert raw data, perform reconstruction/cleanup, and run core preprocessing steps.
2. `a_compcor/`: create CompCor masks/components and estimate noise regressors.
3. `1st_level_analysis/` + `layerfication/`: define ROIs, create layer-informed surfaces/samplings, and compute layer-wise signal estimates. Build and run first-level SPM models for each subject and extract parameter estimates at the layer level.
4. `stats/`: run statistical analyses on layer-wise outputs and generate inference-ready results.
5. `statistical_analysis/`: high-level statistical workflows and reporting scripts.
6. `data/`: input and derived data tables used by the codes found in `statistical_analysis` (behavior, layer estimates, tSNR summaries, mappings).

### Folder purpose (short description)

- `preprocessing/`: utilities and pipelines for data conversion, reconstruction, normalization handling, and pre-fMRI preparation.
- `a_compcor/`: CompCor implementation and helper subfunctions for nuisance regressor construction.
- `1st_level_analysis/`: Contains first level SPM batch to estimate SPM.mat files individually, which can be used later (e.g. the design matrix) to estimate activations. 
- `layerfication/`: Contain the main functions which sample and/or estimate layer level signal from the specified ROI for each individual (using hte preprocessed functional data). Take a look into the functions for more detail. The main function is called BK_layer_sampling_pain_study_pipeline.
- `stats/`: dedicated scripts for downstream layer/statistical testing from modeled outputs.
- `statistical_analysis/`: additional statistical notebooks/scripts for extended analyses.
- `data/`: structured outputs and metadata that support reproducibility and visualization.

## Disclaimer

The code in this repository is not yet organized in its final user-friendly form. However, the full pipeline has been tested for reproducibility multiple times and works reliably. A cleaner and more user-friendly project structure will follow in a future update.



The original repo can be found [here](https://github.com/viktor-pfaffenrot/laminar-fMRI-of-the-perception-of-pain). Many functions are added and the pipeline was fine tuned. The reason of this standalone repository to increase visibility of the project. 
