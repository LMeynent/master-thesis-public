# Human behaviour and infectious diseases
## Using demographically representative survey data as a predictor for the progression of infectious disease outbreaks

This repository contains the code related to my MSc thesis. While the figures between the notebooks in this repository and the thesis may not exactly match, the code used to generate them remains the same.

## Organisation

This repository is organised with the following folders:
 * `figures` is en empty folder in which generated figures are stored
 * `hope_data` contains the pre-processed, publicly available data sources as well as a list of HOPE survey questions
 * `src` contains some useful Python functions shared by multiple analyis notebooks
 * `stan` contains the code necessary to run the hierarchical bayesian model used to model the effective reproduction number

The notebooks at the root of the repository contain most of the analysis presented in the thesis report. We also provide a PDF version of the thesis report. The notebook `Extended HOPE data exploration` contains some additional exploratory analysis of the survey, in a more draft-y format.

## Environment

Conda environment is provided. To set it up, please use:

```Shell
conda env install -f environment.yml
conda activate master-thesis
```

## Data availability

For data privacy reasons, the survey data is not published with the repository. Some external, public data sources are provided in a pre-processed format. 

In addition, the result file of the Makov-Chain Monte-Carlo sampling of our hierarchical Bayesian model is too large to be included in the repository.
