# Unravelling the instability of Mutational Signatures extraction via Archetipal Analysis

## Requisites

- python 3.8.3
- sklearn 0.24.2
- archetypes.py : open source software for archetypal analysis downloadable at https://doi.org/10.25919/5d3958889f7ff
- SigProfilerExtractor : open source software for mutational signature extraction  downloadable at https://github.com/AlexandrovLab/SigProfilerExtractor


## Repository Content
This repository contains the following folders:
- utils: COSMIC signatures profiles and usefull python functions
- notebook: 2 Jupyter notebook to reproduce paper analysis
- synthetic_catalogues: SigsPack synthetic catalogues for each scenario

## Usage
- Archetypal_Analysis.ipynb can be executed indipendently

- Signature_Extraction.ipynb calculate metrics for each scenario after de novo extraction with SigProfilerExtractor. 
  Given the number of simulation and the size of the results to run the notebook is necessary to first perform de novo extraction for each scenario with
  SigProfilerExtractor executing run_SigProfilerExtractor.py. The computation is really expensive and it may takes some times.
 
