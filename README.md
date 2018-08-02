# aswanalysis

Analyze Adaptive Sampling Workflow Results


This package contains routines to co-process different simulation workflows. Say you have a long, single
trajectory from a biomolecular MD simulation and a (separate) set of trajectories from an adaptive
sampling workflow run on the same MD system, then you can use this package to configure and run analyses
that examine and compare all the data. 

### Installation

1. Clone this repo https://github.com/jrossyra/aswanalysis
2. Load conda environment or virtualenvironment container
2. Run `setup.py`

### Analysis Workflow

0. Create Data (outside this package)
1. Dataset Configuration
2. [Create modified version of `analyze.py`]
3. Run `analyze.py` for dataset
