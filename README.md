# Introduction

This code scores the quality of the structure of an RNA sequences.

The main repository contains 3 main script:

* `training.py`, trains the objective function, using interatomic distance distributions that are computed from a dataset of known 3D structures;  
* `plotting.py` plots the scoring profiles, i.e. the score (or estimated Gibbs free energy) as a function of the interatomic distance; 
* `evaluation.py` is an objective function to evaluate predicted structures from the RNA-Puzzles dataset.

### Prerequisite

To limit the storage used, we recommend to use a virtual environment (i.e conda)

```python 

```
