# Introduction

This code scores the quality of the structure of an RNA sequences.

The main repository contains 3 main script:

* `training.py`, trains the objective function, using interatomic distance distributions that are computed from a dataset of known 3D structures;  
* `plotting.py` plots the scoring profiles, i.e. the score (or estimated Gibbs free energy) as a function of the interatomic distance; 
* `scoring.py` is an objective function to evaluate predicted structures from the RNA-Puzzles dataset.
* ``main.ipynb` where you can run the program by step. Here, the 3 steps are detailed.

See bellow the content of our 
```
.
├── README.md
├── data
│   ├── cifs
│   │   ├── ...
│   │   └── RNA.cif
│   ├── pdbs
│   │   ├── decoys
│   │   │   ├── ...
│   │   │   └── PZ7_Chen_1.pdb
│   │   └── native
│   │       ├── ...
│   │       └── 7QR4.pdb
│   ├── plots
│   │   ├── ...
│   │   └── UU.png
│   ├── profiles
│   │   ├── ...
│   │   └── UU.txt
│   └── scores
│       └── scores.csv
├── main.ipynb
├── reference
│   ├── additiononal_info.md
│   └── reference.py
├── src
│   ├── __pycache__
│   │   ├── ...
│   │   └── training.cpython-310.pyc
│   ├── plotting.py
│   ├── scoring.py
│   └── training.py
└── utils
    ├── ...
    ├── __pycache__
    │   ├── ...
    │   └── pdb.cpython-310.pyc
    ├── interpolation.py
    ├── model.py
    ├── pair.py
    └── pdb.py
```

### Prerequisite

To limit the storage used, we recommend to use a virtual environment (i.e conda)

```python 
conda env create -f data/environment.yml
conda activate rna-score
```

### Run a relation
To run the code, activate the :
```bash
python main.py
```
