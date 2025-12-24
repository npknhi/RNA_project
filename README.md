# RNA Structure Scoring Pipeline

This repository provides a complete pipeline for training, visualizing, and applying a knowledge-based statistical potential to evaluate RNA 3D structures.  
The method learns interatomic distance distributions from experimentally solved RNA structures and uses them to score predicted conformations.

---

## Build & Development Status

- **Docker support**: *In progress*  
  A Dockerfile is being prepared to allow fully containerized execution of the entire pipeline.
  
---

## Repository Structure

```
.
├── README.md
├── main.py
├── environment.yml
├── data
│   ├── structures/
│   │   ├── train/
│   │   └── test/
│   ├── plots/
│   ├── profiles/
│   └── scores/
├── reference/
├── src/
│   ├── training.py
│   ├── plotting.py
│   └── scoring.py
└── utils/
    ├── model.py
    ├── pair.py
    ├── rna_extractor.py
    └── interpolation.py
```

---

## Installation

Create and activate the Conda environment:

```bash
conda env create -f environment.yml
conda activate rna-score
```

---

# Running the Pipeline

The full pipeline (training → plotting → scoring) is executed with:

```bash
python main.py
```

All steps run **by default**.

---

## Command-Line Options

### Disable individual steps

```bash
python main.py --no-train     # Skip training
python main.py --no-plot      # Skip plotting
python main.py --no-score     # Skip scoring
```

---

### Use custom input/output directories

```bash
python main.py     --trainset data/structures/train     --profiles data/profiles     --testset data/structures/test     --scores data/scores
```

---

### Override model parameters

```bash
python main.py     --max-distance 25     --position-skip 3     --maximum-score 12     --maximum-score 2
```

---

# Outputs

### **1. Learned Interaction Profiles**

Saved to:

```
data/profiles/*.txt
```

Each file corresponds to a nucleotide pair (e.g., `AU.txt`, `CG.txt`).

---

### **2. Plot Figures**

One plot per base-pair profile:

```
data/plots/
├── AA.png
├── AU.png
├── CG.png
└── ...
```

---

### **3. Scoring Results**

Saved as:

```
data/scores/<timestamp>_scores.csv
```

CSV format:

```
struct_file, score
```

---

# Module Overview

### **`src/training.py`**
Learns distance-based scoring profiles from training structures.

### **`src/plotting.py`**
Generates a separate PNG figure for each base-pair interaction profile.

### **`src/scoring.py`**
Computes an estimated Gibbs free energy–based score for each RNA structure.

### **`main.py`**
Controls the entire pipeline with user-friendly command-line arguments.

### **`utils/`**
Contains core utilities: PDB/CIF parsing, distance computation, interpolation, and model configuration.

# KDE
Kernel Density Estimator (KDE) is a rational and visually pleasant representation of the data distribution. Especially useful for data distributions which are too irregular.

Source: [towardsdatascience.com](https://towardsdatascience.com/kernel-density-estimation-explained-step-by-step-7cc5b5bc4517/)