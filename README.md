# RNA Structure Scoring Pipeline

**Authors**  
- **Phuc Khanh Nhi NGUYEN**  
- **Emilie Jeanne MAYI**  

**Program**  
*Master 2 GENIOMHE-AI*  
*Université d’Évry Paris-Saclay*

---

This repository provides a complete pipeline for training, visualizing, and applying a knowledge-based statistical potential to evaluate RNA 3D structures.  
The method learns interatomic distance distributions from experimentally solved RNA structures and uses them to score predicted conformations.

---

## Build & Development Status

**Docker support**: *In progress*  
  A Dockerfile is being prepared to allow fully containerized execution of the entire pipeline.

---

## Repository Structure

```
.
├── README.md
├── main.py
├── environment.yml
├── data
│ ├── structures/
│ │ ├── train/
│ │ └── test/
│ ├── plots/
│ ├── profiles/
│ ├── scores/
│ └── distances.csv
├── src/
│ ├── training.py
│ ├── plotting.py
│ ├── scoring.py
│ └── KDE.r
└── utils/
│ ├── model.py
│ ├── pair.py
│ ├── rna_extractor.py
│ └── interpolation.py
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

### Atom representation

```bash
python main.py --atom-mode c3prime
python main.py --atom-mode all_atom
```

---

### Use custom input/output directories

```bash
python main.py     --trainset data/structures/train     --profiles data/profiles     --testset data/structures/test     --scores data/scores
```

---

### Override model parameters

```bash
python main.py     --max-distance 25     --position-skip 3     --bin-width 1.0     --maximum-score 2
```

---

# Outputs

### **1. Learned Interaction Profiles**

Saved to:

```
data/profiles/*.txt
data/profiles_all_atom/*.txt
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

## Run KDE (python file)
```bash
python main.py --no-train --no-plot --no-score --make-kde
```

## Run KDE (Rscript)
As of now:
* you implicitly create a KDE plot with python file when creating `distances.csv`.
* **you get this message twice**: Rscript does not open a graphics window like RStudio. R creates a null graphics device if no interactive device is available
```
null device 
          1
```

```bash
# Create the distances file
python main.py --no-train --no-plot --no-score --make-kde
Rscript src/KDE.r data/distances.csv TRUE
```

## Non-log Scoring Formula

In addition to the classical logarithmic potential of mean force, the pipeline implements an alternative **non-logarithmic scoring formulation**, inspired by the total information gain approach.

### Logarithmic formulation (default)

The standard potential of mean force is defined as:

$$
U_{ij}(r) = -\log \left( \frac{f_{ij}(r)}{f_{\text{ref}}(r)} \right)
$$

where:

- $f_{ij}(r)$ denotes the observed distance frequency for nucleotide pair $(i,j)$,
- $f_{\text{ref}}(r)$ denotes the reference distance frequency.

This logarithmic formulation strongly penalizes rare or unfavorable interactions, leading to sharp energy variations when observed frequencies deviate from the reference state.

---
### Non-logarithmic formulatio

As an alternative, a non-logarithmic scoring function is implemented to provide a smoother penalty profile:

$$
U_{ij}(r) =-
\frac{f_{ij}(r) - f_{\text{ref}}(r)}{f_{\text{ref}}(r)}
$$

This expression represents the pairwise contribution of nucleotide pair $(i,j)$ at distance $r$ to the overall score.

The total non-log score of a structure is then obtained by summing over all interacting pairs:

$$
\text{score} = \sum_{i,j} U_{ij}(r)
$$

By avoiding the logarithmic transformation, this formulation reduces the impact of low-frequency events while preserving the same reference state. As a result, it yields a smoother and more robust scoring function for model quality assessment.


---

### Usage

The non-log scoring formulation can be selected at training time using:

```bash
python main.py --score-formula linear

