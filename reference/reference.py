#!/usr/bin/env python3
"""
RNA Structure Scoring Example using rna-tools
Requirements:
    pip install rna-tools biopython
"""

import sys
import os
from rna_tools.tools.rna_calc_rmsd import main as calc_rmsd
from Bio.PDB import PDBParser

def validate_pdb_file(pdb_path):
    """Check if the PDB file exists and is readable."""
    if not os.path.isfile(pdb_path):
        raise FileNotFoundError(f"PDB file not found: {pdb_path}")
    if not pdb_path.lower().endswith(".pdb"):
        raise ValueError("Input file must have a .pdb extension")
    return True

def compute_rmsd(pdb_ref, pdb_model):
    """
    Compute RMSD between reference and model RNA structures.
    Returns RMSD value in Ångströms.
    """
    try:
        # rna-tools RMSD calculation
        sys.argv = ["rna_calc_rmsd", pdb_ref, pdb_model]
        calc_rmsd()
    except Exception as e:
        print(f"Error computing RMSD: {e}")

def count_atoms(pdb_path):
    """Count atoms in the PDB file using Biopython."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("RNA", pdb_path)
    atom_count = sum(1 for _ in structure.get_atoms())
    return atom_count

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python score_rna.py reference.pdb model.pdb")
        sys.exit(1)

    ref_pdb = sys.argv[1]
    model_pdb = sys.argv[2]

    try:
        validate_pdb_file(ref_pdb)
        validate_pdb_file(model_pdb)

        print(f"Reference atoms: {count_atoms(ref_pdb)}")
        print(f"Model atoms: {count_atoms(model_pdb)}")

        print("\n--- RMSD Calculation ---")
        compute_rmsd(ref_pdb, model_pdb)

        print("\nScoring complete.")
    except Exception as err:
        print(f"Error: {err}")
        sys.exit(1)