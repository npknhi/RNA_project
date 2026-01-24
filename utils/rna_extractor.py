# utils/extractor.py

from Bio.PDB import PDBParser, MMCIFParser
import os
import numpy as np

RNA_85_ATOMS = {
    "AC1'", "AC2", "AC2'", "AC3'", "AC4", "AC4'", "AC5", "AC5'",
    "AC6", "AC8", "AN1", "AN3", "AN6", "AN7", "AN9",
    "AO2'", "AO3'", "AO4'", "AO5'", "AOP1", "AOP2", "AP",

    "CC1'", "CC2", "CC2'", "CC3'", "CC4", "CC4'", "CC5", "CC5'",
    "CC6", "CN1", "CN3", "CN4", "CO2", "CO2'", "CO3'", "CO4'",
    "CO5'", "COP1", "COP2", "CP",

    "GC1'", "GC2", "GC2'", "GC3'", "GC4", "GC4'", "GC5", "GC5'",
    "GC6", "GC8", "GN1", "GN2", "GN3", "GN7", "GN9",
    "GO2'", "GO3'", "GO4'", "GO5'", "GO6", "GOP1", "GOP2", "GP",

    "UC1'", "UC2", "UC2'", "UC3'", "UC4", "UC4'", "UC5", "UC5'",
    "UC6", "UN1", "UN3", "UO2", "UO2'", "UO3'", "UO4", "UO4'",
    "UO5'", "UOP1", "UOP2", "UP",
}

def extract_c3_atoms(struct_path):
    """
    Parse a PDB/CIF file using Biopython's PDBParser/MMCIFParser and extract all C3' atoms
    belonging to standard RNA residues (A, U, G, C) from the first model.

    Arguments:

    Returns:
        list of lists: [chain_id, residue_name, x, y, z]
    """
    
    # Automatically choose parser based on file extension 
    ext = os.path.splitext(struct_path)[1].lower()

    if ext in [".cif", ".mmcif"]:
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)

    structure = parser.get_structure("structure", struct_path)

    # Use only the first model
    model = next(structure.get_models())

    atoms = []
    valid_res = {"A", "U", "G", "C"}

    for atom in model.get_atoms():
        if atom.get_name() != "C3'":
            continue

        residue = atom.get_parent()
        resname = residue.get_resname().strip()

        if resname not in valid_res:
            continue

        chain = residue.get_parent()
        x, y, z = atom.coord
        atoms.append([chain.id, resname, float(x), float(y), float(z)])

    return atoms

def extract_all_atom_residues(struct_path):
    """
    Extract all-atom coordinates per residue using the 85 RNA atom types list.
    Returns:
        [(chain_id, resseq, resname, coords_np_array(N,3)), ...]
    """
    ext = os.path.splitext(struct_path)[1].lower()

    if ext in [".cif", ".mmcif"]:
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)

    structure = parser.get_structure("structure", struct_path)
    model = next(structure.get_models())

    valid_res = {"A", "U", "G", "C"}
    residues_out = []

    for chain in model:
        for residue in chain:
            # ignore hetero/water
            if residue.id[0] != " ":
                continue

            resname = residue.get_resname().strip()
            if resname not in valid_res:
                continue

            coords = []
            for atom in residue:
                atom_name = atom.get_name().strip()  # e.g. "C1'", "OP1"
                if f"{resname}{atom_name}" in RNA_85_ATOMS:
                    coords.append(atom.coord.astype(float))

            if coords:
                residues_out.append((chain.id, residue.id[1], resname, np.vstack(coords)))

    return residues_out
