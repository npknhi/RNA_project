from Bio.PDB import PDBParser

def extract_c3_atoms(pdb_path):
    """
    Parse a PDB file using Biopython's PDBParser and extract all C3' atoms
    belonging to standard RNA residues (A, U, G, C) from the first model.

    Returns:
        list of lists: [chain_id, residue_name, x, y, z]
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_path)

    # Use only the first structural model in the PDB file
    model = next(structure.get_models())

    atoms = []
    valid_res = {"A", "U", "G", "C"}

    for atom in model.get_atoms():
        # Select only atoms named C3'
        if atom.get_name() != "C3'":
            continue

        residue = atom.get_parent()
        resname = residue.get_resname()
        if resname not in valid_res:
            continue

        chain = residue.get_parent()
        x, y, z = atom.coord

        atoms.append([chain.id, resname, float(x), float(y), float(z)])

    return atoms
