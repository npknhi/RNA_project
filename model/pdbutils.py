from Bio.PDB import PDBParser

def extract_c3_atoms(pdb_path):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("", pdb_path)

    atoms = []
    for chain in structure.get_chains():
        for residue in chain:
            if "C3'" in residue:
                atoms.append((residue.resname[0], residue["C3'"].coord))

    return atoms
