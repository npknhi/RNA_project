from Bio.PDB import PDBParser, MMCIFParser
import os

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
