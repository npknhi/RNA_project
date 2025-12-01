from Bio.PDB import PDBParser

def extract_c3_atoms(pdb_path):
    """
    PDB parser
    
    Extracts each chain of the first model for atom type "C3'" and its x,y,z coordinates in angstrom.
    
    Returns:
        list of tuples: (chain_identifier, residue_name, x, y, z)
    """
    atoms = []
    model = 0

    with open(pdb_path, "r") as f:
        for row in f:
            if len(row) < 7:
                continue

            atom = row[0:6].strip()
            if atom == "MODEL":
                model += 1
                if model > 1:
                    return atoms

            if atom == "ATOM":
                atom_name = row[12:16].strip()
                if atom_name == "C3'":
                    chain_identifier = row[21:22].strip()
                    residue_name = row[17:20].strip()
                    if residue_name in ["A", "U", "G", "C"]:
                        x = float(row[30:38].strip())
                        y = float(row[38:46].strip())
                        z = float(row[46:54].strip())
                        atoms.append((chain_identifier, residue_name, x, y, z))

    return atoms



# def extract_c3_atoms(pdb_path):
#     """
#     Parse PDB file using Biopython's PDBParser.
#     Extracts chain identifier, residue name, and coordinates of atom "C3'".
    
#     Returns:
#         list of tuples: (chain_identifier, residue_name, x, y, z)
#     """
#     parser = PDBParser(QUIET=True)
#     structure = parser.get_structure("structure", pdb_path)

#     atoms = []
#     # Only take the first model
#     model = next(structure.get_models())
#     for chain in model:
#         for residue in chain:
#             resname = residue.get_resname().strip()
#             if resname in ["A", "U", "G", "C"]:
#                 for atom in residue:
#                     if atom.get_name() == "C3'":
#                         x, y, z = atom.coord
#                         atoms.append((chain.id, resname, x, y, z))
#     return atoms
