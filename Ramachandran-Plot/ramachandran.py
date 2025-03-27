import tkinter as tk
from tkinter import ttk
from Bio.PDB import PDBParser, PPBuilder
import numpy as np
import matplotlib.pyplot as plt
import requests
from io import StringIO

# Amino acids that will be filtered out before viewing the plot
AMINO_ACIDS = [
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS",
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"
]

# Initiaing the GUI
root = tk.Tk()
root.title("Ramachandran Plot GUI")

# Text box for pdb id or path to file
ttk.Label(root, text="PDB ID or File Path:").grid(row=0, column=0, sticky="w", padx=5, pady=5)
pdb_entry = ttk.Entry(root, width=40)
pdb_entry.grid(row=0, column=1, padx=5, pady=5)

# Checkboxes to filter amino acids
ttk.Label(root, text="Select Residues:").grid(row=1, column=0, sticky="nw", padx=5, pady=5)
residue_frame = ttk.Frame(root)
residue_frame.grid(row=1, column=1, padx=5, pady=5, sticky="w")

check_vars = {}
for i, aa in enumerate(AMINO_ACIDS):
    var = tk.BooleanVar(value=True)
    chk = ttk.Checkbutton(residue_frame, text=aa, variable=var)
    chk.grid(row=i//5, column=i%5, sticky="w", padx=2, pady=2)
    check_vars[aa] = var

def analysis():

    pdb_input = pdb_entry.get()
    selected_residues = [aa for aa, var in check_vars.items() if var.get()]
    parser = PDBParser(QUIET=True)

    # Checking if the file is local
    if pdb_input.endswith(".pdb"):
        structure = parser.get_structure("protein", pdb_input)
    # Getting pdb from rcsb if not local
    else:
        url = f"https://files.rcsb.org/download/{pdb_input.upper()}.pdb"
        response = requests.get(url)
        structure = parser.get_structure("protein", StringIO(response.text))


    # finding phi/psi angles for each AA
    phi_psi_angles = []
    for model in structure:
        for chain in model:
            ppb = PPBuilder()
            for pp in ppb.build_peptides(chain):
                phi_psi = pp.get_phi_psi_list()
                for i, (phi, psi) in enumerate(phi_psi):
                    res = pp[i]
                    if phi and psi:
                        if res.get_resname() in selected_residues:
                            phi_psi_angles.append((np.degrees(phi), np.degrees(psi)))

    phi, psi = zip(*phi_psi_angles)
    plt.figure(figsize=(6, 6))
    plt.scatter(phi, psi, color='blue', s=5)
    plt.title('Ramachandran Plot')
    plt.xlim(-180, 180)
    plt.ylim(-180, 180)
    plt.xlabel('Phi (ϕ)')
    plt.ylabel('Psi (ψ)')
    plt.grid(True)
    plt.show()

# Run button in gui
ttk.Button(root, text="Run", command=analysis).grid(row=2, column=0, columnspan=2, pady=10)

root.mainloop()
