import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.distance import squareform
import os

# --- 1. Paths and Configuration ---
base_path = '1000_samples/MutationConformationsBioEmu'

# Native (Wild Type) paths
wt_xtc = f'{base_path}/out_native/samples.xtc'
wt_top = f'{base_path}/out_native/topology.pdb'

# Mutant (P36R) paths
mut_xtc = f'{base_path}/out_mutant_P36R/samples.xtc'
mut_top = f'{base_path}/out_mutant_P36R/topology.pdb'

# --- 2. Processing Function ---
def get_mean_contact_data(xtc, pdb, scheme='ca'):
    traj = md.load(xtc, top=pdb)
    print(f"Loaded trajectory: {traj}")

    # Compute distances for all residue pairs (C-alpha)
    # distances shape: (n_frames, n_pairs)
    distances, _ = md.compute_contacts(traj, contacts='all', scheme=scheme)

    # Mean across frames (axis 0)
    mean_dist_vector = distances.mean(axis=0)

    # Convert vector to square matrix (N_residues x N_residues)
    matrix = squareform(mean_dist_vector)
    return matrix

# --- 3. Matrix Processing ---
print("Processing Native protein (WT)...")
matrix_wt = get_mean_contact_data(wt_xtc, wt_top)

print("\nProcessing Mutant protein (P36R)...")
matrix_mut = get_mean_contact_data(mut_xtc, mut_top)

# --- 4. Difference Calculations ---
cutoff = 0.8  # 0.8 nm = 8 Angstroms

# Contact Maps (Binary)
map_wt = matrix_wt < cutoff
map_mut = matrix_mut < cutoff

# Contact Difference:
# +1: Gained contact in mutant | -1: Lost contact | 0: No change
diff_contact = map_mut.astype(int) - map_wt.astype(int)

# Distance Difference (Continuous)
diff_dist = matrix_mut - matrix_wt

# --- 5. Visualization and Saving ---
fig, ax = plt.subplots(1, 2, figsize=(20, 10))

# Plot 1: Contact Map Difference (Binary)
sns.heatmap(diff_contact, cmap='RdBu_r', center=0, ax=ax[0], square=True)
ax[0].set_title('Contact Map Difference (Binary)\n(Blue: Gained | Red: Lost)')
ax[0].set_xlabel('Residue Index')
ax[0].set_ylabel('Residue Index')

vmax = np.abs(diff_dist).max() * 0.7  # 70% do máximo para dar contraste
# Plot 2: Mean Distance Difference (Continuous)
sns.heatmap(diff_dist, cmap='bwr', center=0, ax=ax[1], square=True, vmin=-vmax, vmax=vmax,)
ax[1].set_title('Mean Distance Difference (nm)\n(Red: Increased | Blue: Decreased)')
ax[1].set_xlabel('Residue Index')
ax[1].set_ylabel('Residue Index')

# Guides for P36R mutation and N-terminal boundary (~res 200)
for a in ax:
    a.axhline(36, color='green', linestyle=':', alpha=0.7, label='P36R Mutation')
    a.axvline(36, color='green', linestyle=':', alpha=0.7)
    a.axhline(200, color='black', lw=1, alpha=0.5, label='N-terminal Limit')
    a.axvline(200, color='black', lw=1, alpha=0.5)

plt.tight_layout()

# Saving with detailed names
filename = "AnnexinA11_P36R_vs_WT_Contact_Distance_Comparison.png"
plt.savefig(filename)
plt.show()
print(f"\nFigure saved at: {filename}")

# --- 6. Diagnostics ---
print(f"\n--- Diagnostic Summary ---")
print(f"Mean Contacts (WT): {np.sum(map_wt)}")
print(f"Mean Contacts (Mutant): {np.sum(map_mut)}")
print(f"Net Change in Contacts: {np.sum(diff_contact)}")
