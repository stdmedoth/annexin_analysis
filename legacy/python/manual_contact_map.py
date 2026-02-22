import mdtraj as md
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

base = '1000_samples/MutationConformationsBioEmu'
xtc_path = f'{base}/out_native/core_aligned_samples.xtc'
topology_path = f'{base}/out_native/core_aligned_topology.pdb'

# 1. Load trajectory
traj = md.load(xtc_path, top=topology_path)
coords = traj.xyz  # Shape: (n_frames, n_atoms, 3)

# 2. Define cutoff (MDTraj is in nm, so 8.0 Angstroms = 0.8 nm)
d_c = 0.8
n_frames, n_residues, _ = coords.shape

# 3. Vectorized calculation for all frames
# We initialize a matrix to accumulate the counts (frequency of contact)
contact_frequency = np.zeros((n_residues, n_residues))

for frame in range(n_frames):
    # Get coordinates for the current frame: (505, 3)
    c = coords[frame]

    # Broadcasting Magic:
    # c[:, np.newaxis, :] has shape (505, 1, 3)
    # c[np.newaxis, :, :] has shape (1, 505, 3)
    # The subtraction results in a (505, 505, 3) matrix of all relative vectors
    diff = c[:, np.newaxis, :] - c[np.newaxis, :, :]

    # Calculate Euclidean distance along the last axis (the 3 coordinates)
    # distances shape: (505, 505)
    distances = np.linalg.norm(diff, axis=2)

    # Logic: If distance < cutoff, add to the frequency matrix
    # Note: I changed > to < to reflect a standard Contact Map
    contact_frequency += (distances < d_c).astype(float)

# 4. Average by number of frames
final_matrix = contact_frequency / n_frames

# Plotting
plt.figure(figsize=(10, 8))
sns.heatmap(final_matrix, cmap='viridis') # Viridis is better for 0 to 1 scales
plt.title('Manual Vectorized Contact Map - Annexin A11')
plt.show()
