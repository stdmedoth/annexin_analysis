import glob
import mdtraj as md
import warnings
import numpy as np

import matplotlib.pyplot as plt

from sklearn.decomposition import PCA

# for bioemu
warnings.filterwarnings("ignore", message="Unlikely unit cell vectors detected")


traj = md.load('MutationConformationsBioEmu/out_native/samples.xtc', top='MutationConformationsBioEmu/out_native/topology.pdb')
traj.xyz = traj.xyz / 10.0

ca_indices = traj.topology.select('name CA and resi 0 to 505')
traj_ca = traj.atom_slice(ca_indices)

core_indices_in_slice = traj_ca.topology.select('resi 199 to 505')
traj_ca.superpose(traj_ca, 0, atom_indices=core_indices_in_slice)
fluctuations = md.rmsf(traj_ca, traj_ca, 0)

print(f"Mean RMSF: {fluctuations.mean()*10} Å")


plt.figure(figsize=(10,5))
#plt.plot(fluctuations, color="#2c3e50", linewidth=1.5)

residue_numbers = range(0, len(fluctuations))
plt.plot(residue_numbers, fluctuations*10, color="#2c3e50", linewidth=1.5)
plt.axvspan(0,199, color='gray', alpha=0.2, label='N-Terminal')
plt.axvspan(199,505, color='blue', alpha=0.1, label='Annexin core')
plt.legend()


plt.title('Conformational Profile: RMSF by Residue - WT Annexin A11', fontsize=14)
plt.xlabel('Residue Number', fontsize=12)
plt.ylabel('RMSF ($\AA$)', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.6)

plt.savefig('conformational_profile_wt.png', dpi=300)
plt.show()



# The MDTraj has shape (n_frames, n_atoms, 3) the PCA needs (n_frame, n_atoms * 3)
n_frames = traj_ca.xyz.shape[0]
n_atoms = traj_ca.xyz.shape[1]
coords = traj_ca.xyz.reshape(n_frames, n_atoms * 3)

pca = PCA(n_components=2)
reduced_coords = pca.fit_transform(coords)


plt.figure(figsize=(8,6))
plt.scatter(reduced_coords[:, 0], reduced_coords[:, 1], c=range(n_frames),cmap='viridis', alpha=0.6)
plt.colorbar(label='Frame index')
plt.title('Conformational Profile: PCA - WT Annexin A11', fontsize=14)
plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)') 
plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)') 
plt.grid(True, alpha=0.3) 
plt.savefig('pca_conformational_space_wt.png') 
plt.show()



coords = traj_ca.xyz.reshape(n_frames, n_atoms * 3)
diff = np.linalg.norm(coords[0] - coords[1])
print(f"Distance between frames 0 and 1: {diff:.6f} nm")
unique_frames = np.unique(coords.round(decimals=4), axis=0)
print(f"Number of real unique strutures: {len(unique_frames)}")
