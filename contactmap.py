import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.distance import squareform


base = '1000_samples/MutationConformationsBioEmu'
# bioemu luis felipe wt samples
#xtc_path = f'{base}/luis_felipe/samples_luisfelipe10k.xtc'
#topology_path = f'{base}/luis_felipe/topology.pdb'

# bioemu calisto wt samples
xtc_path = f'{base}/out_native/core_aligned_samples.xtc'
topology_path = f'{base}/out_native/core_aligned_topology.pdb'

# bioemu calisto P36R samples
xtc_mutation_path = f'{base}/out_mutant_P36R/core_aligned_samples.xtc'
topology_mutation_path = f'{base}/out_mutant_P36R/core_aligned_topology.pdb'


wt_traj = md.load(xtc_path, top=topology_path)
mut_traj = md.load(xtc_mutation_path, top=topology_mutation_path)

wt_distances, wt_residues_pairs = md.compute_contacts(wt_traj, contacts='all', scheme='ca')
mut_distances, mut_residues_pairs = md.compute_contacts(mut_traj, contacts='all', scheme='ca')

wt_mean_distance = wt_distances.mean(axis=0)
mut_mean_distance = mut_distances.mean(axis=0)

wt_matriz_dist = squareform(wt_mean_distance)
mut_matriz_dist = squareform(mut_mean_distance)


cutoff = 0.8
wt_contact_map = wt_matriz_dist < cutoff
mut_contact_map = mut_matriz_dist < cutoff

print(f"Total de contatos Nativa: {np.sum(wt_contact_map)}")
print(f"Total de contatos Mutada: {np.sum(mut_contact_map)}")

diff_map = wt_contact_map.astype(int) - wt_contact_map.astype(int)

plt.figure(figsize=(10, 8))
sns.heatmap(diff_map, cmap='RdBu_r', center=0, cbar=True)

plt.title('Contact Map Comparisson - Annexin A11 WT vs P36R')
plt.savefig('contact_map_WT_P36R_1000samples.png')
plt.xlabel('Residue i')
plt.ylabel('Residue j')
plt.show()


diff_matrix = wt_matriz_dist - mut_matriz_dist

plt.figure(figsize=(10, 8))

# Usamos o inverso para que o que está perto fique "quente" (opcional)
sns.heatmap(diff_matrix, cmap='viridis_r', square=True,
            xticklabels=100, yticklabels=100)
plt.axhline(200, color='white', linestyle='--', alpha=0.6)
plt.axvline(200, color='white', linestyle='--', alpha=0.6)

plt.title('Average Distance Matrix (Comparisson WT vs P36R) ($nm$)')
plt.savefig('average_distance_matrix_WT_R36R_1000samples.png')
plt.show()
