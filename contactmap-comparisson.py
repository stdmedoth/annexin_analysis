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
#xtc_path = f'{base}/out_native/samples.xtc'
#topology_path = f'{base}/out_native/topology.pdb'

# bioemu calisto P36R samples
xtc_path = f'{base}/out_mutant_P36R/samples.xtc'
topology_path = f'{base}/out_mutant_P36R/topology.pdb'


traj = md.load(xtc_path, top=topology_path)

distances, residues_pairs = md.compute_contacts(traj, contacts='all', scheme='ca')

mean_distance = distances.mean(axis=0)


matriz_dist = squareform(mean_distance)

cutoff = 0.8
mapa_contato = matriz_dist < cutoff


plt.figure(figsize=(8, 6))
sns.heatmap(mapa_contato, cmap='Greys', cbar=False, square=True)

plt.title('Contact Map - Annexin A11')
plt.savefig('contact_map_wt_1000samples.png')
plt.xlabel('Residue i')
plt.ylabel('Residue j')
plt.show()


plt.figure(figsize=(10, 8))
# Usamos o inverso para que o que está perto fique "quente" (opcional)
sns.heatmap(matriz_dist, cmap='viridis_r', square=True,
            xticklabels=50, yticklabels=50)

plt.title('Average Distance Matrix ($nm$)')
plt.savefig('average_distance_matrix.png')
plt.show()
