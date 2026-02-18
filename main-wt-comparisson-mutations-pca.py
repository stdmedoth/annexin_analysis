import glob
import mdtraj as md
import warnings
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

# Configurações para bioemu
warnings.filterwarnings("ignore", message="Unlikely unit cell vectors detected")

def get_aligned_coords(path, top_path, ref_traj=None):
    """Carrega, alinha pelo core e retorna as coordenadas (n_frames, n_atoms * 3)."""
    traj = md.load(path, top=top_path)
    traj.xyz = traj.xyz / 10.0 # nm -> Å
    
    ca_indices = traj.topology.select('name CA and resi 0 to 505')
    traj_ca = traj.atom_slice(ca_indices)
    core_indices = traj_ca.topology.select('resi 199 to 505')
    
    # Se não houver referencial, alinha nela mesma (usado para a WT inicial)
    if ref_traj is None:
        traj_ca.superpose(traj_ca, 0, atom_indices=core_indices)
    else:
        # Alinha a mutante seguindo o referencial da WT para consistência no PCA
        traj_ca.superpose(ref_traj, 0, atom_indices=core_indices)
    
    n_frames = traj_ca.xyz.shape[0]
    n_atoms = traj_ca.xyz.shape[1]
    return traj_ca, traj_ca.xyz.reshape(n_frames, n_atoms * 3)

# 1. Processamento da WT (Referencial)
base = 'MutationConformationsBioEmu'
wt_traj, wt_coords = get_aligned_coords(f'{base}/out_native/samples.xtc', f'{base}/out_native/topology.pdb')

# 2. Treinar o PCA apenas na WT
pca = PCA(n_components=2)
wt_reduced = pca.fit_transform(wt_coords)

# 3. Processar Mutação e Projetar no Espaço da WT
# Aqui você pode repetir para outras mutações em um loop se quiser
mut_label = 'P36R'
mut_traj, mut_coords = get_aligned_coords(
    f'{base}/out_mutant_{mut_label}/samples.xtc', 
    f'{base}/out_mutant_{mut_label}/topology.pdb',
    ref_traj=wt_traj # Alinhamento cruzado
)
mut_reduced = pca.transform(mut_coords) # Projeção

# --- VISUALIZAÇÃO ---
plt.figure(figsize=(10, 7))

# Plot WT
plt.scatter(wt_reduced[:, 0], wt_reduced[:, 1], c='blue', alpha=0.4, label='WT (Nativa)', s=10)

# Plot Mutação
plt.scatter(mut_reduced[:, 0], mut_reduced[:, 1], c='red', alpha=0.4, label=f'Mutante {mut_label}', s=10)

# Estética
plt.title(f'Comparação de Espaço Conformacional: WT vs {mut_label}', fontsize=14)
plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)')
plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)')
plt.legend()
plt.grid(True, alpha=0.3)

plt.savefig('pca_comparison_wt_mutant.png', dpi=300)
plt.show()

# Opcional: Print da variância explicada
print(f"Variância Total Explicada (PC1+PC2): {(pca.explained_variance_ratio_.sum())*100:.2f}%")
