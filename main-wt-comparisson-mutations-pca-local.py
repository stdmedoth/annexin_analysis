import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import warnings

# Ignorar avisos de célula unitária do BioEmu
warnings.filterwarnings("ignore")

def get_region_pca_data(path, top_path, align_selection, pca_selection, ref_traj=None):
    """
    Alinha pelo Core, mas retorna coordenadas apenas da região de interesse (N-term).
    """
    # Carrega trajetória
    traj = md.load(path, top=top_path)
    traj.xyz = traj.xyz * 10.0  # Å -> nm (BioEmu output fix)

    # Seleção de átomos para ALINHAMENTO (Sempre o Core estável)
    # Nota: MDTraj é 0-indexed.
    core_indices = traj.topology.select(align_selection)

    # Seleção de átomos para ANÁLISE (A região flexível 0-100)
    pca_indices = traj.topology.select(pca_selection)

    # Cria sub-trajetórias fatiadas para manipulação
    # Precisamos da traj inteira para alinhar, mas depois cortaremos
    traj_ca = traj.atom_slice(traj.topology.select('name CA'))

    # Atualiza índices baseados apenas nos C-alphas agora
    core_indices_ca = traj_ca.topology.select(align_selection)
    pca_indices_ca = traj_ca.topology.select(pca_selection)

    # 1. Superposição (Alinhamento) baseada no CORE
    if ref_traj is None:
        traj_ca.superpose(traj_ca, 0, atom_indices=core_indices_ca)
        ref_out = traj_ca # Guarda essa ref para a mutante usar
    else:
        traj_ca.superpose(ref_traj, 0, atom_indices=core_indices_ca)
        ref_out = None

    # 2. Extração das coordenadas APENAS da região de interesse (0-100)
    # Fatia a trajetória alinhada para pegar só a cauda
    traj_region = traj_ca.atom_slice(pca_indices_ca)

    n_frames = traj_region.n_frames
    n_atoms = traj_region.n_atoms

    # Retorna (trajetória_ref, coordenadas_achatadas)
    return ref_out, traj_region.xyz.reshape(n_frames, n_atoms * 3)

# --- CONFIGURAÇÃO ---
base = '1000_samples/MutationConformationsBioEmu'
# Alinhar pelo Core Rígido
align_sel = 'resi 199 to 505'
# Calcular PCA na Cauda (focando onde o RMSF deu diferença)
pca_sel = 'resi 0 to 50 and name CA'

# 1. Processar WT (Nativa)
print("Processando WT...")
wt_ref, wt_coords = get_region_pca_data(
    f'{base}/out_native/samples.xtc',
    f'{base}/out_native/topology.pdb',
    align_sel,
    pca_sel
)

# 2. Treinar PCA na região N-term da WT
print("Calculando PCA...")
pca = PCA(n_components=2)
wt_reduced = pca.fit_transform(wt_coords)

# 3. Processar Mutante P36R
print("Processando Mutante...")
_, mut_coords = get_region_pca_data(
    f'{base}/out_mutant_P36R/samples.xtc',
    f'{base}/out_mutant_P36R/topology.pdb',
    align_sel,
    pca_sel,
    ref_traj=wt_ref # Usa a WT alinhada como referência
)
mut_reduced = pca.transform(mut_coords)

# --- PLOT ---
plt.figure(figsize=(10, 8))

# Plot WT
plt.scatter(wt_reduced[:, 0], wt_reduced[:, 1], c='blue', alpha=0.4, label='WT (N-Term 35-37)', s=15, edgecolors='none')

# Plot Mutante
plt.scatter(mut_reduced[:, 0], mut_reduced[:, 1], c='red', alpha=0.4, label='P36R (N-Term 35-37)', s=15, edgecolors='none')

plt.title('PCA Local: Dinâmica da Região N-Terminal (35-37)\nAlinhamento pelo Core (199-505)', fontsize=14)
plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)', fontsize=12)
plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)', fontsize=12)
plt.legend(fontsize=12)
plt.grid(True, alpha=0.3, linestyle='--')

plt.tight_layout()
plt.savefig('pca_local_nterm_0_100.png', dpi=300)
plt.show()

print(f"Variância Explicada (Local): {(pca.explained_variance_ratio_.sum())*100:.2f}%")
