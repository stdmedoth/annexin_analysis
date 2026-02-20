import glob
from pathlib import Path
import mdtraj as md
import warnings
import numpy as np
import matplotlib.pyplot as plt

# Configurações para bioemu
warnings.filterwarnings("ignore", message="Unlikely unit cell vectors detected")

def process_trajectory(path, topology_path, global_reference):
    """Carrega a trajetória, realiza o alinhamento e calcula o RMSF."""
    traj = md.load(path, top=topology_path)

    traj = traj[:737]

    traj.xyz = traj.xyz / 10.0

    ca_indices = traj.topology.select('name CA')
    traj_ca = traj.atom_slice(ca_indices)

    # Alinhamento pelo core da Anexina (resi 199 a 505)x
    core_indices = traj_ca.topology.select('resi 199 to 505')
    traj_ca.superpose(global_reference, 0, atom_indices=core_indices)

    file = Path(path)
    traj_ca.save_xtc(f'test_pymol/{file.parent.name}_samples.xtc')
    traj_ca[0].save_pdb(f'test_pymol/{file.parent.name}_topology.pdb')

    # Cálculo do RMSF (multiplicado por 10 para voltar para Å se necessário)
    fluctuations = md.rmsf(traj_ca, global_reference) * 10
    return fluctuations

# Configuração das pastas e labels
base_path = '1000_samples/MutationConformationsBioEmu'
variants = {
    'Filtered': 'out_native',
    'Luis Felipe': 'luis_felipe',
    # Adicione novas mutações aqui: 'MUTAÇÃO': 'PASTA'
}

plt.figure(figsize=(12, 6))
colors = ['#2c3e50', '#e74c3c', '#27ae60', '#f39c12']


pdb_ref_path = f'AF-P50995-F1-model_v6.pdb'
global_reference = md.load(pdb_ref_path)
global_reference.xyz = global_reference.xyz / 10.0
# Seleção de átomos CA (Resíduos 0 a 505 de acordo com seu script original)
ca_indices = global_reference.topology.select('name CA and resi 0 to 505')
global_reference = global_reference.atom_slice(ca_indices)
# Alinhamento pelo core da Anexina (resi 199 a 505)
ref_traj_core_indices = global_reference.topology.select('resi 199 to 505')
global_reference.superpose(global_reference, 0, atom_indices=ref_traj_core_indices)



for i, (label, folder) in enumerate(variants.items()):
    xtc_path = f'{base_path}/{folder}/samples.xtc'
    pdb_path = f'{base_path}/{folder}/topology.pdb'

    try:
        rmsf_data = process_trajectory(xtc_path, pdb_path, global_reference)
        residues = range(len(rmsf_data))

        plt.plot(residues, rmsf_data, label=label, color=colors[i % len(colors)], linewidth=1.5, alpha=0.8)
    except Exception as e:
        print(f"Erro ao processar {label}: {e}")

# Áreas de destaque (Domínios da Anexina A11)
plt.axvspan(0, 199, color='gray', alpha=0.1, label='N-Terminal')
plt.axvspan(199, 505, color='blue', alpha=0.05, label='Annexin Core')

# Estética do Gráfico
plt.title('Comparativo de Perfil Conformacional: RMSF - Anexina A11', fontsize=14)
plt.xlabel('Número do Resíduo', fontsize=12)
plt.ylabel('RMSF ($Å$)', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(loc='upper right')

plt.tight_layout()
plt.savefig('rmsf_comparison.png', dpi=300)
plt.show()
