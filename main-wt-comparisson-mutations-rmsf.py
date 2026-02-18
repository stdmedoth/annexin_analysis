import glob
import mdtraj as md
import warnings
import numpy as np
import matplotlib.pyplot as plt

# Configurações para bioemu
warnings.filterwarnings("ignore", message="Unlikely unit cell vectors detected")

def process_trajectory(path, topology_path):
    """Carrega a trajetória, realiza o alinhamento e calcula o RMSF."""
    traj = md.load(path, top=topology_path)
    # Conversão de nm para Angstrom (BioEmu costuma salvar em nm)
    traj.xyz = traj.xyz / 10.0 
    
    # Seleção de átomos CA (Resíduos 0 a 505 de acordo com seu script original)
    ca_indices = traj.topology.select('name CA and resi 0 to 505')
    traj_ca = traj.atom_slice(ca_indices)
    
    # Alinhamento pelo core da Anexina (resi 199 a 505)
    core_indices = traj_ca.topology.select('resi 199 to 505')
    traj_ca.superpose(traj_ca, 0, atom_indices=core_indices)
    
    # Cálculo do RMSF (multiplicado por 10 para voltar para Å se necessário)
    fluctuations = md.rmsf(traj_ca, traj_ca, 0) * 10
    return fluctuations

# Configuração das pastas e labels
base_path = 'MutationConformationsBioEmu'
variants = {
    'WT': 'out_native',
    'P36R': 'out_mutant_P36R',
    # Adicione novas mutações aqui: 'MUTAÇÃO': 'PASTA'
}

plt.figure(figsize=(12, 6))
colors = ['#2c3e50', '#e74c3c', '#27ae60', '#f39c12']

for i, (label, folder) in enumerate(variants.items()):
    xtc_path = f'{base_path}/{folder}/samples.xtc'
    pdb_path = f'{base_path}/{folder}/topology.pdb'
    
    try:
        print(f"Processando {label}...")
        rmsf_data = process_trajectory(xtc_path, pdb_path)
        residues = range(len(rmsf_data))
        
        plt.plot(residues, rmsf_data, label=label, color=colors[i % len(colors)], linewidth=1.5, alpha=0.8)
        print(f"Média RMSF {label}: {rmsf_data.mean():.4f} Å")
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
plt.savefig('rmsf_comparison_mutations.png', dpi=300)
plt.show()
