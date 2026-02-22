import mdtraj as md
import numpy as np

# Carrega sua trajetória e a topologia
traj = md.load('1000_samples/MutationConformationsBioEmu/out_native/samples.xtc', top='1000_samples/MutationConformationsBioEmu/out_native/topology.pdb')

# Exemplo: Filtrar por raio de giração (evitar modelos muito expandidos/desdobrados)
rg = md.compute_rg(traj)
# Mantém apenas modelos com Rg dentro de um desvio padrão da média
indices_bons = np.where(rg < (np.mean(rg) + np.std(rg)))[0]

# Salva o novo .xtc apenas com as conformações "estáveis"
traj_filtrada = traj[indices_bons]
traj_filtrada.save_xtc('filtered_samples.xtc')
