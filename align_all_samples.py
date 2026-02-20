import mdtraj as md
from pathlib import Path

base = '1000_samples/MutationConformationsBioEmu'

# bioemu luis felipe wt samples
#sample_xtc = f'{base}/luis_felipe/samples_luisfelipe10k.xtc'
#topology = f'{base}/luis_felipe/topology.pdb'

# bioemu calisto wt samples
#sample_xtc = f'{base}/out_native/samples.xtc'
#topology = f'{base}/out_native/topology.pdb'

# bioemu calisto P36R samples
sample_xtc = f'{base}/out_mutant_P36R/samples.xtc'
topology = f'{base}/out_mutant_P36R/topology.pdb'


core_select = '196 to 505' # Annexin A11

traj = md.load(sample_xtc, top=topology)

ca_indices = traj.topology.select("name CA")
traj_ca = traj.atom_slice(ca_indices)

core_indices = traj_ca.topology.select(f'resi {core_select}')
traj_ca.superpose(traj_ca, 0, atom_indices=core_indices)

file = Path(sample_xtc)


traj_ca.save_xtc(f'{file.parent}/core_aligned_samples.xtc')
traj_ca[0].save_pdb(f'{file.parent}/core_aligned_topology.pdb')
