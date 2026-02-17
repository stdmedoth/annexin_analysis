import glob
import mdtraj as md
import warnings
import numpy as np

import matplotlib.pyplot as plt

# for bioemu
warnings.filterwarnings("ignore", message="Unlikely unit cell vectors detected")


traj = md.load('pdbs/samples.xtc', top='pdbs/topology.pdb')
ca_indices = traj.topology.select('name CA')
traj_ca = traj.atom_slice(ca_indices)
traj_ca.superpose(traj_ca, 0)
fluctuations = md.rmsf(traj_ca, traj_ca, 0)

print(fluctuations.mean())


plt.figure(figsize=(10,5))
#plt.plot(fluctuations, color="#2c3e50", linewidth=1.5)

residue_numbers = range(0, len(fluctuations))
plt.plot(residue_numbers, fluctuations, color="#2c3e50", linewidth=1.5)
plt.axvspan(0,199, color='gray', alpha=0.2, label='N-Terminal')
plt.axvspan(199,505, color='blue', alpha=0.1, label='Annexin core')
plt.legend()


plt.title('Conformational Profile: RMSF by Residue - Annexin A11', fontsize=14)
plt.xlabel('Residue Number', fontsize=12)
plt.ylabel('RMSF ($\AA$)', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.6)

plt.savefig('conformational_profile.png', dpi=300)
plt.show()

