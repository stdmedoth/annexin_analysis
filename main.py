import glob
import mdtraj as md
import warnings


import matplotlib.pyplot as plt

# for bioemu
warnings.filterwarnings("ignore", message="Unlikely unit cell vectors detected")

full_ref = md.load('pdbs/AF-P50995-F1-model_v6.pdb')

ca_indices = full_ref.topology.select('name CA')
ref_ca = full_ref.atom_slice(ca_indices)

files = sorted(glob.glob('pdbs/batch_*.pdb'))
target_traj = md.load(files)

target_traj.superpose(ref_ca, frame=0)

#core_indices = target_traj.topology.select("resid 197 to 505")
core_indices = target_traj.topology.select("resid 199 to 504")

core_traj = target_traj.atom_slice(core_indices)
ref_core = ref_ca.atom_slice(core_indices)

core_traj.superpose(ref_core)

rmsf = md.rmsf(core_traj, ref_core)
#rmsf = md.rmsf(core_traj, ref_traj, atom_indices=core_indices)
#rmsd = md.rmsd(target_traj, ref_ca, atom_indices=core_indices)

plt.figure(figsize=(10,5))
#plt.plot(rmsf, color="#2c3e50", linewidth=1.5)
residue_numbers = range(199, 199 + len(rmsf))
plt.plot(residue_numbers, rmsf * 10, color="#2c3e50", linewidth=1.5)
#plt.axvspan(0,199, color='gray', alpha=0.2, label='N-Terminal')
#plt.axvspan(199,505, color='blue', alpha=0.1, label='Annexin core')


plt.title('Conformational Profile: RMSF by Residue - Annexin A11', fontsize=14)
plt.xlabel('Residue Number', fontsize=12)
plt.ylabel('RMSF ($\AA$)', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.6)

plt.savefig('conformational_profile.png', dpi=300)
plt.show()

#mean = rmsd.mean()

#print(rmsf, rmsd)

#print(rmsf)
#print(rmsd, mean)
