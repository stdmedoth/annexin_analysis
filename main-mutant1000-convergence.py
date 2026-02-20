import glob
import mdtraj as md
import warnings
import numpy as np

import matplotlib.pyplot as plt

from sklearn.decomposition import PCA

# for bioemu
warnings.filterwarnings("ignore", message="Unlikely unit cell vectors detected")

samples_qnt = [i for i in range(10, 1010, 10)]
#print(samples_qnt)
#quit()

full_traj = md.load('1000_samples/MutationConformationsBioEmu/out_mutant_P36R/samples.xtc', top='1000_samples/MutationConformationsBioEmu/out_mutant_P36R/topology.pdb')
ca_indices = full_traj.topology.select('name CA and resi 0 to 505')
traj_ca = full_traj.atom_slice(ca_indices)
core_indices = traj_ca.topology.select('resi 199 to 505')

# Reference: RMSF of the entire ensemble
traj_ca.superpose(traj_ca, 0, atom_indices=core_indices)
final_rmsf = md.rmsf(traj_ca, traj_ca, 0)

diffs = []
for n in samples_qnt:
    # Slice the already loaded trajectory
    sub_traj = traj_ca[0:n]
    sub_traj.superpose(sub_traj, 0, atom_indices=core_indices)
    current_rmsf = md.rmsf(sub_traj, sub_traj, 0)

    # Calculate how much the current RMSF profile differs from the final one
    # Using Mean Absolute Error (MAE)
    #error = np.mean(np.abs(current_rmsf - final_rmsf))
    #diffs.append(error)
    diff_norm = np.linalg.norm(current_rmsf - final_rmsf)
    diffs.append(diff_norm)

# Plotting the error decay
plt.plot(samples_qnt, diffs, 'o-')
plt.ylabel("Mean Absolute Difference in RMSF (nm)")
plt.xlabel("Number of Samples")
plt.show();
