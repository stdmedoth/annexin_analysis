import glob
import mdtraj as md


full_ref = md.load('AF-P50995-F1-model_v6.pdb')

ca_indices = full_ref.topology.select('name CA')
ref_ca = full_ref.atom_slice(ca_indices)

files = sorted(glob.glob('batch_*.pdb'))
target_traj = md.load(files)

target_traj.superpose(ref_ca, frame=0)

core_indices = target_traj.topology.select("resid 199 to 504")

#rmsf = md.rmsf(target_traj, ref_ca)
rmsd = md.rmsd(target_traj, ref_ca, atom_indices=core_indices)


#print(rmsf, rmsd)

#print(rmsf)
print(rmsd)
