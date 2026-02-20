import mdtraj as md

traj = md.load('1000_samples/MutationConformationsBioEmu/out_native/samples.xtc', top='1000_samples/MutationConformationsBioEmu/out_native/topology.pdb')

print(traj.xyz.shape)
