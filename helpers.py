import glob
import numpy as np
from pathlib import Path

class Helper:

    def map_1_to_3(self, letter):
        map = {
            'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE',
            'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'K': 'LYS', 'L': 'LEU',
            'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG',
            'S': 'SER', 'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'
        }

        return map[letter]
    

    def convert_npz_to_pdb(self, npz_filename):
        npzdata = np.load(npz_filename)
        
        positions = npzdata['pos'][0]
        sequence = npzdata['sequence'].item()
        
        output_filename = f"{Path(npz_filename).stem}.pdb"

        
        with open(output_filename, "w") as f:
            for i, (position, aa) in enumerate(zip(positions, sequence)):
                #print(i, position, aa)
                #print(self.map_1_to_3(aa))

                x, y, z = position
                res = self.map_1_to_3(aa)

                line = f"ATOM  {i+1:>5}  CA  {res:>3} A{i+1:>4}    {x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00  0.00           C\n"
                f.write(line)

            f.write("END\n")





helper = Helper()
files = glob.glob("MutationConformationsBioEmu/out_native/*.npz")
print(files)
for file in files:
    helper.convert_npz_to_pdb(file)