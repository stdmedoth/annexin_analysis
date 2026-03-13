#!/usr/bin/env python3
"""
Convert PDB Files to XTC Format

Usage:
    python convert_pdbs_to_xtc.py <input_directory> [output_directory]
"""

import sys
import glob

from pathlib import Path
import mdtraj as md

sys.path.insert(0, str(Path(__file__).parent.parent))

from annexin_analysis import (
    AnnexinConfig,
    TrajectoryLoader,
)


def main():
    """Convert NPZ files to PDB format."""
    config = AnnexinConfig()
    
    # Parse arguments
    if len(sys.argv) < 2:
        print("Usage: python convert_pdbs_to_xtc.py <input_directory>")
        print("\nExample:")
        print("  python convert_pdbs_to_xtc.py MutationConformationsBioEmu/out_native/")
        sys.exit(1)
    
    input_dir = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) > 2 else None
    
    print("PCBs to XTC Converter")
    print("=" * 40)
    print(f"Input directory: {input_dir}")
    
    # Initialize converter
    pdbs_files = glob.glob(str(Path(input_dir) / "*.pdb"))

    output_dir = f"{config.output_dir}"            

    
    print("\nConverting files...")
    
    xtc_filename = 'converted_samples.xtc'
    pdb_filename = 'converted_topology.pdb'

    xtc_path = f"{config.output_dir}/{xtc_filename}"
    pdb_path = f"{config.output_dir}/{pdb_filename}"

    traj = md.load(pdbs_files)
    traj.save_xtc(xtc_path)
    traj[0].save_pdb(pdb_path)
    
    print(f"\nConverted {len(pdbs_files)} files")
    print(f"\nXTC Path: {xtc_path}")
    print(f"\nPDB Path: {pdb_path}")
    
    
    print("\nDone!")


if __name__ == "__main__":
    main()
