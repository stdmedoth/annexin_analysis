#!/usr/bin/env python3
"""
Convert NPZ Files to PDB Format

This utility script converts BioEmu NPZ output files to standard
PDB format for visualization in molecular viewers (PyMOL, VMD, etc.).

BioEmu saves conformational samples in NPZ format with position
arrays and sequence data. This script extracts the first frame
from each NPZ file and writes it as a CA-only PDB.

Usage:
    python convert_npz_to_pdb.py <input_directory> [output_directory]
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from annexin_analysis import NPZConverter


def main():
    """Convert NPZ files to PDB format."""
    
    # Parse arguments
    if len(sys.argv) < 2:
        print("Usage: python convert_npz_to_pdb.py <input_directory> [output_directory]")
        print("\nExample:")
        print("  python convert_npz_to_pdb.py MutationConformationsBioEmu/out_native/")
        sys.exit(1)
    
    input_dir = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) > 2 else None
    
    print("NPZ to PDB Converter")
    print("=" * 40)
    print(f"Input directory: {input_dir}")
    print(f"Output directory: {output_dir or 'same as input'}")
    
    # Initialize converter
    converter = NPZConverter(output_dir)
    
    # Convert all files
    print("\nConverting files...")
    output_files = converter.convert_directory(input_dir)
    
    print(f"\nConverted {len(output_files)} files")
    print("\nDone!")


if __name__ == "__main__":
    main()
