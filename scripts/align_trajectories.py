#!/usr/bin/env python3
"""
Align and Save Trajectories

This utility script loads trajectories, aligns them by the Annexin core,
and saves the aligned trajectories for downstream analysis.

Aligned trajectories are useful for:
- Visualization in molecular viewers
- Pre-processing before contact map analysis
- Combining with other trajectories

Output:
- Aligned XTC trajectory
- Topology PDB (first frame)
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from annexin_analysis import (
    AnnexinConfig,
    TrajectoryLoader
)


def main():
    """Align and save trajectories for WT and P36R."""
    
    # Initialize components
    config = AnnexinConfig()
    loader = TrajectoryLoader(config)
    
    # Define variants
    variants = [config.wt_variant, config.p36r_variant]
    
    print("Trajectory Alignment Utility")
    print("=" * 40)
    
    for variant in variants:
        print(f"\nProcessing: {variant.label}")
        print(f"  Input XTC: {variant.xtc_path}")
        
        try:
            # Process variant
            processed = loader.process_variant(variant)
            
            # Determine output directory
            output_dir = Path(variant.xtc_path).parent
            
            # Save aligned trajectory
            xtc_out, pdb_out = loader.save_aligned_trajectory(
                processed,
                output_dir=output_dir,
                prefix="core_aligned"
            )
            
            print(f"  Output XTC: {xtc_out}")
            print(f"  Output PDB: {pdb_out}")
            print(f"  Frames: {processed.n_frames}")
            print(f"  Atoms: {processed.n_atoms}")
            
        except Exception as e:
            print(f"  Error: {e}")
    
    print("\n" + "=" * 40)
    print("Alignment complete!")


if __name__ == "__main__":
    main()
