#!/usr/bin/env python3
"""
Create an representative cluster use KMeans


Output:
- Aligned PDB with representative trajectory
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
    variants = [config.wt_variant]
    
    print("Trajectory Alignment Utility")
    print("=" * 40)
    
    for variant in variants:
        print(f"\nProcessing: {variant.label}")
        print(f"  Input XTC: {variant.xtc_path}")
        
        try:
            # Determine output directory
            output_dir = Path(variant.xtc_path).parent
            
            # Save aligned trajectory
            kmeans, labels = loader.create_representative_cluster(
                variant,
                output_dir=output_dir,
            )
            
            print(f"Representative pdb file created on {output_dir}!")
        except Exception as e:
            print(f"  Error: {e}")
    
    print("\n" + "=" * 40)


if __name__ == "__main__":
    main()
