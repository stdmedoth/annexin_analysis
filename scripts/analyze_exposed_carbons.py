#!/usr/bin/env python3
"""
Find exposed carbons for Annexin A11 Ensemble


Output:
- Array with the exposed carbons
"""

import numpy as np
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from annexin_analysis import (
    AnnexinConfig,
    TrajectoryLoader,
    ConformationalAnalyzer,
    ConformationalVisualizer,
)


def analyze_exposed_carbons():
    """
    Return the carbons array for a specific variant.
    """
    # Initialize components
    config = AnnexinConfig()
    loader = TrajectoryLoader(config)
    analyzer = ConformationalAnalyzer(config)
    visualizer = ConformationalVisualizer(config)


    variant = config.wt_variant
    
    print("Exposed Carbons")
    print("=" * 40)
    print(f"Variant : {variant.xtc_path}")
    
    print(f"Analysis: {variant.label}")
    print("=" * 40)

    # Process trajectory
    print("Loading trajectory...")
    variant_traj = loader.process_variant(variant, extract_ca_atoms=False)

    print(f"  Total frames: {variant_traj.n_frames}")

    # Run convergence analysis
    print("\nFinding exposed residues ...")
    exposed_residues = analyzer.compute_exposed_residues(variant_traj)    
    print(exposed_residues)

    #visualizer.dump_data(exposed_residues, f'exposed_carbons_{variant.name.lower()}.csv')

    return exposed_residues





def main():
    """Find Exposed Carbons..."""

    print("=" * 50)
    
    # Analyze WT
    print("\n")
    exposed_carbons = analyze_exposed_carbons()
    

    print("\nAnalysis complete!")


if __name__ == "__main__":
    main()
