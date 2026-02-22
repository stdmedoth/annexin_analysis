#!/usr/bin/env python3
"""
RMSF Analysis for Wild-Type Annexin A11

This script computes and visualizes the RMSF (Root Mean Square Fluctuation)
profile for the wild-type Annexin A11 protein ensemble.

RMSF measures the flexibility of each residue, helping identify:
- Highly flexible regions (e.g., N-terminal intrinsically disordered region)
- Stable regions (e.g., Annexin core domain)

Output:
- RMSF profile plot with domain annotations
- Mean RMSF statistics
"""

import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from annexin_analysis import (
    AnnexinConfig,
    TrajectoryLoader,
    ConformationalAnalyzer,
    ConformationalVisualizer
)


def main():
    """Run RMSF analysis for Wild-Type Annexin A11."""
    
    # Initialize configuration and components
    config = AnnexinConfig()
    loader = TrajectoryLoader(config)
    analyzer = ConformationalAnalyzer(config)
    visualizer = ConformationalVisualizer(config)
    
    # Get variant configuration
    wt = config.wt_variant
    print(f"Analyzing: {wt.label}")
    print(f"XTC: {wt.xtc_path}")
    print(f"Topology: {wt.topology_path}")
    
    # Process trajectory
    print("\nLoading and processing trajectory...")
    processed = loader.process_variant(wt)
    
    print(f"  Frames: {processed.n_frames}")
    print(f"  Atoms (CA): {processed.n_atoms}")
    
    # Compute RMSF
    print("\nComputing RMSF...")
    rmsf_result = analyzer.compute_rmsf(processed)
    
    print(f"\nResults:")
    print(f"  Mean RMSF: {rmsf_result.mean_rmsf:.4f} Å")
    
    # Count unique structures
    unique_count = analyzer.count_unique_structures(processed)
    print(f"  Unique structures: {unique_count}")
    
    # Visualize
    print("\nGenerating visualization...")
    visualizer.plot_rmsf(
        rmsf_result,
        title="Conformational Profile: RMSF - Wild-Type Annexin A11",
        filename="rmsf_wt.png"
    )
    
    print("\nAnalysis complete!")


if __name__ == "__main__":
    main()
