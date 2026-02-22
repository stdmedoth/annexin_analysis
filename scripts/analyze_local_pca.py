#!/usr/bin/env python3
"""
Local PCA Analysis: N-Terminal Region Dynamics

This script performs PCA focused on a local region (N-terminal)
while aligning by the stable core region.

This approach is useful when:
- The global PCA is dominated by large-scale motions
- You want to focus on specific regions of interest
- The mutation is located in a particular region (e.g., P36 in N-terminal)

The analysis aligns all structures by the stable Annexin core,
then extracts only the N-terminal coordinates for PCA analysis.

Output:
- Local PCA projection of N-terminal dynamics
- Comparison between WT and P36R local conformations
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from annexin_analysis import (
    AnnexinConfig,
    ConformationalAnalyzer,
    ConformationalVisualizer
)
from annexin_analysis.trajectory import RegionalTrajectoryLoader


def main():
    """Run local PCA analysis on N-terminal region."""
    
    # Initialize components
    config = AnnexinConfig()
    loader = RegionalTrajectoryLoader(config)
    analyzer = ConformationalAnalyzer(config)
    visualizer = ConformationalVisualizer(config)
    
    # Define regions
    align_selection = "resi 199 to 505"  # Align by core
    analysis_selection = "resi 0 to 50 and name CA"  # Analyze N-terminal
    
    # Get variants
    wt = config.wt_variant
    p36r = config.p36r_variant
    
    print("Local PCA Analysis: N-Terminal Region")
    print("=" * 40)
    print(f"Alignment region: {align_selection}")
    print(f"Analysis region: {analysis_selection}")
    
    # Process WT with regional extraction
    print("\nProcessing WT...")
    wt_processed = loader.process_for_regional_analysis(
        wt,
        align_selection=align_selection,
        analysis_selection=analysis_selection
    )
    
    print(f"  Frames: {wt_processed.n_frames}")
    print(f"  Atoms in region: {wt_processed.n_atoms}")
    
    # Compute PCA on WT
    print("\nComputing PCA on WT...")
    wt_pca = analyzer.compute_pca(wt_processed)
    
    # Process P36R with cross-alignment
    print("\nProcessing P36R...")
    p36r_processed = loader.process_for_regional_analysis(
        p36r,
        align_selection=align_selection,
        analysis_selection=analysis_selection,
        reference_trajectory=wt_processed.trajectory
    )
    
    # Project P36R onto WT space
    print("Projecting P36R onto WT PCA space...")
    p36r_pca = analyzer.project_to_pca_space(p36r_processed, wt_pca)
    
    # Visualize
    print("\nGenerating visualization...")
    visualizer.plot_pca_comparison(
        wt_pca,
        {"Mutant P36R": p36r_pca},
        wt_label="Wild Type",
        title="Local PCA: N-Terminal Dynamics (0-50)\nAligned by Core (199-505)",
        filename="pca_local_nterminal.png"
    )
    
    # Statistics
    print("\n" + "=" * 40)
    print("Results:")
    print("=" * 40)
    print(f"Variance explained: {wt_pca.total_variance_explained:.1f}%")
    print(f"  PC1: {wt_pca.pc1_variance:.1f}%")
    print(f"  PC2: {wt_pca.pc2_variance:.1f}%")
    
    print("\nAnalysis complete!")


if __name__ == "__main__":
    main()
