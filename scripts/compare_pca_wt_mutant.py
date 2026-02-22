#!/usr/bin/env python3
"""
PCA Comparison: Wild-Type vs Mutant P36R

This script compares the conformational space of wild-type Annexin A11
and the P36R mutant using PCA.

The key methodology:
1. Train PCA on the WT ensemble (defines the conformational space)
2. Project the P36R mutant onto the same space
3. Visualize to show if mutation shifts the conformational landscape

This provides evidence that the P36R mutation alters the accessible
conformational states of Annexin A11.

Output:
- 2D PCA plot with WT and P36R overlaid
- Variance statistics
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from annexin_analysis import (
    AnnexinConfig,
    TrajectoryLoader,
    ConformationalAnalyzer,
    ConformationalVisualizer
)
from annexin_analysis.analysis import ComparativeAnalyzer


def main():
    """Compare PCA projections of WT and P36R mutant."""
    
    # Initialize components
    config = AnnexinConfig()
    comparative = ComparativeAnalyzer(config)
    visualizer = ConformationalVisualizer(config)
    
    # Define variants
    wt = config.wt_variant
    p36r = config.p36r_variant
    
    print("PCA Comparison: WT vs P36R")
    print("=" * 40)
    print(f"Wild Type: {wt.xtc_path}")
    print(f"Mutant:    {p36r.xtc_path}")
    
    # Perform comparative PCA analysis
    print("\nProcessing...")
    wt_pca, mutant_pcas = comparative.compare_pca(wt, [p36r])
    
    # Get P36R result
    p36r_pca = mutant_pcas.get("p36r")
    
    if p36r_pca is None:
        print("Error: Could not process P36R mutant")
        return
    
    # Visualize
    print("\nGenerating comparison plot...")
    visualizer.plot_pca_comparison(
        wt_pca,
        {"Mutant P36R": p36r_pca},
        wt_label="Wild Type",
        title="Conformational Space Comparison: WT vs P36R",
        filename="pca_comparison_wt_p36r.png"
    )
    
    # Statistics
    print("\n" + "=" * 40)
    print("PCA Statistics:")
    print("=" * 40)
    print(f"PC1: {wt_pca.pc1_variance:.1f}% variance")
    print(f"PC2: {wt_pca.pc2_variance:.1f}% variance")
    print(f"Total: {wt_pca.total_variance_explained:.1f}% variance explained")
    
    print(f"\nWT frames:   {wt_pca.n_frames}")
    print(f"P36R frames: {p36r_pca.n_frames}")
    
    print("\nAnalysis complete!")


if __name__ == "__main__":
    main()
