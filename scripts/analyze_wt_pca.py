#!/usr/bin/env python3
"""
PCA Analysis for Wild-Type Annexin A11

This script performs Principal Component Analysis (PCA) on the
conformational ensemble of wild-type Annexin A11.

PCA reduces the high-dimensional coordinate space to principal components,
enabling visualization of the conformational landscape and identification
of dominant motions.

Output:
- 2D PCA projection plot
- Variance explained by each component
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


def main():
    """Run PCA analysis for Wild-Type Annexin A11."""
    
    # Initialize components
    config = AnnexinConfig()
    loader = TrajectoryLoader(config)
    analyzer = ConformationalAnalyzer(config)
    visualizer = ConformationalVisualizer(config)
    
    # Get variant configuration
    wt = config.wt_variant
    print(f"Analyzing: {wt.label}")
    
    # Process trajectory
    print("\nLoading and processing trajectory...")
    processed = loader.process_variant(wt)
    
    print(f"  Frames: {processed.n_frames}")
    print(f"  Coordinate dimensions: {processed.coordinates_flat.shape}")
    
    # Perform PCA
    print("\nPerforming PCA...")
    pca_result = analyzer.compute_pca(processed)
    
    print(f"\nResults:")
    print(f"  PC1 variance: {pca_result.pc1_variance:.1f}%")
    print(f"  PC2 variance: {pca_result.pc2_variance:.1f}%")
    print(f"  Total variance explained: {pca_result.total_variance_explained:.1f}%")
    
    # Visualize
    print("\nGenerating visualization...")
    visualizer.plot_pca(
        pca_result,
        title="Conformational Space: PCA - Wild-Type Annexin A11",
        filename="pca_wt.png"
    )
    
    print("\nAnalysis complete!")


if __name__ == "__main__":
    main()
