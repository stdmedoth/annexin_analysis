#!/usr/bin/env python3
"""
DSSP Comparison: Wild-Type vs Mutant P36R

This script compares the conformational space of wild-type Annexin A11
and the P36R mutant using DSSP.


Output:

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
    loader = TrajectoryLoader(config)
    analyzer = ConformationalAnalyzer(config)
    visualizer = ConformationalVisualizer(config)

    # Define variants
    wt = config.wt_variant
    p36r = config.p36r_variant

    print("DSSP Comparison: WT vs P36R")
    print("=" * 40)
    print(f"Wild Type: {wt.xtc_path}")
    print(f"Mutant:    {p36r.xtc_path}")

    wt_variant = loader.process_variant(wt, extract_ca_atoms=False)
    p36r_variant = loader.process_variant(p36r, extract_ca_atoms=False)

    # Perform comparative PCA analysis
    print("\nProcessing...")
    wt_dssp_result = analyzer.compute_dssp(wt_variant)
    p36r_dssp_result = analyzer.compute_dssp(p36r_variant)

    visualizer.plot_dssp_profile_comparison(
        wt_result=wt_dssp_result,
        mut_result=p36r_dssp_result,
        filename="dssp_comparisson_wt_p36r.png",
        structure_type='Helix_%'
    )


    print("\nAnalysis complete!")


if __name__ == "__main__":
    main()
