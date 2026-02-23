#!/usr/bin/env python3
"""
Convergence Analysis for Annexin A11 Ensemble

This script assesses whether the conformational ensemble has converged
by analyzing how the RMSF profile changes with increasing sample size.

A converged ensemble shows stable RMSF values that don't change
significantly with more samples. This is crucial for validating
that the ensemble adequately samples the conformational space.

Output:
- Convergence plot (RMSF difference vs sample size)
- Convergence status assessment
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


def analyze_convergence(variant_name: str):
    """
    Analyze convergence for a specific variant.

    Args:
        variant_name: 'wt', 'p36r', ...
    """
    # Initialize components
    config = AnnexinConfig()
    loader = TrajectoryLoader(config)
    analyzer = ConformationalAnalyzer(config)
    visualizer = ConformationalVisualizer(config)

    # Get variant
    if variant_name.lower() == "wt":
        variant = config.wt_variant
    elif variant_name.lower() == "wt2":
        variant = config.wt2_variant
    elif variant_name.lower() == "p36r":
        variant = config.p36r_variant
    else:
        raise ValueError(f"Unknown variant: {variant_name}")

    print(f"Convergence Analysis: {variant.label}")
    print("=" * 40)

    # Process trajectory
    print("Loading trajectory...")
    processed = loader.process_variant(variant)

    print(f"  Total frames: {processed.n_frames}")

    # Run convergence analysis
    print("\nAnalyzing convergence...")
    convergence = analyzer.analyze_convergence(processed)

    # Print results
    print("\nResults:")
    print(f"  Samples tested: {len(convergence.sample_sizes)}")
    #print(f"  Final difference: {convergence.rmsf_list[-1]:.6f}")

    status = "CONVERGED" if convergence.is_converged else "NOT CONVERGED"
    print(f"  Status: {status}")

    # Visualize
    print("\nGenerating plot...")
    visualizer.plot_convergence(
        convergence,
        title=f"RMSF Convergence Analysis - {variant.label}",
        filename=f"convergence_{variant_name.lower()}.png"
    )

    return convergence


def main():
    """Run convergence analysis for both WT, P36R..."""

    print("=" * 50)
    print("CONVERGENCE ANALYSIS - ANNEXIN A11")
    print("=" * 50)

    # Analyze WT
    print("\n")
    wt_conv = analyze_convergence("wt")

    # Analyze WT 2 (luis felipe/ geraldo data)
    print("\n")
    wt2_conv = analyze_convergence("wt2")

    # Analyze P36R
    print("\n")
    p36r_conv = analyze_convergence("p36r")

    # Summary
    print("\n" + "=" * 50)
    print("SUMMARY")
    print("=" * 50)
    print(f"WT converged: {wt_conv.is_converged}")
    print(f"WT2 converged: {wt2_conv.is_converged}")
    print(f"P36R converged: {p36r_conv.is_converged}")

    print("\nAnalysis complete!")


if __name__ == "__main__":
    main()
