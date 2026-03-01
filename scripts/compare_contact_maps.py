#!/usr/bin/env python3
"""
Contact Map Comparison: Wild-Type vs Mutant P36R

This script compares the contact maps between wild-type Annexin A11
and the P36R mutant.

Contact maps show inter-residue contacts (residue pairs closer than
a cutoff distance). Comparing maps reveals:
- Contacts gained in the mutant (new interactions)
- Contacts lost in the mutant (broken interactions)
- Distance changes (residues moving closer/farther)

This analysis provides structural evidence of conformational changes
induced by the P36R mutation.

Output:
- Binary contact difference map (gained/lost contacts)
- Distance difference map (continuous)
- Summary statistics
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from annexin_analysis import (
    AnnexinConfig,
    ContactMapAnalyzer,
    ConformationalVisualizer
)


def main():
    """Compare contact maps between WT and P36R mutant."""
    
    # Initialize components
    config = AnnexinConfig()
    contact_analyzer = ContactMapAnalyzer(config)
    visualizer = ConformationalVisualizer(config)
    
    # Define variants
    wt = config.wt_variant
    p36r = config.p36r_variant
    
    print("Contact Map Comparison: WT vs P36R")
    print("=" * 40)
    
    # Perform full comparison
    comparison = contact_analyzer.analyze_variants(wt, p36r)
    
    # Visualize - generates two separate images
    print("\nGenerating comparison plots...")
    visualizer.plot_contact_comparison(
        comparison,
        mutation_position=36,  # P36R mutation at position 36
        title="Wild-Type vs P36R - Annexin A11",
        filename="contact_map_comparison_wt_p36r.png"
    )
    
    print("\nAnalysis complete!")


if __name__ == "__main__":
    main()
