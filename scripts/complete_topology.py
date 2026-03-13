#!/usr/bin/env python3
"""
Convert an simple topology to an more complete


Output:
- topology.pdb more complete
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from annexin_analysis import (
    AnnexinConfig,
    TrajectoryLoader,
    ConformationalVisualizer
)


def main():
    """Convert the bioemu pdb to an more complete with more atoms."""
    
    # Initialize components
    config = AnnexinConfig()
    loader = TrajectoryLoader()
    visualizer = ConformationalVisualizer(config)
    
    # Define variants
    wt = config.wt_variant

    loader.fix_topology(wt)
    
    print("\Refactorization complete!")


if __name__ == "__main__":
    main()
