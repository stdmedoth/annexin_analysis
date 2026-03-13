#!/usr/bin/env python3
"""
Export an trajectory in frames

Output:
- list of frames.pdb
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from annexin_analysis import (
    AnnexinConfig,
    TrajectoryLoader,
)


def main():
    """Convert the bioemu pdb to an more complete with more atoms."""
    
    # Initialize components
    config = AnnexinConfig()
    loader = TrajectoryLoader()
    
    # Define variants
    variant = config.p36r_variant
    output_dir = f"{config.output_dir}/exported_frames"            

    loader.export_frames(
        variant=variant,
        output_dir=output_dir
    )
    
    print(f"\Exportation complete! {output_dir}")


if __name__ == "__main__":
    main()
