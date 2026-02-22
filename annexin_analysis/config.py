"""
Configuration module for Annexin A11 analysis.

Contains protein-specific constants and configuration classes.
"""

from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional
import warnings


# Suppress BioEmu unit cell warnings globally
warnings.filterwarnings("ignore", message="Unlikely unit cell vectors detected")


@dataclass(frozen=True)
class ProteinRegions:
    """
    Defines structural regions of Annexin A11.

    The protein has two main regions:
    - N-terminal (intrinsically disordered region): residues 0-199
    - Annexin core (structured domain): residues 199-505
    """

    N_TERMINAL_START: int = 0
    N_TERMINAL_END: int = 199
    CORE_START: int = 199
    CORE_END: int = 505
    TOTAL_RESIDUES: int = 506

    # For structural alignment (using the stable core)
    ALIGNMENT_SELECTION: str = "resi 199 to 505"

    # CA selection for full protein
    CA_SELECTION: str = "name CA and resi 0 to 505"

    @property
    def core_selection(self) -> str:
        """Returns MDTraj selection string for the core region."""
        return f"resi {self.CORE_START} to {self.CORE_END}"

    @property
    def n_terminal_selection(self) -> str:
        """Returns MDTraj selection string for the N-terminal region."""
        return f"resi {self.N_TERMINAL_START} to {self.N_TERMINAL_END}"


@dataclass
class VariantConfig:
    """Configuration for a protein variant (WT or mutant)."""

    name: str
    label: str
    xtc_path: str
    topology_path: str
    color: str = "#2c3e50"

    def __post_init__(self):
        """Validate that files exist."""
        if not Path(self.xtc_path).exists():
            warnings.warn(f"XTC file not found: {self.xtc_path}")
        if not Path(self.topology_path).exists():
            warnings.warn(f"Topology file not found: {self.topology_path}")


@dataclass
class AnnexinConfig:
    """
    Main configuration class for Annexin A11 analysis.

    Centralizes all paths, settings, and analysis parameters.

    Attributes:
        base_path: Base directory containing trajectory data.
        output_dir: Directory for saving results.
        regions: Protein structural regions definition.
        variants: Dictionary of protein variants to analyze.
    """

    today = datetime.today().strftime('%Y-%m-%d')
    base_path: str = "1000_samples/MutationConformationsBioEmu"
    output_dir: str = f"results/{today}"
    regions: ProteinRegions = field(default_factory=ProteinRegions)

    # Analysis parameters
    angstrom_to_nm: float = 10.0
    contact_cutoff_nm: float = 0.8  # 8 Angstroms in nm
    pca_components: int = 2

    # Convergence analysis
    convergence_step: int = 10
    convergence_max_samples: int = 1000

    # Visualization defaults
    figure_dpi: int = 300
    default_figsize: tuple = (10, 8)

    # Color palette for variants
    colors: List[str] = field(default_factory=lambda: [
        "#2c3e50",  # Dark blue-gray (WT)
        "#e74c3c",  # Red (Mutant 1)
        "#27ae60",  # Green (Mutant 2)
        "#f39c12",  # Orange (Mutant 3)
        "#9b59b6",  # Purple (Mutant 4)
    ])

    def __post_init__(self):
        """Create output directory if it doesn't exist."""
        Path(self.output_dir).mkdir(parents=True, exist_ok=True)

    def get_variant_path(self, variant_folder: str, filename: str) -> str:
        """
        Construct full path for a variant file.

        Args:
            variant_folder: Folder name within base_path (e.g., 'out_native').
            filename: File name (e.g., 'samples.xtc').

        Returns:
            Full path to the file.
        """
        return f"{self.base_path}/{variant_folder}/{filename}"

    def create_variant(
        self,
        name: str,
        label: str,
        folder: str,
        color_index: int = 0
    ) -> VariantConfig:
        """
        Create a VariantConfig for a protein variant.

        Args:
            name: Internal name for the variant.
            label: Display label for plots.
            folder: Folder name within base_path.
            color_index: Index into the color palette.

        Returns:
            Configured VariantConfig instance.
        """
        return VariantConfig(
            name=name,
            label=label,
            xtc_path=self.get_variant_path(folder, "samples.xtc"),
            topology_path=self.get_variant_path(folder, "topology.pdb"),
            color=self.colors[color_index % len(self.colors)]
        )

    @property
    def wt_variant(self) -> VariantConfig:
        """Returns configuration for the Wild Type variant."""
        return self.create_variant(
            name="wt",
            label="Wild Type",
            folder="out_native",
            color_index=0
        )

    @property
    def p36r_variant(self) -> VariantConfig:
        """Returns configuration for the P36R mutant variant."""
        return self.create_variant(
            name="p36r",
            label="Mutant P36R",
            folder="out_mutant_P36R",
            color_index=1
        )


# Default configuration instance
DEFAULT_CONFIG = AnnexinConfig()
