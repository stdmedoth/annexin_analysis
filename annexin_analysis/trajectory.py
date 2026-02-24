"""
Trajectory loading and preprocessing module.

Provides classes for loading, aligning, and preprocessing molecular dynamics
trajectories from various sources (XTC, PDB, NPZ files).
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Tuple, Union
import numpy as np
import mdtraj as md

from .config import AnnexinConfig, ProteinRegions, VariantConfig, DEFAULT_CONFIG


@dataclass
class ProcessedTrajectory:
    """
    Container for processed trajectory data.

    Attributes:
        trajectory: The MDTraj trajectory object (CA atoms only).
        coordinates_flat: Flattened coordinates array (n_frames, n_atoms * 3).
        n_frames: Number of frames in the trajectory.
        n_atoms: Number of atoms (CA atoms).
        variant_config: Configuration of the protein variant.
    """

    trajectory: md.Trajectory
    coordinates_flat: np.ndarray
    n_frames: int
    n_atoms: int
    variant_config: Optional[VariantConfig] = None

    @property
    def coordinates_3d(self) -> np.ndarray:
        """Returns coordinates in 3D shape (n_frames, n_atoms, 3)."""
        return self.trajectory.xyz

    def get_frame_slice(self, start: int, end: int) -> "ProcessedTrajectory":
        """
        Create a new ProcessedTrajectory from a frame slice.

        Args:
            start: Starting frame index.
            end: Ending frame index (exclusive).

        Returns:
            New ProcessedTrajectory with sliced frames.
        """
        sliced_traj = self.trajectory[start:end]
        n_frames = sliced_traj.n_frames
        n_atoms = sliced_traj.n_atoms
        coords_flat = sliced_traj.xyz.reshape(n_frames, n_atoms * 3)

        return ProcessedTrajectory(
            trajectory=sliced_traj,
            coordinates_flat=coords_flat,
            n_frames=n_frames,
            n_atoms=n_atoms,
            variant_config=self.variant_config
        )


class TrajectoryLoader:
    """
    Handles loading and preprocessing of molecular dynamics trajectories.

    This class provides methods to:
    - Load trajectories from XTC files
    - Extract CA atoms
    - Perform structural alignment on the Annexin core
    - Convert coordinate units (Angstrom <-> nm)

    Attributes:
        config: Configuration object with paths and parameters.
        regions: Protein structural regions definition.
    """

    def __init__(self, config: Optional[AnnexinConfig] = None):
        """
        Initialize the TrajectoryLoader.

        Args:
            config: Configuration object. Uses DEFAULT_CONFIG if not provided.
        """
        self.config = config or DEFAULT_CONFIG
        self.regions = self.config.regions

    def load_trajectory(
        self,
        xtc_path: str,
        topology_path: str,
        convert_units: bool = True,
        max_frames: Optional[int] = None
    ) -> md.Trajectory:
        """
        Load a trajectory from XTC file.

        Args:
            xtc_path: Path to the XTC trajectory file.
            topology_path: Path to the topology PDB file.
            convert_units: If True, convert from Angstrom to nm.
            max_frames: Optional limit on number of frames to load.

        Returns:
            Loaded MDTraj trajectory object.
        """
        traj = md.load(xtc_path, top=topology_path)

        if max_frames is not None:
            traj = traj[:max_frames]

        if convert_units:
            traj.xyz = traj.xyz / self.config.angstrom_to_nm

        return traj

    def extract_ca_atoms(
        self,
        trajectory: md.Trajectory,
        selection: Optional[str] = None
    ) -> md.Trajectory:
        """
        Extract CA (alpha carbon) atoms from trajectory.

        Args:
            trajectory: Input trajectory.
            selection: MDTraj selection string. Uses default CA selection if None.

        Returns:
            Trajectory containing only CA atoms.
        """
        if selection is None:
            selection = self.regions.CA_SELECTION

        ca_indices = trajectory.topology.select(selection)
        return trajectory.atom_slice(ca_indices)

    def align_to_core(
        self,
        trajectory: md.Trajectory,
        reference: Optional[md.Trajectory] = None,
        reference_frame: int = 0
    ) -> md.Trajectory:
        """
        Align trajectory to the stable Annexin core region.

        Args:
            trajectory: Trajectory to align (should contain only CA atoms).
            reference: Reference trajectory for alignment. Uses self if None.
            reference_frame: Frame index to use as reference.

        Returns:
            Aligned trajectory.
        """
        core_selection = self.regions.core_selection
        core_indices = trajectory.topology.select(core_selection)

        if reference is None:
            trajectory.superpose(trajectory, reference_frame, atom_indices=core_indices)
        else:
            trajectory.superpose(reference, reference_frame, atom_indices=core_indices)

        return trajectory

    def process_variant(
        self,
        variant: VariantConfig,
        reference_trajectory: Optional[md.Trajectory] = None,
        convert_units: bool = True,
        max_frames: Optional[int] = None,
        extract_ca_atoms: bool = True,
    ) -> ProcessedTrajectory:
        """
        Full processing pipeline for a protein variant.

        Loads trajectory, extracts CA atoms, aligns to core, and prepares
        flattened coordinates for analysis.

        Args:
            variant: Configuration for the protein variant.
            reference_trajectory: Optional reference for cross-alignment.
            convert_units: If True, convert from Angstrom to nm.
            max_frames: Optional limit on number of frames.

        Returns:
            ProcessedTrajectory containing all processed data.
        """
        # Load raw trajectory
        traj = self.load_trajectory(
            variant.xtc_path,
            variant.topology_path,
            convert_units=convert_units,
            max_frames=max_frames
        )

        traj_ca = traj
        # Extract CA atoms
        if extract_ca_atoms == True:
            traj_ca = self.extract_ca_atoms(traj)

        # Align to core
        self.align_to_core(traj_ca, reference=reference_trajectory)

        # Prepare output
        n_frames = traj_ca.n_frames
        n_atoms = traj_ca.n_atoms
        coords_flat = traj_ca.xyz.reshape(n_frames, n_atoms * 3)

        return ProcessedTrajectory(
            trajectory=traj_ca,
            coordinates_flat=coords_flat,
            n_frames=n_frames,
            n_atoms=n_atoms,
            variant_config=variant
        )

    def load_and_align_for_comparison(
        self,
        wt_variant: VariantConfig,
        mutant_variant: VariantConfig,
        convert_units: bool = True
    ) -> Tuple[ProcessedTrajectory, ProcessedTrajectory]:
        """
        Load and align WT and mutant trajectories for comparison.

        The WT trajectory is self-aligned, then the mutant is cross-aligned
        to the WT reference for consistent comparison.

        Args:
            wt_variant: Configuration for Wild Type.
            mutant_variant: Configuration for mutant.
            convert_units: If True, convert from Angstrom to nm.

        Returns:
            Tuple of (wt_processed, mutant_processed) trajectories.
        """
        # Process WT first (self-aligned)
        wt_processed = self.process_variant(
            wt_variant,
            reference_trajectory=None,
            convert_units=convert_units
        )

        # Process mutant with cross-alignment to WT
        mutant_processed = self.process_variant(
            mutant_variant,
            reference_trajectory=wt_processed.trajectory,
            convert_units=convert_units
        )

        return wt_processed, mutant_processed

    def save_aligned_trajectory(
        self,
        processed: ProcessedTrajectory,
        output_dir: Union[str, Path],
        prefix: str = "core_aligned"
    ) -> Tuple[str, str]:
        """
        Save aligned trajectory to XTC and topology to PDB.

        Args:
            processed: Processed trajectory to save.
            output_dir: Directory to save files in.
            prefix: Filename prefix.

        Returns:
            Tuple of (xtc_path, pdb_path) for saved files.
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        xtc_path = output_dir / f"{prefix}_samples.xtc"
        pdb_path = output_dir / f"{prefix}_topology.pdb"

        processed.trajectory.save_xtc(str(xtc_path))
        processed.trajectory[0].save_pdb(str(pdb_path))

        return str(xtc_path), str(pdb_path)

    def load_reference_structure(
        self,
        pdb_path: str,
        convert_units: bool = True
    ) -> md.Trajectory:
        """
        Load a reference PDB structure (e.g., AlphaFold model).

        Args:
            pdb_path: Path to the PDB file.
            convert_units: If True, convert from Angstrom to nm.

        Returns:
            Loaded and processed reference trajectory.
        """
        ref = md.load(pdb_path)

        if convert_units:
            ref.xyz = ref.xyz / self.config.angstrom_to_nm

        # Extract CA atoms
        ref_ca = self.extract_ca_atoms(ref)

        # Self-align
        self.align_to_core(ref_ca)

        return ref_ca


class RegionalTrajectoryLoader(TrajectoryLoader):
    """
    Extended loader for regional (local) analysis.

    Supports loading trajectories and extracting specific regions
    for focused analysis (e.g., N-terminal dynamics).
    """

    def process_for_regional_analysis(
        self,
        variant: VariantConfig,
        align_selection: str,
        analysis_selection: str,
        reference_trajectory: Optional[md.Trajectory] = None,
        convert_units: bool = True
    ) -> ProcessedTrajectory:
        """
        Process trajectory for regional/local analysis.

        Aligns by one region (e.g., core) but extracts coordinates
        from another region (e.g., N-terminal) for analysis.

        Args:
            variant: Configuration for the protein variant.
            align_selection: MDTraj selection for alignment.
            analysis_selection: MDTraj selection for coordinate extraction.
            reference_trajectory: Optional reference for cross-alignment.
            convert_units: If True, convert from Angstrom to nm.

        Returns:
            ProcessedTrajectory with coordinates from analysis region only.
        """
        # Load raw trajectory
        traj = self.load_trajectory(
            variant.xtc_path,
            variant.topology_path,
            convert_units=convert_units
        )

        # Extract all CA atoms first
        traj_ca = self.extract_ca_atoms(traj, selection="name CA")

        # Get indices for alignment and analysis regions
        align_indices = traj_ca.topology.select(align_selection)
        analysis_indices = traj_ca.topology.select(analysis_selection)

        # Align based on alignment region
        if reference_trajectory is None:
            traj_ca.superpose(traj_ca, 0, atom_indices=align_indices)
        else:
            traj_ca.superpose(reference_trajectory, 0, atom_indices=align_indices)

        # Extract only the analysis region
        traj_region = traj_ca.atom_slice(analysis_indices)

        n_frames = traj_region.n_frames
        n_atoms = traj_region.n_atoms
        coords_flat = traj_region.xyz.reshape(n_frames, n_atoms * 3)

        return ProcessedTrajectory(
            trajectory=traj_region,
            coordinates_flat=coords_flat,
            n_frames=n_frames,
            n_atoms=n_atoms,
            variant_config=variant
        )
