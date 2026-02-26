"""
Conformational analysis module.

Provides classes for analyzing protein conformational dynamics including:
- RMSF (Root Mean Square Fluctuation) analysis
- PCA (Principal Component Analysis)
- Convergence analysis
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple
import numpy as np
import pandas as pd
import mdtraj as md
from sklearn.decomposition import PCA

from .config import AnnexinConfig, VariantConfig, DEFAULT_CONFIG
from .trajectory import ProcessedTrajectory, TrajectoryLoader


@dataclass
class DSSPResult:
    """
    Container for DSSP analysis

    Attributes:
        raw_matrix: The matrix with DSSP values,
        percentages: A dict with the percentages of secondaries structures,
        variant_name: Variant Name
    """
    raw_matrix: np.ndarray
    percentages: dict
    variant_name: str

@dataclass
class RMSFResult:
    """
    Container for RMSF analysis results.

    Attributes:
        fluctuations: RMSF values per residue (in Angstroms).
        residue_numbers: Array of residue indices.
        mean_rmsf: Mean RMSF across all residues.
        variant_name: Name of the protein variant.
    """

    fluctuations: np.ndarray
    residue_numbers: np.ndarray
    mean_rmsf: float
    variant_name: str

    def __post_init__(self):
        """Ensure residue numbers array exists."""
        if self.residue_numbers is None:
            self.residue_numbers = np.arange(len(self.fluctuations))


@dataclass
class PCAResult:
    """
    Container for PCA analysis results.

    Attributes:
        reduced_coords: Projected coordinates in PC space.
        explained_variance_ratio: Variance explained by each PC.
        pca_model: Fitted PCA model (for projecting other data).
        n_frames: Number of frames analyzed.
        variant_name: Name of the protein variant.
    """

    reduced_coords: np.ndarray
    explained_variance_ratio: np.ndarray
    pca_model: PCA
    n_frames: int
    variant_name: str

    @property
    def total_variance_explained(self) -> float:
        """Returns total variance explained by all components."""
        return self.explained_variance_ratio.sum() * 100

    @property
    def pc1_variance(self) -> float:
        """Returns variance explained by PC1 (percentage)."""
        return self.explained_variance_ratio[0] * 100

    @property
    def pc2_variance(self) -> float:
        """Returns variance explained by PC2 (percentage)."""
        return self.explained_variance_ratio[1] * 100


@dataclass
class ConvergenceResult:
    """
    Container for convergence analysis results.

    Attributes:
        sample_sizes: Array of sample sizes tested.
        rmsf_list: RMSF List at each sample size.
        final_rmsf: Reference RMSF from full ensemble.
        is_converged: Whether the analysis has converged.
    """

    sample_sizes: np.ndarray
    rmsf_list: np.ndarray
    final_rmsf: np.ndarray
    is_converged: bool = False
    convergence_threshold: float = 0.01



class ConformationalAnalyzer:
    """
    Main class for conformational space analysis.

    Provides methods for RMSF calculation, PCA analysis, and
    convergence assessment of molecular dynamics ensembles.

    Attributes:
        config: Configuration object with analysis parameters.
        loader: TrajectoryLoader for data loading.
    """

    def __init__(self, config: Optional[AnnexinConfig] = None):
        """
        Initialize the ConformationalAnalyzer.

        Args:
            config: Configuration object. Uses DEFAULT_CONFIG if not provided.
        """
        self.config = config or DEFAULT_CONFIG
        self.loader = TrajectoryLoader(self.config)

    def compute_rmsf(
        self,
        processed: ProcessedTrajectory,
        reference_frame: int = 0,
        to_angstroms: bool = True
    ) -> RMSFResult:
        """
        Compute RMSF for a processed trajectory.

        Args:
            processed: ProcessedTrajectory object.
            reference_frame: Frame to use as reference.
            to_angstroms: If True, convert result to Angstroms.

        Returns:
            RMSFResult containing fluctuation data.
        """
        traj = processed.trajectory

        fluctuations = md.rmsf(traj, traj, reference_frame)

        if to_angstroms:
            fluctuations = fluctuations * self.config.angstrom_to_nm

        residue_numbers = np.arange(len(fluctuations))
        mean_rmsf = float(np.mean(fluctuations))

        variant_name = ""
        if processed.variant_config:
            variant_name = processed.variant_config.name

        return RMSFResult(
            fluctuations=fluctuations,
            residue_numbers=residue_numbers,
            mean_rmsf=mean_rmsf,
            variant_name=variant_name
        )

    def compute_rmsf_vs_reference(
        self,
        processed: ProcessedTrajectory,
        reference: md.Trajectory,
        to_angstroms: bool = True
    ) -> RMSFResult:
        """
        Compute RMSF against an external reference structure.

        Args:
            processed: ProcessedTrajectory object.
            reference: Reference trajectory (e.g., AlphaFold model).
            to_angstroms: If True, convert result to Angstroms.

        Returns:
            RMSFResult containing fluctuation data.
        """
        traj = processed.trajectory

        fluctuations = md.rmsf(traj, reference)

        if to_angstroms:
            fluctuations = fluctuations * self.config.angstrom_to_nm

        residue_numbers = np.arange(len(fluctuations))
        mean_rmsf = float(np.mean(fluctuations))

        variant_name = ""
        if processed.variant_config:
            variant_name = processed.variant_config.name

        return RMSFResult(
            fluctuations=fluctuations,
            residue_numbers=residue_numbers,
            mean_rmsf=mean_rmsf,
            variant_name=variant_name
        )

    def compute_pca(
        self,
        processed: ProcessedTrajectory,
        n_components: Optional[int] = None
    ) -> PCAResult:
        """
        Perform PCA on trajectory coordinates.

        Args:
            processed: ProcessedTrajectory object.
            n_components: Number of PCs to compute. Uses config default if None.

        Returns:
            PCAResult containing reduced coordinates and PCA model.
        """
        if n_components is None:
            n_components = self.config.pca_components

        pca = PCA(n_components=n_components)
        reduced_coords = pca.fit_transform(processed.coordinates_flat)

        variant_name = ""
        if processed.variant_config:
            variant_name = processed.variant_config.name

        return PCAResult(
            reduced_coords=reduced_coords,
            explained_variance_ratio=pca.explained_variance_ratio_,
            pca_model=pca,
            n_frames=processed.n_frames,
            variant_name=variant_name
        )

    def project_to_pca_space(
        self,
        processed: ProcessedTrajectory,
        pca_result: PCAResult
    ) -> PCAResult:
        """
        Project trajectory onto existing PCA space.

        Useful for comparing mutants projected onto WT PCA space.

        Args:
            processed: ProcessedTrajectory to project.
            pca_result: Existing PCAResult with fitted model.

        Returns:
            PCAResult with projected coordinates.
        """
        reduced_coords = pca_result.pca_model.transform(processed.coordinates_flat)

        variant_name = ""
        if processed.variant_config:
            variant_name = processed.variant_config.name

        return PCAResult(
            reduced_coords=reduced_coords,
            explained_variance_ratio=pca_result.explained_variance_ratio,
            pca_model=pca_result.pca_model,
            n_frames=processed.n_frames,
            variant_name=variant_name
        )

    def analyze_convergence(
        self,
        processed: ProcessedTrajectory,
        step: Optional[int] = None,
        max_samples: Optional[int] = None
    ) -> ConvergenceResult:
        """
        Analyze RMSF convergence as function of sample size.

        Tests whether the ensemble has converged by computing RMSF
        profiles for increasing subsets of samples and comparing
        to the final (full ensemble) profile.

        Args:
            processed: ProcessedTrajectory object.
            step: Step size for sample sizes. Uses config default if None.
            max_samples: Maximum samples to test. Uses all if None.

        Returns:
            ConvergenceResult with convergence data.
        """
        if step is None:
            step = self.config.convergence_step

        if max_samples is None:
            max_samples = min(processed.n_frames, self.config.convergence_max_samples)

        # Compute reference RMSF from full ensemble
        final_rmsf = md.rmsf(processed.trajectory, processed.trajectory, 0)

        # Test different sample sizes
        sample_sizes = list(range(step, max_samples + 1, step))
        #differences = []
        rmsf_list = []

        for n in sample_sizes:
            sub_processed = processed.get_frame_slice(0, n)

            # Re-align the subset
            core_indices = sub_processed.trajectory.topology.select(
                self.config.regions.core_selection
            )
            sub_processed.trajectory.superpose(
                sub_processed.trajectory, 0, atom_indices=core_indices
            )

            current_rmsf = md.rmsf(
                sub_processed.trajectory, sub_processed.trajectory, 0
            )

            # Calculate L2 norm of difference
            #diff_norm = np.linalg.norm(current_rmsf - final_rmsf)
            #differences.append(diff_norm)

            rmsf_list.append(current_rmsf)

        differences = []
        for rmsf_n in rmsf_list:
            rmse = np.sqrt(np.mean(rmsf_n-final_rmsf)**2)
            differences.append(rmse)

        differences = np.array(differences)

        is_converged = False
        error_threshold = 0.02
        slope_threshold = 1e-4
        if len(differences) >= 5:
            last_points = differences[-5:]
            x_axis = np.arange(5)
            slope, _ = np.polyfit(x_axis, last_points, 1)
            #print(f'slope {abs(slope)}, last_points {last_points[-1]}')
            is_converged = last_points[-1] < error_threshold and (abs(slope) < slope_threshold or last_points[-1] < 1e-5)


        return ConvergenceResult(
            sample_sizes=np.array(sample_sizes),
            rmsf_list=rmsf_list,
            final_rmsf=final_rmsf,
            is_converged=is_converged
        )

    def count_unique_structures(
        self,
        processed: ProcessedTrajectory,
        decimal_precision: int = 4
    ) -> int:
        """
        Count number of unique structures in the ensemble.

        Args:
            processed: ProcessedTrajectory object.
            decimal_precision: Decimal places for rounding coordinates.

        Returns:
            Number of unique structures.
        """
        unique_frames = np.unique(
            processed.coordinates_flat.round(decimals=decimal_precision),
            axis=0
        )
        return len(unique_frames)

    def compute_rg(
        self,
        processed: ProcessedTrajectory,
        reference_frame: int = 0,
        to_angstroms: bool = True
    ) -> float:
        """
            Compute Mean Radius of Gyration
        """

        traj = processed.trajectory

        rg_per_frame = md.compute_rg(traj)

        rg_angstrom = rg_per_frame * 10

        mean_rg = np.mean(rg_angstrom)

        return mean_rg

    def compute_frame_distance(
        self,
        processed: ProcessedTrajectory,
        frame_i: int,
        frame_j: int
    ) -> float:
        """
        Compute distance between two frames.

        Args:
            processed: ProcessedTrajectory object.
            frame_i: First frame index.
            frame_j: Second frame index.

        Returns:
            Euclidean distance between frames in coordinate space.
        """
        coords = processed.coordinates_flat
        return float(np.linalg.norm(coords[frame_i] - coords[frame_j]))



    def compute_dssp(
        self,
        processed: ProcessedTrajectory

    ) -> DSSPResult:
        traj = processed.trajectory
        variant_name = processed.variant_config.name if processed.variant_config else ""

        # 1. Computar DSSP (8 estados)
        dssp_data = md.compute_dssp(traj, simplified=False)

        # Pegar os IDs reais dos resíduos da topologia
        residue_ids = np.array([res.resSeq for res in traj.topology.residues])

        # 2. Cálculo Vetorizado de Porcentagens (Muito mais rápido)
        # Agrupamos por classes físicas:
        is_helix = np.isin(dssp_data, ['H', 'G', 'I'])
        is_strand = np.isin(dssp_data, ['E', 'B'])
        is_coil = ~is_helix & ~is_strand & (dssp_data != 'NA')

        # Médias temporais por resíduo (axis 0 é o tempo)
        stats = {
            res_id: {
                'Helix_%': np.mean(is_helix[:, i]) * 100,
                'Strand_%': np.mean(is_strand[:, i]) * 100,
                'Coil_%': np.mean(is_coil[:, i]) * 100
            }
            for i, res_id in enumerate(residue_ids)
        }

        return DSSPResult(
            raw_matrix=dssp_data,
            percentages=stats,
            variant_name=variant_name
        )

    def compute_exposed_residues(
        self,
        processed: ProcessedTrajectory
    ) -> np.array:
        """
        Calculate the exposed residues for an trajectory

        Args:
            processed: ProcessedTrajectory object.

        Returns:
            An np array
        """

        threshold = 0.2
        traj = processed.trajectory

        core_indices = traj.topology.select(f'resi {self.config.regions.CORE_START} to {self.config.regions.CORE_END}')
        core_traj = traj.atom_slice(core_indices)


        sasa_per_res = md.shrake_rupley(core_traj, mode='residue')
        avg_sasa = np.mean(sasa_per_res, axis=0)

        exposed_local_idx = np.where(avg_sasa > threshold)[0]
        exposed_res_seqs = [core_traj.topology.residue(i).resSeq for i in exposed_local_idx]

        exposed_ca_indices = traj.topology.select(f"name CA and resSeq " + " ".join(map(str, exposed_res_seqs)))

        return exposed_ca_indices

    

    def compute_idr_to_core_mean_distance(
        self,
        processed: ProcessedTrajectory,
        exposed_ca_indices: np.array
    ) -> np.array:
        """
        Computation of the mean distances between the exposed from CORE carbons and the IDR residue

        Args:
            processed: ProcessedTrajectory object.
            exposed_ca_indices: and np.array with the ca indices of the exposed core

        Returns:
            An np.array with the distances
        """

        traj = processed.trajectory

        idr_ca_indices = md.select('name CA and resi {self.config.regions.N_TERMINAL_START} to {self.config.regions.N_TERMINAL_END}')
        

        # create the pairs [(idr_ca_indices, exposed_ca_indices)]
        pairs = np.array([[i,j] for i in idr_ca_indices for j in exposed_ca_indices])
        distances = md.compute_distances(traj, pairs)
        mean_dist_pair = np.mean(distances, axis=0)

        reshape_dist = mean_dist_pair.reshape((len(idr_ca_indices), len(exposed_ca_indices)))
        
        # axis = 0 : vertical
        # axis = 1 : horizontal
        # axis = 1 because we want to colapse the core and know about the IDR
        avg_dist_to_surface = np.mean(reshape_dist, axis=1)

        return avg_dist_to_surface


        




class ComparativeAnalyzer:
    """
    Class for comparative analysis between protein variants.

    Provides methods to compare conformational spaces of
    wild-type and mutant proteins.
    """

    def __init__(self, config: Optional[AnnexinConfig] = None):
        """
        Initialize the ComparativeAnalyzer.

        Args:
            config: Configuration object. Uses DEFAULT_CONFIG if not provided.
        """
        self.config = config or DEFAULT_CONFIG
        self.loader = TrajectoryLoader(self.config)
        self.analyzer = ConformationalAnalyzer(self.config)

    def compare_rmsf(
        self,
        variants: List[VariantConfig]
    ) -> Dict[str, RMSFResult]:
        """
        Compute and compare RMSF profiles for multiple variants.

        Args:
            variants: List of variant configurations to compare.

        Returns:
            Dictionary mapping variant names to RMSFResult objects.
        """
        results = {}

        for variant in variants:
            try:
                processed = self.loader.process_variant(variant)
                rmsf_result = self.analyzer.compute_rmsf(processed)
                results[variant.name] = rmsf_result
                print(f"Processed {variant.label}: Mean RMSF = {rmsf_result.mean_rmsf:.4f} Å")
            except Exception as e:
                print(f"Error processing {variant.label}: {e}")

        return results

    def compare_pca(
        self,
        wt_variant: VariantConfig,
        mutant_variants: List[VariantConfig]
    ) -> Tuple[PCAResult, Dict[str, PCAResult]]:
        """
        Compare PCA projections of mutants onto WT space.

        The WT is used to train the PCA, then mutants are projected
        onto the same space for fair comparison.

        Args:
            wt_variant: Wild-type variant configuration.
            mutant_variants: List of mutant variant configurations.

        Returns:
            Tuple of (wt_pca_result, dict of mutant_name -> PCAResult).
        """
        # Process WT and train PCA
        wt_processed = self.loader.process_variant(wt_variant)
        wt_pca = self.analyzer.compute_pca(wt_processed)

        print(f"WT PCA - Total variance explained: {wt_pca.total_variance_explained:.2f}%")

        # Project mutants onto WT space
        mutant_results = {}

        for variant in mutant_variants:
            try:
                mutant_processed = self.loader.process_variant(
                    variant,
                    reference_trajectory=wt_processed.trajectory
                )
                mutant_pca = self.analyzer.project_to_pca_space(mutant_processed, wt_pca)
                mutant_results[variant.name] = mutant_pca
                print(f"Projected {variant.label} onto WT PCA space")
            except Exception as e:
                print(f"Error processing {variant.label}: {e}")

        return wt_pca, mutant_results
