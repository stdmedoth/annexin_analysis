"""
Contact map analysis module.

Provides classes for computing and comparing contact maps between
protein conformational ensembles.
"""

from dataclasses import dataclass
from typing import Optional, Tuple
import numpy as np
import mdtraj as md
from scipy.spatial.distance import squareform

from .config import AnnexinConfig, VariantConfig, DEFAULT_CONFIG
from .trajectory import TrajectoryLoader


@dataclass
class ContactMapResult:
    """
    Container for contact map analysis results.
    
    Attributes:
        distance_matrix: Mean distance matrix (N x N residues).
        contact_map: Binary contact map (True where distance < cutoff).
        cutoff: Distance cutoff used (in nm).
        n_contacts: Total number of contacts.
        variant_name: Name of the protein variant.
    """
    
    distance_matrix: np.ndarray
    contact_map: np.ndarray
    cutoff: float
    n_contacts: int
    variant_name: str
    
    @property
    def n_residues(self) -> int:
        """Returns the number of residues in the contact map."""
        return self.distance_matrix.shape[0]


@dataclass
class ContactComparisonResult:
    """
    Container for contact map comparison results.
    
    Attributes:
        wt_result: ContactMapResult for wild-type.
        mutant_result: ContactMapResult for mutant.
        contact_difference: Binary difference map (+1 gained, -1 lost, 0 unchanged).
        distance_difference: Continuous distance difference matrix.
        contacts_gained: Number of contacts gained in mutant.
        contacts_lost: Number of contacts lost in mutant.
        net_change: Net change in contact count.
    """
    
    wt_result: ContactMapResult
    mutant_result: ContactMapResult
    contact_difference: np.ndarray
    distance_difference: np.ndarray
    contacts_gained: int
    contacts_lost: int
    net_change: int


class ContactMapAnalyzer:
    """
    Class for contact map computation and comparison.
    
    Computes inter-residue contact maps from MD trajectories
    and provides methods for comparing maps between variants.
    
    Attributes:
        config: Configuration object with analysis parameters.
        loader: TrajectoryLoader for data loading.
    """
    
    def __init__(self, config: Optional[AnnexinConfig] = None):
        """
        Initialize the ContactMapAnalyzer.
        
        Args:
            config: Configuration object. Uses DEFAULT_CONFIG if not provided.
        """
        self.config = config or DEFAULT_CONFIG
        self.loader = TrajectoryLoader(self.config)
    
    def compute_contact_map(
        self,
        xtc_path: str,
        topology_path: str,
        cutoff: Optional[float] = None,
        scheme: str = "ca"
    ) -> ContactMapResult:
        """
        Compute mean contact map from trajectory.
        
        Uses MDTraj's compute_contacts for efficient calculation.
        
        Args:
            xtc_path: Path to trajectory file.
            topology_path: Path to topology PDB file.
            cutoff: Distance cutoff in nm. Uses config default if None.
            scheme: Contact scheme ('ca' for alpha carbons, 'closest', etc.).
            
        Returns:
            ContactMapResult with distance and contact matrices.
        """
        if cutoff is None:
            cutoff = self.config.contact_cutoff_nm
        
        # Load trajectory (don't convert units - MDTraj uses nm)
        traj = md.load(xtc_path, top=topology_path)
        
        # Compute all pairwise distances
        distances, _ = md.compute_contacts(traj, contacts="all", scheme=scheme)
        
        # Average across frames
        mean_distances = distances.mean(axis=0)
        
        # Convert to square matrix
        distance_matrix = squareform(mean_distances)
        
        # Create binary contact map
        contact_map = distance_matrix < cutoff
        n_contacts = int(np.sum(contact_map))
        
        return ContactMapResult(
            distance_matrix=distance_matrix,
            contact_map=contact_map,
            cutoff=cutoff,
            n_contacts=n_contacts,
            variant_name=""
        )
    
    def compute_contact_map_from_variant(
        self,
        variant: VariantConfig,
        cutoff: Optional[float] = None,
        scheme: str = "ca"
    ) -> ContactMapResult:
        """
        Compute contact map for a protein variant.
        
        Args:
            variant: Variant configuration.
            cutoff: Distance cutoff in nm.
            scheme: Contact scheme.
            
        Returns:
            ContactMapResult with variant name set.
        """
        result = self.compute_contact_map(
            variant.xtc_path,
            variant.topology_path,
            cutoff,
            scheme
        )
        result.variant_name = variant.name
        return result
    
    def compute_contact_map_manual(
        self,
        xtc_path: str,
        topology_path: str,
        cutoff: Optional[float] = None
    ) -> ContactMapResult:
        """
        Compute contact map using manual vectorized calculation.
        
        This method computes the contact frequency (probability of contact)
        rather than mean distance. Useful for validation or when MDTraj's
        compute_contacts is not suitable.
        
        Args:
            xtc_path: Path to trajectory file.
            topology_path: Path to topology PDB file.
            cutoff: Distance cutoff in nm.
            
        Returns:
            ContactMapResult where distance_matrix contains contact frequencies.
        """
        if cutoff is None:
            cutoff = self.config.contact_cutoff_nm
        
        # Load trajectory
        traj = md.load(xtc_path, top=topology_path)
        coords = traj.xyz  # Shape: (n_frames, n_atoms, 3)
        
        n_frames, n_residues, _ = coords.shape
        
        # Initialize contact frequency matrix
        contact_frequency = np.zeros((n_residues, n_residues))
        
        for frame_idx in range(n_frames):
            frame_coords = coords[frame_idx]
            
            # Vectorized distance calculation using broadcasting
            # diff[i, j, :] = coord[i] - coord[j]
            diff = frame_coords[:, np.newaxis, :] - frame_coords[np.newaxis, :, :]
            
            # Euclidean distances
            distances = np.linalg.norm(diff, axis=2)
            
            # Accumulate contacts
            contact_frequency += (distances < cutoff).astype(float)
        
        # Average by number of frames (gives contact probability)
        contact_probability = contact_frequency / n_frames
        
        # Binary contact map (contact if probability > 0.5)
        contact_map = contact_probability > 0.5
        n_contacts = int(np.sum(contact_map))
        
        return ContactMapResult(
            distance_matrix=contact_probability,
            contact_map=contact_map,
            cutoff=cutoff,
            n_contacts=n_contacts,
            variant_name=""
        )
    
    def compare_contact_maps(
        self,
        wt_result: ContactMapResult,
        mutant_result: ContactMapResult
    ) -> ContactComparisonResult:
        """
        Compare contact maps between WT and mutant.
        
        Args:
            wt_result: ContactMapResult for wild-type.
            mutant_result: ContactMapResult for mutant.
            
        Returns:
            ContactComparisonResult with difference analysis.
        """
        # Binary difference: +1 = gained, -1 = lost, 0 = unchanged
        contact_diff = (
            mutant_result.contact_map.astype(int) - 
            wt_result.contact_map.astype(int)
        )
        
        # Continuous distance difference
        distance_diff = mutant_result.distance_matrix - wt_result.distance_matrix
        
        # Count changes
        contacts_gained = int(np.sum(contact_diff > 0))
        contacts_lost = int(np.sum(contact_diff < 0))
        net_change = contacts_gained - contacts_lost
        
        return ContactComparisonResult(
            wt_result=wt_result,
            mutant_result=mutant_result,
            contact_difference=contact_diff,
            distance_difference=distance_diff,
            contacts_gained=contacts_gained,
            contacts_lost=contacts_lost,
            net_change=net_change
        )
    
    def analyze_variants(
        self,
        wt_variant: VariantConfig,
        mutant_variant: VariantConfig,
        cutoff: Optional[float] = None
    ) -> ContactComparisonResult:
        """
        Full comparison pipeline for two variants.
        
        Args:
            wt_variant: Wild-type variant configuration.
            mutant_variant: Mutant variant configuration.
            cutoff: Distance cutoff in nm.
            
        Returns:
            ContactComparisonResult with full analysis.
        """
        print(f"Processing {wt_variant.label}...")
        wt_result = self.compute_contact_map_from_variant(wt_variant, cutoff)
        print(f"  Total contacts: {wt_result.n_contacts}")
        
        print(f"Processing {mutant_variant.label}...")
        mutant_result = self.compute_contact_map_from_variant(mutant_variant, cutoff)
        print(f"  Total contacts: {mutant_result.n_contacts}")
        
        comparison = self.compare_contact_maps(wt_result, mutant_result)
        
        print("\n--- Contact Map Comparison Summary ---")
        print(f"Contacts in WT: {wt_result.n_contacts}")
        print(f"Contacts in Mutant: {mutant_result.n_contacts}")
        print(f"Contacts gained: {comparison.contacts_gained}")
        print(f"Contacts lost: {comparison.contacts_lost}")
        print(f"Net change: {comparison.net_change}")
        
        return comparison
