"""
Utility classes and functions.

Provides helper classes for:
- Amino acid conversion (1-letter to 3-letter codes)
- NPZ to PDB file conversion
- File handling utilities
"""

import glob
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional
import numpy as np
import mdtraj as md


class AminoAcidConverter:
    """
    Converts between amino acid representations.
    
    Supports conversion between 1-letter and 3-letter codes
    for standard amino acids.
    """
    
    ONE_TO_THREE: Dict[str, str] = {
        'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE',
        'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'K': 'LYS', 'L': 'LEU',
        'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG',
        'S': 'SER', 'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'
    }
    
    THREE_TO_ONE: Dict[str, str] = {v: k for k, v in ONE_TO_THREE.items()}
    
    @classmethod
    def one_to_three(cls, letter: str) -> str:
        """
        Convert 1-letter amino acid code to 3-letter code.
        
        Args:
            letter: Single letter amino acid code (e.g., 'A').
            
        Returns:
            Three letter code (e.g., 'ALA').
            
        Raises:
            KeyError: If letter is not a valid amino acid code.
        """
        return cls.ONE_TO_THREE[letter.upper()]
    
    @classmethod
    def three_to_one(cls, code: str) -> str:
        """
        Convert 3-letter amino acid code to 1-letter code.
        
        Args:
            code: Three letter amino acid code (e.g., 'ALA').
            
        Returns:
            Single letter code (e.g., 'A').
            
        Raises:
            KeyError: If code is not a valid amino acid code.
        """
        return cls.THREE_TO_ONE[code.upper()]
    
    @classmethod
    def sequence_to_three(cls, sequence: str) -> List[str]:
        """
        Convert a sequence of 1-letter codes to 3-letter codes.
        
        Args:
            sequence: String of 1-letter amino acid codes.
            
        Returns:
            List of 3-letter codes.
        """
        return [cls.one_to_three(aa) for aa in sequence]


class NPZConverter:
    """
    Converts BioEmu NPZ output files to PDB format.
    
    BioEmu saves conformational samples in NPZ format with
    'pos' (positions) and 'sequence' arrays. This class converts
    them to standard PDB format for visualization and analysis.
    """
    
    def __init__(self, output_dir: Optional[str] = None):
        """
        Initialize the NPZConverter.
        
        Args:
            output_dir: Directory for output PDB files.
                       Uses source directory if None.
        """
        self.output_dir = output_dir
        self.aa_converter = AminoAcidConverter()
    
    def convert_npz_to_pdb(
        self,
        npz_path: str,
        output_path: Optional[str] = None,
        frame_index: int = 0
    ) -> str:
        """
        Convert a single NPZ file to PDB format.
        
        Args:
            npz_path: Path to input NPZ file.
            output_path: Path for output PDB file.
                        Auto-generated if None.
            frame_index: Which frame to extract (default: first).
            
        Returns:
            Path to the created PDB file.
        """
        npz_data = np.load(npz_path)
        
        positions = npz_data['pos'][frame_index]
        sequence = npz_data['sequence'].item()
        
        if output_path is None:
            if self.output_dir:
                output_path = Path(self.output_dir) / f"{Path(npz_path).stem}.pdb"
            else:
                output_path = Path(npz_path).with_suffix('.pdb')
        
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        self._write_pdb(positions, sequence, str(output_path))
        
        return str(output_path)
    
    def _write_pdb(
        self,
        positions: np.ndarray,
        sequence: str,
        output_path: str
    ) -> None:
        """
        Write positions and sequence to PDB file.
        
        Creates a CA-only PDB file (alpha carbons only).
        
        Args:
            positions: Array of (N, 3) coordinates.
            sequence: Amino acid sequence string.
            output_path: Output file path.
        """
        with open(output_path, 'w') as f:
            for i, (position, aa) in enumerate(zip(positions, sequence)):
                x, y, z = position
                residue_name = self.aa_converter.one_to_three(aa)
                
                # Standard PDB ATOM record format
                line = (
                    f"ATOM  {i+1:>5}  CA  {residue_name:>3} A{i+1:>4}    "
                    f"{x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00  0.00           C\n"
                )
                f.write(line)
            
            f.write("END\n")
    
    def convert_directory(
        self,
        input_dir: str,
        pattern: str = "*.npz",
        output_dir: Optional[str] = None
    ) -> List[str]:
        """
        Convert all NPZ files in a directory.
        
        Args:
            input_dir: Directory containing NPZ files.
            pattern: Glob pattern for NPZ files.
            output_dir: Directory for output PDB files.
            
        Returns:
            List of paths to created PDB files.
        """
        if output_dir:
            self.output_dir = output_dir
        
        npz_files = glob.glob(str(Path(input_dir) / pattern))
        output_files = []
        
        for npz_file in npz_files:
            try:
                pdb_file = self.convert_npz_to_pdb(npz_file)
                output_files.append(pdb_file)
                print(f"Converted: {npz_file} -> {pdb_file}")
            except Exception as e:
                print(f"Error converting {npz_file}: {e}")
        
        return output_files


class TrajectoryFilter:
    """
    Filters trajectory frames based on structural criteria.
    
    Useful for removing outliers or selecting conformationally
    stable structures from an ensemble.
    """
    
    @staticmethod
    def filter_by_radius_of_gyration(
        trajectory: md.Trajectory,
        n_std: float = 1.0
    ) -> md.Trajectory:
        """
        Filter frames by radius of gyration.
        
        Keeps frames within n standard deviations of the mean Rg.
        Useful for removing overly expanded/unfolded structures.
        
        Args:
            trajectory: Input trajectory.
            n_std: Number of standard deviations for cutoff.
            
        Returns:
            Filtered trajectory.
        """
        rg = md.compute_rg(trajectory)
        
        mean_rg = np.mean(rg)
        std_rg = np.std(rg)
        
        threshold = mean_rg + n_std * std_rg
        good_indices = np.where(rg < threshold)[0]
        
        return trajectory[good_indices]
    
    @staticmethod
    def filter_by_rmsd(
        trajectory: md.Trajectory,
        reference_frame: int = 0,
        max_rmsd: float = 1.0
    ) -> md.Trajectory:
        """
        Filter frames by RMSD from reference.
        
        Args:
            trajectory: Input trajectory.
            reference_frame: Frame index to use as reference.
            max_rmsd: Maximum allowed RMSD (in nm).
            
        Returns:
            Filtered trajectory.
        """
        rmsd = md.rmsd(trajectory, trajectory, reference_frame)
        good_indices = np.where(rmsd < max_rmsd)[0]
        
        return trajectory[good_indices]
    
    @staticmethod
    def save_filtered_trajectory(
        filtered_traj: md.Trajectory,
        output_xtc: str,
        output_pdb: Optional[str] = None
    ) -> None:
        """
        Save filtered trajectory to files.
        
        Args:
            filtered_traj: Filtered trajectory object.
            output_xtc: Path for output XTC file.
            output_pdb: Path for topology PDB (optional).
        """
        filtered_traj.save_xtc(output_xtc)
        
        if output_pdb:
            filtered_traj[0].save_pdb(output_pdb)


class TrajectoryInfo:
    """
    Utility class for displaying trajectory information.
    """
    
    @staticmethod
    def print_info(trajectory: md.Trajectory, label: str = "") -> None:
        """
        Print summary information about a trajectory.
        
        Args:
            trajectory: MDTraj trajectory object.
            label: Optional label for the trajectory.
        """
        if label:
            print(f"\n=== {label} ===")
        
        print(f"Trajectory: {trajectory}")
        print(f"  Frames: {trajectory.n_frames}")
        print(f"  Atoms: {trajectory.n_atoms}")
        print(f"  Coordinates shape: {trajectory.xyz.shape}")
    
    @staticmethod
    def get_shape(xtc_path: str, topology_path: str) -> tuple:
        """
        Get trajectory shape without loading full data.
        
        Args:
            xtc_path: Path to XTC file.
            topology_path: Path to topology PDB.
            
        Returns:
            Tuple of (n_frames, n_atoms, 3).
        """
        traj = md.load(xtc_path, top=topology_path)
        return traj.xyz.shape
