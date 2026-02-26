"""
Visualization module for conformational analysis.

Provides classes for creating publication-ready plots of:
- RMSF profiles
- PCA projections
- Contact maps and comparisons
- Convergence analysis
"""

from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import seaborn as sns
from .config import AnnexinConfig, ProteinRegions, DEFAULT_CONFIG
from .analysis import RMSFResult, PCAResult, ConvergenceResult, DSSPResult
from .contact_map import ContactMapResult, ContactComparisonResult


class ConformationalVisualizer:
    """
    Class for creating visualizations of conformational analysis results.

    Provides methods for generating publication-ready plots with
    consistent styling and formatting.

    Attributes:
        config: Configuration object with visualization parameters.
        regions: Protein structural regions for annotation.
    """

    def __init__(self, config: Optional[AnnexinConfig] = None):
        """
        Initialize the ConformationalVisualizer.

        Args:
            config: Configuration object. Uses DEFAULT_CONFIG if not provided.
        """
        self.config = config or DEFAULT_CONFIG
        self.regions = self.config.regions

        # Set default matplotlib style
        plt.style.use('seaborn-v0_8-whitegrid')

    def _add_domain_annotations(
        self,
        ax: plt.Axes,
        show_legend: bool = True
    ) -> None:
        """
        Add N-terminal and core domain annotations to a plot.

        Args:
            ax: Matplotlib axes object.
            show_legend: Whether to add legend entries.
        """
        label_nterm = "N-Terminal" if show_legend else None
        label_core = "Annexin Core" if show_legend else None

        ax.axvspan(
            self.regions.N_TERMINAL_START,
            self.regions.N_TERMINAL_END,
            color="gray", alpha=0.15, label=label_nterm
        )
        ax.axvspan(
            self.regions.CORE_START,
            self.regions.CORE_END,
            color="blue", alpha=0.05, label=label_core
        )

    def _save_figure(
        self,
        fig: plt.Figure,
        filename: str,
        output_dir: Optional[str] = None
    ) -> str:
        """
        Save figure to file.

        Args:
            fig: Matplotlib figure object.
            filename: Output filename.
            output_dir: Directory to save in. Uses config default if None.

        Returns:
            Full path to saved file.
        """
        if output_dir is None:
            output_dir = self.config.output_dir

        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        filepath = output_path / filename
        fig.savefig(filepath, dpi=self.config.figure_dpi, bbox_inches="tight")
        print(f"Figure saved: {filepath}")

        return str(filepath)

    # ==================== RMSF Visualizations ====================

    def plot_rmsf(
        self,
        rmsf_result: RMSFResult,
        title: Optional[str] = None,
        filename: Optional[str] = None,
        show: bool = True,
        color: str = "#2c3e50"
    ) -> plt.Figure:
        """
        Plot RMSF profile for a single variant.

        Args:
            rmsf_result: RMSFResult object.
            title: Plot title. Auto-generated if None.
            filename: Output filename. Not saved if None.
            show: Whether to display the plot.
            color: Line color.

        Returns:
            Matplotlib figure object.
        """
        fig, ax = plt.subplots(figsize=(10, 5))

        ax.plot(
            rmsf_result.residue_numbers,
            rmsf_result.fluctuations,
            color=color,
            linewidth=1.5
        )

        self._add_domain_annotations(ax)

        if title is None:
            title = f"Conformational Profile: RMSF - {rmsf_result.variant_name}"

        ax.set_title(title, fontsize=14)
        ax.set_xlabel("Residue Number", fontsize=12)
        ax.set_ylabel(r"RMSF ($\AA$)", fontsize=12)
        ax.legend(loc="upper right")
        ax.grid(True, linestyle="--", alpha=0.6)

        plt.tight_layout()

        if filename:
            self._save_figure(fig, filename)

        if show:
            plt.show()

        return fig

    def plot_rmsf_comparison(
        self,
        rmsf_results: Dict[str, RMSFResult],
        colors: Optional[List[str]] = None,
        title: str = "RMSF Comparison - Annexin A11",
        filename: Optional[str] = None,
        show: bool = True
    ) -> plt.Figure:
        """
        Plot RMSF comparison for multiple variants.

        Args:
            rmsf_results: Dictionary mapping variant names to RMSFResult.
            colors: List of colors for each variant. Uses config defaults if None.
            title: Plot title.
            filename: Output filename. Not saved if None.
            show: Whether to display the plot.

        Returns:
            Matplotlib figure object.
        """
        if colors is None:
            colors = self.config.colors

        fig, ax = plt.subplots(figsize=(12, 6))

        for i, (name, result) in enumerate(rmsf_results.items()):
            ax.plot(
                result.residue_numbers,
                result.fluctuations,
                label=name,
                color=colors[i % len(colors)],
                linewidth=1.5,
                alpha=0.8
            )

        self._add_domain_annotations(ax)

        ax.set_title(title, fontsize=14)
        ax.set_xlabel("Residue Number", fontsize=12)
        ax.set_ylabel(r"RMSF ($\AA$)", fontsize=12)
        ax.legend(loc="upper right")
        ax.grid(True, linestyle="--", alpha=0.6)

        plt.tight_layout()

        if filename:
            self._save_figure(fig, filename)

        if show:
            plt.show()

        return fig

    # ==================== PCA Visualizations ====================

    def plot_pca(
        self,
        pca_result: PCAResult,
        title: Optional[str] = None,
        filename: Optional[str] = None,
        show: bool = True,
        cmap: str = "viridis"
    ) -> plt.Figure:
        """
        Plot PCA projection colored by frame index.

        Args:
            pca_result: PCAResult object.
            title: Plot title. Auto-generated if None.
            filename: Output filename. Not saved if None.
            show: Whether to display the plot.
            cmap: Colormap for frame coloring.

        Returns:
            Matplotlib figure object.
        """
        fig, ax = plt.subplots(figsize=(8, 6))

        scatter = ax.scatter(
            pca_result.reduced_coords[:, 0],
            pca_result.reduced_coords[:, 1],
            c=range(pca_result.n_frames),
            cmap=cmap,
            alpha=0.6,
            s=15
        )

        plt.colorbar(scatter, label="Frame Index")

        if title is None:
            title = f"PCA - {pca_result.variant_name}"

        ax.set_title(title, fontsize=14)
        ax.set_xlabel(f"PC1 ({pca_result.pc1_variance:.1f}%)", fontsize=12)
        ax.set_ylabel(f"PC2 ({pca_result.pc2_variance:.1f}%)", fontsize=12)
        ax.grid(True, alpha=0.3)

        plt.tight_layout()

        if filename:
            self._save_figure(fig, filename)

        if show:
            plt.show()

        return fig

    def plot_pca_comparison(
        self,
        wt_pca: PCAResult,
        mutant_pca_results: Dict[str, PCAResult],
        wt_label: str = "Wild Type",
        wt_color: str = "blue",
        mutant_colors: Optional[Dict[str, str]] = None,
        title: str = "PCA Comparison - Conformational Space",
        filename: Optional[str] = None,
        show: bool = True
    ) -> plt.Figure:
        """
        Plot PCA comparison between WT and mutants.

        Args:
            wt_pca: PCAResult for wild-type.
            mutant_pca_results: Dictionary of mutant name -> PCAResult.
            wt_label: Display label for WT.
            wt_color: Color for WT scatter.
            mutant_colors: Dictionary of mutant name -> color.
            title: Plot title.
            filename: Output filename. Not saved if None.
            show: Whether to display the plot.

        Returns:
            Matplotlib figure object.
        """
        fig, ax = plt.subplots(figsize=(10, 7))

        # Plot WT
        ax.scatter(
            wt_pca.reduced_coords[:, 0],
            wt_pca.reduced_coords[:, 1],
            c=wt_color,
            alpha=0.4,
            label=wt_label,
            s=10
        )

        # Plot mutants
        default_colors = ["red", "green", "orange", "purple"]

        for i, (name, result) in enumerate(mutant_pca_results.items()):
            if mutant_colors and name in mutant_colors:
                color = mutant_colors[name]
            else:
                color = default_colors[i % len(default_colors)]

            ax.scatter(
                result.reduced_coords[:, 0],
                result.reduced_coords[:, 1],
                c=color,
                alpha=0.4,
                label=name,
                s=10
            )

        ax.set_title(title, fontsize=14)
        ax.set_xlabel(f"PC1 ({wt_pca.pc1_variance:.1f}%)", fontsize=12)
        ax.set_ylabel(f"PC2 ({wt_pca.pc2_variance:.1f}%)", fontsize=12)
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)

        plt.tight_layout()

        if filename:
            self._save_figure(fig, filename)

        if show:
            plt.show()

        return fig

    # ==================== Contact Map Visualizations ====================

    def plot_contact_map(
        self,
        contact_result: ContactMapResult,
        title: Optional[str] = None,
        filename: Optional[str] = None,
        show: bool = True,
        cmap: str = "viridis"
    ) -> plt.Figure:
        """
        Plot a single contact map.

        Args:
            contact_result: ContactMapResult object.
            title: Plot title. Auto-generated if None.
            filename: Output filename. Not saved if None.
            show: Whether to display the plot.
            cmap: Colormap for heatmap.

        Returns:
            Matplotlib figure object.
        """
        fig, ax = plt.subplots(figsize=(10, 8))

        sns.heatmap(
            contact_result.distance_matrix,
            cmap=cmap,
            square=True,
            xticklabels=100,
            yticklabels=100,
            ax=ax
        )

        if title is None:
            title = f"Contact Map - {contact_result.variant_name}"

        ax.set_title(title, fontsize=14)
        ax.set_xlabel("Residue Index", fontsize=12)
        ax.set_ylabel("Residue Index", fontsize=12)

        plt.tight_layout()

        if filename:
            self._save_figure(fig, filename)

        if show:
            plt.show()

        return fig

    def plot_contact_comparison(
        self,
        comparison: ContactComparisonResult,
        mutation_position: Optional[int] = None,
        title: str = "Contact Map Comparison",
        filename: Optional[str] = None,
        show: bool = True
    ) -> plt.Figure:
        """
        Plot contact map comparison (binary and distance difference).

        Args:
            comparison: ContactComparisonResult object.
            mutation_position: Residue position of mutation (for annotation).
            title: Plot title.
            filename: Output filename. Not saved if None.
            show: Whether to display the plot.

        Returns:
            Matplotlib figure object.
        """
        fig, axes = plt.subplots(1, 2, figsize=(20, 10))

        # Plot 1: Binary contact difference
        sns.heatmap(
            comparison.contact_difference,
            cmap="RdBu_r",
            center=0,
            square=True,
            ax=axes[0]
        )
        axes[0].set_title(
            "Contact Map Difference (Binary)\n(Blue: Gained | Red: Lost)",
            fontsize=12
        )
        axes[0].set_xlabel("Residue Index", fontsize=10)
        axes[0].set_ylabel("Residue Index", fontsize=10)

        # Plot 2: Distance difference (continuous)
        vmax = np.abs(comparison.distance_difference).max() * 0.7
        sns.heatmap(
            comparison.distance_difference,
            cmap="bwr",
            center=0,
            square=True,
            vmin=-vmax,
            vmax=vmax,
            ax=axes[1]
        )
        axes[1].set_title(
            "Mean Distance Difference (nm)\n(Red: Increased | Blue: Decreased)",
            fontsize=12
        )
        axes[1].set_xlabel("Residue Index", fontsize=10)
        axes[1].set_ylabel("Residue Index", fontsize=10)

        # Add annotations for mutation and N-terminal boundary
        for ax in axes:
            if mutation_position is not None:
                ax.axhline(
                    mutation_position, color="green", linestyle=":",
                    alpha=0.7, linewidth=2
                )
                ax.axvline(
                    mutation_position, color="green", linestyle=":",
                    alpha=0.7, linewidth=2
                )

            # N-terminal boundary
            ax.axhline(
                self.regions.N_TERMINAL_END, color="black",
                linewidth=1, alpha=0.5
            )
            ax.axvline(
                self.regions.N_TERMINAL_END, color="black",
                linewidth=1, alpha=0.5
            )

        fig.suptitle(title, fontsize=14, y=1.02)
        plt.tight_layout()

        if filename:
            self._save_figure(fig, filename)

        if show:
            plt.show()

        return fig

    # ==================== Convergence Visualization ====================

    def plot_convergence(
        self,
        convergence_result: ConvergenceResult,
        title: str = "RMSF Convergence Analysis",
        filename: Optional[str] = None,
        show: bool = True
    ) -> plt.Figure:
        """
        Plot convergence analysis results.

        Args:
            convergence_result: ConvergenceResult object.
            title: Plot title.
            filename: Output filename. Not saved if None.
            show: Whether to display the plot.

        Returns:
            Matplotlib figure object.
        """
        fig, ax = plt.subplots(figsize=(10, 6))

        # Define o mapa de cores (ex: 'viridis', 'plasma' ou 'Blues')
        n_samples = len(convergence_result.sample_sizes)
        colors = cm.Blues(np.linspace(0.3, 1, n_samples))

        for key, sample in enumerate(convergence_result.sample_sizes):
            is_last = (key == n_samples - 1)

            ax.plot(
                convergence_result.rmsf_list[key],
                "-",
                color=colors[key],
                linewidth=2.5 if is_last else 1.0,
                alpha=1.0 if is_last else 0.4,
                label=f"Final ({sample})" if is_last else None,
                zorder=10 if is_last else 1
            )

        ax.set_title(title, fontsize=14)
        ax.set_xlabel("Residue Number", fontsize=12)
        ax.set_ylabel("RMSF Difference (L2 norm)", fontsize=12)
        ax.grid(True, linestyle="--", alpha=0.6)

        # Add convergence status annotation
        status = "Converged" if convergence_result.is_converged else "Not Converged"
        ax.annotate(
            f"Status: {status}",
            xy=(0.95, 0.95),
            xycoords="axes fraction",
            fontsize=10,
            ha="right",
            va="top",
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5)
        )

        plt.tight_layout()

        if filename:
            self._save_figure(fig, filename)

        if show:
            plt.show()

        return fig

    def plot_dssp(
        self,
        dssp_result: DSSPResult,
        title: str = "Secondary Structures Analysis",
        filename=None,
        show=True
    ):

        data = dssp_result.raw_matrix

        # 3. Mapeamento categórico rigoroso
        # 0: Coil, 1: Helix (H, G, I), 2: Strand (E, B)
        numeric_map = np.zeros(data.shape)
        numeric_map[np.isin(data, ['H', 'G', 'I'])] = 1
        numeric_map[np.isin(data, ['E', 'B'])] = 2

        fig, ax = plt.subplots(figsize=(12, 7))

        # Cores discretas: Branco (Coil), Azul (Hélice), Vermelho (Strand)
        cmap = mcolors.ListedColormap(['white', '#3498db', '#e74c3c'])
        bounds = [0, 0.5, 1.5, 2.5]
        norm = mcolors.BoundaryNorm(bounds, cmap.N)

        sns.heatmap(numeric_map.T, cmap=cmap, norm=norm, ax=ax, cbar_kws={
            'ticks': [0, 1, 2],
            'label': 'Structure'
        })

        # Ajustar labels do Colorbar
        cbar = ax.collections[0].colorbar
        cbar.set_ticklabels(['Coil/Turn', 'Helix', 'Strand'])

        ax.set_xlabel('Frames', fontsize=12)
        ax.set_ylabel('Residue ID', fontsize=12)
        ax.set_title(f"{title} - {dssp_result.variant_name}", fontsize=14)

        if filename:
            self._save_figure(fig, filename)
        if show:
            plt.show()
        return fig
    

    def plot_core_idr_mean_distances(
        self,
        distances: np.ndarray,
        title: str = "Mean Distances between IDR Residues to CORE",
        filename: Optional[str] = None,
        show: bool = True
    ) -> plt.Figure:
        """
        Plot the mean distance of each IDR residue to the exposed core surface.
        """
        fig, ax = plt.subplots(figsize=(10, 5))

        # Cria o eixo X baseado nos índices dos resíduos do IDR
        res_indices = np.arange(
            self.regions.N_TERMINAL_START,
            self.regions.N_TERMINAL_START + len(distances)
        )

        ax.plot(res_indices, distances, color='#2c3e50', linewidth=2)
        ax.fill_between(res_indices, distances, alpha=0.2, color='#2c3e50')

        ax.set_xlabel('IDR Residue Index', fontsize=12)
        ax.set_ylabel(r'Mean Distance to Core ($nm$)', fontsize=12)
        ax.set_title(title, fontsize=14)
        ax.grid(True, linestyle="--", alpha=0.6)
        
        # Limita o eixo X exatamente à região do N-terminal
        ax.set_xlim(self.regions.N_TERMINAL_START, self.regions.N_TERMINAL_END)

        plt.tight_layout()

        if filename:
            self._save_figure(fig, filename)
        if show:
            plt.show()

        return fig


    def plot_dssp_profile_comparison(
        self,
        wt_result: DSSPResult,
        mut_result: DSSPResult,
        structure_type='Helix_%',
        title: str = "DSSP Profile Comparisson",
        filename: Optional[str] = None,
        show: bool = True
    ) -> plt.Figure:
        """
        Compair a probability of a specific strucure between WT and Mutant.
        structure_type: 'Helix_%', 'Strand_%' ou 'Coil_%'
        """

        res_ids = np.array(list(wt_result.percentages.keys()))
        wt_vals = np.array([v[structure_type] for v in wt_result.percentages.values()])
        mut_vals = np.array([v[structure_type] for v in mut_result.percentages.values()])

        fig, ax = plt.subplots(figsize=(12, 7))

        ax.plot(res_ids, wt_vals, label='WT', color='blue', alpha=0.6)
        ax.plot(res_ids, mut_vals, label=mut_result.variant_name, color='red', alpha=0.6)

        # Destacar a região da mutação (ex: resíduo 36)
        ax.axvline(x=36, color='black', linestyle='--', label='Mutation Site')

        ax.set_ylabel(f'{structure_type} Occupation')
        ax.set_xlabel('Residue ID')
        ax.legend()

        ax.set_title(f"{structure_type.replace('_%', '')} {title} - WT vs {mut_result.variant_name} - Residue {res_ids[0]} to {res_ids[-1]}", fontsize=14)


        if filename:
            self._save_figure(fig, filename)
        if show:
            plt.show()

        return fig
    



