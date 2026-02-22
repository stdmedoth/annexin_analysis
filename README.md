# Annexin A11 Conformational Analysis

**Undergraduate Research Project - University of São Paulo (USP)**

A computational study investigating the conformational profile of Annexin A11 protein, comparing wild-type (WT) and mutant variants to understand the structural impact of disease-associated mutations.

---

## 📋 Project Overview

Annexin A11 is a calcium-dependent phospholipid-binding protein involved in various cellular processes. Mutations in ANXA11 have been linked to amyotrophic lateral sclerosis (ALS) and other neurological disorders. This project aims to characterize the conformational changes induced by specific mutations using molecular dynamics simulations and structural analysis.

### Mutation Under Study

- **P36R** — Proline to Arginine substitution at position 36 (N-terminal region)

---

## 🧬 Data Source

Conformational ensembles were generated using **BioEmu** (Biological Ensemble Modeling Utility), which produces diverse protein conformations stored in `.npz` format containing:
- `pos`: Atomic coordinates (Cα positions)
- `sequence`: Amino acid sequence

---

## 📁 Project Structure

```
protein_analysis/
├── README.md
├── annexin_analysis/              # Main analysis package
│   ├── __init__.py               # Package exports
│   ├── config.py                 # Configuration and constants
│   ├── trajectory.py             # Trajectory loading and preprocessing
│   ├── analysis.py               # RMSF, PCA, convergence analysis
│   ├── contact_map.py            # Contact map computation
│   ├── visualization.py          # Plotting functions
│   └── utils.py                  # Utility classes (NPZ converter, etc.)
│
├── scripts/                       # Ready-to-run analysis scripts
│   ├── analyze_wt_rmsf.py        # RMSF analysis for WT
│   ├── analyze_wt_pca.py         # PCA analysis for WT
│   ├── compare_rmsf_wt_mutant.py # RMSF comparison WT vs P36R
│   ├── compare_pca_wt_mutant.py  # PCA comparison WT vs P36R
│   ├── compare_contact_maps.py   # Contact map comparison
│   ├── analyze_convergence.py    # Convergence analysis
│   ├── analyze_local_pca.py      # Local/regional PCA analysis
│   ├── align_trajectories.py     # Trajectory alignment utility
│   ├── convert_npz_to_pdb.py     # NPZ to PDB converter
│   └── run_full_analysis.py      # Complete analysis pipeline
│
├── results/                       # Output directory for figures
│
└── 1000_samples/                  # Trajectory data
    └── MutationConformationsBioEmu/
        ├── out_native/            # Wild-type conformations
        │   ├── samples.xtc
        │   └── topology.pdb
        └── out_mutant_P36R/       # P36R mutant conformations
            ├── samples.xtc
            └── topology.pdb
```

---

## 🔬 Analysis Methods

| Analysis | Description | Script |
|----------|-------------|--------|
| **RMSF** | Root Mean Square Fluctuation per residue | `analyze_wt_rmsf.py` |
| **PCA** | Principal Component Analysis of conformations | `analyze_wt_pca.py` |
| **RMSF Comparison** | Compare flexibility between WT and mutant | `compare_rmsf_wt_mutant.py` |
| **PCA Comparison** | Compare conformational space sampling | `compare_pca_wt_mutant.py` |
| **Contact Maps** | Inter-residue distance analysis | `compare_contact_maps.py` |
| **Convergence** | Assess ensemble sampling completeness | `analyze_convergence.py` |
| **Local PCA** | Focused analysis on N-terminal region | `analyze_local_pca.py` |

---

## 🛠️ Dependencies

```bash
pip install numpy mdtraj matplotlib seaborn scikit-learn scipy
```

| Package | Purpose |
|---------|---------|
| `numpy` | Numerical computations, NPZ file handling |
| `mdtraj` | Molecular dynamics trajectory analysis |
| `matplotlib` | Visualization and plotting |
| `seaborn` | Statistical data visualization |
| `scikit-learn` | PCA and machine learning |
| `scipy` | Scientific computing |

---

## 🚀 Quick Start

### Running Individual Analyses

```bash
# Navigate to the project directory
cd protein_analysis

# Run RMSF analysis for wild-type
python scripts/analyze_wt_rmsf.py

# Compare WT vs P36R mutant RMSF
python scripts/compare_rmsf_wt_mutant.py

# Compare PCA projections
python scripts/compare_pca_wt_mutant.py

# Compare contact maps
python scripts/compare_contact_maps.py
```

### Running Complete Analysis Pipeline

```bash
python scripts/run_full_analysis.py
```

This runs all analyses and generates a summary report in `results/`.

---

## 📦 Package Usage

The `annexin_analysis` package can be imported for custom analyses:

```python
from annexin_analysis import (
    AnnexinConfig,
    TrajectoryLoader,
    ConformationalAnalyzer,
    ContactMapAnalyzer,
    ConformationalVisualizer
)

# Initialize with default configuration
config = AnnexinConfig()
loader = TrajectoryLoader(config)
analyzer = ConformationalAnalyzer(config)

# Load and process wild-type
wt = config.wt_variant
processed = loader.process_variant(wt)

# Compute RMSF
rmsf_result = analyzer.compute_rmsf(processed)
print(f"Mean RMSF: {rmsf_result.mean_rmsf:.4f} Å")

# Visualize
visualizer = ConformationalVisualizer(config)
visualizer.plot_rmsf(rmsf_result, filename="my_rmsf_plot.png")
```

### Adding New Variants

```python
from annexin_analysis import AnnexinConfig

config = AnnexinConfig()

# Create a new variant configuration
new_variant = config.create_variant(
    name="g154r",
    label="Mutant G154R",
    folder="out_mutant_G154R",
    color_index=2
)
```

---

## 🏗️ Architecture

The package follows object-oriented design principles:

### Core Classes

| Class | Description |
|-------|-------------|
| `AnnexinConfig` | Centralized configuration for paths and parameters |
| `TrajectoryLoader` | Load, preprocess, and align MD trajectories |
| `ConformationalAnalyzer` | RMSF, PCA, and convergence analysis |
| `ComparativeAnalyzer` | Compare multiple variants |
| `ContactMapAnalyzer` | Contact map computation and comparison |
| `ConformationalVisualizer` | Publication-ready plotting |

### Data Classes

| Class | Description |
|-------|-------------|
| `ProcessedTrajectory` | Container for processed trajectory data |
| `RMSFResult` | RMSF analysis results |
| `PCAResult` | PCA analysis results with fitted model |
| `ConvergenceResult` | Convergence analysis results |
| `ContactMapResult` | Contact map data |
| `ContactComparisonResult` | Contact comparison data |

---

## 📁 Protein Information

- **UniProt ID**: P50995
- **Protein**: Annexin A11 (Human)
- **Length**: 505 residues
- **Reference Structure**: AlphaFold predicted structure (v6)

### Domain Architecture

| Region | Residues | Description |
|--------|----------|-------------|
| N-terminal | 1-198 | Intrinsically disordered, proline-rich |
| Annexin Core | 199-505 | Conserved Ca²⁺-binding domains |

---

## 📊 Key Findings

The analysis demonstrates that the P36R mutation affects the conformational space of Annexin A11:

1. **RMSF Changes**: Altered flexibility profile in the N-terminal region
2. **PCA Shifts**: Different conformational sampling patterns
3. **Contact Changes**: Modified inter-residue interactions

---

## 📚 References

1. AlphaFold Protein Structure Database - [https://alphafold.ebi.ac.uk/](https://alphafold.ebi.ac.uk/)
2. BioEmu - Conformational ensemble generation tool
3. MDTraj - [https://mdtraj.org/](https://mdtraj.org/)

---

## 👨‍🔬 Author

**Calistu** - Undergraduate Research  
University of São Paulo (USP)  
7th Semester - Scientific Initiation (IC) Project

---

## 📄 License

This project is for academic research purposes.

---

*Last updated: February 2026*
