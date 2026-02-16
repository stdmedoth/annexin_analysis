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

### Directory Structure

```
protein_analysis/
├── README.md
├── main.py                          # Main analysis script (RMSF calculation)
├── helpers.py                       # Utility functions (NPZ to PDB conversion)
├── pdbs/                            # Converted PDB files for analysis
├── MutationConformationsBioEmu/
│   ├── out_native/                  # Wild-type conformations
│   │   ├── batch_*.npz              # Conformational ensembles
│   │   ├── sequence.fasta           # WT sequence
│   │   └── samples.xtc              # Trajectory file
│   └── out_mutant_P36R/             # P36R mutant conformations
│       ├── batch_*.npz              # Conformational ensembles
│       └── sequence.fasta           # Mutant sequence
```

---

## 🔬 Analysis Methods

### Planned Analyses

| Analysis | Description | Status |
|----------|-------------|--------|
| **RMSF** | Root Mean Square Fluctuation per residue | ✅ Implemented |
| **RMSD** | Root Mean Square Deviation over conformations | 🔄 In progress |
| **PCA** | Principal Component Analysis of conformational space | ⏳ Planned |
| **Contact Maps** | Difference contact maps (WT vs Mutant) | ⏳ Planned |

### Current Implementation

1. **NPZ to PDB Conversion** (`helpers.py`)
   - Converts BioEmu output to PDB format for visualization and analysis

2. **RMSF Analysis** (`main.py`)
   - Uses AlphaFold structure as reference (`AF-P50995-F1-model_v6.pdb`)
   - Focuses on the Annexin core domain (residues 199-504)
   - Generates conformational profile plots

---

## 🛠️ Dependencies

```bash
pip install numpy mdtraj matplotlib
```

| Package | Purpose |
|---------|---------|
| `numpy` | Numerical computations, NPZ file handling |
| `mdtraj` | Molecular dynamics trajectory analysis |
| `matplotlib` | Visualization and plotting |

---

## 🚀 Usage

### 1. Convert NPZ files to PDB

```bash
python helpers.py
```

### 2. Run RMSF Analysis

```bash
python main.py
```

This will generate `conformational_profile.png` showing the RMSF per residue.

---

## 📊 Expected Outputs

- **Conformational Profile Plot**: RMSF values (in Ångströms) plotted against residue number
- **Comparative Analysis**: Side-by-side comparison of WT and mutant flexibility profiles
- **PCA Projections**: Visualization of conformational sampling in reduced dimensions
- **Difference Contact Maps**: Highlighting altered inter-residue interactions

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

## 📚 References

1. AlphaFold Protein Structure Database - [https://alphafold.ebi.ac.uk/](https://alphafold.ebi.ac.uk/)
2. BioEmu - Conformational ensemble generation tool
3. MDTraj - [https://mdtraj.org/](https://mdtraj.org/)

---

## 👨‍🔬 Author

**Undergraduate Research Student**  
University of São Paulo (USP)  
7th Semester - Scientific Initiation (IC) Project

---

## 📄 License

This project is for academic research purposes.

---

*Last updated: February 2026*
