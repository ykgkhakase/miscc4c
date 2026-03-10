# miscc4c

Miscellaneous scripts for computational chemistry for compounds.

CLI utilities for ligand extraction, download, and bond-order assignment.

## Overview

This project provides a subcommand-based CLI:

- `assign-ligand-bond`
- `download-ligand`
- `download-pdbgz`

The recommended entrypoint is the local wrapper script `./miscc4c`, which always uses this project's `.venv/bin/python`.

## Dependencies

Core runtime dependencies:

- `rdkit`
- `requests`

Development/testing:

- `pytest`

Tested versions in this project environment:

- `rdkit==2025.09.6`
- `pytest==9.0.2`

Dependency versions are managed in `pyproject.toml` and installed via `uv sync`.

## Setup

Install `uv` first (Linux/macOS):

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

Then set up this project:

```bash
# from project root
uv sync
chmod +x miscc4c
```

## How to Run the Wrapper Script

### 1) Run with an absolute path

```bash
./miscc4c --help
or 
<PROJECT_ROOT>/miscc4c --help
```

### 2) Put it in $PATH

```bash
e.g.)
ln -sf <PROJECT_ROOT>/miscc4c ~/path/to/your/home/bin/miscc4c
```

## Command Usage

### 1) Download ligand CIF (`download-ligand`)

Download a ligand CIF by CCD ID (for example, `1RA.cif`):

```bash
./miscc4c download-ligand 1RA --out mytestcache/ligands
```

Optional flags:

- `--timeout 30`
- `--no-cache`

### 2) Download PDB gzip (`download-pdbgz`)

Download an entry structure as `.pdb.gz` by PDB ID (for example, `4KIQ.pdb.gz`):

```bash
./miscc4c download-pdbgz 4KIQ --out mytestcache/pdbs
```

Optional flags:

- `--timeout 30`
- `--no-cache`

### 3) Assign ligand bond orders (`assign-ligand-bond`)

Create a ligand-only SDF with assigned bond orders from PDB coordinates + ligand CIF template.

Example using sample PDB data in `public_data/pdb`:

```bash
./miscc4c assign-ligand-bond \
  public_data/pdb/4KIQ.pdb.gz \
  mytestcache/ligands/1RA.cif \
  --resname 1RA \
  --chain A \
  --resnum 401 \
  --out mytestcache/results
```

Output filename is generated as:

```text
<PDBID>_<RESNAME>_<CHAIN>_<RESNUM>.sdf
```

Example:

```text
4KIQ_1RA_A_401.sdf
```

Optional flags:

- `--no-add-neutral-h` (disable hydrogen addition)

## Notes

- `--out` is required for `assign-ligand-bond`.
- If multiple residues match and `--chain` / `--resnum` are not specific enough, the CLI exits with a clear diagnostic message listing candidate residues.
