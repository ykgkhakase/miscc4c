from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import Iterable

from rdkit import Chem

from .assign_bond_order import (
    add_neutral_hydrogens,
    assign_bond_orders_from_template,
    build_template_from_smiles,
    extract_ligand_pdb_block,
    extract_smiles_from_ccd_cif,
    read_text_auto,
    mol_from_pdb_block,
)
from .download_cif import download_ligand_cif
from .download_pdb import download_pdb_gz


def _iter_residue_candidates(pdb_path: str, resname: str) -> list[tuple[str, int]]:
    text = read_text_auto(pdb_path)
    candidates: set[tuple[str, int]] = set()
    target = resname.upper()

    for line in text.splitlines():
        record = line[:6].strip()
        if record not in {"ATOM", "HETATM"}:
            continue
        if line[17:20].strip().upper() != target:
            continue

        altloc = line[16].strip()
        if altloc not in {"", "A"}:
            continue

        chain = line[21].strip()
        resnum_str = line[22:26].strip()
        try:
            resnum = int(resnum_str)
        except ValueError:
            continue
        candidates.add((chain, resnum))

    return sorted(candidates, key=lambda x: (x[0], x[1]))


def _format_candidates(candidates: Iterable[tuple[str, int]]) -> str:
    return ", ".join(f"(chain={c}, resnum={r})" for c, r in candidates)


def _pdb_label_from_path(pdb_path: str) -> str:
    name = Path(pdb_path).name
    if name.endswith(".pdb.gz"):
        return name[: -len(".pdb.gz")]
    if name.endswith(".pdb"):
        return name[: -len(".pdb")]
    if name.endswith(".gz"):
        return name[: -len(".gz")]
    return Path(pdb_path).stem


def run_assign_ligand_bond(args: argparse.Namespace) -> int:
    pdb_path = str(args.pdb_path)
    cif_path = str(args.cif_path)
    out_dir = Path(args.out)
    resname = args.resname.upper()
    chain = args.chain
    resnum = args.resnum

    if not Path(pdb_path).exists():
        print(f"ERROR: PDB file not found: {pdb_path}", file=sys.stderr)
        return 1
    if not Path(cif_path).exists():
        print(f"ERROR: CIF file not found: {cif_path}", file=sys.stderr)
        return 1

    candidates = _iter_residue_candidates(pdb_path, resname)
    if not candidates:
        print(
            f"ERROR: residue not found for resname={resname} in {pdb_path}",
            file=sys.stderr,
        )
        return 1

    filtered = [
        (c, r)
        for c, r in candidates
        if (chain is None or c == chain) and (resnum is None or r == resnum)
    ]
    if not filtered:
        print(
            (
                "ERROR: no residue matched filters. "
                f"resname={resname}, chain={chain}, resnum={resnum}. "
                f"Available: {_format_candidates(candidates)}"
            ),
            file=sys.stderr,
        )
        return 1

    if len(filtered) > 1 and (chain is None or resnum is None):
        print(
            (
                "ERROR: multiple residues matched. "
                "Specify both --chain and --resnum.\n"
                f"Matched: {_format_candidates(filtered)}"
            ),
            file=sys.stderr,
        )
        return 1

    selected_chain, selected_resnum = filtered[0]

    try:
        smiles = extract_smiles_from_ccd_cif(cif_path)
        pdb_block = extract_ligand_pdb_block(
            pdb_path=pdb_path,
            resname=resname,
            chain=selected_chain,
            resnum=selected_resnum,
        )
        pdb_mol = mol_from_pdb_block(pdb_block)
        template = build_template_from_smiles(smiles)
        assigned = assign_bond_orders_from_template(pdb_mol, template)
        result_mol = add_neutral_hydrogens(assigned) if args.add_neutral_h else assigned
    except Exception as exc:
        print(f"ERROR: assign-ligand-bond failed: {exc}", file=sys.stderr)
        return 1

    out_dir.mkdir(parents=True, exist_ok=True)
    pdb_label = _pdb_label_from_path(pdb_path)
    out_file = out_dir / f"{pdb_label}_{resname}_{selected_chain}_{selected_resnum}.sdf"

    writer = Chem.SDWriter(str(out_file))
    if writer is None:
        print(f"ERROR: failed to open output SDF: {out_file}", file=sys.stderr)
        return 1
    writer.write(result_mol)
    writer.close()

    print(f"Saved: {out_file}")
    return 0


def run_download_ligand(args: argparse.Namespace) -> int:
    try:
        path = download_ligand_cif(
            args.ccd_id,
            out_dir=args.out,
            cache=not args.no_cache,
            timeout=args.timeout,
        )
    except Exception as exc:
        print(f"ERROR: download-ligand failed: {exc}", file=sys.stderr)
        return 1
    print(f"Saved: {path}")
    return 0


def run_download_pdbgz(args: argparse.Namespace) -> int:
    try:
        path = download_pdb_gz(
            args.pdb_id,
            out_dir=args.out,
            cache=not args.no_cache,
            timeout=args.timeout,
        )
    except Exception as exc:
        print(f"ERROR: download-pdbgz failed: {exc}", file=sys.stderr)
        return 1
    print(f"Saved: {path}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="miscc4c")
    subparsers = parser.add_subparsers(dest="command", required=True)

    assign = subparsers.add_parser(
        "assign-ligand-bond",
        help="Assign bond orders to a ligand extracted from a PDB structure.",
    )
    assign.add_argument("pdb_path", help="Input PDB file path (.pdb or .pdb.gz).")
    assign.add_argument("cif_path", help="Input ligand CCD CIF file path.")
    assign.add_argument("--resname", required=True, help="Ligand residue name (e.g. 1RA).")
    assign.add_argument("--chain", default=None, help="Chain ID.")
    assign.add_argument("--resnum", type=int, default=None, help="Residue number.")
    assign.add_argument("--out", required=True, help="Output directory path.")
    assign.add_argument(
        "--add-neutral-h",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Add neutral hydrogens after assigning bond orders (default: true).",
    )
    assign.set_defaults(func=run_assign_ligand_bond)

    dl_lig = subparsers.add_parser(
        "download-ligand",
        help="Download ligand CIF by CCD ID.",
    )
    dl_lig.add_argument("ccd_id", help="CCD ID (e.g. 1RA).")
    dl_lig.add_argument("--out", default="ligand_cache", help="Output directory path.")
    dl_lig.add_argument("--timeout", type=int, default=20, help="Request timeout in seconds.")
    dl_lig.add_argument("--no-cache", action="store_true", help="Disable cache reuse.")
    dl_lig.set_defaults(func=run_download_ligand)

    dl_pdb = subparsers.add_parser(
        "download-pdbgz",
        help="Download entry PDB as .pdb.gz by PDB ID.",
    )
    dl_pdb.add_argument("pdb_id", help="PDB ID (e.g. 4KIQ).")
    dl_pdb.add_argument("--out", default="pdb_cache", help="Output directory path.")
    dl_pdb.add_argument("--timeout", type=int, default=20, help="Request timeout in seconds.")
    dl_pdb.add_argument("--no-cache", action="store_true", help="Disable cache reuse.")
    dl_pdb.set_defaults(func=run_download_pdbgz)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)
