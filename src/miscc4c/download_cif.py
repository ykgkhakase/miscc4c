from __future__ import annotations

import gzip
from pathlib import Path
from typing import Iterable

import requests
import re


PDBJ_LIGAND_URL = "https://pdbj.org/rest/cc/{ccd_id}"
RCSB_LIGAND_URL = "https://files.rcsb.org/ligands/download/{ccd_id}.cif"

RCSB_ENTRY_CIF_URL = "https://files.rcsb.org/download/{pdb_id}.cif"


def _download_bytes(url: str, timeout: int = 20) -> bytes:
    r = requests.get(url, timeout=timeout)
    r.raise_for_status()
    return r.content


def download_ligand_cif(
    ccd_id: str,
    out_dir: str | Path = "ligand_cache",
    cache: bool = True,
    timeout: int = 20,
) -> Path:
    """
    CCD ID の ligand CIF を取得する

    優先順位
    1. PDBj
    2. RCSB PDB
    """

    ccd_id = ccd_id.upper()

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    out_file = out_dir / f"{ccd_id}.cif"

    if cache and out_file.exists() and out_file.stat().st_size > 0:
        return out_file

    # ---- PDBj first ----
    url_pdbj = PDBJ_LIGAND_URL.format(ccd_id=ccd_id)

    try:
        data = _download_bytes(url_pdbj, timeout)
        out_file.write_bytes(data)
        return out_file
    except Exception:
        pass

    # ---- fallback: RCSB ----
    url_rcsb = RCSB_LIGAND_URL.format(ccd_id=ccd_id)

    data = _download_bytes(url_rcsb, timeout)
    out_file.write_bytes(data)

    return out_file


def download_ligand_cif_many(
    ccd_ids: Iterable[str],
    out_dir: str | Path = "ligand_cache",
    cache: bool = True,
    timeout: int = 20,
    skip_errors: bool = False,
):
    """
    複数 CCD を一括ダウンロード
    """

    results = {}

    for cid in ccd_ids:

        try:
            path = download_ligand_cif(
                cid,
                out_dir=out_dir,
                cache=cache,
                timeout=timeout,
            )
            results[cid] = path

        except Exception:

            if skip_errors:
                results[cid] = None
            else:
                raise

    return results


def download_pdb_mmcif(
    pdb_id: str,
    cache_dir: str | Path = "pdb_cache",
    cache: bool = True,
    timeout: int = 20,
) -> Path:

    pdb_id = pdb_id.upper()

    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)

    out_file = cache_dir / f"{pdb_id}.cif"

    if cache and out_file.exists():
        return out_file

    url = RCSB_ENTRY_CIF_URL.format(pdb_id=pdb_id)

    data = _download_bytes(url, timeout)

    out_file.write_bytes(data)

    return out_file

def _split_cif_tokens(line: str) -> list[str]:
    pattern = r"""(?:'[^']*'|"[^"]*"|\S+)"""
    tokens = re.findall(pattern, line.strip())
    out = []
    for t in tokens:
        if len(t) >= 2 and t[0] == t[-1] and t[0] in {"'", '"'}:
            out.append(t[1:-1])
        else:
            out.append(t)
    return out

def extract_ligand_ccd_ids_from_mmcif_text(
    text: str,
    exclude_common_solvents: bool = True,
) -> list[str]:
    """
    entry mmCIF の _chem_comp loop から ligand CCD ID を抽出する。
    """
    lines = text.splitlines()

    for i, line in enumerate(lines):
        if line.strip() != "loop_":
            continue

        headers = []
        j = i + 1
        while j < len(lines) and lines[j].strip().startswith("_"):
            headers.append(lines[j].strip())
            j += 1

        if "_chem_comp.id" not in headers or "_chem_comp.type" not in headers:
            continue

        col_idx = {h: k for k, h in enumerate(headers)}
        id_idx = col_idx["_chem_comp.id"]
        type_idx = col_idx["_chem_comp.type"]

        ligand_ids = []
        while j < len(lines):
            s = lines[j].strip()

            if not s:
                j += 1
                continue

            if s == "#":
                break

            if s.startswith("loop_") or s.startswith("_"):
                break

            row = _split_cif_tokens(s)
            if len(row) > max(id_idx, type_idx):
                comp_id = row[id_idx].upper()
                comp_type = row[type_idx].lower()

                if ("non-polymer" in comp_type) or ("branched" in comp_type):
                    ligand_ids.append(comp_id)

            j += 1

        seen = set()
        ligand_ids = [x for x in ligand_ids if not (x in seen or seen.add(x))]

        if exclude_common_solvents:
            common_exclude = {
                "HOH", "DOD", "WAT",
                "NA", "K", "CL", "MG", "CA", "ZN", "MN", "CU", "CO", "FE",
                "SO4", "PO4", "GOL", "EDO", "PEG", "ACT", "ACY", "FMT", "EOH",
            }
            ligand_ids = [x for x in ligand_ids if x not in common_exclude]

        return ligand_ids

    raise ValueError("mmCIF 中に _chem_comp loop が見つかりませんでした。")


def extract_ligand_ccd_ids_from_pdb_id(
    pdb_id: str,
    cache_dir: str | Path = "pdb_cache",
    cache: bool = True,
    timeout: int = 20,
    exclude_common_solvents: bool = True,
) -> list[str]:
    cif_path = download_pdb_mmcif(
        pdb_id=pdb_id,
        cache_dir=cache_dir,
        cache=cache,
        timeout=timeout,
    )
    text = cif_path.read_text(encoding="utf-8")
    return extract_ligand_ccd_ids_from_mmcif_text(
        text,
        exclude_common_solvents=exclude_common_solvents,
    )

def download_all_ligand_cifs_for_pdb(
    pdb_id: str,
    ligand_out_dir: str | Path = "ligand_cache",
    entry_cache_dir: str | Path = "pdb_cache",
    cache: bool = True,
    timeout: int = 20,
    exclude_common_solvents: bool = True,
    skip_errors: bool = False,
):
    ccd_ids = extract_ligand_ccd_ids_from_pdb_id(
        pdb_id,
        cache_dir=entry_cache_dir,
        cache=cache,
        timeout=timeout,
        exclude_common_solvents=exclude_common_solvents,
    )

    return download_ligand_cif_many(
        ccd_ids,
        out_dir=ligand_out_dir,
        cache=cache,
        timeout=timeout,
        skip_errors=skip_errors,
    )