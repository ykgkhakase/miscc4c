from __future__ import annotations

from pathlib import Path

import requests


PDBJ_ENTRY_PDB_GZ_URL = "https://pdbj.org/rest/newweb/fetch/file?type=pdb&id={pdb_id}&format=gz"
RCSB_ENTRY_PDB_GZ_URL = "https://files.rcsb.org/download/{pdb_id}.pdb.gz"


def _download_bytes(url: str, timeout: int = 20) -> bytes:
    r = requests.get(url, timeout=timeout)
    r.raise_for_status()
    return r.content


def download_pdb_gz(
    pdb_id: str,
    out_dir: str | Path = "pdb_cache",
    cache: bool = True,
    timeout: int = 20,
) -> Path:
    """
    PDB ID のエントリを .pdb.gz 形式で取得する。
    """
    pdb_id = pdb_id.upper()

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    out_file = out_dir / f"{pdb_id}.pdb.gz"

    if cache and out_file.exists() and out_file.stat().st_size > 0:
        return out_file

    url_pdbj = PDBJ_ENTRY_PDB_GZ_URL.format(pdb_id=pdb_id)
    try:
        data = _download_bytes(url_pdbj, timeout)
        out_file.write_bytes(data)
        return out_file
    except Exception:
        pass

    url_rcsb = RCSB_ENTRY_PDB_GZ_URL.format(pdb_id=pdb_id)
    data = _download_bytes(url_rcsb, timeout)
    out_file.write_bytes(data)

    return out_file
