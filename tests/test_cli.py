from pathlib import Path

import pytest

from miscc4c import cli


@pytest.fixture(scope="session")
def cache_dir():
    return Path("mytestcache")


def test_download_ligand_subcommand(monkeypatch, cache_dir):
    out_dir = cache_dir / "cli_ligands"
    calls = {}

    def fake_download_ligand_cif(ccd_id, out_dir, cache, timeout):
        calls["ccd_id"] = ccd_id
        calls["out_dir"] = out_dir
        calls["cache"] = cache
        calls["timeout"] = timeout
        p = Path(out_dir) / f"{ccd_id.upper()}.cif"
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text("dummy", encoding="utf-8")
        return p

    monkeypatch.setattr(cli, "download_ligand_cif", fake_download_ligand_cif)

    rc = cli.main(["download-ligand", "1RA", "--out", str(out_dir), "--timeout", "11"])
    assert rc == 0
    assert calls == {
        "ccd_id": "1RA",
        "out_dir": str(out_dir),
        "cache": True,
        "timeout": 11,
    }


def test_download_pdbgz_subcommand(monkeypatch, cache_dir):
    out_dir = cache_dir / "cli_pdbs"
    calls = {}

    def fake_download_pdb_gz(pdb_id, out_dir, cache, timeout):
        calls["pdb_id"] = pdb_id
        calls["out_dir"] = out_dir
        calls["cache"] = cache
        calls["timeout"] = timeout
        p = Path(out_dir) / f"{pdb_id.upper()}.pdb.gz"
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(b"dummy")
        return p

    monkeypatch.setattr(cli, "download_pdb_gz", fake_download_pdb_gz)

    rc = cli.main(["download-pdbgz", "4KIQ", "--out", str(out_dir), "--no-cache"])
    assert rc == 0
    assert calls == {
        "pdb_id": "4KIQ",
        "out_dir": str(out_dir),
        "cache": False,
        "timeout": 20,
    }
