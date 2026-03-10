import gzip
from pathlib import Path

import pytest

from miscc4c import download_pdb
from miscc4c.download_pdb import download_pdb_gz


@pytest.fixture(scope="session")
def cache_dir():
    return Path("mytestcache")


def test_download_pdb_gz_saves_gzip(cache_dir):
    out_dir = cache_dir / "pdbs_real"
    path = download_pdb_gz("4KIQ", out_dir=out_dir, cache=False)

    assert path.exists()
    assert path.name == "4KIQ.pdb.gz"

    with gzip.open(path, "rt", encoding="utf-8") as f:
        head = f.read(200)

    assert "HEADER" in head or "ATOM" in head or "HETATM" in head


def test_download_pdb_gz_cache(cache_dir):
    out_dir = cache_dir / "pdbs_real"
    path1 = download_pdb_gz("4KIQ", out_dir=out_dir, cache=False)
    path2 = download_pdb_gz("4KIQ", out_dir=out_dir)

    assert path1 == path2
    assert path1.exists()


def test_download_pdb_gz_prefers_pdbj(monkeypatch, cache_dir):
    out_dir = cache_dir / "pdbs_mock"
    called_urls = []

    def fake_download(url: str, timeout: int = 20):
        called_urls.append(url)
        return b"dummy-gz-bytes"

    monkeypatch.setattr(download_pdb, "_download_bytes", fake_download)

    path = download_pdb_gz("4KIQ", out_dir=out_dir, cache=False)

    assert path.exists()
    assert called_urls == [download_pdb.PDBJ_ENTRY_PDB_GZ_URL.format(pdb_id="4KIQ")]


def test_download_pdb_gz_fallback_to_rcsb_when_pdbj_fails(monkeypatch, cache_dir):
    out_dir = cache_dir / "pdbs_mock"
    called_urls = []

    def fake_download(url: str, timeout: int = 20):
        called_urls.append(url)
        if "pdbj.org" in url:
            raise RuntimeError("pdbj failed")
        return b"dummy-gz-bytes"

    monkeypatch.setattr(download_pdb, "_download_bytes", fake_download)

    path = download_pdb_gz("4KIQ", out_dir=out_dir, cache=False)

    assert path.exists()
    assert called_urls == [
        download_pdb.PDBJ_ENTRY_PDB_GZ_URL.format(pdb_id="4KIQ"),
        download_pdb.RCSB_ENTRY_PDB_GZ_URL.format(pdb_id="4KIQ"),
    ]
