import pytest
from pathlib import Path
from miscc4c import download_cif

# 将来的にここへケースを追加していく
TEST_CASES = [
    ("4KIQ", "1RA"),
    ("2CVD", "HQL")
    # ("1ABC", "LIG"),
    # ("2XYZ", "ATP"),
]

CASE_IDS = [f"{pdb_id}-{ligand_id}" for pdb_id, ligand_id in TEST_CASES]

#@pytest.fixture
#def cache_dir(tmp_path):
#    """pytest の一時ディレクトリをキャッシュとして使用"""
#    return tmp_path

@pytest.fixture(scope="session")
def cache_dir():
    return Path("mytestcache")


@pytest.mark.parametrize("pdb_id, ligand_id", TEST_CASES, ids=CASE_IDS)
def test_extract_ligand_ccd_ids_from_pdb(pdb_id, ligand_id, cache_dir):
    """PDB ID から CCD ID を抽出できる"""
    ccd_ids = download_cif.extract_ligand_ccd_ids_from_pdb_id(
        pdb_id,
        cache_dir=cache_dir / "entries",
    )

    assert isinstance(ccd_ids, list)
    assert ligand_id in ccd_ids


@pytest.mark.parametrize("pdb_id, ligand_id", TEST_CASES, ids=CASE_IDS)
def test_download_ligand_cif_from_expected_ligand_id(pdb_id, ligand_id, cache_dir):
    """期待される CCD ID の CIF を取得できる"""
    path = download_cif.download_ligand_cif(
        ligand_id,
        out_dir=cache_dir / "ligands",
    )

    assert path.exists()

    text = path.read_text(encoding="utf-8")

    assert "data_" in text
    assert "_chem_comp.id" in text
    assert ligand_id in text


@pytest.mark.parametrize("pdb_id, ligand_id", TEST_CASES, ids=CASE_IDS)
def test_cache_functionality(pdb_id, ligand_id, cache_dir):
    """キャッシュが機能する"""
    out_dir = cache_dir / "ligands"

    path1 = download_cif.download_ligand_cif(
        ligand_id,
        out_dir=out_dir,
    )

    path2 = download_cif.download_ligand_cif(
        ligand_id,
        out_dir=out_dir,
    )

    assert path1 == path2
    assert path1.exists()


@pytest.mark.parametrize("pdb_id, ligand_id", TEST_CASES)
def test_download_all_ligands_for_pdb(pdb_id, ligand_id, cache_dir):
    """PDB ID から全 ligand を取得できる"""
    results = download_cif.download_all_ligand_cifs_for_pdb(
        pdb_id,
        ligand_out_dir=cache_dir / "ligands",
        entry_cache_dir=cache_dir / "entries",
    )

    assert isinstance(results, dict)
    assert ligand_id in results
    assert results[ligand_id] is not None
    assert results[ligand_id].exists()