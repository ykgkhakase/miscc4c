from pathlib import Path

import pytest
from rdkit import Chem

from miscc4c.assign_bond_order import (
    add_neutral_hydrogens,
    assign_bond_orders_from_template,
    build_template_from_smiles,
    extract_ligand_pdb_block,
    extract_smiles_from_ccd_cif,
    mol_from_pdb_block,
)


TEST_CASES = [
    ("4KIQ", "1RA", "public_data/pdb/4KIQ.pdb.gz", "public_data/pdb/1RA.cif", "A", 401),
    ("2CVD", "HQL", "public_data/pdb/2CVD.pdb.gz", "public_data/pdb/HQL.cif", "A", 2201),
]

CASE_IDS = [f"{pdb_id}-{ligand_id}" for pdb_id, ligand_id, *_ in TEST_CASES]


@pytest.mark.parametrize(
    "pdb_id, ligand_id, pdb_rel_path, cif_rel_path, chain, resnum",
    TEST_CASES,
    ids=CASE_IDS,
)
def test_assign_bond_order_pipeline(
    pdb_id: str,
    ligand_id: str,
    pdb_rel_path: str,
    cif_rel_path: str,
    chain: str,
    resnum: int,
):
    repo_root = Path(__file__).resolve().parents[1]
    pdb_path = repo_root / pdb_rel_path
    cif_path = repo_root / cif_rel_path

    if not pdb_path.exists():
        pytest.skip(f"PDB file not found: {pdb_path}")
    if not cif_path.exists():
        pytest.skip(f"CIF file not found: {cif_path}")

    smiles = extract_smiles_from_ccd_cif(str(cif_path))
    assert isinstance(smiles, str)
    assert smiles

    pdb_block = extract_ligand_pdb_block(
        pdb_path=str(pdb_path),
        resname=ligand_id,
        chain=chain,
        resnum=resnum,
    )
    assert "HETATM" in pdb_block or "ATOM" in pdb_block

    pdb_mol = mol_from_pdb_block(pdb_block)
    template = build_template_from_smiles(smiles)

    assigned = assign_bond_orders_from_template(pdb_mol, template)
    assigned_h = add_neutral_hydrogens(assigned)

    assigned_no_h = Chem.RemoveHs(assigned)
    template_no_h = Chem.RemoveHs(template)

    pdb_smiles = Chem.MolToSmiles(pdb_mol, canonical=True)
    assigned_smiles = Chem.MolToSmiles(assigned_no_h, canonical=True)
    template_smiles = Chem.MolToSmiles(template_no_h, canonical=True)

    print(f"[{pdb_id}-{ligand_id}] PDBbased SMILES: {pdb_smiles}")
    print(f"[{pdb_id}-{ligand_id}] assigned SMILES: {assigned_smiles}")
    print(f"[{pdb_id}-{ligand_id}] template SMILES: {template_smiles}")
    print(
        f"[{pdb_id}-{ligand_id}] atoms/bonds assigned={assigned_no_h.GetNumAtoms()}/{assigned_no_h.GetNumBonds()} "
        f"template={template_no_h.GetNumAtoms()}/{template_no_h.GetNumBonds()}"
    )

    assert assigned_no_h.GetNumAtoms() == template_no_h.GetNumAtoms()
    assert assigned_no_h.GetNumBonds() == template_no_h.GetNumBonds()
    assert assigned_h.GetNumAtoms() >= assigned.GetNumAtoms()
