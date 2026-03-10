from __future__ import annotations

import gzip
import shlex
from typing import List, Optional

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

def read_text_auto(path: str) -> str:
    if path.endswith(".gz"):
        with gzip.open(path, "rt", encoding="utf-8") as f:
            return f.read()
    else:
        with open(path, "r", encoding="utf-8") as f:
            return f.read()


def extract_ligand_pdb_block(
    pdb_path: str,
    resname: str,
    chain: Optional[str] = None,
    resnum: Optional[int] = None,
) -> str:
    text = read_text_auto(pdb_path)
    selected_lines: List[str] = []

    for line in text.splitlines():
        record = line[:6].strip()
        if record not in {"ATOM", "HETATM"}:
            continue

        line_resname = line[17:20].strip()
        line_chain = line[21].strip()
        line_resnum_str = line[22:26].strip()
        line_altloc = line[16].strip()

        if line_resname != resname:
            continue
        if chain is not None and line_chain != chain:
            continue
        if resnum is not None:
            try:
                line_resnum = int(line_resnum_str)
            except ValueError:
                continue
            if line_resnum != resnum:
                continue

        # altLoc は空白または A のみ採用
        if line_altloc not in {"", "A"}:
            continue

        selected_lines.append(line)

    if not selected_lines:
        raise ValueError(
            f"リガンドが見つかりません: resname={resname}, chain={chain}, resnum={resnum}"
        )

    return "\n".join(selected_lines) + "\nEND\n"


def mol_from_pdb_block(pdb_block: str) -> Chem.Mol:
    mol = Chem.MolFromPDBBlock(
        pdb_block,
        sanitize=False,
        removeHs=False,
        proximityBonding=True,
    )
    if mol is None:
        raise ValueError("PDB block から分子を作成できませんでした。")
    return mol


def extract_smiles_from_ccd_cif(cif_path: str, prefer_canonical: bool = True) -> str:
    """
    wwPDB chemical component CIF から SMILES を抽出する。
    優先順位:
      1. OpenEye SMILES_CANONICAL
      2. CACTVS SMILES_CANONICAL
      3. OpenEye SMILES
      4. CACTVS SMILES
      5. その他の SMILES
    """
    text = read_text_auto(cif_path)
    lines = text.splitlines()

    in_loop = False
    headers = []
    data_rows = []

    for i, line in enumerate(lines):
        s = line.strip()
        if not s:
            continue

        if s == "loop_":
            in_loop = True
            headers = []
            data_rows = []
            continue

        if in_loop and s.startswith("_"):
            headers.append(s)
            continue

        if in_loop and headers:
            # 必要な descriptor loop か確認
            if (
                "_pdbx_chem_comp_descriptor.comp_id" in headers
                and "_pdbx_chem_comp_descriptor.type" in headers
                and "_pdbx_chem_comp_descriptor.program" in headers
                and "_pdbx_chem_comp_descriptor.descriptor" in headers
            ):
                if s.startswith("#"):
                    break
                data_rows.append(shlex.split(s))
            elif s.startswith("#"):
                in_loop = False
                headers = []
                data_rows = []

    if not headers or not data_rows:
        raise ValueError("CIF 中に descriptor loop が見つかりませんでした。")

    col_idx = {h: i for i, h in enumerate(headers)}

    candidates = []
    for row in data_rows:
        try:
            dtype = row[col_idx["_pdbx_chem_comp_descriptor.type"]]
            program = row[col_idx["_pdbx_chem_comp_descriptor.program"]]
            desc = row[col_idx["_pdbx_chem_comp_descriptor.descriptor"]]
        except (IndexError, KeyError):
            continue

        if "SMILES" not in dtype:
            continue

        score = 0
        if prefer_canonical and dtype == "SMILES_CANONICAL":
            score += 100
        if "OpenEye" in program:
            score += 20
        elif "CACTVS" in program:
            score += 10
        elif dtype == "SMILES":
            score += 1

        candidates.append((score, desc, dtype, program))

    if not candidates:
        raise ValueError("CIF 中に SMILES descriptor が見つかりませんでした。")

    candidates.sort(reverse=True, key=lambda x: x[0])
    return candidates[0][1]


def neutralize_molecule(mol: Chem.Mol) -> Chem.Mol:
    """
    簡便な中性化。
    厳密なプロトン化状態決定ではなく、初期化用の処理とみなす。
    """
    pattern_replacements = (
        ("[n+;H]", "n"),
        ("[N+;!H0]", "N"),
        ("[$([O-]);!$([O-][#7])]", "O"),
        ("[S-;X1]", "S"),
        ("[$([N-;X2]S(=O)=O)]", "N"),
        ("[$([N-;X2][C,N]=C)]", "N"),
        ("[n-]", "[nH]"),
        ("[$([S-]=O)]", "S"),
        ("[$([N-]C=O)]", "N"),
    )

    out = Chem.Mol(mol)
    for smarts, repl in pattern_replacements:
        patt = Chem.MolFromSmarts(smarts)
        repl_mol = Chem.MolFromSmiles(repl)
        while out.HasSubstructMatch(patt):
            rms = AllChem.ReplaceSubstructs(out, patt, repl_mol, replaceAll=False)
            out = rms[0]

    Chem.SanitizeMol(out)
    return out


def build_template_from_smiles(smiles: str) -> Chem.Mol:
    template = Chem.MolFromSmiles(smiles)
    if template is None:
        raise ValueError(f"SMILES を解釈できませんでした: {smiles}")
    template = neutralize_molecule(template)
    Chem.SanitizeMol(template)
    return template


def assign_bond_orders_from_template(
    pdb_mol: Chem.Mol,
    template_mol: Chem.Mol,
) -> Chem.Mol:
    pdb_no_h = Chem.RemoveHs(pdb_mol)
    template_no_h = Chem.RemoveHs(template_mol)

    assigned = AllChem.AssignBondOrdersFromTemplate(template_no_h, pdb_no_h)
    Chem.SanitizeMol(assigned)
    return assigned


def add_neutral_hydrogens(mol: Chem.Mol) -> Chem.Mol:
    Chem.SanitizeMol(mol)
    mol_h = Chem.AddHs(mol, addCoords=True)
    return mol_h
