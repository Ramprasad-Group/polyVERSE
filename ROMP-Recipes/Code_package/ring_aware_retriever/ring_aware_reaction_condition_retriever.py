#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import json
import math
import re
from collections import Counter
from pathlib import Path

import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
FG_SMARTS = {
    'alkene': '[CX3]=[CX3]',
    'alkyne': '[CX2]#[CX2]',
    'alcohol': '[OX2H]',
    'phenol': 'c[OX2H]',
    'amine': '[NX3;H2,H1;!$(NC=O)]',
    'amide': '[NX3][CX3](=[OX1])',
    'ester': '[CX3](=[OX1])[OX2][#6]',
    'ether': '[OD2]([#6])[#6]',
    'thioether': '[SD2]([#6])[#6]',
    'thiol': '[SX2H]',
    'carbonyl': '[CX3]=[OX1]',
    'carboxylic_acid': '[CX3](=O)[OX2H1]',
    'nitrile': '[CX2]#N',
    'halide': '[F,Cl,Br,I]',
    'sulfone': '[SX4](=[OX1])(=[OX1])([#6])[#6]',
    'sulfoxide': '[SX3](=[OX1])([#6])[#6]'
}
ROMP_ALKENE_SMARTS = '[C;R]=[C;R]'
GENERIC_ALKENE_SMARTS = '[CX3]=[CX3]'


def load_table(path: str) -> pd.DataFrame:
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix == '.csv':
        return pd.read_csv(path)
    if suffix in ['.xlsx', '.xls']:
        return pd.read_excel(path)
    if suffix == '.parquet':
        return pd.read_parquet(path)
    raise ValueError(f'Unsupported file format: {suffix}')


def normalize_text(x):
    if pd.isna(x):
        return ''
    return str(x).strip()


def to_float_safe(x):
    try:
        return float(x)
    except Exception:
        return math.nan


def mol_from_smiles(smiles: str):
    if not smiles or pd.isna(smiles):
        return None
    try:
        return Chem.MolFromSmiles(str(smiles))
    except Exception:
        return None


_MORGAN_GENERATORS = {}


def morgan_fp_from_smiles(smiles: str, radius=2, nbits=2048):
    mol = mol_from_smiles(smiles)
    if mol is None:
        return None
    key = (int(radius), int(nbits))
    generator = _MORGAN_GENERATORS.get(key)
    if generator is None:
        generator = rdFingerprintGenerator.GetMorganGenerator(radius=key[0], fpSize=key[1])
        _MORGAN_GENERATORS[key] = generator
    return generator.GetFingerprint(mol)


def tanimoto(fp1, fp2):
    if fp1 is None or fp2 is None:
        return None
    return DataStructs.TanimotoSimilarity(fp1, fp2)


def has_functional_group(smiles: str, fg_name: str) -> bool:
    mol = mol_from_smiles(smiles)
    if mol is None:
        return False
    smarts = FG_SMARTS.get(fg_name)
    if smarts is None:
        raise ValueError(f'Unknown functional group: {fg_name}')
    patt = Chem.MolFromSmarts(smarts)
    if patt is None:
        return False
    return mol.HasSubstructMatch(patt)


def infer_family_from_mechanism(mechanism_text: str) -> str:
    s = normalize_text(mechanism_text).lower()
    if not s:
        return ''
    if s == 'romp':
        return 'ROMP'
    if 'rop' in s:
        return 'ROP'
    if 'radical' in s:
        return 'Radical'
    return mechanism_text


def first_existing(df: pd.DataFrame, candidates, default=None):
    lower_map = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand in df.columns:
            return cand
        if cand.lower() in lower_map:
            return lower_map[cand.lower()]
    return default


def smart_split(text: str):
    s = normalize_text(text)
    if not s:
        return []
    s = s.replace('\n', ' ').replace('\r', ' ')
    s = re.sub('\\s+', ' ', s).strip()
    for sep in [';', '/', '|']:
        parts = [p.strip() for p in s.split(sep)]
        if len(parts) > 1:
            return parts
    return [s]


def choose_best_split(values_text: str, categories_text: str):
    cats = smart_split(categories_text)
    raw = normalize_text(values_text).replace('\n', ' ').replace('\r', ' ')
    raw = re.sub('\\s+', ' ', raw).strip()
    if not raw:
        return (cats, [])
    candidate_splits = []
    for sep in [';', '/', '|']:
        parts = [p.strip() for p in raw.split(sep)]
        candidate_splits.append((sep, parts))
    if cats:
        exact = [parts for _, parts in candidate_splits if len(parts) == len(cats)]
        if exact:
            return (cats, exact[0])
        best = min(candidate_splits, key=lambda x: abs(len(x[1]) - max(len(cats), 1)))
        return (cats, best[1])
    return (cats, smart_split(values_text))


def parse_ingredients_row(ingredients, categories):
    cats, vals = choose_best_split(ingredients, categories)
    parsed = {}
    n = min(len(cats), len(vals))
    for i in range(n):
        c = normalize_text(cats[i]).lower()
        v = normalize_text(vals[i])
        if not c:
            continue
        parsed[c] = v
    return parsed


def extract_temperature(texts):
    patterns = [
        '(-?\\d+(?:\\.\\d+)?)\\s*(?:°C|ºC|oC)',
        '(-?\\d+(?:\\.\\d+)?)\\s*(?:deg\\s*C)',
        '(?:at\\s*)?room temperature',
        '\\bR\\.T\\.\\b',
        '\\bRT\\b'
    ]
    hits = []
    joined = ' ; '.join([normalize_text(t) for t in texts if normalize_text(t)])
    if not joined:
        return ('', math.nan)
    for m in re.finditer(patterns[0], joined, flags=re.I):
        hits.append(float(m.group(1)))
    for m in re.finditer(patterns[1], joined, flags=re.I):
        hits.append(float(m.group(1)))
    if hits:
        return (
            '; '.join((str(x).rstrip('0').rstrip('.') if isinstance(
                x,
                float
            ) else str(x) for x in hits)),
            hits[0]
        )
    if re.search(
        patterns[2],
        joined,
        flags=re.I
    ) or re.search(patterns[3], joined, flags=re.I) or re.search(patterns[4], joined, flags=re.I):
        return ('R.T.', 25.0)
    return ('', math.nan)


def find_reactive_alkene_bond(mol):
    if mol is None:
        return None
    patt_ring = Chem.MolFromSmarts(ROMP_ALKENE_SMARTS)
    if patt_ring is not None:
        matches = mol.GetSubstructMatches(patt_ring)
        if matches:
            for a1, a2 in matches:
                bond = mol.GetBondBetweenAtoms(a1, a2)
                if bond is not None and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    return bond.GetIdx()
    patt_generic = Chem.MolFromSmarts(GENERIC_ALKENE_SMARTS)
    if patt_generic is not None:
        matches = mol.GetSubstructMatches(patt_generic)
        if matches:
            for a1, a2 in matches:
                bond = mol.GetBondBetweenAtoms(a1, a2)
                if bond is not None and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    return bond.GetIdx()
    return None


def get_atoms_within_hops(mol, center_atom_ids, max_hops=2):
    visited = set(center_atom_ids)
    frontier = set(center_atom_ids)
    for _ in range(max_hops):
        next_frontier = set()
        for aid in frontier:
            atom = mol.GetAtomWithIdx(aid)
            for nbr in atom.GetNeighbors():
                nid = nbr.GetIdx()
                if nid not in visited:
                    visited.add(nid)
                    next_frontier.add(nid)
        frontier = next_frontier
        if not frontier:
            break
    return sorted(visited)


def build_local_env_smiles(smiles, hops=2):
    mol = mol_from_smiles(smiles)
    if mol is None:
        return None
    bond_idx = find_reactive_alkene_bond(mol)
    if bond_idx is None:
        return None
    bond = mol.GetBondWithIdx(bond_idx)
    a1 = bond.GetBeginAtomIdx()
    a2 = bond.GetEndAtomIdx()
    atom_ids = get_atoms_within_hops(mol, [a1, a2], max_hops=hops)
    try:
        env = Chem.MolFragmentToSmiles(
            mol,
            atomsToUse=atom_ids,
            canonical=True,
            isomericSmiles=True
        )
        return env
    except Exception:
        return None


def local_env_fp_from_smiles(smiles, hops=2, radius=2, nbits=2048):
    local_smiles = build_local_env_smiles(smiles, hops=hops)
    if not local_smiles:
        return None
    return morgan_fp_from_smiles(local_smiles, radius=radius, nbits=nbits)


def classify_ring_topology(mol, target_ring_atoms):
    ring_info = mol.GetRingInfo()
    atom_rings = [set(r) for r in ring_info.AtomRings()]
    target_set = set(target_ring_atoms)
    overlaps = []
    for r in atom_rings:
        if r == target_set:
            continue
        inter = len(target_set & r)
        if inter > 0:
            overlaps.append(inter)
    if not overlaps:
        return 'simple'
    n_ge2 = sum((x >= 2 for x in overlaps))
    n_eq1 = sum((x == 1 for x in overlaps))
    if n_ge2 == 1 and n_eq1 == 0:
        return 'fused'
    if n_ge2 == 0 and n_eq1 >= 1:
        return 'spiro'
    if n_ge2 >= 2 or (n_ge2 >= 1 and n_eq1 >= 1):
        return 'bridged'
    return 'simple'


def get_ring_info_for_bond(mol, bond_idx):
    if mol is None or bond_idx is None:
        return {
            'ring_atom_ids': None,
            'ring_sizes': [],
            'min_ring_size': None,
            'ring_topology': 'none'
        }
    ring_info = mol.GetRingInfo()
    atom_rings = list(ring_info.AtomRings())
    bond_rings = list(ring_info.BondRings())
    matched_ring_atom_ids = []
    ring_sizes = []
    for a_ring, b_ring in zip(atom_rings, bond_rings):
        if bond_idx in b_ring:
            matched_ring_atom_ids.append(tuple(sorted(a_ring)))
            ring_sizes.append(len(a_ring))
    if not matched_ring_atom_ids:
        return {
            'ring_atom_ids': None,
            'ring_sizes': [],
            'min_ring_size': None,
            'ring_topology': 'none'
        }
    min_size = min(ring_sizes)
    target_ring = None
    for ring_atoms, size in zip(matched_ring_atom_ids, ring_sizes):
        if size == min_size:
            target_ring = ring_atoms
            break
    topology = classify_ring_topology(mol, target_ring)
    return {
        'ring_atom_ids': list(target_ring),
        'ring_sizes': ring_sizes,
        'min_ring_size': min_size,
        'ring_topology': topology
    }


def build_reactive_ring_smiles(smiles):
    mol = mol_from_smiles(smiles)
    if mol is None:
        return None
    bond_idx = find_reactive_alkene_bond(mol)
    if bond_idx is None:
        return None
    ring_data = get_ring_info_for_bond(mol, bond_idx)
    atom_ids = ring_data['ring_atom_ids']
    if not atom_ids:
        return None
    try:
        frag = Chem.MolFragmentToSmiles(
            mol,
            atomsToUse=atom_ids,
            canonical=True,
            isomericSmiles=True
        )
        return frag
    except Exception:
        return None


def reactive_ring_fp_from_smiles(smiles, radius=2, nbits=2048):
    ring_smiles = build_reactive_ring_smiles(smiles)
    if not ring_smiles:
        return None
    return morgan_fp_from_smiles(ring_smiles, radius=radius, nbits=nbits)


def ring_size_score(size1, size2):
    if size1 is None or size2 is None:
        return 0.0
    diff = abs(size1 - size2)
    if diff == 0:
        return 1.0
    elif diff == 1:
        return 0.6
    elif diff == 2:
        return 0.3
    else:
        return 0.0


def ring_topology_score(t1, t2):
    if not t1 or not t2 or t1 == 'none' or (t2 == 'none'):
        return 0.0
    if t1 == t2:
        return 1.0
    close_pairs = {('fused', 'bridged'), ('bridged', 'fused')}
    if (t1, t2) in close_pairs:
        return 0.4
    return 0.0


def ring_meta_score(query_ring_size, cand_ring_size, query_topology, cand_topology):
    s1 = ring_size_score(query_ring_size, cand_ring_size)
    s2 = ring_topology_score(query_topology, cand_topology)
    return 0.5 * s1 + 0.5 * s2


def combined_similarity_advanced(
    global_fp1,
    global_fp2,
    local_fp1,
    local_fp2,
    ring_fp1,
    ring_fp2,
    ring_size1,
    ring_size2,
    ring_topology1,
    ring_topology2,
    w_global=0.2,
    w_local=0.3,
    w_ring=0.3,
    w_ring_meta=0.2


):
    g = tanimoto(global_fp1, global_fp2)
    l = tanimoto(local_fp1, local_fp2)
    r = tanimoto(ring_fp1, ring_fp2)
    rm = ring_meta_score(ring_size1, ring_size2, ring_topology1, ring_topology2)
    parts = []
    weights = []
    if g is not None:
        parts.append(g)
        weights.append(w_global)
    if l is not None:
        parts.append(l)
        weights.append(w_local)
    if r is not None:
        parts.append(r)
        weights.append(w_ring)
    if ring_size1 is not None and ring_size2 is not None:
        parts.append(rm)
        weights.append(w_ring_meta)
    if not parts:
        return None
    weighted_sum = sum((p * w for p, w in zip(parts, weights)))
    total_w = sum(weights)
    return weighted_sum / total_w if total_w > 0 else None


def canonicalize_input_table(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    smiles_col = first_existing(
        out,
        ['monomer_smiles', 'Monomer_smiles', 'Monomer_SMILES', 'SMILES']
    )
    if smiles_col is None:
        raise ValueError('No monomer SMILES column was found (for example, monomer_smiles or Monomer_smiles).')
    if smiles_col != 'monomer_smiles':
        out['monomer_smiles'] = out[smiles_col]
    poly_col = first_existing(out, ['polymer_psmiles', 'Polymer_smiles', 'polymer_smiles'])
    if poly_col and poly_col != 'polymer_psmiles':
        out['polymer_psmiles'] = out[poly_col]
    summary_sources = [c for c in [
        first_existing(out, ['text_summary']),
        first_existing(out, ['Procedure']),
        first_existing(out, ['Conditions']),
        first_existing(out, ['Yield'])
    ] if c]
    if 'text_summary' not in out.columns:
        if summary_sources:
            out['text_summary'] = out[summary_sources].fillna('').astype(str).agg(
                ' | '.join,
                axis=1
            )
        else:
            out['text_summary'] = ''
    mech_col = first_existing(out, ['mechanism', 'Mechanism'])
    if mech_col is None:
        out['mechanism'] = 'ROMP'
    elif mech_col != 'mechanism':
        out['mechanism'] = out[mech_col]
    fam_col = first_existing(out, ['family', 'Family'])
    if fam_col is None:
        out['family'] = out['mechanism'].apply(infer_family_from_mechanism)
    elif fam_col != 'family':
        out['family'] = out[fam_col]
    for src, dst in [
        (first_existing(out, ['DOI']), 'DOI'),
        (first_existing(out, ['Monomer_name']), 'Monomer_name'),
        (first_existing(out, ['Polymer_name']), 'Polymer_name'),
        (first_existing(out, ['Ingredients']), 'Ingredients'),
        (first_existing(out, ['Ingredients_category']), 'Ingredients_category'),
        (first_existing(out, ['Ingredients_quantity']), 'Ingredients_quantity'),
        (first_existing(out, ['Conditions']), 'Conditions'),
        (first_existing(out, ['Yield']), 'Yield'),
        (first_existing(out, ['Procedure']), 'Procedure')
    ]:
        if src and src != dst:
            out[dst] = out[src]
    parsed_rows = []
    for _, row in out.iterrows():
        parsed = parse_ingredients_row(
            row.get('Ingredients', ''),
            row.get('Ingredients_category', '')
        )
        parsed_rows.append(parsed)
    parsed_df = pd.DataFrame(parsed_rows)
    for key in parsed_df.columns:
        colname = f'_parsed_{key}'
        out[colname] = parsed_df[key]

    def prefer_existing_or_parsed(
        target,
        candidates_existing=None,
        candidates_parsed=None,
        fallback=''
    ):
        if candidates_existing is None:
            candidates_existing = []
        if candidates_parsed is None:
            candidates_parsed = []
        if target in out.columns:
            return
        for c in candidates_existing:
            src = first_existing(out, [c])
            if src is not None:
                out[target] = out[src]
                return
        for c in candidates_parsed:
            pcol = f'_parsed_{c.lower()}'
            if pcol in out.columns:
                out[target] = out[pcol]
                return
        out[target] = fallback
    prefer_existing_or_parsed('Catalyst', ['Catalyst'], ['catalyst', 'initiator'])
    prefer_existing_or_parsed('Initiator', ['Initiator'], ['initiator', 'catalyst'])
    prefer_existing_or_parsed('Solvent', ['Solvent'], ['solvent'])
    prefer_existing_or_parsed('Atmosphere', ['Atmosphere'], ['atmosphere'])
    prefer_existing_or_parsed('Quencher', ['Quencher'], ['quencher'])
    out['Solvent'] = out['Solvent'].apply(normalize_text)
    if 'Temperature' not in out.columns:
        temp_pairs = out.apply(
            lambda r: extract_temperature([
                r.get('Conditions', ''),
                r.get('Yield', ''),
                r.get('Procedure', '')
            ]),
            axis=1
        )
        out['Temperature'] = [tp[0] for tp in temp_pairs]
        out['_temp_numeric_from_template'] = [tp[1] for tp in temp_pairs]
    if 'n_rows' not in out.columns:
        out['n_rows'] = 1
    if 'index_list' not in out.columns:
        out['index_list'] = out.index.astype(str)
    drop_cols = [c for c in out.columns if str(c).startswith('Unnamed:')]
    if drop_cols:
        out = out.drop(columns=drop_cols)
    return out


class ReactionRetriever:

    def __init__(
        self,
        df: pd.DataFrame,
        smiles_col='monomer_smiles',
        family_col='family',
        mechanism_col='mechanism',
        catalyst_col='Catalyst',
        initiator_col='Initiator',
        solvent_col='Solvent',
        temp_col='Temperature',
        summary_col='text_summary',
        polymer_col='polymer_psmiles',
        atmosphere_col='Atmosphere',
        quencher_col='Quencher',
        nrows_col='n_rows',
        index_list_col='index_list',
        radius=2,
        nbits=2048,
        local_hops=2,
        w_global=0.15,
        w_local=0.3,
        w_ring=0.35,
        w_ring_meta=0.2
    ):
        self.df = canonicalize_input_table(df)
        self.smiles_col = smiles_col
        self.family_col = family_col
        self.mechanism_col = mechanism_col
        self.catalyst_col = catalyst_col
        self.initiator_col = initiator_col
        self.solvent_col = solvent_col
        self.temp_col = temp_col
        self.summary_col = summary_col
        self.polymer_col = polymer_col
        self.atmosphere_col = atmosphere_col
        self.quencher_col = quencher_col
        self.nrows_col = nrows_col
        self.index_list_col = index_list_col
        self.radius = radius
        self.nbits = nbits
        self.local_hops = local_hops
        self.w_global = w_global
        self.w_local = w_local
        self.w_ring = float(w_ring)
        self.w_ring_meta = float(w_ring_meta)
        self._validate_similarity_weights()
        self._prepare()

    def _validate_similarity_weights(self):
        weights = [self.w_global, self.w_local, self.w_ring, self.w_ring_meta]
        if any((w < 0 for w in weights)):
            raise ValueError('Similarity weights must be non-negative.')
        if not math.isclose(sum(weights), 1.0, rel_tol=0.0, abs_tol=1e-08):
            raise ValueError(f'Similarity weights must sum to 1.0; received {sum(weights):.8f}.')

    def _extract_ring_meta(self, smiles):
        mol = mol_from_smiles(smiles)
        if mol is None:
            return {'ring_size': None, 'ring_topology': 'none'}
        bond_idx = find_reactive_alkene_bond(mol)
        if bond_idx is None:
            return {'ring_size': None, 'ring_topology': 'none'}
        info = get_ring_info_for_bond(mol, bond_idx)
        return {'ring_size': info['min_ring_size'], 'ring_topology': info['ring_topology']}

    def _prepare(self):
        if self.family_col not in self.df.columns:
            if self.mechanism_col in self.df.columns:
                self.df[self.family_col] = self.df[self.mechanism_col].apply(infer_family_from_mechanism)
            else:
                self.df[self.family_col] = ''
        text_cols = [
            self.smiles_col,
            self.family_col,
            self.mechanism_col,
            self.catalyst_col,
            self.initiator_col,
            self.solvent_col,
            self.summary_col,
            self.polymer_col,
            self.atmosphere_col,
            self.quencher_col,
            self.index_list_col,
            'DOI',
            'Monomer_name',
            'Polymer_name',
            'Conditions',
            'Yield',
            'Procedure',
            'Ingredients',
            'Ingredients_category'
        ]
        for col in text_cols:
            if col in self.df.columns:
                self.df[col] = self.df[col].apply(normalize_text)
        if self.nrows_col in self.df.columns:
            self.df[self.nrows_col] = pd.to_numeric(self.df[self.nrows_col], errors='coerce')
        if self.temp_col in self.df.columns:
            self.df['_temp_text'] = self.df[self.temp_col].apply(normalize_text)
            self.df['_temp_numeric'] = self.df.get(
                '_temp_numeric_from_template',
                self.df[self.temp_col].apply(to_float_safe)
            )
        else:
            self.df['_temp_text'] = ''
            self.df['_temp_numeric'] = math.nan
        if self.smiles_col not in self.df.columns:
            raise ValueError(f'Required column is missing: {self.smiles_col}')
        self.df['_mol'] = self.df[self.smiles_col].apply(mol_from_smiles)
        self.df['_fp'] = self.df[self.smiles_col].apply(lambda s: morgan_fp_from_smiles(
            s,
            radius=self.radius,
            nbits=self.nbits
        ))
        self.df['_local_env_smiles'] = self.df[self.smiles_col].apply(lambda s: build_local_env_smiles(
            s,
            hops=self.local_hops
        ))
        self.df['_local_fp'] = self.df[self.smiles_col].apply(lambda s: local_env_fp_from_smiles(
            s,
            hops=self.local_hops,
            radius=self.radius,
            nbits=self.nbits
        ))
        self.df['_ring_smiles'] = self.df[self.smiles_col].apply(build_reactive_ring_smiles)
        self.df['_ring_fp'] = self.df[self.smiles_col].apply(lambda s: reactive_ring_fp_from_smiles(
            s,
            radius=self.radius,
            nbits=self.nbits
        ))
        ring_meta = self.df[self.smiles_col].apply(self._extract_ring_meta)
        self.df['_ring_size'] = [x['ring_size'] for x in ring_meta]
        self.df['_ring_topology'] = [x['ring_topology'] for x in ring_meta]
        self.df['_valid_smiles'] = self.df['_fp'].notna()
        self.df['_has_local_alkene_env'] = self.df['_local_fp'].notna()
        self.df['_has_ring_alkene'] = self.df['_ring_size'].notna()
        for fg in FG_SMARTS:
            self.df[f'_fg_{fg}'] = self.df[self.smiles_col].apply(lambda s: has_functional_group(
                s,
                fg
            ))

    def filter_candidates(
        self,
        family=None,
        mechanism=None,
        required_fgs=None,
        require_ring_alkene=False
    ):
        sub = self.df.copy()
        if family and self.family_col in sub.columns:
            sub = sub[sub[self.family_col].str.lower() == family.strip().lower()]
        if mechanism and self.mechanism_col in sub.columns:
            sub = sub[sub[self.mechanism_col].str.lower() == mechanism.strip().lower()]
        if required_fgs:
            for fg in required_fgs:
                col = f'_fg_{fg}'
                if col not in sub.columns:
                    raise ValueError(f'Unknown functional group: {fg}')
                sub = sub[sub[col] == True]
        sub = sub[sub['_valid_smiles'] == True]
        if require_ring_alkene:
            sub = sub[sub['_has_ring_alkene'] == True]
        return sub

    def retrieve(
        self,
        query_smiles: str,
        topk=5,
        family=None,
        mechanism=None,
        required_fgs=None,
        min_similarity=0.6,
        require_ring_alkene=False
    ):
        qfp = morgan_fp_from_smiles(query_smiles, radius=self.radius, nbits=self.nbits)
        if qfp is None:
            raise ValueError('query_smiles could not be parsed as a valid molecule.')
        q_local_fp = local_env_fp_from_smiles(
            query_smiles,
            hops=self.local_hops,
            radius=self.radius,
            nbits=self.nbits
        )
        q_ring_fp = reactive_ring_fp_from_smiles(
            query_smiles,
            radius=self.radius,
            nbits=self.nbits
        )
        q_mol = mol_from_smiles(query_smiles)
        q_bond_idx = find_reactive_alkene_bond(q_mol)
        if q_bond_idx is not None:
            q_ring_info = get_ring_info_for_bond(q_mol, q_bond_idx)
            q_ring_size = q_ring_info['min_ring_size']
            q_ring_topology = q_ring_info['ring_topology']
        else:
            q_ring_size = None
            q_ring_topology = 'none'
        if require_ring_alkene and q_ring_size is None:
            raise ValueError('The query does not contain a ring-contained alkene, but require_ring_alkene is enabled.')
        sub = self.filter_candidates(
            family=family,
            mechanism=mechanism,
            required_fgs=required_fgs,
            require_ring_alkene=require_ring_alkene
        ).copy()
        if len(sub) == 0:
            return {
                'query': {
                    'query_smiles': query_smiles,
                    'family_filter': family,
                    'mechanism_filter': mechanism,
                    'required_fgs': required_fgs or [],
                    'require_ring_alkene': require_ring_alkene
                },
                'top_cases': [],
                'summary': {
                    'n_candidates_after_filter': 0,
                    'message': 'No candidate reactions remained after filtering.'
                }
            }
        sub['_global_similarity'] = sub['_fp'].apply(lambda fp: tanimoto(qfp, fp))
        sub['_local_similarity'] = sub['_local_fp'].apply(lambda fp: tanimoto(q_local_fp, fp))
        sub['_ring_similarity'] = sub['_ring_fp'].apply(lambda fp: tanimoto(q_ring_fp, fp))
        sub['_ring_meta_score'] = sub.apply(
            lambda r: ring_meta_score(
                q_ring_size,
                r['_ring_size'],
                q_ring_topology,
                r['_ring_topology']
            ),
            axis=1
        )
        sub['_similarity'] = sub.apply(
            lambda r: combined_similarity_advanced(
                qfp,
                r['_fp'],
                q_local_fp,
                r['_local_fp'],
                q_ring_fp,
                r['_ring_fp'],
                q_ring_size,
                r['_ring_size'],
                q_ring_topology,
                r['_ring_topology'],
                w_global=self.w_global,
                w_local=self.w_local,
                w_ring=self.w_ring,
                w_ring_meta=self.w_ring_meta
            ),
            axis=1
        )
        sub = sub[sub['_similarity'].notna()]
        sub = sub[sub['_similarity'] >= min_similarity]
        sub = sub.sort_values('_similarity', ascending=False).head(topk).copy()
        summary = self._make_summary(sub)
        summary['minimum_similarity'] = float(min_similarity)
        if sub.empty:
            summary['message'] = 'No precedent met the minimum weighted-similarity threshold.'
        return {
            'query': {
                'query_smiles': query_smiles,
                'family_filter': family,
                'mechanism_filter': mechanism,
                'required_fgs': required_fgs or [],
                'require_ring_alkene': require_ring_alkene,
                'minimum_similarity': float(min_similarity),
                'local_hops': self.local_hops,
                'w_global': self.w_global,
                'w_local': self.w_local,
                'w_ring': self.w_ring,
                'w_ring_meta': self.w_ring_meta,
                'query_ring_size': q_ring_size,
                'query_ring_topology': q_ring_topology
            },
            'top_cases': self._format_top_cases(sub),
            'summary': summary
        }

    def _format_top_cases(self, sub: pd.DataFrame):
        cols_keep = [
            'DOI',
            'Monomer_name',
            'Polymer_name',
            self.polymer_col,
            self.smiles_col,
            self.family_col,
            self.mechanism_col,
            self.atmosphere_col,
            self.quencher_col,
            self.catalyst_col,
            self.initiator_col,
            self.solvent_col,
            self.temp_col,
            self.nrows_col,
            self.index_list_col,
            'Conditions',
            'Yield',
            self.summary_col,
            'Procedure',
            'Ingredients',
            'Ingredients_category',
            '_local_env_smiles',
            '_ring_smiles',
            '_ring_size',
            '_ring_topology',
            '_global_similarity',
            '_local_similarity',
            '_ring_similarity',
            '_ring_meta_score',
            '_similarity'
        ]
        cols_keep = [c for c in cols_keep if c in sub.columns]
        records = []
        for i, row in sub[cols_keep].iterrows():
            rec = {
                'row_id': int(i),
                'global_similarity': round(
                    float(row['_global_similarity']),
                    4
                ) if pd.notna(row['_global_similarity']) else None,
                'local_similarity': round(
                    float(row['_local_similarity']),
                    4
                ) if pd.notna(row['_local_similarity']) else None,
                'ring_similarity': round(
                    float(row['_ring_similarity']),
                    4
                ) if pd.notna(row['_ring_similarity']) else None,
                'ring_meta_score': round(
                    float(row['_ring_meta_score']),
                    4
                ) if pd.notna(row['_ring_meta_score']) else None,
                'similarity': round(
                    float(row['_similarity']),
                    4
                ) if pd.notna(row['_similarity']) else None
            }
            for col in cols_keep:
                if col in [
                    '_global_similarity',
                    '_local_similarity',
                    '_ring_similarity',
                    '_ring_meta_score',
                    '_similarity'
                ]:
                    continue
                if col in row.index:
                    val = row[col]
                    rec[col] = None if pd.isna(val) else val
            records.append(rec)
        return records

    def _count_nonempty(self, series: pd.Series, topn=10):
        vals = []
        for x in series:
            x = normalize_text(x)
            if x:
                vals.append(x)
        return Counter(vals).most_common(topn)

    def _make_summary(self, sub: pd.DataFrame):
        family_counts = self._count_nonempty(sub[self.family_col]) if self.family_col in sub.columns else []
        mechanism_counts = self._count_nonempty(sub[self.mechanism_col]) if self.mechanism_col in sub.columns else []
        catalyst_counts = self._count_nonempty(sub[self.catalyst_col]) if self.catalyst_col in sub.columns else []
        initiator_counts = self._count_nonempty(sub[self.initiator_col]) if self.initiator_col in sub.columns else []
        solvent_counts = self._count_nonempty(sub[self.solvent_col]) if self.solvent_col in sub.columns else []
        atmosphere_counts = self._count_nonempty(sub[self.atmosphere_col]) if self.atmosphere_col in sub.columns else []
        quencher_counts = self._count_nonempty(sub[self.quencher_col]) if self.quencher_col in sub.columns else []
        temp_text_counts = self._count_nonempty(sub['_temp_text']) if '_temp_text' in sub.columns else []
        ring_topology_counts = self._count_nonempty(sub['_ring_topology']) if '_ring_topology' in sub.columns else []
        temp_stats = {}
        if '_temp_numeric' in sub.columns:
            temp_valid = pd.to_numeric(sub['_temp_numeric'], errors='coerce').dropna()
            if len(temp_valid) > 0:
                temp_stats = {
                    'count': int(len(temp_valid)),
                    'mean': round(float(temp_valid.mean()), 2),
                    'median': round(float(temp_valid.median()), 2),
                    'min': round(float(temp_valid.min()), 2),
                    'max': round(float(temp_valid.max()), 2)
                }
        similarity_stats = {}
        if '_similarity' in sub.columns and len(sub) > 0:
            simv = sub['_similarity'].dropna()
            if len(simv) > 0:
                similarity_stats = {
                    'count': int(len(simv)),
                    'mean': round(float(simv.mean()), 4),
                    'min': round(float(simv.min()), 4),
                    'max': round(float(simv.max()), 4)
                }
        global_similarity_stats = {}
        if '_global_similarity' in sub.columns and len(sub) > 0:
            simv = pd.to_numeric(sub['_global_similarity'], errors='coerce').dropna()
            if len(simv) > 0:
                global_similarity_stats = {
                    'count': int(len(simv)),
                    'mean': round(float(simv.mean()), 4),
                    'min': round(float(simv.min()), 4),
                    'max': round(float(simv.max()), 4)
                }
        local_similarity_stats = {}
        if '_local_similarity' in sub.columns and len(sub) > 0:
            simv = pd.to_numeric(sub['_local_similarity'], errors='coerce').dropna()
            if len(simv) > 0:
                local_similarity_stats = {
                    'count': int(len(simv)),
                    'mean': round(float(simv.mean()), 4),
                    'min': round(float(simv.min()), 4),
                    'max': round(float(simv.max()), 4)
                }
        ring_similarity_stats = {}
        if '_ring_similarity' in sub.columns and len(sub) > 0:
            simv = pd.to_numeric(sub['_ring_similarity'], errors='coerce').dropna()
            if len(simv) > 0:
                ring_similarity_stats = {
                    'count': int(len(simv)),
                    'mean': round(float(simv.mean()), 4),
                    'min': round(float(simv.min()), 4),
                    'max': round(float(simv.max()), 4)
                }
        ring_meta_stats = {}
        if '_ring_meta_score' in sub.columns and len(sub) > 0:
            simv = pd.to_numeric(sub['_ring_meta_score'], errors='coerce').dropna()
            if len(simv) > 0:
                ring_meta_stats = {
                    'count': int(len(simv)),
                    'mean': round(float(simv.mean()), 4),
                    'min': round(float(simv.min()), 4),
                    'max': round(float(simv.max()), 4)
                }
        ring_size_stats = {}
        if '_ring_size' in sub.columns:
            rv = pd.to_numeric(sub['_ring_size'], errors='coerce').dropna()
            if len(rv) > 0:
                ring_size_stats = {
                    'count': int(len(rv)),
                    'mean': round(float(rv.mean()), 2),
                    'median': round(float(rv.median()), 2),
                    'min': round(float(rv.min()), 2),
                    'max': round(float(rv.max()), 2)
                }
        nrows_stats = {}
        if self.nrows_col in sub.columns:
            nr = pd.to_numeric(sub[self.nrows_col], errors='coerce').dropna()
            if len(nr) > 0:
                nrows_stats = {
                    'count': int(len(nr)),
                    'sum': round(float(nr.sum()), 2),
                    'mean': round(float(nr.mean()), 2),
                    'median': round(float(nr.median()), 2),
                    'min': round(float(nr.min()), 2),
                    'max': round(float(nr.max()), 2)
                }
        aggregated_priors = {
            'top_families': family_counts[:5],
            'top_mechanisms': mechanism_counts[:5],
            'top_atmospheres': atmosphere_counts[:5],
            'top_quenchers': quencher_counts[:5],
            'top_catalysts': catalyst_counts[:5],
            'top_initiators': initiator_counts[:5],
            'top_solvents': solvent_counts[:5],
            'top_temperatures_text': temp_text_counts[:5],
            'top_ring_topologies': ring_topology_counts[:5],
            'temperature_numeric_stats': temp_stats,
            'ring_size_stats': ring_size_stats,
            'global_similarity_stats': global_similarity_stats,
            'local_similarity_stats': local_similarity_stats,
            'ring_similarity_stats': ring_similarity_stats,
            'ring_meta_stats': ring_meta_stats,
            'similarity_stats': similarity_stats,
            'n_rows_stats': nrows_stats
        }
        return {
            'n_retrieved': int(len(sub)),
            'family_distribution': family_counts,
            'mechanism_distribution': mechanism_counts,
            'atmosphere_distribution': atmosphere_counts,
            'quencher_distribution': quencher_counts,
            'catalyst_distribution': catalyst_counts,
            'initiator_distribution': initiator_counts,
            'solvent_distribution': solvent_counts,
            'temperature_distribution_text': temp_text_counts,
            'ring_topology_distribution': ring_topology_counts,
            'temperature_numeric_stats': temp_stats,
            'ring_size_stats': ring_size_stats,
            'global_similarity_stats': global_similarity_stats,
            'local_similarity_stats': local_similarity_stats,
            'ring_similarity_stats': ring_similarity_stats,
            'ring_meta_stats': ring_meta_stats,
            'similarity_stats': similarity_stats,
            'n_rows_stats': nrows_stats,
            'aggregated_priors': aggregated_priors
        }


def save_results_excel(result_dict, output_path):
    top_cases = pd.DataFrame(result_dict['top_cases'])
    summary_rows = []
    for k, v in result_dict['summary'].items():
        summary_rows.append({
            'item': k,
            'value': json.dumps(v, ensure_ascii=False) if isinstance(v, (list, dict)) else v
        })
    summary_df = pd.DataFrame(summary_rows)
    with pd.ExcelWriter(output_path, engine='openpyxl') as writer:
        top_cases.to_excel(writer, index=False, sheet_name='top_cases')
        summary_df.to_excel(writer, index=False, sheet_name='summary')


DEFAULT_CONFIG = {
    'topk': 5,
    'family': None,
    'mechanism': 'ROMP',
    'required_fgs': None,
    'min_similarity': 0.6,
    'local_hops': 2,
    'w_global': 0.15,
    'w_local': 0.3,
    'w_ring': 0.35,
    'w_ring_meta': 0.2,
    'require_ring_alkene': True
}


def load_config(path):
    config = dict(DEFAULT_CONFIG)
    if path is None:
        return config
    with open(path, 'r', encoding='utf-8') as handle:
        user_config = json.load(handle)
    unknown = sorted(set(user_config) - set(DEFAULT_CONFIG))
    if unknown:
        raise ValueError(f'Unknown configuration keys: {unknown}')
    config.update(user_config)
    return config


def build_parser():
    parser = argparse.ArgumentParser(description='Retrieve ring-aware ROMP reaction-condition precedents.')
    parser.add_argument(
        '--input',
        required=True,
        help='Reaction table in CSV, XLSX, XLS, or Parquet format.'
    )
    parser.add_argument('--query-smiles', required=True, help='Query monomer SMILES.')
    parser.add_argument('--config', default=None, help='Optional JSON configuration file.')
    parser.add_argument(
        '--topk',
        type=int,
        default=None,
        help='Number of ranked precedents to return.'
    )
    parser.add_argument('--family', default=None, help='Optional polymerization-family filter.')
    parser.add_argument(
        '--mechanism',
        default=None,
        help='Optional mechanism filter; paper default is ROMP.'
    )
    parser.add_argument(
        '--fg',
        nargs='*',
        default=None,
        help='Optional required functional groups, for example alkene thioether.'
    )
    parser.add_argument(
        '--min-similarity',
        type=float,
        default=None,
        help='Minimum weighted similarity; paper default is 0.6.'
    )
    parser.add_argument(
        '--local-hops',
        type=int,
        default=None,
        help='Graph radius in bonds around the reactive alkene.'
    )
    parser.add_argument(
        '--w-global',
        type=float,
        default=None,
        help='Whole-molecule similarity weight.'
    )
    parser.add_argument(
        '--w-local',
        type=float,
        default=None,
        help='Reactive-alkene-local similarity weight.'
    )
    parser.add_argument(
        '--w-ring',
        type=float,
        default=None,
        help='Alkene-containing ring-substructure weight.'
    )
    parser.add_argument(
        '--w-ring-meta',
        type=float,
        default=None,
        help='Ring-size/topology metadata weight.'
    )
    parser.add_argument(
        '--require-ring-alkene',
        action=argparse.BooleanOptionalAction,
        default=None,
        help='Require a ring-contained alkene in candidate records.'
    )
    parser.add_argument('--output', default=None, help='Optional XLSX output path.')
    parser.add_argument('--json-output', default=None, help='Optional JSON output path.')
    return parser


def resolve_parameter(cli_value, config, key):
    return config[key] if cli_value is None else cli_value


def main():
    args = build_parser().parse_args()
    config = load_config(args.config)
    topk = resolve_parameter(args.topk, config, 'topk')
    family = resolve_parameter(args.family, config, 'family')
    mechanism = resolve_parameter(args.mechanism, config, 'mechanism')
    required_fgs = args.fg if args.fg is not None else config['required_fgs']
    min_similarity = resolve_parameter(args.min_similarity, config, 'min_similarity')
    local_hops = resolve_parameter(args.local_hops, config, 'local_hops')
    w_global = resolve_parameter(args.w_global, config, 'w_global')
    w_local = resolve_parameter(args.w_local, config, 'w_local')
    w_ring = resolve_parameter(args.w_ring, config, 'w_ring')
    w_ring_meta = resolve_parameter(args.w_ring_meta, config, 'w_ring_meta')
    require_ring_alkene = resolve_parameter(
        args.require_ring_alkene,
        config,
        'require_ring_alkene'
    )
    if topk < 1:
        raise ValueError('topk must be at least 1.')
    if not 0.0 <= min_similarity <= 1.0:
        raise ValueError('min_similarity must be between 0 and 1.')
    df = load_table(args.input)
    retriever = ReactionRetriever(
        df,
        local_hops=local_hops,
        w_global=w_global,
        w_local=w_local,
        w_ring=w_ring,
        w_ring_meta=w_ring_meta
    )
    result = retriever.retrieve(
        query_smiles=args.query_smiles,
        topk=topk,
        family=family,
        mechanism=mechanism,
        required_fgs=required_fgs,
        min_similarity=min_similarity,
        require_ring_alkene=require_ring_alkene
    )
    payload = json.dumps(result, ensure_ascii=False, indent=2)
    print(payload)
    if args.json_output:
        Path(args.json_output).write_text(payload + '\n', encoding='utf-8')
        print(f'[Saved] {args.json_output}')
    if args.output:
        save_results_excel(result, args.output)
        print(f'[Saved] {args.output}')


if __name__ == '__main__':
    main()
