#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import json
from dataclasses import dataclass, asdict
from typing import List, Dict, Any, Tuple, Optional

import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, rdmolops


@dataclass
class Reason:
    rule_id: str
    severity: str
    title: str
    reason: str
    suggestion: str


SMARTS_LIBRARY = {
    'free_amine': [('[NX3;H2,H1,H0;!$(N[C,S,P]=O)]', 0)],
    'aromatic_n': [('[n;X2]', 0)],
    'thioether': [('[#6]-[SX2]-[#6]', 1)],
    'thiol': [('[SH]', 0)],
    'phosphine': [('[PX3]', 0)],
    'polar_oxygen_donor': [
        ('[OX1]=[CX3]', 0),
        ('[SX3](=[OX1])([#6])[#6]', 1),
        ('[SX4](=[OX1])(=[OX1])([#6])[#6]', 1),
        ('[SX4](=[OX1])(=[OX1])([#6])[#6]', 2)
    ]
}


def safe_mol_from_smiles(smiles: str):
    if smiles is None or (isinstance(smiles, float) and pd.isna(smiles)):
        return None
    smiles = str(smiles).strip()
    if not smiles:
        return None
    try:
        mol = Chem.MolFromSmiles(smiles)
        return mol
    except Exception:
        return None


def compile_smarts_list(smarts_list):
    pats = []
    for item in smarts_list:
        try:
            smarts, donor_pos = item
            q = Chem.MolFromSmarts(smarts)
            if q is not None:
                pats.append((q, donor_pos))
        except Exception:
            pass
    return pats


COMPILED_SMARTS = {k: compile_smarts_list(v) for k, v in SMARTS_LIBRARY.items()}


def match_labeled_atoms(mol, pattern_list) -> List[int]:
    atom_ids = set()
    for patt, donor_pos in pattern_list:
        try:
            ms = mol.GetSubstructMatches(patt, uniquify=True)
            for m in ms:
                if donor_pos < len(m):
                    atom_ids.add(m[donor_pos])
        except Exception:
            pass
    return sorted(atom_ids)


def get_ring_double_bonds(mol) -> List[Dict[str, Any]]:
    ring_double_bonds = []
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and bond.IsInRing():
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()
            atom1 = mol.GetAtomWithIdx(a1)
            atom2 = mol.GetAtomWithIdx(a2)
            if atom1.GetAtomicNum() != 6 or atom2.GetAtomicNum() != 6:
                continue
            ring_sizes = []
            for ring in atom_rings:
                if a1 in ring and a2 in ring:
                    ring_sizes.append(len(ring))
            ring_double_bonds.append({
                'bond_idx': bond.GetIdx(),
                'atoms': (a1, a2),
                'ring_sizes': sorted(set(ring_sizes)),
                'min_ring_size': min(ring_sizes) if ring_sizes else None
            })
    return ring_double_bonds


def count_extra_alkenes(mol) -> Tuple[int, List[int]]:
    n = 0
    extra_bonds = []
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                if not bond.IsInRing():
                    n += 1
                    extra_bonds.append(bond.GetIdx())
    return (n, extra_bonds)


def get_bridgehead_info(mol) -> Dict[str, Any]:
    n_bridgeheads = rdMolDescriptors.CalcNumBridgeheadAtoms(mol)
    return {'num_bridgeheads': int(n_bridgeheads), 'has_bridgehead': bool(n_bridgeheads > 0)}


def get_donor_atoms(mol) -> Dict[str, List[int]]:
    free_amine_atoms = match_labeled_atoms(mol, COMPILED_SMARTS['free_amine'])
    aromatic_n_atoms = match_labeled_atoms(mol, COMPILED_SMARTS['aromatic_n'])
    thioether_atoms = match_labeled_atoms(mol, COMPILED_SMARTS['thioether'])
    thiol_atoms = match_labeled_atoms(mol, COMPILED_SMARTS['thiol'])
    phosphine_atoms = match_labeled_atoms(mol, COMPILED_SMARTS['phosphine'])
    polar_oxygen_donor_atoms = match_labeled_atoms(mol, COMPILED_SMARTS['polar_oxygen_donor'])
    return {
        'free_amine_atoms': free_amine_atoms,
        'aromatic_n_atoms': aromatic_n_atoms,
        'thioether_atoms': thioether_atoms,
        'thiol_atoms': thiol_atoms,
        'phosphine_atoms': phosphine_atoms,
        'polar_oxygen_donor_atoms': polar_oxygen_donor_atoms
    }


def min_distance_to_alkene(
    mol,
    donor_idxs: List[int],
    alkene_atoms: Tuple[int, int]


) -> Optional[int]:
    if not donor_idxs or alkene_atoms is None:
        return None
    a, b = alkene_atoms
    best = None
    for d in donor_idxs:
        if d == a or d == b:
            cur = 0
        else:
            try:
                da = len(rdmolops.GetShortestPath(mol, d, a)) - 1
                db = len(rdmolops.GetShortestPath(mol, d, b)) - 1
                cur = min(da, db)
            except Exception:
                continue
        if best is None or cur < best:
            best = cur
    return best


def estimate_alkene_substitution_score(mol, alkene_atoms: Tuple[int, int]) -> int:
    a, b = alkene_atoms
    score = 0
    for idx_self, idx_other in [(a, b), (b, a)]:
        atom = mol.GetAtomWithIdx(idx_self)
        for nbr in atom.GetNeighbors():
            nidx = nbr.GetIdx()
            if nidx == idx_other:
                continue
            if nbr.GetAtomicNum() > 1:
                score += 1
                heavy_nbrs = sum((1 for x in nbr.GetNeighbors() if x.GetAtomicNum() > 1))
                if heavy_nbrs >= 3:
                    score += 1
                if nbr.GetIsAromatic():
                    score += 1
    return score


def choose_primary_ring_alkene(ring_double_bonds: List[Dict[str, Any]]):
    if not ring_double_bonds:
        return None
    ranked = sorted(
        ring_double_bonds,
        key=lambda x: (999 if x['min_ring_size'] is None else x['min_ring_size'], x['bond_idx'])
    )
    return ranked[0]


def extract_features(smiles: str) -> Dict[str, Any]:
    mol = safe_mol_from_smiles(smiles)
    if mol is None:
        return {'valid': False, 'smiles': smiles, 'error': 'SMILES could not be parsed'}
    ring_double_bonds = get_ring_double_bonds(mol)
    extra_alkene_count, extra_alkene_bonds = count_extra_alkenes(mol)
    bridge_info = get_bridgehead_info(mol)
    donor_info = get_donor_atoms(mol)
    primary_alkene = choose_primary_ring_alkene(ring_double_bonds)
    strong_donor_atoms = sorted(set(donor_info['free_amine_atoms'] + donor_info['aromatic_n_atoms'] + donor_info['thioether_atoms'] + donor_info['thiol_atoms'] + donor_info['phosphine_atoms']))
    if primary_alkene is not None:
        alkene_atoms = primary_alkene['atoms']
        min_ring_size = primary_alkene['min_ring_size']
        alkene_sub_score = estimate_alkene_substitution_score(mol, alkene_atoms)
        min_strong_donor_distance = min_distance_to_alkene(mol, strong_donor_atoms, alkene_atoms)
        min_polar_oxygen_distance = min_distance_to_alkene(
            mol,
            donor_info['polar_oxygen_donor_atoms'],
            alkene_atoms
        )
    else:
        alkene_atoms = None
        min_ring_size = None
        alkene_sub_score = None
        min_strong_donor_distance = None
        min_polar_oxygen_distance = None
    features = {
        'valid': True,
        'smiles': smiles,
        'has_ring_alkene': len(ring_double_bonds) > 0,
        'num_ring_alkenes': len(ring_double_bonds),
        'ring_alkene_bonds': ring_double_bonds,
        'primary_alkene_atoms': alkene_atoms,
        'min_ring_size': min_ring_size,
        'num_bridgeheads': bridge_info['num_bridgeheads'],
        'has_bridgehead': bridge_info['has_bridgehead'],
        'num_extra_alkenes': extra_alkene_count,
        'extra_alkene_bonds': extra_alkene_bonds,
        'free_amine_atoms': donor_info['free_amine_atoms'],
        'aromatic_n_atoms': donor_info['aromatic_n_atoms'],
        'thioether_atoms': donor_info['thioether_atoms'],
        'thiol_atoms': donor_info['thiol_atoms'],
        'phosphine_atoms': donor_info['phosphine_atoms'],
        'polar_oxygen_donor_atoms': donor_info['polar_oxygen_donor_atoms'],
        'has_free_amine': len(donor_info['free_amine_atoms']) > 0,
        'has_aromatic_n': len(donor_info['aromatic_n_atoms']) > 0,
        'has_thioether': len(donor_info['thioether_atoms']) > 0,
        'has_thiol': len(donor_info['thiol_atoms']) > 0,
        'has_phosphine': len(donor_info['phosphine_atoms']) > 0,
        'has_polar_oxygen_donor': len(donor_info['polar_oxygen_donor_atoms']) > 0,
        'min_strong_donor_distance': min_strong_donor_distance,
        'min_polar_oxygen_distance': min_polar_oxygen_distance,
        'alkene_substitution_score': alkene_sub_score
    }
    return features


def apply_rules(features: Dict[str, Any]) -> List[Reason]:
    reasons: List[Reason] = []
    if not features.get('valid', False):
        reasons.append(Reason(
            rule_id='INVALID_SMILES',
            severity='high',
            title='SMILES parsing failed',
            reason='The input could not be parsed by RDKit, so no structure-based heuristic analysis could be performed.',
            suggestion='Please check whether the SMILES string is valid or standardize it before analysis.'
        ))
        return reasons
    has_ring_alkene = features['has_ring_alkene']
    min_ring_size = features['min_ring_size']
    has_bridgehead = features['has_bridgehead']
    num_extra_alkenes = features['num_extra_alkenes']
    has_free_amine = features['has_free_amine']
    has_aromatic_n = features['has_aromatic_n']
    has_thioether = features['has_thioether']
    has_thiol = features['has_thiol']
    has_phosphine = features['has_phosphine']
    has_polar_oxygen_donor = features['has_polar_oxygen_donor']
    min_strong_donor_distance = features['min_strong_donor_distance']
    min_polar_oxygen_distance = features['min_polar_oxygen_distance']
    alkene_sub_score = features['alkene_substitution_score']
    if not has_ring_alkene:
        reasons.append(Reason(
            rule_id='I_1_NO_RING_ALKENE',
            severity='high',
            title='Lack of a typical ROMP-active motif',
            reason='No ring-contained C=C bond was detected, so the structure does not contain the usual cyclic olefin motif associated with conventional ROMP monomers.',
            suggestion='Consider redesigning the polymerizable site as a ring alkene, preferably a strained cycloalkene or bicyclic olefin scaffold.'
        ))
    if has_ring_alkene and min_ring_size is not None and (min_ring_size == 6) and (not has_bridgehead):
        reasons.append(Reason(
            rule_id='II_1_SIX_MEMBERED_LOW_DRIVE',
            severity='medium',
            title='Ring-opening driving force may be limited',
            reason='The target alkene is embedded in a non-bridged six-membered ring. For some such scaffolds, the thermodynamic driving force for ring opening can be limited; this is a qualitative caution rather than a prediction of experimental failure.',
            suggestion='Consider evaluating ring-opening thermodynamics or modifying the scaffold through ring-size adjustment, conformational constraint, or incorporation of a strained bicyclic motif.'
        ))
    if has_ring_alkene and alkene_sub_score is not None and (alkene_sub_score >= 4):
        reasons.append(Reason(
            rule_id='III_1_ALKENE_CONGESTION',
            severity='medium',
            title='The target alkene may be sterically congested',
            reason='The local substitution pattern around the target ring alkene is relatively crowded, which may hinder approach of the metal carbene and slow initiation or propagation.',
            suggestion='Consider reducing steric bulk near the alkene or moving functional groups farther away using a flexible spacer.'
        ))
    if has_free_amine:
        reasons.append(Reason(
            rule_id='IV_1_FREE_AMINE',
            severity='high',
            title='A free amine may interfere with the catalyst',
            reason='A free-amine-like site was detected. Free amines can interact unfavorably with metathesis catalysts and therefore represent an important compatibility risk.',
            suggestion='Consider protecting the amine as an amide, carbamate, or Boc-protected amine, or evaluating a protonated/salt form if appropriate.'
        ))
    if has_aromatic_n:
        reasons.append(Reason(
            rule_id='IV_2_AROMATIC_N',
            severity='medium',
            title='An aromatic nitrogen donor is present',
            reason='An aromatic nitrogen donor motif was detected. Such sites may compete for coordination with the catalyst center and reduce catalyst compatibility.',
            suggestion='Consider lowering the basicity of that nitrogen, relocating it to a more distant side chain, or inserting a spacer between that donor and the polymerizable alkene.'
        ))
    if has_thioether or has_thiol or has_phosphine:
        groups = []
        if has_thioether:
            groups.append('thioether')
        if has_thiol:
            groups.append('thiol')
        if has_phosphine:
            groups.append('phosphine')
        grp_text = ', '.join(groups)
        reasons.append(Reason(
            rule_id='IV_3_STRONG_DONOR',
            severity='high',
            title='Strong donor functionality is present',
            reason=f'A potentially strong donor group such as {grp_text} was detected. Such groups may compete for coordination with the catalyst center and increase the risk of catalyst deactivation or inhibition.',
            suggestion='Consider protecting, removing, or relocating these donor groups farther from the target alkene.'
        ))
    if has_ring_alkene and has_polar_oxygen_donor and (min_polar_oxygen_distance is not None):
        if min_polar_oxygen_distance <= 2:
            reasons.append(Reason(
                rule_id='IV_4_POLAR_OXYGEN_TOO_CLOSE_HIGH',
                severity='high',
                title='A nearby polar oxygen donor may interfere with monomer insertion',
                reason='A polar oxygen donor, such as a carbonyl, sulfoxide, or sulfone oxygen, was detected very close to the target ring alkene. Coordination of this site may perturb the local reactive environment and interfere with monomer insertion.',
                suggestion='Consider increasing the separation between the oxygen-containing group and the target alkene, for example by introducing a -(CH2)n- spacer or relocating that group to a more remote side chain.'
            ))
        elif min_polar_oxygen_distance == 3:
            reasons.append(Reason(
                rule_id='IV_4_POLAR_OXYGEN_TOO_CLOSE_MED',
                severity='medium',
                title='A nearby polar oxygen donor may affect the local insertion environment',
                reason='A polar oxygen donor, such as a carbonyl, sulfoxide, or sulfone oxygen, lies at short graph distance from the target ring alkene and may influence the local coordination environment during monomer insertion.',
                suggestion='Consider moving the oxygen-containing group farther away from the polymerizable alkene.'
            ))
    if has_ring_alkene and min_strong_donor_distance is not None:
        if min_strong_donor_distance <= 2:
            reasons.append(Reason(
                rule_id='IV_5_STRONG_DONOR_TOO_CLOSE_HIGH',
                severity='high',
                title='A strong donor site is too close to the target alkene',
                reason='A strong donor site is located very close to the target ring alkene, which may increase local coordination interference and negatively affect catalyst approach or productive insertion.',
                suggestion='Consider inserting a spacer between the strong donor group and the alkene, or relocating that donor to a more remote side chain.'
            ))
        elif min_strong_donor_distance == 3:
            reasons.append(Reason(
                rule_id='IV_5_STRONG_DONOR_TOO_CLOSE_MED',
                severity='medium',
                title='A strong donor site is relatively close to the target alkene',
                reason='A strong donor group lies at short graph distance from the target ring alkene and may perturb the local reaction environment.',
                suggestion='Consider increasing the separation between the strong donor site and the polymerizable alkene.'
            ))
    if num_extra_alkenes == 1:
        reasons.append(Reason(
            rule_id='V_1_EXTRA_ALKENE',
            severity='medium',
            title='An additional alkene is present',
            reason='Besides the target ring alkene, one extra non-ring C=C bond was detected. This may introduce competing metathesis, chain transfer, or other side-reaction pathways.',
            suggestion='Consider deleting, masking, or installing that extra alkene after polymerization to simplify the ROMP reaction network.'
        ))
    elif num_extra_alkenes >= 2:
        reasons.append(Reason(
            rule_id='V_2_MULTI_ALKENE',
            severity='high',
            title='This is a polyene-like system',
            reason='Besides the target ring alkene, multiple additional non-ring C=C bonds were detected. Such structures often broaden the side-reaction window and make clean, well-controlled ROMP less likely.',
            suggestion='Consider simplifying the structure to a single ROMP-active alkene and introducing the other alkenes in a post-polymerization modification step.'
        ))
    low_drive_flag = any((r.rule_id == 'II_1_SIX_MEMBERED_LOW_DRIVE' for r in reasons))
    extra_alkene_flag = num_extra_alkenes >= 1
    if low_drive_flag and extra_alkene_flag:
        reasons.append(Reason(
            rule_id='V_3_LOW_DRIVE_PLUS_EXTRA_ALKENE',
            severity='high',
            title='Limited driving force and side-reaction risk are both present',
            reason='The structure combines a potentially weak ring-opening driving force with one or more additional alkene sites, which may increase the likelihood of competing metathesis or secondary reactivity.',
            suggestion='It is usually better to address the scaffold driving-force issue first, and then decide whether the additional alkene should be retained.'
        ))
    return reasons


def summarize_risk(reasons: List[Reason]) -> str:
    if not reasons:
        return 'No obvious structural red flags were detected by this heuristic rule set (this does not guarantee experimental ROMP success).'
    high = sum((r.severity == 'high' for r in reasons))
    medium = sum((r.severity == 'medium' for r in reasons))
    if high >= 2:
        return 'This structure appears to be a lower-priority ROMP candidate.'
    elif high == 1 or medium >= 3:
        return 'This structure may be less favorable for conventional ROMP.'
    elif medium >= 1:
        return 'This structure shows several ROMP-relevant cautionary features.'
    else:
        return 'Only minor structural cautions were identified.'


def format_reason_list(reasons: List[Reason]) -> str:
    if not reasons:
        return 'No heuristic rules were triggered.'
    lines = []
    for i, r in enumerate(reasons, 1):
        lines.append(f'{i}. {r.title}: {r.reason} Suggestion: {r.suggestion}')
    return '\n'.join(lines)


def format_suggestions(reasons: List[Reason]) -> str:
    suggs = []
    seen = set()
    for r in reasons:
        if r.suggestion not in seen:
            suggs.append(r.suggestion)
            seen.add(r.suggestion)
    return '; '.join(suggs) if suggs else ''


def analyze_smiles(smiles: str) -> Dict[str, Any]:
    features = extract_features(smiles)
    reasons = apply_rules(features)
    overall = summarize_risk(reasons)
    result = {
        'tool_version': '2.1',
        'interpretation_scope': 'Qualitative structural cautions; not an experimental verdict.',
        'smiles': smiles,
        'valid': features.get('valid', False),
        'overall_interpretation': overall,
        'num_reasons': len(reasons),
        'has_ring_alkene': features.get('has_ring_alkene'),
        'num_ring_alkenes': features.get('num_ring_alkenes'),
        'min_ring_size': features.get('min_ring_size'),
        'has_bridgehead': features.get('has_bridgehead'),
        'num_bridgeheads': features.get('num_bridgeheads'),
        'num_extra_alkenes': features.get('num_extra_alkenes'),
        'has_free_amine': features.get('has_free_amine'),
        'has_aromatic_n': features.get('has_aromatic_n'),
        'has_thioether': features.get('has_thioether'),
        'has_thiol': features.get('has_thiol'),
        'has_phosphine': features.get('has_phosphine'),
        'has_polar_oxygen_donor': features.get('has_polar_oxygen_donor'),
        'min_strong_donor_distance': features.get('min_strong_donor_distance'),
        'min_polar_oxygen_distance': features.get('min_polar_oxygen_distance'),
        'alkene_substitution_score': features.get('alkene_substitution_score'),
        'rule_ids': '; '.join([r.rule_id for r in reasons]),
        'severity_list': '; '.join([r.severity for r in reasons]),
        'titles': '; '.join([r.title for r in reasons]),
        'explanations': format_reason_list(reasons),
        'suggestions': format_suggestions(reasons),
        'features_json': json.dumps(features, ensure_ascii=False),
        'reasons_json': json.dumps([asdict(r) for r in reasons], ensure_ascii=False)
    }
    return result


def read_table(path: str) -> pd.DataFrame:
    if path.lower().endswith('.csv'):
        return pd.read_csv(path)
    elif path.lower().endswith(('.xlsx', '.xls')):
        return pd.read_excel(path)
    else:
        raise ValueError('Only CSV / Excel input is supported')


def write_table(df: pd.DataFrame, path: str):
    if path.lower().endswith('.csv'):
        df.to_csv(path, index=False, encoding='utf-8-sig')
    elif path.lower().endswith(('.xlsx', '.xls')):
        df.to_excel(path, index=False)
    else:
        raise ValueError('Only CSV / Excel output is supported')


def build_argparser():
    p = argparse.ArgumentParser(description='ROMP heuristic explainer v2.1')
    p.add_argument('--smiles', type=str, default=None, help='single SMILES')
    p.add_argument('--input', type=str, default=None, help='input CSV/Excel')
    p.add_argument('--smiles-col', type=str, default='smiles', help='SMILES column name')
    p.add_argument('--output', type=str, default=None, help='output CSV/Excel')
    return p


def main():
    args = build_argparser().parse_args()
    if args.smiles is None and args.input is None:
        raise SystemExit('Please provide either --smiles or --input')
    if args.smiles is not None:
        result = analyze_smiles(args.smiles)
        print('=' * 80)
        print('ROMP Heuristic Analysis Report')
        print('=' * 80)
        print(f"SMILES: {result['smiles']}")
        print(f"Overall interpretation: {result['overall_interpretation']}")
        print(f"Number of triggered rules: {result['num_reasons']}")
        print('-' * 80)
        print('Rationales:')
        print(result['explanations'])
        print('-' * 80)
        print('Suggestions:')
        print(result['suggestions'])
        print('-' * 80)
        print('Structural features:')
        for k in [
            'has_ring_alkene',
            'num_ring_alkenes',
            'min_ring_size',
            'num_extra_alkenes',
            'has_free_amine',
            'has_aromatic_n',
            'has_thioether',
            'has_thiol',
            'has_phosphine',
            'has_polar_oxygen_donor',
            'min_strong_donor_distance',
            'min_polar_oxygen_distance',
            'alkene_substitution_score'
        ]:
            print(f'{k}: {result[k]}')
    if args.input is not None:
        df = read_table(args.input)
        if args.smiles_col not in df.columns:
            raise ValueError(f'SMILES column not found: {args.smiles_col}')
        records = []
        for _, row in df.iterrows():
            smi = row[args.smiles_col]
            res = analyze_smiles(smi)
            out_row = dict(row)
            out_row.update(res)
            records.append(out_row)
        out_df = pd.DataFrame(records)
        if args.output:
            write_table(out_df, args.output)
            print(f'Output written to: {args.output}')
        else:
            print(out_df.head().to_string(index=False))


if __name__ == '__main__':
    main()
