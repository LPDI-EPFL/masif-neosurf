from pathlib import Path
from functools import partial

import numpy as np
from tqdm import tqdm
from Bio.PDB import PDBParser, Selection
from Bio.PDB.Polypeptide import PPBuilder, is_aa
from Bio.PDB.StructureBuilder import StructureBuilder
from Bio.PDB.SASA import ShrakeRupley
import pandas as pd


def parse_score_file(pdb_file):
    score_file = pdb_file.with_suffix('.score')

    with open(score_file, 'r') as f:
        content = {x.split(': ')[0]: x.split(': ')[1] for x in f.read().split(', ')}

    return content


def get_score(pdb_file, score_type='score'):
    scores = parse_score_file(pdb_file)
    return float(scores[score_type])


def rmsd(coords1, coords2):
    assert coords1.shape == coords2.shape
    diff = coords1 - coords2
    return np.sqrt(np.mean(np.sum(diff * diff, axis=-1)))


def compute_rmsd(ref_file, moved_file):

    def atom_filter(atom):
        return atom.element != 'H' and is_aa(atom.parent)

    fixed = Selection.unfold_entities(PDBParser(QUIET = True).get_structure('' , ref_file)[0], 'A')
    moving = Selection.unfold_entities(PDBParser(QUIET = True).get_structure('' , moved_file)[0], 'A')

    fixed = [a for a in fixed if atom_filter(a)]
    moving = [a for a in moving if atom_filter(a)]

    fixed_coord = np.stack([a.get_coord() for a in fixed])
    moving_coord = np.stack([a.get_coord() for a in moving])

    return rmsd(fixed_coord, moving_coord)


def compute_irmsd(ref_file, moved_file, ref_binder_file, iface_all_atom_thresh=5.0):

    def atom_filter(atom):
        return atom.element != 'H' and is_aa(atom.parent)

    fixed = Selection.unfold_entities(PDBParser(QUIET = True).get_structure('' , ref_file)[0], 'A')
    moving = Selection.unfold_entities(PDBParser(QUIET = True).get_structure('' , moved_file)[0], 'A')
    binder = Selection.unfold_entities(PDBParser(QUIET = True).get_structure('' , ref_binder_file)[0], 'A')

    fixed = [a for a in fixed if atom_filter(a)]
    moving = [a for a in moving if atom_filter(a)]
    binder = [a for a in binder if atom_filter(a)]

    fixed_coord = np.stack([a.get_coord() for a in fixed])
    moving_coord = np.stack([a.get_coord() for a in moving])
    binder_coord = np.stack([a.get_coord() for a in binder])

    dist = np.sqrt(np.sum((fixed_coord[:, None, :] - binder_coord[None, :, :])**2, axis=-1))
    iface = np.any(dist < iface_all_atom_thresh, axis=1)

    return rmsd(fixed_coord[iface], moving_coord[iface])


def parse_results(result_dir, binder_dir="without_ligand/data_preparation/01-benchmark_pdbs", benchmark_definition="benchmark_pdbs.txt"):
    results = {}

    score_fn = 'score'  # NN score
    # score_fn = 'desc_dist_score'

    correct_partner = {}
    with open(benchmark_definition, 'r') as f:
        for line in f.readlines():
            line = line.strip()

            if line.startswith("#"):
                continue

            pdb_id, chain_1, chain_2, drug, sdf_name = line.split(',')

            correct_partner[f"{pdb_id}_{chain_1}"] = f"{pdb_id}_{chain_2}"
            correct_partner[f"{pdb_id}_{chain_2}"] = f"{pdb_id}_{chain_1}"

    target_dirs = [d for d in result_dir.iterdir() if d.is_dir()]
    for target_dir in tqdm(target_dirs):
        if not target_dir.is_dir():
            continue

        target = target_dir.name
        results[target] = []

        for site_dir in target_dir.glob("site_*"):

            for chain_match in site_dir.iterdir():
                if not chain_match.is_dir():
                    continue

                matches = list(chain_match.glob("[!.]*.pdb"))

                if len(matches) < 1:
                    continue

                matches = sorted(matches, key=partial(get_score, score_type=score_fn), reverse=True)

                best_match = matches[0]

                # get score
                best_match_score = get_score(best_match, score_type=score_fn)

                match_info = {'match': f"{site_dir.name}_{best_match.name}", 'score': best_match_score}

                # compute RMSD
                if chain_match.name == correct_partner[target]:  # correct match
                    match_info['RMSD'] = compute_rmsd(Path(binder_dir, correct_partner[target] + '.pdb'), best_match)
                    match_info['iRMSD'] = compute_irmsd(Path(binder_dir, correct_partner[target] + '.pdb'),
                                                        best_match, Path(binder_dir, target + '.pdb'))

                # append to results
                results[target].append(match_info)

    return results


def get_topk(results, k_values):
    found_cumulated = []
    solved = set()
    for k in k_values:
        if k == 0:
            found_cumulated.append(0)
            continue

        for pdbid, result_list in results.items():
            if len(result_list) <= k:
                continue

            sorted_results = sorted(result_list, key=lambda x: x['score'], reverse=True)
            kth_result = sorted_results[k - 1]
            if 'RMSD' in kth_result:  # correct match (determined above)
                if kth_result['iRMSD'] < 5.0:
                    solved.add(pdbid)
                else:
                    print(f"Found match for {pdbid} but with RMSD {kth_result['RMSD']}")

        found_cumulated.append(len(solved))

    print("Not solved:", set(results.keys()) - solved)

    return found_cumulated


def extract_chains(input_struct, chains, keep_hetatm=False):
    builder = StructureBuilder()
    builder.init_structure('complex')
    builder.init_model(0)
    for c in chains:
        builder.init_chain(c)
    out_struct = builder.get_structure()[0]

    for c in chains:
        for res in input_struct[c].get_residues():
            if is_aa(res) or keep_hetatm:
                _res = res.copy()
                out_struct[c].add(_res)

    return out_struct


def compute_sasa_values(pdb_file, target_chain, binder_chain, ligand_name, ligand_chain):
    pdb_file = Path(pdb_file)
    target_chain = ''.join(sorted(set(target_chain) | set(ligand_chain)))

    out = {'name': pdb_file.name}

    whole_complex = PDBParser(QUIET=True).get_structure('', pdb_file)[0]
    target = extract_chains(whole_complex, target_chain, keep_hetatm=True)
    binder = extract_chains(whole_complex, binder_chain, keep_hetatm=False)

    # Target
    ShrakeRupley().compute(target, level="R")  # compute SASA
    out['target_sasa'] = sum([r.sasa for r in target.get_residues()])

    ligand = [res for res in target[ligand_chain].get_residues() if res.get_resname() == ligand_name]
    assert len(ligand) == 1
    ligand = ligand[0]
    out['ligand_sasa'] = ligand.sasa

    # Binder
    ShrakeRupley().compute(binder, level="R")  # compute SASA
    out['binder_sasa'] = sum([r.sasa for r in binder.get_residues()])

    # Complex
    ShrakeRupley().compute(whole_complex, level="R")
    out['complex_sasa'] = sum([r.sasa for r in whole_complex.get_residues()])
    target_residues_in_complex = [res for c in target_chain for res in whole_complex[c].get_residues()]
    out['target_complex_sasa'] = sum([r.sasa for r in target_residues_in_complex])

    out['target_buried_sasa'] = out['target_sasa'] - out['target_complex_sasa']

    ligand_in_complex = [res for res in whole_complex[ligand_chain].get_residues() if res.get_resname() == ligand_name]
    assert len(ligand_in_complex) == 1
    ligand_in_complex = ligand_in_complex[0]

    out['ligand_complex_sasa'] = ligand_in_complex.sasa
    out['ligand_buried_sasa'] = out['ligand_sasa'] - out['ligand_complex_sasa']
    out['ligand_sasa_per_target_buried_sasa'] = out['ligand_buried_sasa'] / out['target_buried_sasa']

    return out


def parse_results_pd(result_dir):
    results = []

    target_dirs = [d for d in Path(result_dir).iterdir() if d.is_dir()]
    for target_dir in tqdm(target_dirs):
        if not target_dir.is_dir():
            continue

        target = target_dir.name

        for site_dir in target_dir.glob("site_*"):

            for chain_match in site_dir.iterdir():
                if not chain_match.is_dir():
                    continue

                matches = list(chain_match.glob("[!.]*.pdb"))

                if len(matches) < 1:
                    continue

                for match in matches:
                    summary = parse_score_file(match)
                    match_info = {
                        'target': target,
                        'target_path': Path(target_dir.resolve(), f'{target}.pdb'),
                        'target_site': site_dir.name,
                        'matched_protein': summary['name'],
                        'matched_patch_id': int(summary['point id']),
                        'matched_protein_path': match.resolve(),
                        'score': float(summary['score']),  # NN score
                        'desc_dist_score': float(summary['desc_dist_score']),
                        'clashing_ca': int(summary['clashing_ca']),
                        'clashing_heavy': int(summary['clashing_heavy']),
                    }

                    # append to results
                    results.append(match_info)

    return pd.DataFrame(results)
