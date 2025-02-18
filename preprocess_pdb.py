import sys
import os
import shutil
from pathlib import Path
from argparse import ArgumentParser
import tempfile
import time

import numpy as np
import pymesh

# import MaSIF modules
masif_neosurf_root = Path(__file__).resolve().parent
masif_root = Path(masif_neosurf_root, 'masif')
sys.path.append(str(Path(masif_root, 'source').resolve()))
from masif.source.default_config.masif_opts import masif_opts
from masif.source.triangulation.computeMSMS import computeMSMS
from masif.source.triangulation.fixmesh import fix_mesh
from masif.source.triangulation.ligand_utils import extract_ligand
from masif.source.input_output.extractPDB import extractPDB
from masif.source.input_output.save_ply import save_ply
from masif.source.input_output.protonate import protonate
from masif.source.triangulation.computeHydrophobicity import computeHydrophobicity
from masif.source.triangulation.computeCharges import computeCharges, assignChargesToNewMesh
from masif.source.triangulation.computeAPBS import computeAPBS
from masif.source.triangulation.compute_normal import compute_normal
from masif.source.masif_modules.read_data_from_surface import read_data_from_surface, compute_shape_complementarity
from masif.source.masif_modules.MaSIF_site import MaSIF_site
from masif.source.masif_modules.train_masif_site import run_masif_site
from masif.source.masif_modules.MaSIF_ppi_search import MaSIF_ppi_search
from masif.source.masif_modules.train_ppi_search import compute_val_test_desc
from masif.data.masif_site.nn_models.all_feat_3l.custom_params import custom_params as masif_site_params
from masif.data.masif_ppi_search.nn_models.sc05.all_feat.custom_params import custom_params as masif_search_params


def extract_het_dict(ligand_code, in_file, out_file):

    with open(in_file, "r") as f:
        blocks = f.read().split("\n\n")

    # find the relevant entry
    block = [x for x in blocks if x.startswith(f"RESIDUE   {ligand_code}")]
    assert len(block) == 1, f"Need to find exactly one more for ligand code {ligand_code}, found {len(block)}"
    block = block[0]

    # replace original name with three-letter code
    block = block.replace(ligand_code, ligand_code[:3])

    with open(out_file, "w") as f:
        f.write(block)


def extract_and_triangulate(pdb_filename, name_chain, outdir, tmp_dir, ligand_name_chain=None, sdf_file=None, mol2_patch=None):
    # Process inputs
    pdb_id, chain_ids1 = name_chain.split("_")
    if ligand_name_chain is None:
        ligand_code, ligand_chain = None, None
    else:
        ligand_code, ligand_chain = ligand_name_chain.split('_')
    ligand_tla = ligand_code[:3]

    # Output locations
    pdb_chain_dir = Path(outdir, masif_opts['pdb_chain_dir'])
    ply_chain_dir = Path(outdir, masif_opts['ply_chain_dir'])
    pdb_chain_dir.mkdir(parents=True, exist_ok=True)
    ply_chain_dir.mkdir(parents=True, exist_ok=True)
    outfile_stem = Path(pdb_id + "_" + chain_ids1)
    outfile_pdb = Path(pdb_chain_dir, outfile_stem.with_suffix(".pdb"))
    outfile_ply = Path(ply_chain_dir, outfile_stem.with_suffix(".ply"))
    tmp_file_base = Path(tmp_dir, pdb_id + "_" + chain_ids1)

    # Protonate structure.
    protonated_file = Path(tmp_dir, pdb_id + "_protonated.pdb")
    if len(ligand_code) > 3:
        # if the ligand code is too long, reduce can't find the correct entry based on the abbreviated three-letter code in the pdb file
        # so we hack the hetero atom dictionary instead
        with tempfile.NamedTemporaryFile() as tmp_het_dict:
            default_het_dict = os.environ.get('REDUCE_HET_DICT')
            extract_het_dict(ligand_code, default_het_dict, tmp_het_dict.name)
            protonate(pdb_filename, protonated_file, het_dict=tmp_het_dict.name)
    else:
        protonate(pdb_filename, protonated_file)
    pdb_filename = protonated_file

    # Extract chains of interest.
    extractPDB(pdb_filename, str(outfile_pdb), chain_ids1, ligand_tla, ligand_chain)

    # Compute MSMS of surface w/hydrogens,
    vertices1, faces1, normals1, names1, areas1 = computeMSMS(outfile_pdb, protonate=True, ligand_code=ligand_tla)

    # Get and RDKit molecule object
    mol2_file, rdmol = None, None
    if ligand_code is not None and ligand_chain is not None:
        mol2_file = str(Path(tmp_dir, f"{ligand_code}_{ligand_chain}.mol2").resolve())
        rdmol = extract_ligand(outfile_pdb, ligand_code, ligand_chain, mol2_file, sdf_template=sdf_file, patched_mol2_file=mol2_patch)

    # Compute "charged" vertices
    if masif_opts['use_hbond']:
        vertex_hbond = computeCharges(str(outfile_pdb.with_suffix("")), vertices1, names1, ligand_tla, rdmol)

    # For each surface residue, assign the hydrophobicity of its amino acid. 
    if masif_opts['use_hphob']:
        vertex_hphobicity = computeHydrophobicity(names1, ligand_tla, rdmol)

    # If protonate = false, recompute MSMS of surface, but without hydrogens (set radius of hydrogens to 0).
    vertices2 = vertices1
    faces2 = faces1

    # Fix the mesh.
    mesh = pymesh.form_mesh(vertices2, faces2)
    print(f"Fixing mesh...")
    regular_mesh = fix_mesh(mesh, masif_opts['mesh_res'])
    print(f"Fixmesh done!")

    # Compute the normals
    vertex_normal = compute_normal(regular_mesh.vertices, regular_mesh.faces)
    # Assign charges on new vertices based on charges of old vertices (nearest
    # neighbor)

    if masif_opts['use_hbond']:
        vertex_hbond = assignChargesToNewMesh(regular_mesh.vertices, vertices1, vertex_hbond, masif_opts)

    if masif_opts['use_hphob']:
        vertex_hphobicity = assignChargesToNewMesh(regular_mesh.vertices, vertices1, vertex_hphobicity, masif_opts)

    if masif_opts['use_apbs']:
        print(f"Computing APBS...")
        shutil.copy(outfile_pdb, tmp_file_base.with_suffix(".pdb"))
        vertex_charges = computeAPBS(regular_mesh.vertices, str(outfile_pdb), str(tmp_file_base), mol2_file)
        print(f"APBS done!")

    # Convert to ply and save.
    save_ply(str(outfile_ply), regular_mesh.vertices, regular_mesh.faces, normals=vertex_normal, charges=vertex_charges, normalize_charges=True, hbond=vertex_hbond, hphob=vertex_hphobicity, iface=np.zeros(len(regular_mesh.vertices)))


def masif_precompute(ppi_pair_list, masif_app, output_root):
    assert masif_app in {'masif_ppi_search', 'masif_site'}

    if masif_app == 'masif_ppi_search': 
        params = masif_opts['ppi_search']
    elif masif_app == 'masif_site':
        params = masif_opts['site']
        params['ply_chain_dir'] = masif_opts['ply_chain_dir']
    elif masif_app == 'masif_ligand':
        params = masif_opts['ligand']

    ply_file_template = str(Path(output_root, masif_opts['ply_file_template']))

    np.random.seed(0)
    print('Reading data from input ply surface files.')
    for ppi_pair_id in ppi_pair_list:

        my_precomp_dir = Path(output_root, params['masif_precomputation_dir'], ppi_pair_id)
        my_precomp_dir.mkdir(exist_ok=True, parents=True)
        
        # Read directly from the ply file.
        fields = ppi_pair_id.split('_')
        ply_file = {}
        ply_file['p1'] = ply_file_template.format(fields[0], fields[1])

        if len (fields) == 2 or fields[2] == '':
            pids = ['p1']
        else:
            ply_file['p2'] = ply_file_template.format(fields[0], fields[2])
            pids = ['p1', 'p2']
            
        # Compute shape complementarity between the two proteins. 
        rho = {}
        neigh_indices = {}
        mask = {}
        input_feat = {}
        theta = {}
        iface_labels = {}
        verts = {}

        for pid in pids:
            input_feat[pid], rho[pid], theta[pid], mask[pid], neigh_indices[pid], iface_labels[pid], verts[pid] = read_data_from_surface(ply_file[pid], params)

        if len(pids) > 1 and masif_app == 'masif_ppi_search':
            start_time = time.time()
            p1_sc_labels, p2_sc_labels = compute_shape_complementarity(ply_file['p1'], ply_file['p2'], neigh_indices['p1'],neigh_indices['p2'], rho['p1'], rho['p2'], mask['p1'], mask['p2'], params)
            np.save(Path(my_precomp_dir, 'p1_sc_labels'), p1_sc_labels)
            np.save(Path(my_precomp_dir, 'p2_sc_labels'), p2_sc_labels)
            end_time = time.time()
            print("Computing shape complementarity took {:.2f}".format(end_time-start_time))

        # Save data only if everything went well. 
        for pid in pids: 
            np.save(Path(my_precomp_dir, pid + '_rho_wrt_center'), rho[pid])
            np.save(Path(my_precomp_dir, pid + '_theta_wrt_center'), theta[pid])
            np.save(Path(my_precomp_dir, pid + '_input_feat'), input_feat[pid])
            np.save(Path(my_precomp_dir, pid + '_mask'), mask[pid])
            np.save(Path(my_precomp_dir, pid + '_list_indices'), neigh_indices[pid])
            np.save(Path(my_precomp_dir, pid + '_iface_labels'), iface_labels[pid])
            # Save x, y, z
            np.save(Path(my_precomp_dir, pid + '_X.npy'), verts[pid][:,0])
            np.save(Path(my_precomp_dir, pid + '_Y.npy'), verts[pid][:,1])
            np.save(Path(my_precomp_dir, pid + '_Z.npy'), verts[pid][:,2])


def mask_input_feat(input_feat, mask):
        """Apply mask to input_feat"""
        mymask = np.where(np.array(mask) == 0.0)[0]
        return np.delete(input_feat, mymask, axis=2)


def predict_binding_sites(ppi_pair_ids, output_root):

    params = masif_opts["site"]
    for key in masif_site_params:
        print("Setting {} to {} ".format(key, masif_site_params[key]))
        params[key] = masif_site_params[key]

    # Shape precomputation dir.
    parent_in_dir = Path(output_root, params["masif_precomputation_dir"])
    model_base = Path(masif_root, 'data', 'masif_site', params["model_dir"], 'model')

    # Build the neural network model
    learning_obj = MaSIF_site(
        params["max_distance"],
        n_thetas=4,
        n_rhos=3,
        n_rotations=4,
        idx_gpu="/gpu:0",
        feat_mask=params["feat_mask"],
        n_conv_layers=params["n_conv_layers"],
    )
    print("Restoring model from: " + str(model_base))
    learning_obj.saver.restore(learning_obj.session, str(model_base))
    out_pred_dir = Path(output_root, params["out_pred_dir"])
    out_surf_dir = Path(output_root, params["out_surf_dir"])
    out_pred_dir.mkdir(exist_ok=True, parents=True)
    out_surf_dir.mkdir(exist_ok=True, parents=True)

    for ppi_pair_id in ppi_pair_ids:
        print(ppi_pair_id)
        in_dir = Path(parent_in_dir, ppi_pair_id)

        fields = ppi_pair_id.split('_')
        if len(fields) < 2:
            continue
        pdbid = ppi_pair_id.split("_")[0]
        chain1 = ppi_pair_id.split("_")[1]
        pids = ["p1"]
        chains = [chain1]
        if len(fields) == 3 and fields[2] != "":
            chain2 = fields[2]
            pids = ["p1", "p2"]
            chains = [chain1, chain2]

        for ix, pid in enumerate(pids):
            pdb_chain_id = pdbid + "_" + chains[ix]

            print("Evaluating {}".format(pdb_chain_id))

            rho_wrt_center = np.load(Path(in_dir, pid + "_rho_wrt_center.npy"))
            theta_wrt_center = np.load(Path(in_dir, pid + "_theta_wrt_center.npy"))
            input_feat = np.load(Path(in_dir, pid + "_input_feat.npy"))
            input_feat = mask_input_feat(input_feat, params["feat_mask"])
            mask = np.load(Path(in_dir, pid + "_mask.npy"))
            indices = np.load(Path(in_dir, pid + "_list_indices.npy"), encoding="latin1", allow_pickle=True)

            print("Total number of patches:{} \n".format(len(mask)))

            tic = time.time()
            scores = run_masif_site(
                params,
                learning_obj,
                rho_wrt_center,
                theta_wrt_center,
                input_feat,
                mask,
                indices,
            )
            toc = time.time()
            print(
                "Total number of patches for which scores were computed: {}\n".format(
                    len(scores[0])
                )
            )
            print("GPU time (real time, not actual GPU time): {:.3f}s".format(toc-tic))
            np.save(
                Path(out_pred_dir, "pred_" + pdbid + "_" + chains[ix] + ".npy"),
                scores,
            )

            # Save as ply file
            ply_file = Path(output_root, masif_opts["ply_file_template"].format(pdbid, chains[ix]))
            mymesh = pymesh.load_mesh(str(ply_file))
            mymesh.remove_attribute("vertex_iface")
            mymesh.add_attribute("iface")
            mymesh.set_attribute("iface", scores[0])
            mymesh.remove_attribute("vertex_x")
            mymesh.remove_attribute("vertex_y")
            mymesh.remove_attribute("vertex_z")
            mymesh.remove_attribute("face_vertex_indices")

            out_surf_file = Path(out_surf_dir, pdbid + "_" + chains[ix] + ".ply")
            pymesh.save_mesh(
                str(out_surf_file),
                mymesh,
                *mymesh.get_attribute_names(),
                use_float=True,
                ascii=True
            )
            print(f"Successfully saved file {out_surf_file}")


def compute_descriptors(ppi_list, output_root):
    params = masif_opts["ppi_search"]

    for key in masif_search_params:
        print("Setting {} to {} ".format(key, masif_search_params[key]))
        params[key] = masif_search_params[key]

    # Read the positive first
    parent_in_dir = Path(output_root, params["masif_precomputation_dir"])
    model_base = Path(masif_root, 'data', 'masif_ppi_search', params["model_dir"], 'model')

    np.random.seed(0)

    #   Load existing network.
    print("Reading pre-trained network")
    learning_obj = MaSIF_ppi_search(
        params["max_distance"],
        n_thetas=16,
        n_rhos=5,
        n_rotations=16,
        idx_gpu="/gpu:0",
        feat_mask=params["feat_mask"],
    )
    learning_obj.saver.restore(learning_obj.session, str(model_base))

    desc_dir = Path(output_root, params["desc_dir"])
    desc_dir.mkdir(exist_ok=True, parents=True)

    for count, ppi_pair_id in enumerate(ppi_list):

        in_dir = Path(parent_in_dir, ppi_pair_id)
        print(ppi_pair_id)

        def _compute_descs(pid):
            tic = time.time()
            rho_wrt_center = np.load(Path(in_dir, pid + "_rho_wrt_center.npy"))
            theta_wrt_center = np.load(Path(in_dir, pid + "_theta_wrt_center.npy"))
            input_feat = np.load(Path(in_dir, pid + "_input_feat.npy"))
            input_feat = mask_input_feat(input_feat, params["feat_mask"])
            mask = np.load(Path(in_dir, pid + "_mask.npy"))
            idx = np.array(range(len(rho_wrt_center)))
            print("Data loading time: {:.2f}s".format(time.time() - tic))
            tic = time.time()
            desc_str = compute_val_test_desc(
                learning_obj,
                idx,
                rho_wrt_center,
                theta_wrt_center,
                input_feat,
                mask,
                batch_size=1000,
                flip=False,
            )
            desc_flip = compute_val_test_desc(
                learning_obj,
                idx,
                rho_wrt_center,
                theta_wrt_center,
                input_feat,
                mask,
                batch_size=1000,
                flip=True,
            )
            print("Running time: {:.2f}s".format(time.time() - tic))
            return desc_str, desc_flip

        out_desc_dir = Path(desc_dir, ppi_pair_id)
        out_desc_dir.mkdir(exist_ok=True)

        desc1_str, desc1_flip = _compute_descs("p1")
        np.save(Path(out_desc_dir, "p1_desc_straight.npy"), desc1_str)
        np.save(Path(out_desc_dir, "p1_desc_flipped.npy"), desc1_flip)

        # second protein specified
        if len(ppi_pair_id.split("_")) > 2:
            desc2_str, desc2_flip = _compute_descs("p2")
            np.save(Path(out_desc_dir, "p2_desc_straight.npy"), desc2_str)
            np.save(Path(out_desc_dir, "p2_desc_flipped.npy"), desc2_flip)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("pdbfile", type=Path, help="Input PDB file.")
    parser.add_argument("name_chain", type=str, help="Protein name and chain, e.g. '1A7X_A'.")
    parser.add_argument("-o", "--outdir", type=Path, help="Folder in which processed files will be saved.")
    parser.add_argument("-l", "--ligand", type=str, default=None,
                        help="Three letter code and chain of a ligand, e.g. 'FKA_B'.")
    parser.add_argument("-s", "--sdf", type=str, default=None, 
                        help="Optional SDF file used to infer the ligand bond types.")
    parser.add_argument("-m", "--mol2", type=Path, default=None, 
                        help="Optional custom mol2 file. Should not be necessary in most cases.")
    parser.add_argument("--tmp_dir", type=Path, default=None,
                        help="Directory where temporary files will be saved. Provide a path if you would like to inspect these files for debugging.")
    args = parser.parse_args()

    if args.tmp_dir is not None:
        args.tmp_dir.mkdir(exist_ok=True, parents=True)

    with tempfile.TemporaryDirectory() as _tmp_dir:
        tmp_dir = args.tmp_dir or _tmp_dir
        extract_and_triangulate(args.pdbfile, args.name_chain, args.outdir, tmp_dir, ligand_name_chain=args.ligand, sdf_file=args.sdf, mol2_patch=args.mol2)
        masif_precompute([args.name_chain], "masif_site", args.outdir)
        masif_precompute([args.name_chain], "masif_ppi_search", args.outdir)
        predict_binding_sites([args.name_chain], args.outdir)
        compute_descriptors([args.name_chain], args.outdir)
