import os
import sys
import shutil
from pathlib import Path
from argparse import ArgumentParser
from packaging import version

import pymesh
from scipy.spatial import cKDTree
import numpy as np
from Bio.PDB import PDBParser, PDBIO
import open3d as o3d

# import MaSIF modules
masif_neosurf_dir = Path(__file__).resolve().parent
sys.path.append(str(Path(masif_neosurf_dir, 'masif', 'source').resolve()))
sys.path.append(str(Path(masif_neosurf_dir, 'masif_seed_search', 'source').resolve()))
from masif.source.default_config.masif_opts import masif_opts
from masif_seed_search.source.alignment_evaluation_nn import AlignmentEvaluationNN
from masif_seed_search.source.alignment_utils import get_patch_coords, load_protein_pcd, get_patch_geo, get_target_vix, match_descriptors, align_protein, compute_nn_score


def parse_transform_arg(raw_transform):
    """
    Parse a rigid transform from either:
    - path to .npy file containing a (4, 4) matrix (or flat 16-vector)
    - comma-separated 16-value string in row-major order
    """
    path_candidate = Path(raw_transform)
    if str(path_candidate).endswith(".npy"):
        if not path_candidate.exists():
            raise ValueError(f"Transform file does not exist: {path_candidate}")
        transform = np.load(path_candidate)
        transform_source = f"file:{path_candidate}"
    else:
        parts = [x.strip() for x in raw_transform.split(",") if x.strip()]
        if len(parts) != 16:
            raise ValueError(
                "Inline transform must contain exactly 16 comma-separated numeric values."
            )
        try:
            transform = np.asarray([float(x) for x in parts], dtype=float)
        except ValueError as exc:
            raise ValueError(
                "Inline transform contains non-numeric values."
            ) from exc
        transform_source = "inline"

    transform = np.asarray(transform, dtype=float)
    if transform.shape == (16,):
        transform = transform.reshape(4, 4)
    elif transform.shape != (4, 4):
        raise ValueError(
            f"Transform must have shape (4, 4) or flat 16 values, got {transform.shape}."
        )
    if not np.all(np.isfinite(transform)):
        raise ValueError("Transform contains non-finite values.")
    return transform, transform_source


def apply_transform_to_point_cloud(pcd, transform):
    """Apply a rigid transform to point-cloud points and normals."""
    rotation = transform[:3, :3]
    translation = transform[:3, 3]

    points = np.asarray(pcd.points)
    normals = np.asarray(pcd.normals)

    transformed_points = points @ rotation.T + translation
    transformed_normals = normals @ rotation.T
    normal_norm = np.linalg.norm(transformed_normals, axis=1, keepdims=True)
    normal_norm[normal_norm < 1e-12] = 1.0
    transformed_normals = transformed_normals / normal_norm

    pcd.points = o3d.utility.Vector3dVector(transformed_points)
    pcd.normals = o3d.utility.Vector3dVector(transformed_normals)

# Helper function for writing debug PDB files to visualize the patches used for scoring
def _write_debug_transformed_binder_pdb(source_name, source_paths, debug_dir, transform=None):
    """Write binder PDB with optional rigid transform applied."""
    pdb_parser = PDBParser(QUIET=True)
    binder_pdb_path = Path(source_paths["pdb_dir"]) / f"{source_name}.pdb"
    if not binder_pdb_path.exists():
        print(f"Warning: binder PDB not found for debug export: {binder_pdb_path}")
        return

    binder_struct = pdb_parser.get_structure(str(binder_pdb_path), str(binder_pdb_path))
    if transform is not None:
        rotation = transform[:3, :3]
        translation = transform[:3, 3]
        for atom in binder_struct.get_atoms():
            atom.set_coord(atom.get_coord() @ rotation.T + translation)

    out_pdb = Path(debug_dir) / "binder_transformed.pdb"
    pdb_writer = PDBIO()
    pdb_writer.set_structure(binder_struct)
    pdb_writer.save(str(out_pdb))


def score_complex(
    target_name,
    target_vertices,
    source_name,
    target_pcd,
    target_coord,
    target_desc,
    target_pcd_tree,
    target_iface,
    source_paths,
    params,
    nn_score,
    flip_target_normals=True,
):
    if params.get("site_debug", None) is not None:
        site_debug_dir = Path(params["site_debug"])
        site_debug_dir.mkdir(parents=True, exist_ok=True)
        _write_debug_transformed_binder_pdb(
            source_name=source_name,
            source_paths=source_paths,
            debug_dir=site_debug_dir,
            transform=params.get("score_binder_transform", None),
        )

    # Go through every selected site
    for site_ix, target_vix in enumerate(target_vertices):

        # Get the geodesic patch and descriptor patch for each target patch
        target_patch, target_patch_descs, target_patch_idx = get_patch_geo(
            target_pcd, target_coord, target_vix, target_desc,
            flip_normals=flip_target_normals,
            outward_shift=params['surface_outward_shift']
        )
        if params.get("site_debug", None) is not None:
            site_debug_dir = Path(params["site_debug"])
            site_debug_dir.mkdir(parents=True, exist_ok=True)
            out_query_patch = site_debug_dir / f"site_vix_{target_vix}.vert"
            with out_query_patch.open("w+") as out_patch:
                for point in target_patch.points:
                    out_patch.write('{}, {}, {}\n'.format(point[0], point[1], point[2]))

        # Make a ckdtree with the target vertices.
        target_ckdtree = cKDTree(target_patch.points)

        # Load binder
        pid = 'p1'
        if pid == 'p1':
            chain = source_name.split('_')[1]
            chain_number = 1
        else:
            chain = source_name.split('_')[2]
            chain_number = 2
        source_pcd, source_desc, source_iface = load_protein_pcd(source_name, chain_number, source_paths, flipped_features=False, read_mesh=False)
        if params.get("score_binder_transform", None) is not None:
            apply_transform_to_point_cloud(source_pcd, params["score_binder_transform"])


        # Find closest point on binder
        target_vix_coord = np.asarray(target_pcd.points)[target_vix]
        source_points = np.asarray(source_pcd.points)
        dists = np.linalg.norm(target_vix_coord[None, :] - source_points, axis=-1)
        source_vix = np.argmin(dists)

        # Compute descriptor distance
        desc_dist = np.sqrt(np.sum(np.square(source_desc[0, source_vix] - target_desc[0, target_vix]), axis=-1))

        # Compute NN and descriptor distance scores
        source_coord = get_patch_coords(params['seed_precomp_dir'], source_name, pid, cv=[source_vix])
        source_patch, source_patch_descs, source_patch_idx = get_patch_geo(source_pcd, source_coord, source_vix, source_desc, outward_shift=params['surface_outward_shift'])
        if params.get("site_debug", None) is not None:
            site_debug_dir = Path(params["site_debug"])
            site_debug_dir.mkdir(parents=True, exist_ok=True)
            out_binder_patch = site_debug_dir / f"binder_vix_{source_vix}.vert"
            with out_binder_patch.open("w+") as out_patch:
                for point in source_patch.points:
                    out_patch.write('{}, {}, {}\n'.format(point[0], point[1], point[2]))
        d_vi_at, _= target_pcd_tree.query(np.asarray(source_patch.points), k=1)

        align_scores, _ = compute_nn_score(
            target_patch, source_patch, None, target_patch_descs,
            source_patch_descs, target_ckdtree, nn_score, d_vi_at, 1.0
        )

        results = {
            "query_name": target_name,
            "query_site": site_ix,
            "query_vix": target_vix,
            "binder_name": source_name,
            "binder_vix": source_vix,
            "distance_between_center_points": dists[source_vix],
            "query_iface_score": target_iface[target_vix],
            "binder_iface_score": source_iface[source_vix],
            "descriptor_distance": desc_dist,
            "nn_score": align_scores[0],
            "desc_dist_score": align_scores[1],
            "mean_desc_dist_score": align_scores[2],
        }
        print(', '.join(f"{k}: {v}" for k, v in results.items()))


def masif_search(params):

    # Decide whether to run complementarity or similarity search
    if params['similarity_mode']:
        params['surface_outward_shift'] = 0.0
        params['allowed_CA_clashes'] = float('inf')
        params['allowed_heavy_atom_clashes'] = float('inf')
        params['nn_score_cutoff'] = -1.0  # neural network score does not work here
        flip_target_features = False
        flip_target_normals = False
        print("Running surface similarity search. (experimental feature)")
    else:
        flip_target_features = True
        flip_target_normals = True
        print("Running surface complementarity search.")

    # # Load target patches.
    target_ppi_pair_id = params['target_name']
    target_pid = 'p1'
    target_chain_ix = 1

    seed_ppi_pair_ids = params['seed_ppi_pair_ids']

    # Initialize two neural networks - one that does not account for atomic clashes (initial filter) and one with clashes.
    nn_score_atomic = AlignmentEvaluationNN(params['nn_score_atomic_fn'], selected_features=[0,1,2,3], max_npoints=params['max_npoints'])
    nn_score_atomic.restore_model()

    # Go through every 12A patch in the target protein -- get a sorted least in order of the highest iface mean in the patch
    target_ply_fn = os.path.join(params['target_ply_iface_dir'], target_ppi_pair_id + '.ply')
    mymesh = pymesh.load_mesh(target_ply_fn)
    iface = mymesh.get_attribute('vertex_iface')
    target_coord = get_patch_coords(params['target_precomp_dir'], target_ppi_pair_id, target_pid)


    # Define target and source paths (for interface scores, descriptors, ply files)
    target_paths = {}
    target_paths['surf_dir'] = params['target_surf_dir']
    target_paths['iface_dir'] = params['target_iface_dir']
    target_paths['desc_dir'] = params['target_desc_dir']

    source_paths = {}
    source_paths['surf_dir'] = params['seed_surf_dir']
    source_paths['iface_dir'] = params['seed_iface_dir']
    source_paths['desc_dir'] = params['seed_desc_dir']
    source_paths['pdb_dir'] = params['seed_pdb_dir']

    # Load the target point cloud, descriptors, interface and mesh.
    target_pcd, target_desc, target_iface, target_mesh = load_protein_pcd(target_ppi_pair_id, target_chain_ix, target_paths, flipped_features=flip_target_features, read_mesh=True)

    # Open the pdb structure of the target, load into point clouds for fast access.
    parser = PDBParser()
    target_pdb_path = os.path.join(params['target_pdb_dir'], '{}.pdb'.format(target_ppi_pair_id))
    target_struct = parser.get_structure(target_pdb_path, target_pdb_path)
    target_atom_coords = [atom.get_coord() for atom in target_struct.get_atoms() if not atom.get_name().startswith('H') ]
    # Create kdtree search trees (for fast comparision).
    target_pcd_tree = cKDTree(np.array(target_atom_coords))

    # NOTE: target CA coords don't seem to be used
    # target_ca_coords = [atom.get_coord() for atom in target_struct.get_atoms() if atom.get_id() == 'CA']
    target_ca_pcd_tree = None # cKDTree(np.array(target_ca_coords))

    # If a specific residue is selected, then go after that residue
    if 'site_vix' in params and params['site_vix'] is not None:
        assert params['num_sites'] is None or params['num_sites'] == len(params['site_vix'])
        target_vertices = params['site_vix']

    else:
        assert params['num_sites'] is not None, "Please specify for how many sites the search should be performed."

        if 'target_residue' in params:
            target_chain = params['target_residue']['chain']
            # Use the tuple for biopython: (' ', resid, ' ')
            target_resid = [x.id for x in target_struct[0][target_chain].get_residues() if x.id[1] == params['target_residue']['resid']]
            assert len(target_resid) == 1, print(f"Target residue ID not unique: {target_resid}")
            target_resid = target_resid[0]
            print(f"Using residue: {target_resid}")

            target_cutoff = params['target_cutoff']
            if 'target_atom' in params:
                target_atom_id = params['target_atom']['atom_id']
                coord = target_struct[0][target_chain][target_resid][target_atom_id].get_coord()
                # find atom indices close to the target.
                dists = np.sqrt(np.sum(np.square(mymesh.vertices - coord), axis=1))
                neigh_indices = np.where(dists < target_cutoff)[0]

            else:
                residue_coords = np.stack([a.get_coord() for a in target_struct[0][target_chain][target_resid].get_atoms() if a.element != 'H'])
                dists = np.sqrt(np.sum(np.square(np.asarray(mymesh.vertices).reshape(-1, 1, 3) - residue_coords.reshape(1, -1, 3)), axis=-1))
                neigh_indices = np.where(np.any(dists < target_cutoff, axis=-1))[0]

            # Get a target vertex for every target site.
            target_vertices = get_target_vix(target_coord, iface, num_sites=params['num_sites'], selected_vertices=neigh_indices)

        elif 'target_point' in params:
            coord = np.array(params['target_point']['coord'])
            target_cutoff = params['target_cutoff']

            # find atom indices close to the target.
            dists = np.sqrt(np.sum(np.square(mymesh.vertices - coord), axis=1))
            neigh_indices = np.where(dists < target_cutoff)[0]
            # Get a target vertex for every target site.
            target_vertices = get_target_vix(target_coord, iface, num_sites=params['num_sites'], selected_vertices=neigh_indices)
        else:
            # Get a target vertex for every target site.
            target_vertices = get_target_vix(target_coord, iface, num_sites=params['num_sites'])

    if params.get("score_binder_transform", None) is not None:
        det_rotation = np.linalg.det(params["score_binder_transform"][:3, :3])
        tnorm = np.linalg.norm(params["score_binder_transform"][:3, 3])
        print(
            "Using --transform from {} (det(R)={:.6f}, |t|={:.3f})".format(
                params.get("score_binder_transform_source", "unknown"),
                det_rotation,
                tnorm,
            )
        )

    if params.get('score_binder', None) is not None:
        score_complex(
            target_name=target_ppi_pair_id,
            target_vertices=target_vertices,
            source_name=params['score_binder'],
            target_pcd=target_pcd,
            target_coord=target_coord,
            target_desc=target_desc,
            target_pcd_tree=target_pcd_tree,
            target_iface=target_iface,
            source_paths=source_paths,
            params=params,
            nn_score=nn_score_atomic,
            flip_target_normals=True,
        )
        return  # descriptor matching and alignment are not necessary for scoring

    outdir = os.path.join(params['out_dir'], params['target_name'])
    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)

    # Copy the pdb structure and the ply file of the target
    shutil.copy(target_pdb_path, outdir)
    shutil.copy(target_ply_fn, outdir)

    # Write out the targeted patches
    with open(os.path.join(outdir, 'selected_sites.vert'), 'w+') as out_patch:
        for site_vix in target_vertices:
            point = target_pcd.points[site_vix]
            out_patch.write('{}, {}, {}\n'.format(point[0], point[1], point[2]))


    # Go through every selected site
    for site_ix, site_vix in enumerate(target_vertices):
        site_outdir = os.path.join(outdir, 'site_{}'.format(site_ix))
        if not os.path.exists(site_outdir):
            os.makedirs(site_outdir, exist_ok=True)

        # Get the geodesic patch and descriptor patch for each target patch
        target_patch, target_patch_descs, target_patch_idx = get_patch_geo(
            target_pcd, target_coord, site_vix, target_desc,
            flip_normals=flip_target_normals,
            outward_shift=params['surface_outward_shift'])

        # Make a ckdtree with the target vertices.
        target_ckdtree = cKDTree(target_patch.points)

        # Write out the patch itself.
        with open(site_outdir + '/target.vert', 'w+') as out_patch:
            for point in target_patch.points:
                out_patch.write('{}, {}, {}\n'.format(point[0], point[1], point[2]))

        with open(site_outdir + '/target_info.txt', 'w') as out_info:
            out_info.write(f'name: {target_ppi_pair_id}, site: {site_ix}, vix: {site_vix}')

        # Match the top descriptors in the database based on descriptor distance.
        print('Starting to match target descriptor to descriptors from {} proteins; this may take a while.'.format(len(seed_ppi_pair_ids)))
        matched_dict, scores_dict = match_descriptors(seed_ppi_pair_ids, ['p1', 'p2'], target_desc[0][site_vix], params, return_scores=True)

        if len(matched_dict.keys())==0:
            continue

        # for (matched_name, _), matched_inds in matched_dict.items():
        #     matched_outdir = os.path.join(site_outdir, matched_name)
        #     os.makedirs(matched_outdir, exist_ok=True)
        #     with open(matched_outdir + '/first_stage_matches.txt', 'w') as f:
        #         f.write(' '.join(map(str, matched_inds)))

        print(" ")
        print("Second stage of MaSIF seed search: each matched descriptor is aligned and scored; this may take a while..")
        count_matched_fragments = 0
        for ix, name in enumerate(matched_dict.keys()):
            try:
                align_protein(
                    name, \
                    target_patch, \
                    target_patch_descs, \
                    target_ckdtree, \
                    target_ca_pcd_tree, \
                    target_pcd_tree, \
                    source_paths, \
                    matched_dict,\
                    nn_score_atomic, \
                    site_outdir, \
                    params, \
                    first_stage_scores=scores_dict[name],
                    n_retry_alignment=params['n_retry_alignment'],
                )
            except Exception as e:
                print(f"Error '{e}' while trying to align {name}.")
            if (ix + 1) % 1000 == 0:
                print('So far, MaSIF has aligned {} fragments from {} proteins.'.format(count_matched_fragments, ix + 1))
            count_matched_fragments += len(matched_dict[name])

    print("Done!")


if __name__ == "__main__":

    parser = ArgumentParser()
    parser.add_argument("--target", dest="target_name", type=str)
    parser.add_argument("--target_dir", dest="masif_target_root", type=Path)
    parser.add_argument("--database", dest="top_seed_dir", type=Path)
    parser.add_argument("--out_dir", type=Path, default=None)
    parser.add_argument("--resid", dest="target_resid", type=int, default=None)
    parser.add_argument("--chain", dest="target_chain", type=str, default=None)
    parser.add_argument("--atom_id", dest="target_atom_id", type=str, default=None)
    parser.add_argument("--coord", dest="target_coord", type=float, nargs="+", default=None)

    parser.add_argument("--subset", dest="database_subset", type=Path, default=None)
    parser.add_argument("--cutoff", dest="target_cutoff", type=float, default=10.0)
    parser.add_argument("--num_sites", type=int, default=None)
    parser.add_argument("--site_vix", type=int, nargs='+', default=None)
    parser.add_argument("--site_vix_file", type=Path, default=None)
    parser.add_argument("--desc_dist_cutoff", type=float, default=2.0, help="Recommended values: [1.5-2.0] (lower is stricter)")
    parser.add_argument("--iface_cutoff", type=float, default=0.75, help="Recommended values: [0.75-0.95] range (higher is stricter)")
    parser.add_argument("--nn_score_cutoff", type=float, default=0.8, help="# Recommended values: [0.8-0.95] (higher is stricter)")
    parser.add_argument("--desc_dist_score_cutoff", type=float, default=0.0, help="# Recommended values: [0.0-20.0] (higher is stricter)")
    parser.add_argument("--allowed_CA_clashes", type=int, default=0)
    parser.add_argument("--allowed_heavy_atom_clashes", type=int, default=5)
    parser.add_argument("--ransac_iter", type=int, default=100000)  # Note: previous default was 2000
    parser.add_argument("--n_retry_alignment", type=int, default=1)
    parser.add_argument("--sim", dest="similarity_mode", action="store_true")
    parser.add_argument("--score_binder", type=str, default=None, help="Specify a name of a processed protein to score it without alignment.")
    parser.add_argument(
        "--site_debug",
        type=Path,
        default=None,
        help="Write score-only query/binder patch .vert files into this directory.",
    )
    parser.add_argument(
        "--transform",
        type=str,
        default=None,
        help=(
            "Rigid transform for --score_binder. "
            "Accepted formats: path to .npy file (4x4 or 16 values) "
            "or inline comma-separated 16-value row-major string."
        ),
    )
    parser.add_argument("--random_seed", type=int, default=None)
    args = parser.parse_args()
    if args.site_debug is not None and args.score_binder is None:
        parser.error("--site_debug is only supported together with --score_binder.")
    if args.transform is not None and args.score_binder is None:
        parser.error("--transform is only supported together with --score_binder.")
    if args.transform is not None:
        try:
            transform, transform_source = parse_transform_arg(args.transform)
        except ValueError as exc:
            parser.error(str(exc))
        args.score_binder_transform = transform
        args.score_binder_transform_source = transform_source
    else:
        args.score_binder_transform = None
        args.score_binder_transform_source = None

    if args.random_seed is not None:
        assert version.parse('0.14.1') <= version.parse(o3d.__version__) <= version.parse('0.15.2'), "Random seed not supported by all Open3D versions"
        np.random.seed(args.random_seed)
        args.maybe_seed = {"seed": args.random_seed}

    # Definition of the target patch
    if args.target_coord is not None:
        args.target_point = {'coord': args.target_coord}
    elif args.target_resid is not None and args.target_chain is not None:
        args.target_residue = {'resid': args.target_resid, 'chain': args.target_chain}
        if args.target_atom_id is not None:
            args.target_atom = {'atom_id': args.target_atom_id}

    if args.site_vix_file is not None and args.site_vix is None:
        target_vix = [int(line) for line in args.site_vix_file.read_text().strip().split('\n')]
        args.site_vix = target_vix

    # Database locations
    args.seed_surf_dir = os.path.join(args.top_seed_dir, masif_opts['ply_chain_dir'])
    args.seed_iface_dir = os.path.join(args.top_seed_dir, masif_opts['site']['out_pred_dir'])
    args.seed_ply_iface_dir = os.path.join(args.top_seed_dir, masif_opts['site']['out_surf_dir'])
    args.seed_pdb_dir = os.path.join(args.top_seed_dir, masif_opts['pdb_chain_dir'])
    args.seed_desc_dir = os.path.join(args.top_seed_dir, masif_opts['ppi_search']['desc_dir'])

    # 9 A
    # args.seed_precomp_dir = os.path.join(args.top_seed_dir,masif_opts['site']['masif_precomputation_dir'])
    # 12 A
    args.seed_precomp_dir = os.path.join(args.top_seed_dir, masif_opts['ppi_search']['masif_precomputation_dir'])

    # Database subset
    if args.database_subset is None:
        # use all available protein IDs in the database
        args.seed_ppi_pair_ids = np.array(os.listdir(args.seed_desc_dir))
    else:
        # use a predefined subset of protein IDs
        with open(args.database_subset, "r") as f:
            args.seed_ppi_pair_ids = [x.rstrip() for x in f.readlines()]

    # Target locations
    args.top_target_dir = os.path.join(args.masif_target_root)
    args.target_surf_dir = os.path.join(args.top_target_dir, masif_opts['ply_chain_dir'])
    args.target_iface_dir = os.path.join(args.masif_target_root, masif_opts['site']['out_pred_dir'])
    args.target_ply_iface_dir = os.path.join(args.masif_target_root, masif_opts['site']['out_surf_dir'])
    args.target_pdb_dir = os.path.join(args.top_target_dir, masif_opts['pdb_chain_dir'])
    args.target_desc_dir = os.path.join(args.top_target_dir, masif_opts['ppi_search']['desc_dir'])
    args.target_desc_dir = os.path.join(args.top_target_dir, masif_opts['ppi_search']['desc_dir'])

    # 9 A
    # args.target_precomp_dir = os.path.join(args.top_target_dir,masif_opts['site']['masif_precomputation_dir'])
    # 12 A
    args.target_precomp_dir = os.path.join(args.top_target_dir, masif_opts['ppi_search']['masif_precomputation_dir'])

    # Some hard-coded parameters
    # Neural network scores.
    args.nn_score_atomic_fn = os.path.join(masif_neosurf_dir, "masif_seed_search/data/scoring_nn/models_std/weights_12A_0129")
    args.max_npoints = 200

    # Ransac parameters
    if version.parse(o3d.__version__) > version.parse('0.11.0'):
        args.ransac_convergence_kwargs = {'confidence': 0.999}
    # Ransac radius - should not be changed.
    args.ransac_radius = 1.5
    # How much to expand the surface for alignment.
    args.surface_outward_shift = 0.25

    masif_search(vars(args))
