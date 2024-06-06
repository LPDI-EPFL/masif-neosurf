from geometry.open3d_import import *
import scipy.linalg
import networkx as nx 
from scipy.spatial import cKDTree
import copy 
from IPython.core.debugger import set_trace
import numpy as np
from simple_mesh import Simple_mesh
from pathlib import Path
import scipy.sparse as spio
import pymesh
from Bio.PDB import PDBParser, PDBIO, Selection
import tensorflow as tf
from tensorflow import keras
import os

def rand_rotation_matrix(deflection=1.0, randnums=None):
    """
    Creates a random rotation matrix. used to randomize initial pose. 
    deflection: the magnitude of the rotation. For 0, no rotation; for 1, competely random
    rotation. Small deflection => small perturbation.
    randnums: 3 random numbers in the range [0, 1]. If `None`, they will be auto-generated.
    """
    # from http://www.realtimerendering.com/resources/GraphicsGems/gemsiii/rand_rotation.c

    if randnums is None:
        randnums = np.random.uniform(size=(3,))

    theta, phi, z = randnums

    theta = theta * 2.0 * deflection * np.pi  # Rotation about the pole (Z).
    phi = phi * 2.0 * np.pi  # For direction of pole deflection.
    z = z * 2.0 * deflection  # For magnitude of pole deflection.

    # Compute a vector V used for distributing points over the sphere
    # via the reflection I - V Transpose(V).  This formulation of V
    # will guarantee that if x[1] and x[2] are uniformly distributed,
    # the reflected points will be uniform on the sphere.  Note that V
    # has length sqrt(2) to eliminate the 2 in the Householder matrix.

    r = np.sqrt(z)
    Vx, Vy, Vz = V = (np.sin(phi) * r, np.cos(phi) * r, np.sqrt(2.0 - z))

    st = np.sin(theta)
    ct = np.cos(theta)

    R = np.array(((ct, st, 0), (-st, ct, 0), (0, 0, 1)))

    # Construct the rotation matrix  ( V Transpose(V) - I ) R.

    M = (np.outer(V, V) - np.eye(3)).dot(R)
    return M

def get_center_and_random_rotate(pcd):
    """
        Get the center of a point cloud and randomly rotate it.
        pcd: the point cloud.
    """
    pts = pcd.points
    mean_pt = np.mean(pts, axis=0)
    # pts = pts - mean_pt
    rand_mat = rand_rotation_matrix()
    # pts = Vector3dVector(np.dot(pts,rand_mat))
    transform = np.vstack([rand_mat.T, -mean_pt]).T
    # transform = np.vstack([np.diag([1,1,1]),-mean_pt]).T
    transform = np.vstack([transform, [0, 0, 0, 1]])
    return transform

# Apply transformation to source struct
def random_rotate_source_struct(source_struct, transformation):
    # Randomly rotate the source struct 
    structure_atoms = [atom for atom in source_struct.get_atoms()]
    structure_coords = np.array([atom.get_coord() for atom in structure_atoms])
    structure_coord_pcd = PointCloud()
    structure_coord_pcd.points = Vector3dVector(structure_coords)
    structure_coord_pcd_notTransformed = copy.deepcopy(structure_coord_pcd)
    structure_coord_pcd.transform(transformation)
    # - This is a bit of a bad programming practice: shouldn't be a for loop..
    for ix, v in enumerate(structure_coord_pcd.points):
        structure_atoms[ix].set_coord(v)


def load_protein_pcd(full_pdb_id, chain_number, paths, flipped_features=False, read_mesh=True):
    """
    Returns the protein surface point cloud, list of descriptors, interface predictions
    and patch coordinates given full pdb name (str), chain number (int) and if descriptors
    should be flipped or not (bool).
    """

    full_pdb_id_split = full_pdb_id.split('_')
    pdb_id = full_pdb_id_split[0]
    chain_ids = full_pdb_id_split[1:]
    pdb_chain = chain_ids[chain_number - 1]

    # Read the protein surface. 
    pdb_pcd = read_point_cloud(
        str(Path(paths['surf_dir']) /
            '{}_{}.ply'.format(
            pdb_id,
            pdb_chain)))

    if read_mesh == True:
        # Read both the vertices and their connection (i.e. the mesh)
        pdb_mesh = read_triangle_mesh(
            str(Path(paths['surf_dir']) /
            '{}_{}.ply'.format(
            pdb_id,
            pdb_chain))) 

    # Read the masif-site predictions.
    pdb_iface = np.squeeze(np.load(
        Path(paths['iface_dir']) /
        'pred_{}_{}.npy'.format(
            pdb_id,
            pdb_chain)))

    if flipped_features:
        pdb_desc_sc05 = np.load(
            Path(paths['desc_dir']) /
            full_pdb_id /
            'p{}_desc_flipped.npy'.format(chain_number))
    else:
        pdb_desc_sc05 = np.load(
            Path(paths['desc_dir']) /
            full_pdb_id /
            'p{}_desc_straight.npy'.format(chain_number))


    pdb_desc = np.array([pdb_desc_sc05,])

    if read_mesh:
        return pdb_pcd, pdb_desc, pdb_iface, pdb_mesh
    else:
        return pdb_pcd, pdb_desc, pdb_iface 

def get_patch_geo(
        pcd,
        patch_coords,
        center,
        descriptors,
        outward_shift=0.25,
        flip_normals=False):
    """
    Returns a patch from a point cloud pcd with center point center (int),
    based on geodesic distances from patch coords and corresponding Feature descriptors.
    """
    patch_idxs = patch_coords[center]
    patch_pts = np.asarray(pcd.points)[patch_idxs, :]
    patch_nrmls = np.asarray(pcd.normals)[patch_idxs, :]
    patch_pts = patch_pts + outward_shift * patch_nrmls
    if flip_normals:
        patch_nrmls = -patch_nrmls

    patch = PointCloud()
    patch.points = Vector3dVector(patch_pts)
    patch.normals = Vector3dVector(patch_nrmls)
    patch_descs = [Feature(), Feature(), Feature()]
    patch_descs[0].data = descriptors[0,patch_idxs, :].T
    return patch, patch_descs, patch_idxs


def get_patch_coords(top_dir, pdb, pid, cv=None):
    """ 
    Load precomputed patch coordinates.
    """
    if cv is None:
        patch_coords = np.load(os.path.join(
            top_dir, pdb, pid+'_list_indices.npy'), allow_pickle=True)
        cv = np.arange(0, patch_coords.shape[0])
    else:
        patch_coords = np.load(os.path.join(
            top_dir, pdb, pid+'_list_indices.npy'), allow_pickle=True)[cv]
    patch_coords = {key: patch_coords[ii] for ii, key in enumerate(cv)}
    return patch_coords 


# Get the center of the interface based on the ground truth: find the most shape complementary patch. 
def geodists(verts, faces):
    # Graph 
    G=nx.Graph()
    n = len(verts)
    G.add_nodes_from(np.arange(n))

    # Get edges
    f = np.array(faces, dtype = int)
    rowi = np.concatenate([f[:,0], f[:,0], f[:,1], f[:,1], f[:,2], f[:,2]], axis = 0)
    rowj = np.concatenate([f[:,1], f[:,2], f[:,0], f[:,2], f[:,0], f[:,1]], axis = 0)
    edges = np.stack([rowi, rowj]).T
    
    # Get weights 
    edgew = verts[rowi] - verts[rowj]
    edgew = scipy.linalg.norm(edgew, axis=1)
    wedges = np.stack([rowi, rowj, edgew]).T
    G.add_weighted_edges_from(wedges)
    dists = nx.all_pairs_dijkstra_path_length(G, cutoff=12.0)
    d2 = {}
    for key_tuple in dists:
        d2[key_tuple[0]] = key_tuple[1]
    return d2

def get_target_vix_sc(target_mesh, binder_mesh):
    binder_vert_ckd = cKDTree(binder_mesh.vertices)
    dists, ckd_results = binder_vert_ckd.query(target_mesh.vertices)
    n1 = np.array([target_mesh.get_attribute('vertex_nx'), target_mesh.get_attribute('vertex_ny'), target_mesh.get_attribute('vertex_nz')]).T
    n2 = np.array([binder_mesh.get_attribute('vertex_nx'), binder_mesh.get_attribute('vertex_ny'), binder_mesh.get_attribute('vertex_nz')]).T
    v1 = target_mesh.vertices
    v2 = binder_mesh.vertices

    geodists_target = geodists(target_mesh.vertices, target_mesh.faces)
    w = 0.5
    num_rings = 10
    radius =12
    scales = np.arange(0, radius, radius/10)
    n = len(target_mesh.vertices)
    v1_sc_25 = np.zeros((n, 10))
    v1_sc_50 = np.zeros((n, 10))

    for vix in range(len(target_mesh.vertices)):
        if dists[vix] <= 1.5:
            # Compute sc
            p1ix = np.array([int(x) for x in geodists_target[vix].keys()])
            patch_dists = [geodists_target[vix][x] for x in p1ix]
            p2ix = ckd_results[p1ix]

            # First v1->v2
            comp1 = [np.dot(n1[p1ix[x]], -n2[p2ix[x]]) for x in range(len(p1ix))]
            comp1 = np.multiply(comp1, np.exp(-w * np.square(dists[p1ix])))
            # Use 10 rings such that each ring has equal weight in shape complementarity
            comp_rings1_25 = np.zeros(num_rings)
            comp_rings1_50 = np.zeros(num_rings)

            for ring in range(num_rings-1):
                scale = scales[ring]
                members = [int(x) for x in range(len(p1ix)) if patch_dists[x] >= scales[ring] and \
                                                                        patch_dists[x] < scales[ring+1]]
                if len(members) == 0:
                    comp_rings1_25[ring] = 0.0
                    comp_rings1_50[ring] = 0.0
                else:
                    comp_rings1_25[ring] = np.percentile(comp1[members], 25)
                    comp_rings1_50[ring] = np.percentile(comp1[members], 50)
            
            v1_sc_25[vix] = comp_rings1_25
            v1_sc_50[vix] = comp_rings1_50

    #center_point = np.argmax(np.median(np.nan_to_num(v1_sc_25), axis=1))
    center_point = np.argmax(np.median(np.nan_to_num(v1_sc_50), axis=1))
    return center_point

# Get the vertex indices of target sites.
def get_target_vix(pc, iface, num_sites=1, selected_vertices=None):
    target_vertices = []
    iface = np.copy(iface)
    if selected_vertices is None:
        selected_vertices = np.arange(len(iface))
    for site_ix in range(num_sites):
        iface_patch_vals = []
        # Go through each patch.
        best_ix = -1
        best_val = float("-inf")
        best_neigh = []
        for ii in selected_vertices:
            neigh = pc[ii]
            val = np.mean(iface[neigh])
            if val > best_val:
                best_ix = ii
                best_val = val
                best_neigh = neigh

        # Now that a site has been identified, clear the iface values so that the same site is not picked again.
        iface[best_neigh] *= 0.5
        target_vertices.append(best_ix)

    return target_vertices

def test_alignments(transformation, source_structure, target_pcd_tree, interface_dist=10.0):
    """
        test_alignments - evaluate quality of the alignments.
    """
    structure_coords = np.array([atom.get_coord() for atom in source_structure.get_atoms()])
    structure_coord_pcd = PointCloud()
    structure_coord_pcd.points = Vector3dVector(structure_coords)
    structure_coord_pcd.transform(transformation)
        
    d_nn_interface, i_nn_interface = target_pcd_tree.query(np.asarray(structure_coord_pcd.points),k=1,distance_upper_bound=interface_dist)
    interface_atoms = np.where(d_nn_interface<=interface_dist)[0]
    rmsd = np.sqrt(np.mean(np.square(np.linalg.norm(structure_coords[interface_atoms,:]-np.asarray(structure_coord_pcd.points)[interface_atoms,:],axis=1))))
    return rmsd

def match_descriptors(directory_list, pids, target_desc, params):
    """ 
        Match descriptors to the target descriptor.
    """

    all_matched_names = []
    all_matched_vix = []
    all_matched_desc_dist = []
    count_proteins = 0
    count_descriptors = 0

    for ppi_pair_id in directory_list:
        if 'seed_pdb_list' in params and ppi_pair_id not in params['seed_pdb_list']:
            continue

        if '.npy' in ppi_pair_id or '.txt' in ppi_pair_id:
            continue

        mydescdir = os.path.join(params['seed_desc_dir'], ppi_pair_id)
        for pid in pids:
            try:
                fields = ppi_pair_id.split('_')
                if pid == 'p1':
                    pdb_chain_id = fields[0]+'_'+fields[1]
                elif pid == 'p2':
                    pdb_chain_id = fields[0]+'_'+fields[2]
                iface = np.load(params['seed_iface_dir']+'/pred_'+pdb_chain_id+'.npy')[0]
                descs = np.load(mydescdir+'/'+pid+'_desc_straight.npy')
            except:
                continue
            #print(pdb_chain_id)
            name = (ppi_pair_id, pid)
            count_proteins += 1
            count_descriptors += len(iface)

            diff = np.sqrt(np.sum(np.square(descs - target_desc), axis=1))

            true_iface = np.where(iface > params['iface_cutoff'])[0]
            near_points = np.where(diff < params['desc_dist_cutoff'])[0]

            selected = np.intersect1d(true_iface, near_points)
            if len(selected > 0):
                all_matched_names.append([name]*len(selected))
                all_matched_vix.append(selected)
                all_matched_desc_dist.append(diff[selected])
                #print('Matched {}'.format(ppi_pair_id))
                #print('Scores: {} {}'.format(iface[selected], diff[selected]))
            if count_proteins % 1000 == 0:
                print('Compared target with {} \'fragments\' from {} proteins'.format(count_descriptors, count_proteins))


    # Flatten the matches and convert to a dictionary 
    try:
        matched_names = np.concatenate(all_matched_names, axis=0)
        matched_vix = np.concatenate(all_matched_vix, axis=0)
        print('Iterated over {} fragments from {} proteins; matched {} based on descriptor similarity.'.format(count_descriptors, count_proteins, len(matched_vix)))
        matched_desc_dist = np.concatenate(all_matched_desc_dist, axis=0)
    except:
        print("matched no descriptors")
        return {}
    matched_dict = {}
    for name_ix, name in enumerate(matched_names):
        name = (name[0], name[1])
        if name not in matched_dict:
            matched_dict[name] = []
        matched_dict[name].append(matched_vix[name_ix])
            
    return matched_dict

def count_clashes(transformation, source_surface_vertices, source_structure, \
        target_ca_pcd_tree, target_pcd_tree, radius=2.0, clashing_ca_thresh=1.0, clashing_thresh=5.0):
    """
        count_clashes and align if the threshold passes. 
    """

    structure_ca_coords = np.array([atom.get_coord() for atom in source_structure.get_atoms() if atom.get_id() == 'CA'])
    structure_ca_coord_pcd = PointCloud()
    structure_ca_coord_pcd.points = Vector3dVector(structure_ca_coords)
    structure_ca_coord_pcd_notTransformed = copy.deepcopy(structure_ca_coord_pcd)
    structure_ca_coord_pcd.transform(transformation)
    d_nn_ca, i_nn_ca = target_pcd_tree.query(np.asarray(structure_ca_coord_pcd.points),k=1,distance_upper_bound=radius)
    clashing_ca = np.sum(d_nn_ca<=radius)
    if clashing_ca > clashing_ca_thresh:
        return clashing_ca, 0.0

    structure_atoms = [atom for atom in source_structure.get_atoms() if not atom.get_name().startswith('H')]
    structure_coords = np.array([atom.get_coord() for atom in structure_atoms])
    structure_coord_pcd = PointCloud()
    structure_coord_pcd.points = Vector3dVector(structure_coords)
    structure_coord_pcd_notTransformed = copy.deepcopy(structure_coord_pcd)
    structure_coord_pcd.transform(transformation)
    # Invoke cKDTree with (structure_coord_pcd)
    d_nn, i_nn = target_pcd_tree.query(np.asarray(structure_coord_pcd.points),k=1,distance_upper_bound=radius)
    clashing = np.sum(d_nn<=radius)

    if clashing > clashing_thresh:
        return clashing_ca, clashing 

    # Transform the input structure. - This is a bit of a bad programming practice: shouldn't be a for loop..
    for ix, v in enumerate(structure_coord_pcd.points):
        structure_atoms[ix].set_coord(v)

    return clashing_ca,clashing


def multidock(source_pcd, source_patch_coords, source_descs, 
            cand_pts, target_pcd, target_descs,
               params):
    ransac_radius=params['ransac_radius'] 
    ransac_iter=params['ransac_iter']
    all_results = []
    all_source_patch = []
    all_source_scores = []
    all_source_desc = []
    all_source_idx = []
    for pt in cand_pts:
        source_patch, source_patch_descs, source_patch_idx = get_patch_geo(
            source_pcd, source_patch_coords, pt, source_descs, outward_shift=params['surface_outward_shift'])
           
        result = registration_ransac_based_on_feature_matching(
            source_patch, target_pcd, source_patch_descs[0], target_descs[0],
            ransac_radius,
            TransformationEstimationPointToPoint(False), 3,
            [CorrespondenceCheckerBasedOnEdgeLength(0.9),
             CorrespondenceCheckerBasedOnDistance(1.0),
             CorrespondenceCheckerBasedOnNormal(np.pi/2)],
            RANSACConvergenceCriteria(ransac_iter, 500)
        )
        ransac_transformation = result.transformation 
        
        # TODO: there is a potential bug here in benchmark cases only. If a random rotation is not previously applied, 
        # there exists a possibility that masif-search/masif-site will get the right patches, but ransac will fail to 
        # get a transformation. In these rare cases, icp will start from the ground truth. To get around this, make
        # sure a rotation is applied beforehand (for benchmarks at least)
        result_icp = registration_icp(source_patch, target_pcd,
                    1.0, result.transformation, TransformationEstimationPointToPlane())

        source_patch.transform(result_icp.transformation)
        all_results.append(result_icp)
        all_source_patch.append(source_patch)
        all_source_desc.append(source_patch_descs)
        all_source_idx.append(source_patch_idx)

    return all_results, all_source_patch, all_source_desc, all_source_idx 


def align_and_save(out_filename_base, patch, source_structure, 
                   point_importance, 
                   ):
    
    io = PDBIO()
    # Remove hydrogens.
    for atom in Selection.unfold_entities(source_structure, 'A'):
        if atom.get_name().startswith('H'):
            parent = atom.get_parent()
            parent.detach_child(atom.get_id())
    io.set_structure(source_structure)
    io.save(out_filename_base+'.pdb')
    # Save patch - not necessary for benchmarking purposes. Uncomment if wanted.
    #mesh = Simple_mesh(np.asarray(patch.points))
    #mesh.set_attribute('vertex_charge', point_importance)
    #mesh.save_mesh(out_filename_base+'_patch.ply')
    

# Compute different types of scores: 
# -- Compute the neural network score
def compute_nn_score(target_pcd, source_pcd, corr, 
        target_desc, source_desc, 
        target_ckdtree, nn_score, 
        vertex_to_atom_dist, nn_score_cutoff=1.0):


    # Compute nn scores 
    # Compute all points correspondences and distances for nn
    d, r = target_ckdtree.query(source_pcd.points)

    # r: for every point in source, what is its correspondence in target
    # feat0 - distance between poitnts
    distance = d # distances between points 
    distance[distance<0.5] = 0.5
    distance = 1.0/distance

    # feat1: descriptors
    desc_dist = target_desc[0].data[:,r] - source_desc[0].data[:,:]
    desc_dist = np.sqrt(np.sum(np.square(desc_dist.T), axis=1))
    desc_dist = 1.0/desc_dist
    # desc dist score: sum over those within 1.5A
    neigh = np.where(d < 1.5)[0]
    desc_dist_score = np.sum(np.square(desc_dist[neigh]))

    # feat3: normal dot product
    n1 = np.asarray(source_pcd.normals)
    n2 = np.asarray(target_pcd.normals)[r]
    normal_dp = np.multiply(n1, n2).sum(1)

    # feat9 - vertex to atom distance
    vert_atom_dist = vertex_to_atom_dist
    vert_atom_dist[vert_atom_dist < 0.25] = 0.25
    vert_atom_dist = 1.0/vert_atom_dist

    features = np.vstack([distance, desc_dist, normal_dp, vert_atom_dist]).T

    nn_score_pred, point_importance = \
                    nn_score.eval_model(
#                            features, nn_score_cutoff
                            features, 0.9 
                            )

    ret = (np.array([nn_score_pred[0][0], desc_dist_score]).T, point_importance)
    return ret 

def align_protein(name, \
        target_patch, \
        target_patch_descs, \
        target_ckdtree, \
        target_ca_pcd_tree, \
        target_pcd_tree, \
        source_paths, \
        matched_dict, \
        nn_score, \
        site_outdir, \
        params):
    ppi_pair_id = name[0]
    pid = name[1]
    pdb = ppi_pair_id.split('_')[0]

    # PDB Parser
    parser = PDBParser()

    if pid == 'p1':
        chain = ppi_pair_id.split('_')[1]
        chain_number = 1
    else: 
        chain = ppi_pair_id.split('_')[2]
        chain_number = 2
        
    # Load source ply file, coords, and descriptors.
    source_pcd, source_desc, source_iface = load_protein_pcd(ppi_pair_id, chain_number, source_paths, flipped_features=False, read_mesh=False)
    
    # Randomly rotate the source ply file, and store the random transformation matrix (for benchmark purposes only)
    random_transformation = get_center_and_random_rotate(source_pcd)
    source_pcd.transform(random_transformation)
    

    # Get coordinates for all matched vertices.
    source_vix = matched_dict[name]

    source_coord =  get_patch_coords(params['seed_precomp_dir'], ppi_pair_id, pid, cv=source_vix)
    
    # Perform all alignments to target. 
    all_results, all_source_patch, all_source_patch_desc, all_source_idx = multidock(
            source_pcd, source_coord, source_desc,
            source_vix, target_patch, target_patch_descs, 
            params
            ) 

    # Score the results using a 'lightweight' scoring function.
    all_source_scores = [None]*len(source_vix)
    for viii in range(len(source_vix)):
        if len(all_results[viii].correspondence_set)/float(len(np.asarray(all_source_patch[viii].points))) < 0.3:
            # Ignore those with poor fitness.
            all_source_scores[viii] = ([0,0], np.zeros_like(all_source_idx[viii]))
        else:
            # Compute the distance between every source_surface_vertices and every target vertex.
            d_vi_at, _= target_pcd_tree.query(np.asarray(all_source_patch[viii].points), k=1)

            all_source_scores[viii] = compute_nn_score(target_patch, all_source_patch[viii], np.asarray(all_results[viii].correspondence_set),
                 target_patch_descs, all_source_patch_desc[viii], \
                target_ckdtree, nn_score, d_vi_at, 1.0) # Ignore point importance for speed.
            if len(np.asarray(all_results[viii].correspondence_set)) <= 2.0:
                if all_source_scores[viii][0][0] > 0.9: 
                    print('Error in masif_seed_search; check scoring.')
                    return

    # All_point_importance: impact of each vertex on the score. 
    all_point_importance = [x[1] for x in all_source_scores]
    all_source_scores = [x[0] for x in all_source_scores]
    scores = np.asarray(all_source_scores)
    
    # Filter anything below neural network score cutoff
    top_scorers = np.where(scores[:,0] > params['nn_score_cutoff'])[0]
    
    if len(top_scorers) > 0:
        source_outdir = os.path.join(site_outdir, '{}'.format(ppi_pair_id))
        if not os.path.exists(source_outdir):
            os.makedirs(source_outdir)

        # Load source structure 
        source_struct = parser.get_structure('{}_{}'.format(pdb,chain), os.path.join(params['seed_pdb_dir'],'{}_{}.pdb'.format(pdb,chain)))
        # Use the preexisting random rotation matrix that was applied to the patch.
        random_rotate_source_struct(source_struct, random_transformation)


        # Perform the transformation on the atoms
        for j in top_scorers:
            res = all_results[j]

            source_structure_cp = copy.deepcopy(source_struct)
            
            # Count clashes and if threshold passes, align the pdb.
            clashing_ca, clashing_total = count_clashes(res.transformation, np.asarray(all_source_patch[j].points), source_structure_cp, \
                target_ca_pcd_tree, target_pcd_tree, clashing_ca_thresh=params['allowed_CA_clashes'], clashing_thresh=params['allowed_heavy_atom_clashes'])

            # Check if the number of clashes exceeds the number allowed. 
            if clashing_ca <= params['allowed_CA_clashes'] and clashing_total <= params['allowed_heavy_atom_clashes']:

                # Align and save the pdb + patch 
                out_fn = source_outdir+'/{}_{}_{}'.format(pdb, chain, j)
                align_and_save(out_fn, all_source_patch[j], source_structure_cp, \
                                            point_importance=all_point_importance[j])

                # Recompute the score with clashes.
                print('Selected fragment: {} fragment_id: {} score: {:.4f} desc_dist_score: {:.4f} clashes(CA): {} clashes(total):{}\n'.format(j , ppi_pair_id, scores[j][0], scores[j][1], clashing_ca, clashing_total))

                # Output the score for convenience. 
                out_score = open(out_fn+'.score', 'w+')
                out_score.write('name: {}, point id: {}, score: {:.4f}, clashing_ca: {}, clashing_heavy: {}, desc_dist_score: {}\n'.format(ppi_pair_id, j, scores[j][0], clashing_ca,clashing_total, scores[j][1]))
                out_score.close()
                
                print ('Surface alignment is deactivated - aligning only pdbs. ')
                continue

                # Align and save the ply file for convenience.     
                mesh = Simple_mesh()
                mesh.load_mesh(os.path.join(params['seed_surf_dir'],'{}.ply'.format(pdb+'_'+chain)))
                
                source_pcd_copy = copy.deepcopy(source_pcd)
                source_pcd_copy.transform(res.transformation)
                out_vertices = np.asarray(source_pcd_copy.points)
                out_normals = np.asarray(source_pcd_copy.normals)
                mesh.set_attribute('vertex_x', out_vertices[:,0])
                mesh.set_attribute('vertex_y', out_vertices[:,1])
                mesh.set_attribute('vertex_z', out_vertices[:,2])
                mesh.set_attribute('vertex_nx', out_normals[:,0])
                mesh.set_attribute('vertex_ny', out_normals[:,1])
                mesh.set_attribute('vertex_nz', out_normals[:,2])
                mesh.vertices = out_vertices
                mesh.set_attribute('vertex_iface', source_iface) 
                mesh.save_mesh(out_fn+'.ply')
                #save_ply(out_fn+'.ply', out_vertices, mesh.faces, out_normals, charges=mesh.get_attribute('vertex_charge'))
