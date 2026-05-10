import numpy as np
import networkx as nx
from pymatgen.core import Element
from pymatgen.core.structure import Structure
from pymatgen.util import coord as pmg_coord


def _split_hydrogens(ma_indices, struct, graph):
    c_index = [i for i in ma_indices if struct[i].specie.symbol == "C"][0]
    n_index = [i for i in ma_indices if struct[i].specie.symbol == "N"][0]
    h_indices = [i for i in ma_indices if struct[i].specie.symbol == "H"]

    h_inds_c, h_inds_n = [], []
    for h in h_indices:
        nbrs = set(graph.neighbors(h))
        if c_index in nbrs:
            h_inds_c.append(h)
        if n_index in nbrs:
            h_inds_n.append(h)

    dist_c = {h: struct.get_distance(c_index, h) for h in h_indices}
    dist_n = {h: struct.get_distance(n_index, h) for h in h_indices}

    if len(h_inds_c) == 3 and len(h_inds_n) == 3 and set(h_inds_c).isdisjoint(h_inds_n):
        return sorted(h_inds_c), sorted(h_inds_n)

    h_inds_c, h_inds_n = [], []
    for h in sorted(h_indices, key=lambda x: abs(dist_c[x] - dist_n[x])):
        prefer_c = dist_c[h] <= dist_n[h]
        if prefer_c and len(h_inds_c) < 3:
            h_inds_c.append(h)
        elif not prefer_c and len(h_inds_n) < 3:
            h_inds_n.append(h)
        elif len(h_inds_c) < 3:
            h_inds_c.append(h)
        else:
            h_inds_n.append(h)

    h_inds_c = h_inds_c[:3]
    remaining = [h for h in h_indices if h not in h_inds_c]
    h_inds_n = remaining[:3]
    return sorted(h_inds_c), sorted(h_inds_n)


def _shortest_vec(lattice, frac_coords, i, j):
    return pmg_coord.pbc_shortest_vectors(lattice, frac_coords[i], frac_coords[j])[0, 0, :]


def _pbc_distance(lattice, frac_coords, i, j):
    return np.linalg.norm(_shortest_vec(lattice, frac_coords, i, j))


def _closest_along_direction(lattice, frac_coords, origin_idx, dirvec, shift_val, indices):
    if len(indices) == 0:
        return None
    normdir = dirvec / np.linalg.norm(dirvec)
    origin_cart = lattice.get_cartesian_coords(frac_coords[origin_idx])
    target = origin_cart + normdir * shift_val
    target_frac = lattice.get_fractional_coords(target)
    best = None
    for idx in indices:
        vec = pmg_coord.pbc_shortest_vectors(lattice, target_frac, frac_coords[idx])[0, 0, :]
        dist = np.linalg.norm(vec)
        if best is None or dist < best[0]:
            best = (dist, idx)
    return None if best is None else best[1]


def _build_br_axes(lattice, frac_coords, origin_idx, axis_vecs, br_indices, br_neighbor_candidates):
    br_axis, br_axis_neg = [], []
    for dirvec in axis_vecs:
        pos_idx = _closest_along_direction(
            lattice, frac_coords, origin_idx, dirvec, np.linalg.norm(dirvec) / 2.0, br_indices
        )
        if pos_idx is not None:
            br_axis.append(int(pos_idx))

        neg_idx = _closest_along_direction(
            lattice, frac_coords, origin_idx, -dirvec, np.linalg.norm(dirvec) / 2.0, br_neighbor_candidates
        )
        if neg_idx is not None:
            br_axis_neg.append(int(neg_idx))
    return br_axis, br_axis_neg


def _get_pb_ortho_axes(
    lattice,
    frac_coords,
    origin_idx,
    axis_vecs,
    pb_axis,
    pb_indices_second,
    pb_indices_fourth,
    shift_val,
    shift_val_pc=None,
):
    reqvecs = (
        (np.array([0.0, 1.0, -1.0]), np.array([0.0, 1.0, 1.0])),
        (np.array([-1.0, 0.0, 1.0]), np.array([1.0, 0.0, 1.0])),
        (np.array([1.0, 1.0, 0.0]), np.array([-1.0, 1.0, 0.0])),
    )

    ortho_axes = []
    for i in range(3):
        base = pb_axis[1][i]
        dir1 = (
            reqvecs[i][0][0] * axis_vecs[0]
            + reqvecs[i][0][1] * axis_vecs[1]
            + reqvecs[i][0][2] * axis_vecs[2]
        )
        dir2 = (
            reqvecs[i][1][0] * axis_vecs[0]
            + reqvecs[i][1][1] * axis_vecs[1]
            + reqvecs[i][1][2] * axis_vecs[2]
        )

        idx1 = _closest_along_direction(
            lattice, frac_coords, origin_idx, dir1, shift_val, pb_indices_second
        )
        idx2 = _closest_along_direction(
            lattice, frac_coords, origin_idx, dir2, shift_val, pb_indices_second
        )

        axis0 = base
        if shift_val_pc is not None and len(pb_indices_fourth) > 0:
            axis0_candidate = _closest_along_direction(
                lattice,
                frac_coords,
                origin_idx,
                axis_vecs[i],
                shift_val_pc,
                pb_indices_fourth,
            )
            if axis0_candidate is not None:
                axis0 = axis0_candidate

        axis_list = [int(axis0)]
        if idx1 is not None:
            axis_list.append(int(idx1))
        if idx2 is not None:
            axis_list.append(int(idx2))
        ortho_axes.append(axis_list)

    return (int(origin_idx), ortho_axes)


def _get_pb_inds_orig_env(lattice, frac_coords, c_ind, n_ind, pb_indices, cutoff=6.0):
    c_pos = lattice.get_cartesian_coords(frac_coords[c_ind])
    n_pos = lattice.get_cartesian_coords(frac_coords[n_ind])
    com_frac = lattice.get_fractional_coords((c_pos + n_pos) / 2.0)

    pb_vecs = []
    pb_dists = []
    for idx in pb_indices:
        vec = pmg_coord.pbc_shortest_vectors(lattice, com_frac, frac_coords[idx])[0, 0, :]
        pb_vecs.append(vec)
        pb_dists.append(np.linalg.norm(vec))
    pb_vecs = np.array(pb_vecs)
    pb_dists = np.array(pb_dists)

    pb_env_local = np.where(pb_dists < cutoff)[0]
    if len(pb_env_local) == 0:
        raise ValueError("No Pb neighbors found within cutoff")

    neg_candidates = [idx for idx in pb_env_local if np.all(pb_vecs[idx] < 0.0)]
    if len(neg_candidates) == 0:
        origin_local = pb_env_local[0]
    else:
        origin_local = neg_candidates[int(np.argmin(pb_vecs[neg_candidates][:, 0]))]

    pb_env_local = [i for i in pb_env_local if i != origin_local]
    return pb_indices[origin_local], [pb_indices[i] for i in pb_env_local]


def _build_chn_graph(struct):
    species = struct.symbol_set

    def _covalent(sym):
        radius = getattr(Element(sym), "average_covalent_radius", None)
        return float(radius) if radius is not None else 0.7

    cov_r = {sym: _covalent(sym) for sym in species}

    def bond_cutoff(s1, s2):
        return 1.2 * (cov_r[s1] + cov_r[s2])

    chn_indices = [i for i, site in enumerate(struct) if site.specie.symbol in {"C", "N", "H"}]
    chn_set = set(chn_indices)
    max_cut = max(bond_cutoff(a, b) for a in ("C", "N", "H") for b in ("C", "N", "H"))
    graph = nx.Graph()
    graph.add_nodes_from(range(len(struct)))
    for i in chn_indices:
        for nn in struct.get_neighbors(struct[i], max_cut):
            j = nn.index
            if j <= i or j not in chn_set:
                continue
            cutoff = bond_cutoff(struct[i].specie.symbol, struct[j].specie.symbol)
            if nn.nn_distance <= cutoff:
                graph.add_edge(i, j)
    return graph


def get_unitcell_indexdata(
    system,
    cage_cutoff=6.5,
    pb_cutoff=6.5,
    br_orig_cutoff=4.5,
    exclude_br_axis=False,
    pb_orig_shift_val=5.926,
    ortho_axes=False,
):
    """Extract MAPbBr3 unit-cell index data from an ASE Atoms reference."""
    struct = Structure(
        lattice=system.cell,
        species=system.get_chemical_symbols(),
        coords=system.get_positions(),
        coords_are_cartesian=True,
    )
    graph = _build_chn_graph(struct)
    all_ma_indices = [
        list(comp)
        for comp in nx.connected_components(graph)
        if len(comp) == 8 and "C" in [struct[i].specie.symbol for i in comp]
    ]
    if not all_ma_indices:
        raise ValueError("No methylammonium molecules found in the structure.")

    pb_indices = [i for i, site in enumerate(struct) if site.specie.symbol == "Pb"]
    br_indices = [i for i, site in enumerate(struct) if site.specie.symbol == "Br"]
    lattice = struct.lattice
    frac_coords = struct.frac_coords

    unitcell_indexdata = []
    for ma_indices in all_ma_indices:
        c_index = [i for i in ma_indices if struct[i].specie.symbol == "C"][0]
        n_index = [i for i in ma_indices if struct[i].specie.symbol == "N"][0]
        h_indices_c, h_indices_n = _split_hydrogens(ma_indices, struct, graph)

        cage_pb_indices = [
            i for i in pb_indices if _pbc_distance(lattice, frac_coords, c_index, i) <= cage_cutoff
        ]
        cage_br_indices = [
            i for i in br_indices if _pbc_distance(lattice, frac_coords, c_index, i) <= cage_cutoff
        ]
        if len(cage_pb_indices) != 8:
            continue

        origin_pb, pb_env_inds = _get_pb_inds_orig_env(
            lattice, frac_coords, c_index, n_index, pb_indices, cutoff=pb_cutoff
        )
        axis_pbs = []
        for dirvec in lattice.matrix:
            cand = _closest_along_direction(
                lattice, frac_coords, origin_pb, dirvec, pb_orig_shift_val, pb_env_inds
            )
            if cand is not None:
                axis_pbs.append(int(cand))
        pb_axis = (origin_pb, axis_pbs)
        axis_vecs = [_shortest_vec(lattice, frac_coords, origin_pb, pb_idx) for pb_idx in axis_pbs]

        unitcell_data = {
            "c_ind": c_index,
            "n_ind": n_index,
            "h_inds_c": h_indices_c,
            "h_inds_n": h_indices_n,
            "pb_inds": cage_pb_indices,
            "br_inds": cage_br_indices,
            "pb_axis": pb_axis,
        }

        if not exclude_br_axis:
            br_neigh = [
                i for i in br_indices if _pbc_distance(lattice, frac_coords, origin_pb, i) <= br_orig_cutoff
            ]
            br_axis, br_axis_neg = _build_br_axes(
                lattice, frac_coords, origin_pb, axis_vecs, br_indices, br_neigh
            )
            unitcell_data["br_axis"] = br_axis
            unitcell_data["br_axis_neg"] = br_axis_neg

        if ortho_axes:
            pb_dists_from_origin = [
                np.linalg.norm(_shortest_vec(lattice, frac_coords, origin_pb, idx))
                for idx in pb_indices
            ]
            pb_inds_second = [
                idx
                for idx, dist in zip(pb_indices, pb_dists_from_origin)
                if cage_cutoff < dist <= np.sqrt(2) * cage_cutoff
            ]
            pb_inds_fourth = [
                idx
                for idx, dist in zip(pb_indices, pb_dists_from_origin)
                if np.sqrt(3) * cage_cutoff < dist <= 2 * cage_cutoff
            ]
            unitcell_data["pb_ortho_axes_small"] = _get_pb_ortho_axes(
                lattice,
                frac_coords,
                origin_pb,
                axis_vecs,
                pb_axis,
                pb_inds_second,
                [],
                shift_val=np.sqrt(2) * pb_orig_shift_val,
                shift_val_pc=None,
            )
            unitcell_data["pb_ortho_axes"] = _get_pb_ortho_axes(
                lattice,
                frac_coords,
                origin_pb,
                axis_vecs,
                pb_axis,
                pb_inds_second,
                pb_inds_fourth,
                shift_val=np.sqrt(2) * pb_orig_shift_val,
                shift_val_pc=2 * pb_orig_shift_val,
            )

        unitcell_indexdata.append(unitcell_data)

    return unitcell_indexdata


def ma_indices_from_unitcell_data(udata):
    return [udata["c_ind"], *udata["h_inds_c"], udata["n_ind"], *udata["h_inds_n"]]


def host_indices_from_unitcell_data(udata):
    return [udata["pb_axis"][0], *udata["br_axis"]]


def unit_order_from_unitcell_data(unitcell_data, unitcell_index=0):
    udata = unitcell_data[unitcell_index]
    return [*ma_indices_from_unitcell_data(udata), *host_indices_from_unitcell_data(udata)]


def unit_orders_from_unitcell_data(unitcell_data):
    return [unit_order_from_unitcell_data(unitcell_data, i) for i in range(len(unitcell_data))]


def unit_order_from_reference(reference, unitcell_index=0, **kwargs):
    return unit_order_from_unitcell_data(
        get_unitcell_indexdata(reference, **kwargs),
        unitcell_index=unitcell_index,
    )
