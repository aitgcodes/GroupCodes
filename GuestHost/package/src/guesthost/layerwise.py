import numpy as np


dirs_perp = [(1, 2), (2, 0), (0, 1)]


def _default_cellshape(lattice):
    return (lattice.nx, lattice.ny, lattice.nz)


def _cellshape_permute(values, order, cellshape):
    arr = np.array(values, dtype=object).reshape(cellshape)
    return np.transpose(arr, order)


def _index_array(lattice, cellshape=None, order=None):
    if cellshape is None:
        cellshape = _default_cellshape(lattice)
    if order is None:
        order = [0, 1, 2]
    indices = np.empty(cellshape, dtype=object)
    for ind in np.ndindex(cellshape):
        indices[ind] = ind
    return np.transpose(indices, order)


def compute_all_ucells(func, lattice, **kwargs):
    return [func(lattice, ind, **kwargs) for ind in np.ndindex(lattice.nx, lattice.ny, lattice.nz)]


def compute_all_ucells_matrix(func, lattice, order=None, cellshape=None, **kwargs):
    if cellshape is None:
        cellshape = _default_cellshape(lattice)
    ind_arr = _index_array(lattice, cellshape=cellshape, order=order)
    out = np.empty_like(ind_arr, dtype=object)
    for idx in np.ndindex(ind_arr.shape):
        out[idx] = np.array(func(lattice, ind_arr[idx], **kwargs))
    return out


def compute_all_ucells_data(func, lattice, order=None, cellshape=None, **kwargs):
    if cellshape is None:
        cellshape = _default_cellshape(lattice)
    ind_arr = _index_array(lattice, cellshape=cellshape, order=order)
    out = np.empty_like(ind_arr, dtype=object)
    for idx in np.ndindex(ind_arr.shape):
        out[idx] = func(lattice, ind_arr[idx], **kwargs)
    return out


def layerwise_OmegaLat(
    lattice,
    dir_coup,
    cellshape=None,
    k_vec=None,
    orig_cell=(0, 0),
    alongaxis=False,
    swap_xz=False,
):
    if cellshape is None:
        cellshape = _default_cellshape(lattice)
    if k_vec is None:
        k_vec = 2 * np.pi * np.array([0.5, 0.5])

    br_dirs = dirs_perp[dir_coup]
    order = [dir_coup, br_dirs[0], br_dirs[1]]
    ind_arr = _index_array(lattice, cellshape=cellshape, order=order)
    if swap_xz:
        cellshape = cellshape[::-1]
    layershape = (cellshape[br_dirs[0]], cellshape[br_dirs[1]])

    omega_layers = []
    for i in range(cellshape[dir_coup]):
        ind_mat = ind_arr[i, :, :]
        omega = np.zeros(2, dtype=float)
        n1, n2 = layershape
        for j in range(n1):
            for k in range(n2):
                R = np.array([k, j]) + np.array(orig_cell)
                U = lattice.U_r(R, k_vec, MA=False)
                if alongaxis:
                    eta = lattice.ucell_eta_lat_alongaxis(ind_mat[j, k], dir_coup)
                else:
                    eta = lattice.ucell_eta_lat(ind_mat[j, k], dir_coup)
                omega += U @ np.array(eta)
        omega_layers.append(omega / (n1 * n2))
    return omega_layers


def layerwise_eta_Lat(
    lattice,
    dir_coup,
    cellshape=None,
    k_vec=None,
    V=None,
    orig_cell=(0, 0),
):
    if cellshape is None:
        cellshape = _default_cellshape(lattice)
    if k_vec is None:
        k_vec = 2 * np.pi * np.array([0.5, 0.5])
    if V is None:
        V = 0.5 * np.array([1.0, -1.0])

    br_dirs = dirs_perp[dir_coup]
    order = [dir_coup, br_dirs[0], br_dirs[1]]
    ind_arr = _index_array(lattice, cellshape=cellshape, order=order)
    layershape = (cellshape[br_dirs[0]], cellshape[br_dirs[1]])

    eta_layers = []
    for i in range(cellshape[dir_coup]):
        ind_mat = ind_arr[i, :, :]
        n1, n2 = layershape
        layer_vals = np.zeros((n1, n2))
        for j in range(n1):
            for k in range(n2):
                R = np.array([k, j]) + np.array(orig_cell)
                U = lattice.U_r(R, k_vec, MA=False)
                eta = lattice.ucell_eta_lat(ind_mat[j, k], dir_coup)
                layer_vals[j, k] = V @ (U @ np.array(eta))
        eta_layers.append(layer_vals)
    return eta_layers


def global_layerwise_xi(
    lattice,
    dir_coup,
    q=0,
    cellshape=None,
    orig_cell=(0, 0),
    alongaxis=False,
    swap_xz=False,
):
    omega_layers = layerwise_OmegaLat(
        lattice, dir_coup, cellshape, None, orig_cell, alongaxis, swap_xz
    )
    ph = np.exp(1j * q)
    return [
        float(np.real((ph ** i) * 0.5 * (om[0] - om[1])))
        for i, om in enumerate(omega_layers)
    ]


def global_layerwise_eta_alloctahedra(lattice, dir_coup, cellshape=None, alongaxis=False):
    if cellshape is None:
        cellshape = _default_cellshape(lattice)
    fn = (
        lambda lat, ind, dir_coup: lat.ucell_eta_lat_alongaxis(ind, dir_coup)
        if alongaxis
        else lat.ucell_eta_lat(ind, dir_coup)
    )
    eta_all = compute_all_ucells_data(
        fn,
        lattice,
        order=[dir_coup, *dirs_perp[dir_coup]],
        cellshape=cellshape,
        dir_coup=dir_coup,
    )
    dims = [2, 1, 0][dir_coup]
    return [np.take(eta_all, i, axis=dims) for i in range(eta_all.shape[dims])]


def global_layerwise_xi_alloctahedra(lattice, dir_coup, cellshape=None, alongaxis=False):
    eta_layers = global_layerwise_eta_alloctahedra(
        lattice, dir_coup, cellshape=cellshape, alongaxis=alongaxis
    )
    results = []
    for lmat in eta_layers:
        d1 = np.vectorize(lambda x: x[0])(lmat)
        d2 = np.vectorize(lambda x: x[1])(lmat)
        doct = (d1 - d2) / 2.0
        sign = np.fromfunction(lambda i, j: (-1) ** (i + j), doct.shape)
        results.append(sign * doct)
    return results


def global_layerwise_xi_new(lattice, dir_coup, cellshape=None, alongaxis=False):
    xi_layers = global_layerwise_xi_alloctahedra(
        lattice, dir_coup, cellshape=cellshape, alongaxis=alongaxis
    )
    return [float(np.mean(xi)) for xi in xi_layers]


def layerwise_OmegaMA(
    lattice,
    dir,
    cellshape=None,
    k_vec=None,
    change_phi_domain=True,
    rot_MA=True,
    phi_domain=(-45, 135),
    eps=None,
    norm=False,
    orig_cell=(0, 0),
):
    if cellshape is None:
        cellshape = _default_cellshape(lattice)
    if k_vec is None:
        k_vec = 2 * np.pi * np.array([0.5, 0.5])

    perp_dirs = dirs_perp[dir]
    order = [dir, perp_dirs[0], perp_dirs[1]]
    ind_arr = _index_array(lattice, cellshape=cellshape, order=order)
    layershape = (cellshape[perp_dirs[0]], cellshape[perp_dirs[1]])

    omega_layers = []
    for i in range(cellshape[dir]):
        ind_mat = ind_arr[i, :, :]
        omega = np.zeros(2, dtype=float)
        n1, n2 = layershape
        for j in range(n1):
            for k in range(n2):
                R = np.array([k - 1, j - 1]) + np.array(orig_cell)
                U = lattice.U_r(R, k_vec, MA=rot_MA)
                eta = lattice.ucell_eta_MA(
                    ind_mat[j, k],
                    dir,
                    change_phi_domain=change_phi_domain,
                    phi_domain=phi_domain,
                    eps=eps,
                    norm=norm,
                )
                omega += U @ np.array(eta)
        omega_layers.append(omega / (n1 * n2))
    return omega_layers


def layerwise_eta_MA(
    lattice,
    dir,
    cellshape=None,
    k_vec=None,
    change_phi_domain=True,
    rot_MA=True,
    V=None,
    phi_domain=(-45, 135),
    eps=None,
    orig_cell=(0, 0),
):
    if cellshape is None:
        cellshape = _default_cellshape(lattice)
    if k_vec is None:
        k_vec = 2 * np.pi * np.array([0.5, 0.5])
    if V is None:
        V = np.array([1.0, 0.0])

    perp_dirs = dirs_perp[dir]
    order = [dir, perp_dirs[0], perp_dirs[1]]
    ind_arr = _index_array(lattice, cellshape=cellshape, order=order)
    layershape = (cellshape[perp_dirs[0]], cellshape[perp_dirs[1]])

    eta_layers = []
    for i in range(cellshape[dir]):
        ind_mat = ind_arr[i, :, :]
        n1, n2 = layershape
        layer_vals = np.zeros((n1, n2))
        for j in range(n1):
            for k in range(n2):
                R = np.array([k - 1, j - 1]) + np.array(orig_cell)
                U = lattice.U_r(R, k_vec, MA=rot_MA)
                eta = lattice.ucell_eta_MA(
                    ind_mat[j, k],
                    dir,
                    change_phi_domain=change_phi_domain,
                    phi_domain=phi_domain,
                    eps=eps,
                )
                omega = U @ np.array(eta)
                layer_vals[j, k] = np.dot(V, omega) if V is not None else omega
        eta_layers.append(layer_vals)
    return eta_layers


def global_xi_MA(
    lattice,
    dir,
    q=0,
    cellshape=None,
    k_vec=None,
    change_phi_domain=True,
    rot_MA=True,
    V=None,
    phi_domain=(-45, 135),
    eps=None,
    orig_cell=(0, 0),
):
    if V is None:
        V = np.array([1.0, 0.0])
    omega_layers = layerwise_OmegaMA(
        lattice,
        dir,
        cellshape,
        k_vec,
        change_phi_domain,
        rot_MA,
        phi_domain,
        eps,
        False,
        orig_cell,
    )
    ph = np.exp(1j * q)
    total = 0.0
    for i, omega in enumerate(omega_layers):
        total += np.real((ph ** i) * np.dot(omega, V))
    return total / len(omega_layers)


def layerwise_OmegaCoupling(
    lattice,
    dir,
    cellshape=None,
    k_vec=None,
    change_phi_domain=True,
    rot_MA=True,
    phi_domain=(-45, 135),
    eps=None,
    norm=False,
    orig_cell=(0, 0),
):
    if cellshape is None:
        cellshape = _default_cellshape(lattice)
    if k_vec is None:
        k_vec = 2 * np.pi * np.array([0.5, 0.5])

    perp_dirs = dirs_perp[dir]
    order = [dir, perp_dirs[0], perp_dirs[1]]
    ind_arr = _index_array(lattice, cellshape=cellshape, order=order)
    layershape = (cellshape[perp_dirs[0]], cellshape[perp_dirs[1]])

    omega_layers = []
    for i in range(cellshape[dir]):
        ind_mat = ind_arr[i, :, :]
        omega = np.zeros(2, dtype=float)
        n1, n2 = layershape
        for j in range(n1):
            for k in range(n2):
                R = np.array([k - 1, j - 1]) + np.array(orig_cell)
                U_MA = lattice.U_r(R, k_vec, MA=rot_MA)
                eta_MA = lattice.ucell_eta_MA(
                    ind_mat[j, k],
                    dir,
                    change_phi_domain=change_phi_domain,
                    phi_domain=phi_domain,
                    eps=eps,
                    norm=norm,
                )
                omega_MA = U_MA @ np.array(eta_MA)

                U_lat = lattice.U_r(R, k_vec, MA=False)
                eta_lat = lattice.ucell_eta_lat(ind_mat[j, k], dir)
                omega_lat = U_lat @ np.array(eta_lat)

                omega += omega_lat * omega_MA
        omega_layers.append(omega / (n1 * n2))
    return omega_layers


def global_xi_coupling(
    lattice,
    dir,
    q=0,
    cellshape=None,
    k_vec=None,
    change_phi_domain=True,
    rot_MA=True,
    V=None,
    phi_domain=(-45, 135),
    eps=None,
    norm=False,
    orig_cell=(0, 0),
):
    if V is None:
        V = np.array([1.0, 0.0])
    omega_layers = layerwise_OmegaCoupling(
        lattice,
        dir,
        cellshape,
        k_vec,
        change_phi_domain,
        rot_MA,
        phi_domain,
        eps,
        norm,
        orig_cell,
    )
    ph = np.exp(1j * q)
    total = 0.0
    for i, omega in enumerate(omega_layers):
        total += np.real((ph ** i) * np.dot(omega, V))
    return total / len(omega_layers)


def global_volume_local(lattice):
    vols = [lattice.ucell_localcellparameters(ind)["vol"] for ind in np.ndindex(lattice.nx, lattice.ny, lattice.nz)]
    return float(np.mean(vols))


def global_xi(
    lattice,
    dir_coup,
    q=0,
    cellshape=None,
    V=None,
    orig_cell=(0, 0),
    alongaxis=False,
    swap_xz=False,
):
    if cellshape is None:
        cellshape = _default_cellshape(lattice)
    if V is None:
        V = 0.5 * np.array([1.0, -1.0])
    omega_layers = layerwise_OmegaLat(
        lattice, dir_coup, cellshape, None, orig_cell, alongaxis, swap_xz
    )
    nlayers = cellshape[::-1][dir_coup] if swap_xz else cellshape[dir_coup]
    ph = np.exp(1j * q)
    return float(sum(np.real((ph ** i) * (om @ V)) for i, om in enumerate(omega_layers)) / nlayers)


def global_S(lattice, dir_coup, cellshape=None, lags=None, demean=False):
    if cellshape is None:
        cellshape = _default_cellshape(lattice)
    if lags is None:
        lags = range(0, cellshape[dir_coup])
    omega_layers = layerwise_OmegaLat(lattice, dir_coup, cellshape)
    nlayers = cellshape[dir_coup]
    xi_layers = np.array([0.5 * (om[0] - om[1]) for om in omega_layers], dtype=float)
    if demean:
        xi_layers = xi_layers - xi_layers.mean()
    ac = []
    for lag in lags:
        n = nlayers - lag
        if n <= 0:
            ac.append(0.0)
            continue
        ac.append(float(np.sum(xi_layers[:n] * xi_layers[lag:lag + n]) / n))
    return ac


def compute_layerwise(obj, func, **kwargs):
    from guesthost.lattice import LatticeTrajectory

    if isinstance(obj, LatticeTrajectory):
        return [func(lattice, **kwargs) for lattice in obj.lattices]
    return func(obj, **kwargs)
