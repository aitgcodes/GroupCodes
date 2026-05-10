import numpy as np

from guesthost.layerwise import dirs_perp


def _wrap180(x):
    x = np.asarray(x, dtype=float)
    x = np.where(x <= -180.0, x + 360.0, x)
    x = np.where(x >= 180.0, x - 360.0, x)
    return x


def _polar_op_from_phi(phi_arr, dir_coup):
    phi_prev = np.roll(phi_arr, shift=-1, axis=dir_coup)
    dphi = _wrap180(phi_prev - phi_arr)

    z = np.exp(1j * np.pi * dphi / 180.0)

    axis_len = phi_arr.shape[dir_coup]
    phase_shape = [1, 1, 1]
    phase_shape[dir_coup] = axis_len
    phase_pi = np.exp(-1j * np.pi * np.arange(axis_len)).reshape(phase_shape)
    phase_pi2 = np.exp(-1j * (np.pi / 2.0) * np.arange(axis_len)).reshape(phase_shape)

    z_pi = z * phase_pi
    z_pi_by2 = z * phase_pi2

    perp = dirs_perp[dir_coup]
    z_chain = z.mean(axis=perp)
    z_chain_pi = z_pi.mean(axis=perp)
    z_chain_pi_by2 = z_pi_by2.mean(axis=perp)

    return {
        "dphi": dphi,
        "z": z,
        "z_pi": z_pi,
        "z_pi_by2": z_pi_by2,
        "z_chain": z_chain,
        "z_chain_pi": z_chain_pi,
        "z_chain_pi_by2": z_chain_pi_by2,
        "S_chain": np.abs(z_chain) ** 2,
        "S_chain_pi": np.abs(z_chain_pi) ** 2,
        "S_chain_pi_by2": np.abs(z_chain_pi_by2) ** 2,
    }


def polar_order_parameter(lattice, dir_coup, cellshape=None):
    if cellshape is None:
        cellshape = (lattice.nx, lattice.ny, lattice.nz)

    phi_arr = np.empty(cellshape, dtype=float)
    for ind in np.ndindex(cellshape):
        _, phi = lattice.ucell_theta_phi(ind, dir=dir_coup)
        phi_arr[ind] = phi

    return _polar_op_from_phi(phi_arr, dir_coup)
