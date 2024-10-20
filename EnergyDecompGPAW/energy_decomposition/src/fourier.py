import numpy as np

from gpaw.tddft.units import eV_to_au


def fourier(time_t, a_tX, omega_w, folding='Gauss', width=0.08, sign=1.0):
    if width is None:
        folding = None

    if folding is None:
        def envelope(t):
            return 1.0
    else:
        width = width * eV_to_au
        if folding == 'Gauss':
            # f(w) = Nw exp(-w^2/2/sigma^2)
            # f(t) = exp(-t^2 sigma^2/2)
            def envelope(t):
                return np.exp(-0.5 * width**2 * t**2)
        elif folding == 'Lorentz':
            # f(w) = Nw eta / [omega^2 + (eta)^2]
            # f(t) = exp(-t eta)
            def envelope(t):
                return np.exp(-width * t)
        else:
            raise RuntimeError('unknown folding "' + folding + '"')

    dt_t = np.insert(time_t[1:] - time_t[:-1], 0, 0.0)
    f_wt = np.exp(sign * 1j * np.outer(omega_w, time_t))
    a_Xt = np.rollaxis(a_tX, 0, len(a_tX.shape))
    a_wX = np.tensordot(f_wt, dt_t * envelope(time_t) * a_Xt,
                        axes=(1, len(a_Xt.shape) - 1))
    return a_wX


def inversefourier(omega_w, a_wX, time_t, folding='Gauss', width=0.08):
    # This function assumes that time-domain quantity is real
    a_tX = fourier(omega_w, a_wX / (2 * np.pi), time_t, folding, width,
                   sign=-1)
    return 2 * a_tX.real
