import numpy as np
from scipy.interpolate import interp1d
from gpaw.tddft.units import as_to_au, fs_to_au, eV_to_au
import fourier


def do_tX(f, a_tX, *args, **kwargs):
    # Flatten the X dimensions
    Nt = a_tX.shape[0]
    N_x = a_tX.shape[1:]
    Np = np.prod(N_x)
    a_tp = a_tX.reshape((Nt, Np))

    # Result array
    r_pY = []

    # Loop over each p
    for p in range(Np):
        r_Y = f(a_tp[:, p], *args, **kwargs)
        r_pY.append(r_Y)

    # Reshape p to X
    r_Yp = np.moveaxis(np.array(r_pY), 0, -1)
    N_y = r_Yp.shape[:-1]
    r_YX = r_Yp.reshape(N_y + N_x)
    return r_YX


class Convolver(object):

    def __init__(self, pulse):
        self.pulse = pulse

    def convolve(self, time, units='au'):
        # Convert time to ndarray
        time_t = np.array(time)
        if units == 'fs':
            time_t *= fs_to_au
        elif units == 'as':
            time_t *= as_to_au
        elif units != 'au':
            raise RuntimeError('unknown units')

        if time_t.ndim == 0:
            # If argument is a single number,
            # return value without time axis
            return self._convolve(time_t[np.newaxis])[0]

        return self._convolve(time_t)

    def _convolve(self, time_t):
        raise NotImplementedError()


class FourierConvolver(Convolver):
    """
    Convolve using Fourier transforms.

    Use given frequencies in Fourier transform.
    """
    def __init__(self, pulse, omega_w, data_wX):
        Convolver.__init__(self, pulse)

        # Impulse response data
        self.omega_w = omega_w
        d_wX = data_wX

        # Pulse in frequency space
        pulse_w = self.pulse.fourier(omega_w)

        # Convolution product
        d_Xw = np.moveaxis(d_wX, 0, -1)
        pd_Xw = d_Xw * pulse_w
        self.pd_wX = np.moveaxis(pd_Xw, -1, 0)

    def _convolve(self, time_t):
        pd_tX = fourier.inversefourier(self.omega_w, self.pd_wX, time_t, None)
        return pd_tX


class TimeConvolver(Convolver):
    """
    Convolve using time-domain integration

    """

    def __init__(self, pulse, time_t, data_tX):
        Convolver.__init__(self, pulse)

        # Impulse response data
        self.time_t = time_t
        self.d_tX = data_tX

        # Pulse in time
        pulse_t = self.pulse.strength(time_t)

        # Convolution integral
        dt_t = time_t[1:] - time_t[:-1]
        #if len(np.unique(np.around(dt_t, 6))) != 1:
        #    raise RuntimeError('multiple different time steps not supported')
        dt = time_t[1] - time_t[0]
        pd_TX = do_tX(np.convolve, self.d_tX, pulse_t) * dt

        # Convolution has data for points that do not fully overlap
        # Construct the corresponding time vector
        time_T = np.arange(pd_TX.shape[0]) * dt

        # Interpolator
        self.interp_pd_tX = interp1d(time_T, pd_TX, axis=0,
                                     kind='linear',
                                     bounds_error=False, fill_value=0.0)

    def _convolve(self, time_t):
        # Interpolate data to given time
        pd_tX = self.interp_pd_tX(time_t)
        return pd_tX

    def fourier(self, omega_w=None):
        if omega_w is None:
            omega_w = np.arange(0, 10 + 1e-3, 0.01) * eV_to_au

        d_wX = fourier.fourier(self.time_t, self.d_tX, omega_w, None)
        return FourierConvolver(self.pulse, omega_w, d_wX)
