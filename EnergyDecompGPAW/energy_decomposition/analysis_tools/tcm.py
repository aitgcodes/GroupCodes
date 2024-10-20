import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LogNorm


def generate_gridspec(**kwargs):
    from matplotlib.gridspec import GridSpec
    width = 0.84
    bottom = 0.12
    left = 0.12
    return GridSpec(2, 2, width_ratios=[3, 1], height_ratios=[1, 3],
                    bottom=bottom, top=bottom + width,
                    left=left, right=left + width,
                    **kwargs)

class TCM(object):

    def __init__(self, energy_o, energy_u, fermilevel):
        self.energy_o = energy_o
        self.energy_u = energy_u
        self.fermilevel = fermilevel 

        self.base_o = np.zeros_like(energy_o)
        self.base_u = np.zeros_like(energy_u)

    def __getattr__(self, attr):
        # Generate axis only when needed
        #if attr in ['ax_occ_dos', 'ax_unocc_dos', 'ax_tcm']:
        if attr in ['ax_tcm']:
            gs = generate_gridspec(hspace=0.05, wspace=0.05)
            #self.ax_occ_dos = plt.subplot(gs[0])
            #self.ax_unocc_dos = plt.subplot(gs[3])
            self.ax_tcm = plt.subplot(gs[2])
            return getattr(self, attr)
        if attr in ['ax_spec']:
            gs = generate_gridspec(hspace=0.8, wspace=0.8)
            self.ax_spec = plt.subplot(gs[1])
            return getattr(self, attr)
        if attr in ['ax_cbar']:
            self.ax_cbar = plt.axes((0.15, 0.6, 0.02, 0.1))
            return getattr(self, attr)
        raise AttributeError('%s object has no attribute %s' %
                             (repr(self.__class__.__name__), repr(attr)))

    def plot_TCM(self, tcm_ou, vmax='80%', vmin='symmetrize', cmap='seismic',
                 log=False, colorbar=True, lw=None):
        if lw is None:
            lw = mpl.rcParams['lines.linewidth']
        energy_o = self.energy_o
        energy_u = self.energy_u
        fermilevel = self.fermilevel

        tcmmax = np.max(np.absolute(tcm_ou))
        print('tcmmax', tcmmax)

        # Plot TCM
        ax = self.ax_tcm
        plt.sca(ax)
        plt.cla()
        if isinstance(vmax, str):
            assert vmax[-1] == '%'
            tcmmax = np.max(np.absolute(tcm_ou))
            vmax = tcmmax * float(vmax[:-1]) / 100.0
        if vmin == 'symmetrize':
            vmin = -vmax
        if tcm_ou.dtype == complex:
            linecolor = 'w'
            from matplotlib.colors import hsv_to_rgb

            def transform_to_hsv(z, rmin, rmax, hue_start=90):
                amp = np.absolute(z)  # **2
                amp = np.where(amp < rmin, rmin, amp)
                amp = np.where(amp > rmax, rmax, amp)
                ph = np.angle(z, deg=1) + hue_start
                h = (ph % 360) / 360
                s = 1.85 * np.ones_like(h)
                v = (amp - rmin) / (rmax - rmin)
                return hsv_to_rgb(np.dstack((h, s, v)))

            img = transform_to_hsv(tcm_ou.T, 0, vmax)
            plt.imshow(img, origin='lower',
                       extent=[energy_o[0], energy_o[-1],
                               energy_u[0], energy_u[-1]],
                       interpolation='bilinear',
                       )
        else:
            linecolor = 'k'
            if cmap == 'magma':
                linecolor = 'w'
            if log:
                norm = LogNorm(vmin=vmin, vmax=vmax)
            else:
                norm = Normalize(vmin=vmin, vmax=vmax)
            plt.pcolormesh(energy_o, energy_u, tcm_ou.T,
                           cmap=cmap, rasterized=True, norm=norm,
                           )
        if colorbar:
            ax = self.ax_cbar
            ax.clear()
            cb = plt.colorbar(cax=ax)
            cb.outline.set_edgecolor(linecolor)
            ax.tick_params(axis='both', colors=linecolor)
            # ax.yaxis.label.set_color(linecolor)
            # ax.xaxis.label.set_color(linecolor)
        ax = self.ax_tcm
        plt.sca(ax)
        plt.axhline(fermilevel, c=linecolor, lw=lw)
        plt.axvline(fermilevel, c=linecolor, lw=lw)

        ax.tick_params(axis='both', which='major', pad=2)
        plt.xlabel(r'Hole energy (eV)', labelpad=0)
        plt.ylabel(r'Electron energy (eV)', labelpad=0)
        plt.xlim(np.take(energy_o, (0, -1)))
        plt.ylim(np.take(energy_u, (0, -1)))

    def plot_spectrum(self):
        raise NotImplementedError()

    def plot_TCM_diagonal(self, energy, **kwargs):
        x_o = np.take(self.energy_o, (0, -1))
        self.ax_tcm.plot(x_o, x_o + energy, **kwargs)

    def set_title(self, *args, **kwargs):
        self.ax_occ_dos.set_title(*args, **kwargs)


class TCMPlotter(TCM):

    def __init__(self, ksd, energy_o, energy_u, sigma, zero_fermilevel=True):
        eig_n, fermilevel = ksd.get_eig_n(zero_fermilevel)
        TCM.__init__(self, energy_o, energy_u, fermilevel)
        self.ksd = ksd
        self.sigma = sigma
        self.eig_n = eig_n

    def plot_TCM(self, tcm_ou, **kwargs):
        
        TCM.plot_TCM(self, tcm_ou, **kwargs)
