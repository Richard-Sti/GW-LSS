# Copyright (C) 2023 Richard Stiskalek
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
"""
Plot event posterior localisation samples.
"""
import matplotlib.pyplot as plt
import numpy
import scienceplots  # noqa

import plt_utils

try:
    import gwlss
except ImportError:
    import sys
    sys.path.append("../")
    import gwlss


def plot_event_localisation(event):
    """
    Plot event posterior localisation samples: ra-dec and redshift histogram.

    Parameters
    ----------
    event : str
        Event name.
    """
    paths = gwlss.Paths(gwlss.paths_glamdring)
    samples = paths.load_event(event)

    with plt.style.context(plt_utils.mplstyle):
        plt.figure()

        plt.hist(samples["redshift"], bins="auto")
        try:
            plt.axvline(gwlss.EM_counterpart[event]["redshift"], c="red",
                        label="EM counterpart")
            plt.legend()
        except KeyError:
            pass
        plt.xlabel(r"$z$")
        plt.ylabel("Counts")

        plt.tight_layout()
        for ext in plt_utils.ext:
            fout = f"../plots/{event}_redshift.{ext}"
            print(f"Saving to `{fout}`.")
            plt.savefig(fout, dpi=plt_utils.dpi)
        plt.close()

    with plt.style.context(plt_utils.mplstyle):
        plt.figure()
        plt.scatter(numpy.rad2deg(samples['ra']),
                    numpy.rad2deg(samples['dec']), s=0.01)
        try:
            plt.scatter(gwlss.EM_counterpart[event]["RA"],
                        gwlss.EM_counterpart[event]["dec"], c="red",
                        label="EM counterpart", s=5, marker="x")
            plt.legend()
        except KeyError:
            pass
        plt.xlabel("RA [deg]")
        plt.ylabel("dec [deg]")

        plt.tight_layout()
        for ext in plt_utils.ext:
            fout = f"../plots/{event}_sky.{ext}"
            print(f"Saving to `{fout}`.")
            plt.savefig(fout, dpi=plt_utils.dpi)
        plt.close()


def plot_radec_rotated(event, nrot):
    """
    Plot a randomly rotated event's posterior RA-dec samples.

    Parameters
    ----------
    event : str
        Event name.
    nrot : int
        Number of random rotations.
    """
    paths = gwlss.Paths(gwlss.paths_glamdring)
    samples = paths.load_event(event)
    ra = samples["ra"][:]
    dec = samples["dec"][:]

    with plt.style.context("science"):
        plt.figure()
        plt.scatter(ra, dec, s=0.05, c="black", rasterized=True)
        for __ in range(nrot):
            rotmat = gwlss.rand_rotation_matrix()
            ra_rot, dec_rot = gwlss.rotate_radec(ra, dec, rotmat)
            plt.scatter(ra_rot, dec_rot, s=0.05, rasterized=True, alpha=0.5)

        plt.xlim(0, 2 * numpy.pi)
        plt.ylim(-numpy.pi / 2, numpy.pi / 2)
        plt.xlabel("RA [rad]")
        plt.ylabel("dec [rad]")

        plt.tight_layout()
        for ext in plt_utils.ext:
            fout = f"../plots/{event}_rotate_sky.{ext}"
            print(f"Saving to `{fout}`.")
            plt.savefig(fout, dpi=plt_utils.dpi)
        plt.close()


def plot_field(event, kind, nsims, grid, plot_rand=False, smooth_scale=None):
    r"""
    Plot a CSiBORG field evaluated at locations of the event's posterior
    localisation samples.

    Parameters
    ----------
    event : str
        Event name.
    kind : str
        Field kind.
    nsim : list of int
        Simulation indices.
    grid : int
        Grid size.
    plot_rand : bool, optional
        Whether to plot random rotations.
    smooth_scale : float, optional
        Smoothing scale in :math:`\mathrm{Mpc}/h`.
    """
    paths = gwlss.Paths(gwlss.paths_glamdring)

    with plt.style.context("science"):
        plt.figure()

        for nsim in nsims:
            f = paths.evaluated_field(event, kind, nsim, grid,
                                      is_rand=False, smooth_scale=smooth_scale)
            data = numpy.load(f)
            bins = numpy.linspace(data.min(), data.max() + 2, 50)
            plt.hist(data, bins=bins, density=1, histtype="step")

        if plot_rand:
            f = paths.evaluated_field(event, kind, nsim, grid,
                                      is_rand=True)
            data = numpy.load(f)
            for i in range(35):
                plt.hist(data[i, :], bins=bins, density=1, histtype="step",
                         label="Random" if i == 0 else None, ls="dotted")
            plt.legend()

        plt.xlabel(r"$\rho / \langle \rho \rangle$")
        plt.ylabel("Normalized counts")

        plt.tight_layout()
        for ext in plt_utils.ext:
            fout = f"../plots/{event}_{kind}.{ext}"
            print(f"Saving to `{fout}`.")
            plt.savefig(fout, dpi=plt_utils.dpi)
        plt.close()


if __name__ == "__main__":
    if True:
        grid = 512
        nsims = plt_utils.get_nsims([-1])
        plot_field("GW170817", "overdensity", nsims, grid)
