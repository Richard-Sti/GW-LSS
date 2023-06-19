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
Plot GW170817 posterior localisation samples.
"""
import matplotlib.pyplot as plt
import numpy
import scienceplots  # noqa

import utils

try:
    import gwlss
except ImportError:
    import sys
    sys.path.append("../")
    import gwlss


def plot_gw170817():
    """
    Plot GW170817 posterior localisation samples: ra-dec and redshift
    histogram.
    """
    paths = gwlss.Paths(gwlss.paths_glamdring)
    samples = gwlss.open_gw170817()

    with plt.style.context(utils.mplstyle):
        plt.figure()

        plt.hist(samples["redshift"], bins="auto")
        plt.axvline(gwlss.gw170817["redshift"], c="red",
                    label="EM counterpart")
        plt.legend()
        plt.xlabel(r"$z$")
        plt.ylabel("Counts")
        plt.tight_layout()

        for ext in utils.ext:
            fout = f"../plots/GW170817_distance.{ext}"
            print(f"Saving to `{fout}`.")
            plt.savefig(fout, dpi=utils.dpi)
        plt.close()

    with plt.style.context(utils.mplstyle):
        plt.figure()
        plt.scatter(numpy.rad2deg(samples['ra']),
                    numpy.rad2deg(samples['dec']), s=0.01)
        plt.scatter(gwlss.gw170817["RA"], gwlss.gw170817["dec"],
                    c="red", label="EM counterpart",
                    s=5, marker="x")
        plt.legend()
        plt.xlabel("RA [deg]")
        plt.ylabel("dec [deg]")
        plt.tight_layout()

        for ext in utils.ext:
            fout = f"../plots/GW170817_sky.{ext}"
            print(f"Saving to `{fout}`.")
            plt.savefig(fout, dpi=utils.dpi)
        plt.close()


def plot_gw170817_rotated(nrot):
    """
    Plot randomly rotated GW170817 posterior ra-dec samples.

    Parameters
    ----------
    nrot : int
        Number of random rotations.
    """
    samples = gwlss.open_gw170817()
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
        for ext in utils.ext:
            fout = f"../plots/GW170817_rotate_sky.{ext}"
            print(f"Saving to `{fout}`.")
            plt.savefig(fout, dpi=utils.dpi)
        plt.close()


def plot_gw170817_field(kind, nsim, smooth_scale=None):
    r"""
    Plot GW170817-evaluated CSiBORG field.

    Parameters
    ----------
    kind : str
        Field kind.
    nsim : int
        Simulation index.
    smooth_scale : float, optional
        Smoothing scale in :math:`\mathrm{Mpc}/h`.
    """
    paths = gwlss.Paths(gwlss.paths_glamdring)

    with plt.style.context("science"):
        plt.figure()

        f = paths.evaluated_field("GW170817", "density", nsim, 256,
                                  is_rand=False, smooth_scale=smooth_scale)
        data = numpy.load(f)
        bins = numpy.linspace(data.min(), data.max() + 2, 50)
        plt.hist(data, bins=bins, density=1, histtype="step",
                 label="GW170817")

        f = paths.evaluated_field("GW170817", "density", nsim, 256,
                                  is_rand=True)
        data = numpy.load(f)
        for i in range(35):
            plt.hist(data[i, :], bins=bins, density=1, histtype="step",
                     label="Random" if i == 0 else None, ls="dotted")

        plt.xlabel(r"$\rho / \langle \rho \rangle$")
        plt.ylabel("Normalized counts")
        plt.legend()

        # plt.yscale("log")
        plt.tight_layout()
        for ext in utils.ext:
            fout = f"../plots/GW170817_{kind}.{ext}"
            print(f"Saving to `{fout}`.")
            plt.savefig(fout, dpi=utils.dpi)
        plt.close()


if __name__ == "__main__":
    if False:
        plot_gw170817()

    if False:
        plot_gw170817_rotated(50)

    if True:
        plot_gw170817_field("overdensity", 7444)
