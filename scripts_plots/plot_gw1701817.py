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
import scienceplots  # noqa
from h5py import File

import utils


def plot_gw170817():
    data = File("/mnt/extraspace/rstiskalek/GWLSS/H1L1V1-EXTRACT_POSTERIOR_GW170817-1187008600-400.hdf", 'r')  # noqa
    samples = data["samples"]

    with plt.style.context(utils.mplstyle):
        plt.figure()
        plt.hist(samples["distance"], bins="auto")
        plt.xlabel(r"$D_{\rm L} ~ [\mathrm{Mpc}]$")
        plt.ylabel("Counts")
        plt.tight_layout()

        for ext in utils.ext:
            fout = f"../plots/GW170817_distance.{ext}"
            print(f"Saving to `{fout}`.")
            plt.savefig(fout, dpi=utils.dpi)
        plt.close()

    with plt.style.context(['science', 'notebook']):
        plt.figure()
        plt.scatter(samples['ra'], samples['dec'], s=0.01)
        plt.scatter([3.44616], [-0.408084], c="red", label="EM counterpart",
                    s=5, marker="x")
        plt.legend()
        plt.xlabel("RA")
        plt.ylabel("dec")
        plt.tight_layout()

        for ext in utils.ext:
            fout = f"../plots/GW170817_sky.{ext}"
            print(f"Saving to `{fout}`.")
            plt.savefig(fout, dpi=utils.dpi)
        plt.close()


if __name__ == "__main__":
    plot_gw170817()
