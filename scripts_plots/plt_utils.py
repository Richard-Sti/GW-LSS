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

import csiborgtools

dpi = 600
mplstyle = ["science"]
ext = ["png"]


def str2label(label):
    """
    Convert string to LaTeX label.

    Parameters
    ----------
    label : str
        String to convert.

    Returns
    -------
    converted : str
        Converted LaTeX label. If no conversion is found, the input string is
        returned.
    """
    data = {"density": r"$\rho$",
            "overdensity": r"$\delta$",
            "potential": r"$\Phi$"}
    return data.get(label, label)


def str2label_units(label):
    """
    Convert string to LaTeX label with units.

    Parameters
    ----------
    label : str
        String to convert.

    Returns
    -------
    converted : str
        Converted LaTeX label with units. If no conversion is found, the input
        label is returned.
    """
    data = {},
    return data.get(label, label)


def get_nsims(nsims):
    """
    Get CSiBORG simulations to use from the command line arguments.

    Parameters
    ----------
    args : argparse.Namespace
        Command line arguments.

    Returns
    -------
    nsims : list of int
    """
    if nsims[0] == -1:
        csiborg_paths = csiborgtools.read.Paths(**csiborgtools.paths_glamdring)
        return list(csiborg_paths.get_ics("csiborg"))
    else:
        return args.nsims