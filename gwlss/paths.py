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

from os import makedirs
from os.path import isdir, join
from warnings import warn

from h5py import File


class Paths:
    """
    Object for handling paths to data for the GW-LSS project.

    Parameters
    ----------
    paths_dict : dict
        Dictionary of basic paths.
    """
    def __init__(self, paths_dict):
        assert isinstance(paths_dict, dict)
        self._paths_dict = paths_dict

    def load_event(self, event):
        r"""
        Load a specific event

        Parameters
        ----------
        event : str
            Event name.

        Returns
        -------
        path : str
        """
        if event == "GW170817":
            return File(self["GW170817_darkPE"], 'r')["samples"]
        else:
            raise KeyError(f"Event `{event}` not found.")

    def evaluated_field(self, event, kind, nsim, grid, MAS="PCS",
                        is_rand=False, in_rsp=True, smooth_scale=None):
        r"""
        Paths to the files containing the evaluated density fields in CSiBORG
        of a given event.

        Parameters
        ----------
        event : str
            Event name.
        kind : str
            Field type.
        nsim : int
            Simulation index.
        grid : int
            Grid size.
        MAS : str, optional
           Mass-assignment scheme.
        is_rand : bool, optional
            Whether the event is randomly rotated.
        in_rsp : bool, optional
            Whether the calculation is performed in redshift space.
        smooth_scale : float, optional
            Smoothing scale in :math:`\mathrm{Mpc}/h`

        Returns
        -------
        path : str
        """
        fdir = join(self["dumpdir"], "evaluated")
        if not isdir(fdir):
            makedirs(fdir)
            warn(f"Created directory `{fdir}`.", UserWarning, stacklevel=1)
        if is_rand:
            event = "rand_" + event
        if in_rsp:
            kind = kind + "_rsp"
        fname = f"{event}_{kind}_{MAS}_{str(nsim).zfill(5)}_grid{grid}.npz"
        if smooth_scale is not None and smooth_scale > 0:
            smooth_scale = float(smooth_scale)
            fname = fname.replace(".npz", f"smooth{smooth_scale}.npz")
        return join(fdir, fname)

    def __getitem__(self, name):
        try:
            return self._paths_dict[name]
        except KeyError:
            raise AttributeError(f"Path `{name}` not found in `paths_dict`.")
