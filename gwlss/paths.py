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

from os.path import join

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
        str
        """
        if event == "GW170817":
            return File(self["GW170817_darkPE"], 'r')["samples"]
        else:
            raise KeyError(f"Event `{event}` not found.")

    def evaluated_field(self, event, kind, nsim, grid, MAS="PCS",
                        in_rsp=True):
        r"""
        Paths to the file containing the evaluated field in CSiBORG for a
        given event.

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
        in_rsp : bool, optional
            Whether the calculation is performed in redshift space.

        Returns
        -------
        str
        """
        kind = kind + "_rsp" if in_rsp else kind
        fdir = self["maindir"]
        if event == "GW170817":
            fname = f"{kind}_{MAS}_{grid}_{nsim}_H1L1V1-EXTRACT_POSTERIOR_GW170817-1187008600-400.npz"  # noqa
        else:
            raise KeyError(f"Event `{event}` not found.")

        return join(fdir, fname)

    def __getitem__(self, name):
        try:
            return self._paths_dict[name]
        except KeyError:
            raise AttributeError(f"Path `{name}` not found in `paths_dict`.")
