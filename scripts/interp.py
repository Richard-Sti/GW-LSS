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

from argparse import ArgumentParser
import csiborgtools
import numpy
from tqdm import trange
from mpi4py import MPI

from taskmaster import work_delegation

try:
    import gwlss
except ImportError:
    import sys
    sys.path.append("../")
    import gwlss


def load_field(kind, nsim, grid, MAS, in_rsp=False, smooth_scale=None,
               save_info=True):
    r"""
    Load a single CSiBORG field.

    Parameters
    ----------
    kind : str
        Field kind.
    nsim : int
        Simulation index.
    grid : int
        Grid size.
    MAS : str
        Mass assignment scheme.
    in_rsp : bool, optional
        Whether to load the field in redshift space.
    smooth_scale : float, optional
        Smoothing scale in :math:`\mathrm{Mpc} / h`.

    Returns
    -------
    field : 3-dimensional array
        Field evaluated on a 3D grid.
    info : dict
        Dictionary with information about the field.
    """
    smooth_scale = 0 if smooth_scale is None else smooth_scale
    paths = csiborgtools.read.Paths(**csiborgtools.paths_glamdring)
    # If the field is an overdensity, we need to load the density field.
    if kind == "overdensity":
        isoverdensity = True
        kind = "density"

    # For density field we smoothen it here instead of loading smoothed.
    if kind == "density" and smooth_scale > 0:
        field = numpy.load(paths.field(kind, MAS, grid, nsim, in_rsp=in_rsp))
        field = csiborgtools.field.smoothen_field(field, smooth_scale,
                                                  boxsize=677.7)
    else:
        field = numpy.load(paths.field(kind, MAS, grid, nsim, in_rsp=in_rsp,
                                       smooth_scale=smooth_scale))
    if isoverdensity:
        field /= field.mean()
        field -= 1

    info = {"mean": numpy.mean(field),
            "std": numpy.std(field)}

    return field, info


def evaluate_event(event, kind, nsim, grid, smooth_scale=None, to_save=True,
                   nrot=None, verbose=True, seed=None):
    r"""
    Evaluate a CSiBORG field at the positions of a given event.

    Parameters
    ----------
    event : str
        Event name.
    kind : str
        Field kind.
    nsim : int
        Simulation index.
    grid : int
        Grid size.
    smooth_scale : float, optional
        Smoothing scale in :math:`\mathrm{Mpc} / h`.
    to_save : bool, optional
        Whether to save the result.
    nrot : int, optional
        Number of random rotations.
    verbose : bool, optional
        Verbose flag.
    seed : int, optional
        Random seed for `nrot`.

    Returns
    -------
    val : 1-dimensional array
        Field values at the event positions.
    """
    # First load the event
    paths = gwlss.Paths(gwlss.paths_glamdring)
    samples = paths.load_event(event)
    ra0 = samples["ra"][:]
    dec0 = samples["dec"][:]
    dist = gwlss.cosmo_csiborg.comoving_distance(samples["redshift"][:]).value

    # Load the CSiBORG field
    csiborg_paths = csiborgtools.read.Paths(**csiborgtools.paths_glamdring)
    box = csiborgtools.read.CSiBORGBox(
        max(csiborg_paths.get_snapshots(nsim)), nsim, csiborg_paths)
    field, info = load_field(kind, nsim, grid, "PCS", in_rsp=True,
                             smooth_scale=smooth_scale)

    # Create the position array
    pos = numpy.vstack([dist, ra0, dec0]).T
    # Either evaluate it straight away or rotate the event randomly.
    if nrot is None:
        val = csiborgtools.field.evaluate_sky(field, pos=pos, box=box,
                                              isdeg=False)
    else:
        assert isinstance(nrot, int)
        val = numpy.full((nrot, ra0.size), numpy.nan)
        for i in trange(nrot) if verbose else range(nrot):
            rotmat = gwlss.rand_rotation_matrix(seed=seed)
            ra_rot, dec_rot = gwlss.rotate_radec(ra0, dec0, rotmat)
            pos[:, 1] = ra_rot
            pos[:, 2] = dec_rot
            val[i, :] = csiborgtools.field.evaluate_sky(field, pos=pos,
                                                        box=box, isdeg=False)

    if to_save:
        fout = paths.evaluated_field("GW170817", kind, nsim, grid,
                                     is_rand=nrot is not None,
                                     smooth_scale=smooth_scale)
        if verbose:
            print(f"Saving output to `{fout}`.", flush=True)
        numpy.savez(fout, values=val, mean=info["mean"], std=info["std"])

    return val


def get_nsims(args):
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
    if args.nsims[0] == -1:
        csiborg_paths = csiborgtools.read.Paths(**csiborgtools.paths_glamdring)
        return list(csiborg_paths.get_ics("csiborg"))
    else:
        return args.nsims


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--event", type=str, help="Event name.")
    parser.add_argument("--kind", type=str, help="Field kind.")
    parser.add_argument("--smooth_scale", type=float, default=0.)
    parser.add_argument("--nsims", type=int, nargs="+", default=None,
                        help="IC realisations. `-1` for all simulations.")
    parser.add_argument("--grid", type=int, help="Grid size.")
    args = parser.parse_args()
    COMM = MPI.COMM_WORLD

    if True:
        def main(nsim):
            return evaluate_event(
                args.event, args.kind, nsim, smooth_scale=args.smooth_scale,
                grid=args.grid, nrot=None)
        nsims = get_nsims(args)
        work_delegation(main, nsims, COMM, master_verbose=True)
