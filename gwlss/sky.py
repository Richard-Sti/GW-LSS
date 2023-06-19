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
Sky coordinates related functions.
"""
import healpy
import numpy


def rand_rotation_matrix(deflection=1.0, seed=None):
    """
    Create a random rotation matrix. Taken from [1].

    Parameters
    ----------
    deflection : float, optional
        The magnitude of the rotation. For 0, no rotation; for 1, competely
        random rotation.
    seed : int, optional
        Random seed.

    Returns
    -------
    rotmat : 3-dimensional array

    References
    ----------
    [1] http://blog.lostinmyterminal.com/python/2015/05/12/random-rotation-matrix.html  # noqa
    """
    gen = numpy.random.RandomState(seed)
    randnums = gen.uniform(size=(3,))

    theta, phi, z = randnums

    theta = theta * 2.0 * deflection * numpy.pi  # Rotation about the pole (Z)
    phi = phi * 2.0*numpy.pi                     # Direction of pole deflection
    z = z * 2.0 * deflection                     # Magnitude of pole deflection
    # Compute a vector V used for distributing points over the sphere
    # via the reflection I - V Transpose(V).  This formulation of V
    # will guarantee that if x[1] and x[2] are uniformly distributed,
    # the reflected points will be uniform on the sphere.  Note that V
    # has length sqrt(2) to eliminate the 2 in the Householder matrix.
    r = numpy.sqrt(z)
    Vx, Vy, Vz = V = (numpy.sin(phi) * r,
                      numpy.cos(phi) * r,
                      numpy.sqrt(2.0 - z))
    st = numpy.sin(theta)
    ct = numpy.cos(theta)

    R = numpy.array(((ct, st, 0), (-st, ct, 0), (0, 0, 1)))
    # Construct the rotation matrix  ( V Transpose(V) - I ) R.
    M = (numpy.outer(V, V) - numpy.eye(3)).dot(R)
    return M


def rotate_radec(ra, dec, rotmat):
    """
    Rotate a set of `ra`, `dec` by a rotation matrix. Assume RA is in
    `[0, 2pi]` and dec is in `[-pi/2, pi/2]`.

    Parameters
    ----------
    ra : 1-dimensional array
        Right ascension.
    dec : 1-dimensional array
        Declination.
    rotmat : 3-dimensional array
        Rotation matrix.

    Returns
    -------
    ra_rot : array-like
        Rotated right ascension.
    dec_rot : array-like
        Rotated declination.
    """
    dec_rot, ra_rot = healpy.rotator.rotateDirection(
        rotmat, ra - numpy.pi, numpy.pi / 2 - dec)

    dec_rot -= numpy.pi / 2
    ra_rot += numpy.pi
    return ra_rot, dec_rot
