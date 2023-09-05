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

from .importance_sampler import bias_from_density           # noqa
from .sky import rand_rotation_matrix, rotate_radec         # noqa
from .paths import Paths                                    # noqa
from .utils import load_field                               # noqa

from astropy.cosmology import FlatLambdaCDM


gw170817 = {"RA": 197.45042353,
            "dec": -23.38192721,
            "redshift": 0.0099}

EM_counterpart = {"GW170817": gw170817,
                  }

paths_glamdring = {
    "GW170817_darkPE" : "/mnt/extraspace/rstiskalek/GWLSS/H1L1V1-EXTRACT_POSTERIOR_GW170817-1187008600-400.hdf",  # noqa
    "dumpdir": "/mnt/extraspace/rstiskalek/GWLSS/dumpdir",
    "maindir": "/mnt/extraspace/rstiskalek/GWLSS/",
    }

cosmo_csiborg = FlatLambdaCDM(H0=70.5, Om0=0.307, Ob0=0.04825, Tcmb0=2.728,
                              Neff=3.046, m_nu=0.0, name='csiborg')
