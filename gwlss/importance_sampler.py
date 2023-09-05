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
Importance sampling related functions.
"""

import numpy as np


def bias_from_density(density, bias_type, bias_params):
    """
    Calculate bias from density using a given bias model.

    Parameters
    ----------
    density : float
        Density in h^2 Msun / kpc^3.
    bias_type : str
        Bias model type. Must be one of `simple_bias` or `sigmoid_bias`.
    bias_params : list
        Bias model parameters. Must be a list of floats.

    Returns
    -------
    bias : float
    """
    rho_m = 88.11787915  # matter density in h^2 Msun / kpc^3

    if bias_type == "simple_bias":
        assert len(bias_params) == 1, "Simple bias only has one parameter."
        beta = bias_params[0]
        return (density / rho_m)**beta
    elif bias_type == "sigmoid_bias":
        assert len(bias_params) == 2, "Sigmoid bias has two parameters."
        a_t, t = bias_params
        assert t > 0, "Sigmoid bias parameter `t` must be positive."
        x = (np.log(density / rho_m) - a_t) / t
        return 1.0 / (1.0 + np.exp(-x))
    else:
        raise ValueError(f"Unrecognised bias type `{bias_type}`.")
