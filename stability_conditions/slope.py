r"""
Modules for computations about slope stability.

EXAMPLES::

<Lots and lots of examples>

.. TODO::
    - Add examples.
"""

# ****************************************************************************
#       Copyright (C) 2021 Benjamin Schmidt <schmbe@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
# ****************************************************************************

# noinspection PyUnresolvedReferences
import sage.all
# noinspection PyUnresolvedReferences
from sage.functions.other import imag_part, real_part
# noinspection PyUnresolvedReferences
from sage.symbolic.constants import I
# noinspection PyUnresolvedReferences
from sage.rings.infinity import infinity


def delta(v, w=None):
    r"""
    Computes the Bogomolov discriminant :math:`\Delta(v, w)`.

    TESTS::

        sage: from stability_conditions import *

        sage: v = Element([1, -1, 1/2, 5, 3, 4])
        sage: slope.delta(v, v)
        0
        sage: slope.delta(v)
        0

        sage: v = Element([1, 0, -10])
        sage: slope.delta(v)
        20

        sage: v = Element([2, -1, -1/2])
        sage: w = Element([1, 1, -1/2])
        sage: slope.delta(v, w)
        1/2

        sage: var('r, c, d', domain = RR)
        (r, c, d)
        sage: v = Element([r, c, d])
        sage: bool(slope.delta(v) == c^2 - 2*r*d)
        True
    """
    if w is None:
        return delta(v, v)
    else:
        return v[1]*w[1] - v[0]*w[2] - v[2]*w[0]


def mu(v):
    r"""
    Computes the Mumford slope of `v`.

    TESTS::

        sage: from stability_conditions import *

        sage: slope.mu(Element([1, -1]))
        -1
        sage: slope.mu(Element([3, -1, -1/2, -1/6]))
        -1/3
        sage: slope.mu(Element([-1, 5]))
        -5
        sage: slope.mu(Element([0, 1]))
        +Infinity

        sage: var('r, c', domain = RR)
        (r, c)
        sage: slope.mu(Element([r, c]))
        c/r
    """
    central_charge = z(v)
    if imag_part(central_charge) == 0:
        return infinity
    return -real_part(central_charge)/imag_part(central_charge)


def z(v):
    r"""
    Computes the central charge for slope stability.

    TESTS::

        sage: from stability_conditions import *

        sage: slope.z(Element([4, -2, 3, 4, 8]))
        4*I + 2

        sage: slope.z(Element([0, 4]))
        -4

        sage: slope.z(Element([-5, 1/2, 4, 3/2]))
        -5*I - 1/2
    """
    return -v[1] + I*v[0]
