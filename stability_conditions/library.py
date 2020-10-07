r"""
Collections of short functions required in stability computations.
"""

# ****************************************************************************
#       Copyright (C) 2020 Benjamin Schmidt <schmbe@gmail.com>
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
from sage.arith.misc import xgcd
from sage.functions.other import floor
from sage.rings.all import Integer


def previous_farey(x, n):
    r"""
    Previous element in the Farey sequence for the fraction x in `F_n`.

    Here, `F_n` is the Farey sequence of all fractions with denominator
    smaller than or equal to n. Check
    https://en.wikipedia.org/wiki/Farey_sequence for more details on the
    Farey sequence.

    INPUT:

    - ``x`` -- rational number with denominator at most `n`
    - ``n`` -- positive integer

    EXAMPLES::

        sage: from stability_conditions.all import *

        sage: library.previous_farey(0, 1)
        -1

        sage: library.previous_farey(1/2, 2)
        0

        sage: library.previous_farey(2, 1)
        1

        sage: library.previous_farey(-2/5, 6)
        -1/2

        sage: library.previous_farey(-8/7, 11)
        -7/6
    """
    a = Integer(x.numerator())
    b = Integer(x.denominator())
    d, u, v = xgcd(a, -b)
    # Need to find 0 < u <= n with u maximal
    # so that v/u is actually the Farey neighbor in F_n.
    k = floor((n - u)/b)
    return (v + k*a)/(u + k*b)
