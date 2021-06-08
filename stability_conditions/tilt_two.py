r"""
Module for computations in second tilt stability as conjecturally described
in [BMT14].

EXAMPLES::

<Lots and lots of examples>

.. TODO::
    - Add examples.
    - Add additional parameters for Z from Bayer-Macri-Stellari and adjust Q.
    - Add some different parameters of k to tests for Q.
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
from sage.calculus.functional import expand
# noinspection PyUnresolvedReferences
from sage.functions.other import imag_part, real_part
# noinspection PyUnresolvedReferences
from sage.symbolic.all import i
# noinspection PyUnresolvedReferences
from sage.rings.all import Integer
# noinspection PyUnresolvedReferences
from sage.rings.infinity import infinity

from .slope import delta
from .varieties.variety import ch


def mu(v, a, b, s):
    r"""
        Computes the slope for second tilt stability.

        TESTS::

            sage: from stability_conditions import *

            sage: v = Element((1, 0, 0, 0))
            sage: tilt_two.mu(v, 1, 0, 0)
            0
            sage: tilt_two.mu(v, 1, -1, 1)
            +Infinity
            sage: tilt_two.mu(v, 2, 1, 4)
            -11

            sage: v = Element((3, -5, 1/2, 1/6, 1/24))
            sage: tilt_two.mu(v, 5, 6, 3/2)
            1515/94

            sage: var('a, b, s, r, c, d, e', domain = RR)
            (a, b, s, r, c, d, e)
            sage: p = tilt_two.mu(Element([r, c, d, e]), a, b, s)
            sage: output = 6*a^2*b*r*s + a^2*b*r - b^3*r - 6*a^2*c*s - a^2*c
            sage: output += 3*b^2*c - 6*b*d + 6*e
            sage: output /= a^2*r - b^2*r + 2*b*c - 2*d
            sage: output /= 3
            sage: bool(p == -output)
            True
    """
    central_charge = z(v, a, b, s)
    if imag_part(central_charge) == 0:
        return infinity
    return -real_part(central_charge) / imag_part(central_charge)


def q(v, w=None, a=0, b=0, k=Integer(1)):
    r"""
    Computes the quadratic form :math:`Q_{a, b, k}(v, w)`.

    INPUT:

    - ``v`` -- Element of the numerical Chow ring
    - ``w`` -- Element of the numerical Chow ring or None. If `w` is None,
      then it assumes `w = v`. (default: `None`)
    - ``a`` -- positive real number (default: `0`)
    - ``b`` -- arbitrary real number (default: `0`)
    - ``k`` -- real number with :math:`k \geq 1` (default: `1`)

    TESTS::

        sage: from stability_conditions import *

        sage: tilt.q(Element([1, 0, -4, 9]), a=0, b=-3)
        -26

        sage: v = Element([2, -1, -1/2, 5/6])
        sage: w = Element([3, -1, -1/2, -1/6])
        sage: tilt.q(v, w, a=4, b=5)
        174
        sage: tilt.q(v, w, a=4, b=5)
        174

        sage: var('s, t', domain = RR)
        (s, t)
        sage: p = tilt.q(Element([3, -2, 0, 2/3]), a=s + 1, b=t - 1)
        sage: bool(p == 4*s^2 + 4*t^2 + 8*s + 4*t + 4)
        True
    """
    if w is None:
        return q(v, v, a, b, k)
    v_tw = ch(v, b)
    w_tw = ch(w, b)
    return (k * a ** 2 * delta(v, w) + 4 * v_tw[2] * w_tw[2]
            - 3 * v_tw[1] * w_tw[3] - 3 * v_tw[3] * w_tw[1])


def wall(v, w, a, b, s):
    r"""
    Returns an equation whose solutions describe the wall between `v` and `w`.

    INPUT:

    - ``v`` - Element
    - ``w`` - Element
    - ``a`` - symbolic expression
    - ``b`` - symbolic expression
    - ``s`` - symbolic expression

    TESTS::

        sage: from stability_conditions import *
        sage: var('a, b, s, r, c, d, e, t, x, y, z', domain = RR)
        (a, b, s, r, c, d, e, t, x, y, z)

        sage: v = Element((1, 0, -4, 8))
        sage: w = Element((1, -1, 1/2, -7/6))
        sage: p = tilt_two.wall(v, w, a, b, s)
        sage: output = -1/2*a^4*s - 1/2*a^2*b^2*s - 1/12*a^4 + 1/6*a^2*b^2
        sage: output += - 1/12*b^4 - 9/2*a^2*b*s + 3/2*a^2*b - 3/2*b^3
        sage: output += - 4*a^2*s + 47/12*a^2 - 79/12*b^2 - 8*b + 2/3
        sage: bool(p == output)
        True

        sage: v = Element((r, c, d, e))
        sage: w = Element((t, x, y, z))
        sage: p = tilt_two.wall(v, w, a, b, s)
        sage: output = -1/2*a^4*c*s*t - 1/2*a^2*b^2*c*s*t + 1/2*a^4*r*s*x
        sage: output += 1/2*a^2*b^2*r*s*x - 1/12*a^4*c*t + 1/6*a^2*b^2*c*t
        sage: output += - 1/12*b^4*c*t + a^2*b*d*s*t + 1/12*a^4*r*x
        sage: output += - 1/6*a^2*b^2*r*x + 1/12*b^4*r*x - a^2*b*r*s*y
        sage: output += - 1/3*a^2*b*d*t + 1/3*b^3*d*t - a^2*d*s*x
        sage: output += 1/3*a^2*b*r*y - 1/3*b^3*r*y + a^2*c*s*y + 1/2*a^2*e*t
        sage: output += - 1/2*b^2*e*t - 1/6*a^2*d*x - 1/2*b^2*d*x
        sage: output += 1/6*a^2*c*y + 1/2*b^2*c*y - 1/2*a^2*r*z + 1/2*b^2*r*z
        sage: output += b*e*x - b*c*z - e*y + d*z
        sage: bool(p == output)
        True
    """
    charge1 = z(v, a, b, s)
    charge2 = z(w, a, b, s)
    return expand((real_part(charge1) * imag_part(charge2)
                   - imag_part(charge1) * real_part(charge2)))


def z(v, a, b, s):
    r"""
    Computes the central charge for second tilt stability.

    TESTS::

        sage: from stability_conditions import *

        sage: v = Element((1, 0, 0, 0))
        sage: tilt_two.z(v, 1, 0, 0)
        -1/2*I
        sage: tilt_two.z(v, 1, -1, 1)
        1
        sage: tilt_two.z(v, 2, 1, 4)
        -3/2*I - 33/2

        sage: v = Element((3, -5, 1/2, 1/6, 1/24))
        sage: tilt_two.z(v, 5, 6, 3/2)
        47*I - 1515/2

        sage: var('a, b, s, r, c, d, e', domain = RR)
        (a, b, s, r, c, d, e)
        sage: p = tilt_two.z(Element([r, c, d, e]), a, b, s)
        sage: output = -a^2*b*r*s - 1/6*a^2*b*r + 1/6*b^3*r + a^2*c*s
        sage: output += 1/6*a^2*c - 1/2*b^2*c - 1/2*I*a^2*r + 1/2*I*b^2*r
        sage: output += - I*b*c + b*d + I*d - e
        sage: bool(p ==  output)
        True
    """
    twisted_ch = ch(v, b)
    return (-twisted_ch[3] + (Integer(1)/Integer(6) + s)*a**2*twisted_ch[1]
            + i*(twisted_ch[2] - a ** 2 / 2 * twisted_ch[0]))
