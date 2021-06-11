r"""
Module for computations in tilt stability

This module contains various functionality for computations in tilt
stability as introduced in [Bri08] and further generalized in [AB13] and
[BMT14].

EXAMPLES::

<Lots and lots of examples>

.. TODO::
    - Add examples.
    - Document algorithms for wall computations: explain max rank calc., etc.
    - Document attributes of TiltWall.
    - Check speed issues on wall computation for (0, 3, -1/2) on p(2)
    - Check wall computations with other custom bounds classes.
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
from sage.arith.misc import gcd
# noinspection PyUnresolvedReferences
from sage.functions.other import floor, imag_part, real_part, sqrt
# noinspection PyUnresolvedReferences
from sage.rings.infinity import infinity
# noinspection PyUnresolvedReferences
from sage.rings.all import Integer
# noinspection PyUnresolvedReferences
from sage.rings.rational_field import QQ
# noinspection PyUnresolvedReferences
from sage.structure.all import SageObject
# noinspection PyUnresolvedReferences
from sage.symbolic.all import i

from .library import previous_farey
from .slope import delta, mu
from .varieties.bounds import Bounds
from .varieties.variety import ch, Element


class TiltWall(SageObject):
    r"""
    A wall in tilt stability.

    Walls are either semicircles centered along the :math:`\beta`-axis with
    center :math:`\beta = s` and radius squared `p`, vertical rays given by
    :math:`\beta = s`, or empty.


    INPUT:

    - ``s`` -- center if the wall is a semicircle, constant :math:`\beta`
      value if the wall is vertical, or arbitrary if the wall is
      empty. (default: None)
    - ``p`` -- radius squared of a semicircle, None if the wall is vertical,
      or arbitrary if the wall is empty. (default: `None`)
    - ``vertical`` -- boolean describing whether the wall is vertical.
      (default: `False`)
    - ``empty`` -- boolean describing whether the wall is empty.
      (default: `False`)


    TESTS::

        sage: from stability_conditions import *

        sage: tilt.TiltWall(s=None, p=None, vertical=False, empty=False)
        Traceback (most recent call last):
        ...
        ValueError: Semicircular walls must have radius.

        sage: tilt.TiltWall(s=None, p=1, vertical=False, empty=False)
        Traceback (most recent call last):
        ...
        ValueError: Semicircular walls must have center.

        sage: tilt.TiltWall(s=1, p=1, vertical=True, empty=False)
        Traceback (most recent call last):
        ...
        ValueError: Vertical walls have no radius.

        sage: tilt.TiltWall(s=None, p=None, vertical=True, empty=False)
        Traceback (most recent call last):
        ...
        ValueError: Vertical walls must have s not None.

        sage: tilt.TiltWall(s=1, p=None, vertical=True, empty=True)
        Traceback (most recent call last):
        ...
        ValueError: Vertical walls cannot be empty.

        sage: tilt.TiltWall(s=0, p=-1, vertical=False, empty=False)
        TiltWall(empty=True)

        sage: tilt.TiltWall(s=0, p=1).is_empty
        False

        sage: tilt.TiltWall(s=-4/5, p=-1/100).is_empty
        True

        sage: tilt.TiltWall(empty=True).is_empty
        True

        sage: tilt.TiltWall(s=9, vertical=True).is_empty
        False

        sage: tilt.TiltWall(s=0, p=1).is_vertical
        False

        sage: tilt.TiltWall(empty=True).is_vertical
        False

        sage: tilt.TiltWall(s=7, vertical=True).is_vertical
        True
        """

    def __init__(self, s=None, p=None, vertical=False, empty=False):
        if not vertical and not empty:
            if p is None:
                raise ValueError("Semicircular walls must have radius.")
            elif s is None:
                raise ValueError("Semicircular walls must have center.")

        if vertical:
            if p is not None:
                raise ValueError("Vertical walls have no radius.")
            elif s is None:
                raise ValueError("Vertical walls must have s not None.")
            elif empty:
                raise ValueError("Vertical walls cannot be empty.")

        self.s = s
        self.p = p
        self.is_vertical = vertical
        self.is_empty = empty

        if p is not None and p < 0:
            self.is_empty = True

    def __eq__(self, other):
        r"""
        Walls are identical if they consist of exactly the same points.

        TESTS::

            sage: from stability_conditions import *

            sage: wall1 = tilt.TiltWall(s=1, vertical=True)
            sage: wall2 = tilt.TiltWall(s=1, vertical=True)
            sage: wall3 = tilt.TiltWall(s=-1, vertical=True)
            sage: wall1 == wall2
            True
            sage: wall1 == wall3
            False

            sage: wall4 = tilt.TiltWall(empty=True)
            sage: wall5 = tilt.TiltWall(empty=True)
            sage: wall1 == wall4
            False
            sage: wall4 == wall5
            True

            sage: wall6 = tilt.TiltWall(s=1, p=4)
            sage: wall7 = tilt.TiltWall(s=1, p=4)
            sage: wall1 == wall6
            False
            sage: wall4 == wall6
            False
            sage: wall6 == wall7
            True

            sage: wall8 = tilt.TiltWall(s = -1, p = 10)
            sage: wall6 == wall8
            False
        """
        if self.is_empty:
            return other.is_empty
        elif self.is_vertical:
            return other.is_vertical and self.s == other.s
        else:
            return self.p == other.p and self.s == other.s

    def __ge__(self, other):
        r"""
        Implemented as other < self or self == other.

        TESTS::

            sage: from stability_conditions import *

            sage: wall1 = tilt.TiltWall(s=1, vertical=True)
            sage: wall2 = tilt.TiltWall(s=1, vertical=True)
            sage: wall3 = tilt.TiltWall(s=-1, vertical=True)
            sage: wall1 >= wall2
            True
            sage: wall1 >= wall3
            False

            sage: wall4 = tilt.TiltWall(empty=True)
            sage: wall5 = tilt.TiltWall(s=0, p=-1)
            sage: wall4 >= wall1
            False
            sage: wall1 >= wall4
            True
            sage: wall4 >= wall5
            True

            sage: wall6 = tilt.TiltWall(s=0, p=9)
            sage: wall7 = tilt.TiltWall(s=2, p=1)
            sage: wall8 = tilt.TiltWall(s=-2, p=1/4)
            sage: wall9 = tilt.TiltWall(s=3, p=1)
            sage: wall6 >= wall1
            False
            sage: wall7 >= wall6
            False
            sage: wall8 >= wall6
            False
            sage: wall4 >= wall9
            False
            sage: wall9 >= wall6
            False
            sage: wall7 >= wall8
            False
        """
        return other < self or self == other

    def __gt__(self, other):
        r"""
        Implemented as other < self.

        TESTS::

            sage: from stability_conditions import *

            sage: wall1 = tilt.TiltWall(s=1, vertical=True)
            sage: wall2 = tilt.TiltWall(s=1, vertical=True)
            sage: wall3 = tilt.TiltWall(s=-1, vertical=True)
            sage: wall1 > wall2
            False
            sage: wall1 > wall3
            False

            sage: wall4 = tilt.TiltWall(empty=True)
            sage: wall5 = tilt.TiltWall(s=0, p=-1)
            sage: wall4 > wall1
            False
            sage: wall1 > wall4
            True
            sage: wall4 > wall5
            False

            sage: wall6 = tilt.TiltWall(s=0, p=9)
            sage: wall7 = tilt.TiltWall(s=2, p=1)
            sage: wall8 = tilt.TiltWall(s=-2, p=1/4)
            sage: wall9 = tilt.TiltWall(s=3, p=1)
            sage: wall6 > wall1
            False
            sage: wall7 > wall6
            False
            sage: wall8 > wall6
            False
            sage: wall4 > wall9
            False
            sage: wall9 > wall6
            False
            sage: wall7 > wall8
            False
        """
        return other < self

    def __hash__(self):
        r"""
        Returns a hash value for this object.

        TESTS::

            sage: from stability_conditions import *

            sage: wall1 = tilt.TiltWall(s=1, vertical=True)
            sage: wall2 = tilt.TiltWall(s=1, vertical=True)
            sage: wall3 = tilt.TiltWall(s=-1, vertical=True)
            sage: hash(wall1) == hash(wall2)
            True
            sage: hash(wall1) == hash(wall3)
            False

            sage: wall4 = tilt.TiltWall(empty=True)
            sage: wall5 = tilt.TiltWall(empty=True)
            sage: hash(wall1) == hash(wall4)
            False
            sage: hash(wall4) == hash(wall5)
            True

            sage: wall6 = tilt.TiltWall(s=1, p=4)
            sage: wall7 = tilt.TiltWall(s=1, p=4)
            sage: hash(wall1) == hash(wall6)
            False
            sage: hash(wall4) == hash(wall6)
            False
            sage: hash(wall6) == hash(wall7)
            True

            sage: wall8 = tilt.TiltWall(s = -1, p = 10)
            sage: hash(wall6) == hash(wall8)
            False
        """
        return hash((self.s, self.p, self.is_vertical, self.is_empty))

    def __le__(self, other):
        r"""
        Implemented as self < other or self == other.

        TESTS::

            sage: from stability_conditions import *

            sage: wall1 = tilt.TiltWall(s=1, vertical=True)
            sage: wall2 = tilt.TiltWall(s=1, vertical=True)
            sage: wall3 = tilt.TiltWall(s=-1, vertical=True)
            sage: wall1 <= wall2
            True
            sage: wall1 <= wall3
            False

            sage: wall4 = tilt.TiltWall(empty=True)
            sage: wall5 = tilt.TiltWall(s=0, p=-1)
            sage: wall4 <= wall1
            True
            sage: wall1 <= wall4
            False
            sage: wall4 <= wall5
            True

            sage: wall6 = tilt.TiltWall(s=0, p=9)
            sage: wall7 = tilt.TiltWall(s=2, p=1)
            sage: wall8 = tilt.TiltWall(s=-2, p=1/4)
            sage: wall9 = tilt.TiltWall(s=3, p=1)
            sage: wall6 <= wall1
            False
            sage: wall7 <= wall6
            True
            sage: wall8 <= wall6
            True
            sage: wall4 <= wall9
            True
            sage: wall9 <= wall6
            False
            sage: wall7 <= wall8
            False
        """
        return self < other or self == other

    # noinspection DuplicatedCode
    def __lt__(self, other):
        r"""
        Checks whether this wall is strictly contained in the wall `other`.

        Semicircular walls are compared by inclusion in the enclosed
        semidisk of the larger wall. Vertical walls are considered to be
        larger than all semicircular walls that they do not intersect for
        :math:`\alpha > 0`. Empty walls are considered smaller than all
        non-empty walls.

        TESTS::

            sage: from stability_conditions import *

            sage: wall1 = tilt.TiltWall(s=1, vertical=True)
            sage: wall2 = tilt.TiltWall(s=1, vertical=True)
            sage: wall3 = tilt.TiltWall(s=-1, vertical=True)
            sage: wall1 < wall2
            False
            sage: wall1 < wall3
            False

            sage: wall4 = tilt.TiltWall(empty=True)
            sage: wall5 = tilt.TiltWall(s=0, p=-1)
            sage: wall4 < wall1
            True
            sage: wall1 < wall4
            False
            sage: wall4 < wall5
            False

            sage: wall6 = tilt.TiltWall(s=0, p=9)
            sage: wall7 = tilt.TiltWall(s=2, p=1)
            sage: wall8 = tilt.TiltWall(s=-2, p=1/4)
            sage: wall9 = tilt.TiltWall(s=3, p=1)
            sage: wall6 < wall1
            False
            sage: wall7 < wall6
            True
            sage: wall8 < wall6
            True
            sage: wall4 < wall9
            True
            sage: wall9 < wall6
            False
            sage: wall7 < wall8
            False
        """
        if self.is_empty:
            return not other.is_empty
        elif other.is_empty:
            return False
        elif self.is_vertical:
            return False
        elif other.is_vertical:
            if self.s == other.s:
                return False
            elif self.s > other.s:
                # Need to check whether self.s - sqrt(self.p) >= other.s
                return self.p <= (self.s - other.s) ** 2
            elif self.s < other.s:
                # Need to check whether self.s + sqrt(self.p) <= other.s
                return self.p <= (other.s - self.s) ** 2
        elif self.s == other.s:
            return self.p < other.p
        elif self.s > other.s:
            # In this case self.s + sqrt(self.p) < other.s + sqrt(other.p) is
            # the only thing to check, since it implies
            # self.s - sqrt(self.p) > other.s - sqrt(other.p).
            # Some basic math on inequalities allows us to compare without
            # computing square roots as follows:
            disc = (other.p - self.s ** 2 - other.s ** 2 - self.p +
                    2 * self.s * other.s) / (2 * self.s - 2 * other.s)
            if disc < 0:
                return False
            else:
                return self.p <= disc ** 2
        elif self.s < other.s:
            # In this case self.s - sqrt(self.p) > other.s - sqrt(other.p) is
            # the only thing to check, since it implies
            # self.s + sqrt(self.p) < other.s + sqrt(other.p).
            # Some basic math on inequalities allows us to compare without
            # computing square roots as follows:
            disc = (other.p - self.s ** 2 - other.s ** 2 - self.p +
                    2 * self.s * other.s) / (2 * other.s - 2 * self.s)
            if disc < 0:
                return False
            else:
                return self.p <= disc ** 2

    def __ne__(self, other):
        r"""
        Implemented as the opposite of self == other.

        TESTS::

            sage: from stability_conditions import *

            sage: wall1 = tilt.TiltWall(s=1, vertical=True)
            sage: wall2 = tilt.TiltWall(s=1, vertical=True)
            sage: wall3 = tilt.TiltWall(s=-1, vertical=True)
            sage: wall1 != wall2
            False
            sage: wall1 != wall3
            True

            sage: wall4 = tilt.TiltWall(empty=True)
            sage: wall5 = tilt.TiltWall(empty=True)
            sage: wall1 != wall4
            True
            sage: wall4 != wall5
            False

            sage: wall6 = tilt.TiltWall(s=1, p=4)
            sage: wall7 = tilt.TiltWall(s=1, p=4)
            sage: wall1 != wall6
            True
            sage: wall4 != wall6
            True
            sage: wall6 != wall7
            False

            sage: wall8 = tilt.TiltWall(s = -1, p = 10)
            sage: wall6 != wall8
            True
        """
        return not self == other

    def _repr_(self):
        r"""
        Returns the object as a string in the syntax that would be used to
        create it newly.

        TESTS::

            sage: from stability_conditions import *

            sage: tilt.TiltWall(s=-1/2, p=None, vertical=True, empty=False)
            TiltWall(s=-1/2, vertical=True)

            sage: tilt.TiltWall(s=None, p=None, vertical=False, empty=True)
            TiltWall(empty=True)

            sage: tilt.TiltWall(s=1, p=1, vertical=False, empty=False)
            TiltWall(s=1, p=1)

            sage: tilt.TiltWall(s=1, p=1, vertical=False, empty=True)
            TiltWall(empty=True)
        """
        if self.is_empty:
            return "TiltWall(empty=True)"
        elif self.is_vertical:
            return "TiltWall(s=%s, vertical=True)" % str(self.s)
        else:
            return "TiltWall(s=%s, p=%s)" % (str(self.s), str(self.p))

    def _latex_(self):
        r"""
        Returns a latex expression for the defining equation of the wall.

        TESTS::

            sage: from stability_conditions import *

            sage: latex(tilt.TiltWall(s=-1/2, p=None, vertical=True, \
                                      empty=False))
            \beta = -1/2

            sage: latex(tilt.TiltWall(s=None, p=None, vertical=False, \
                                      empty=True))
            \emptyset

            sage: latex(tilt.TiltWall(s=1, p=1, vertical=False, empty=False))
            \alpha^2 + \left(\beta - 1\right)^2 = 1

            sage: latex(tilt.TiltWall(s=-1, p=3, vertical=False, empty=False))
            \alpha^2 + \left(\beta + 1\right)^2 = 3

            sage: latex(tilt.TiltWall(s=0, p=2, vertical=False, empty=False))
            \alpha^2 + \beta^2 = 2

            sage: latex(tilt.TiltWall(s=1, p=1, vertical=False, empty=True))
            \emptyset
        """
        if self.is_empty:
            return "\\emptyset"
        elif self.is_vertical:
            return "\\beta = %s" % str(self.s)
        elif self.s == 0:
            return "\\alpha^2 + \\beta^2 = %s" % str(self.p)
        elif self.s < 0:
            return ("\\alpha^2 + \\left(\\beta + %s\\right)^2 = %s"
                    % (str(-self.s), str(self.p)))
        else:  # elif self.s > 0
            return ("\\alpha^2 + \\left(\\beta - %s\\right)^2 = %s"
                    % (str(self.s), str(self.p)))

    def equation(self, a, b):
        r"""
        Substitutes :math:`a = \alpha` and :math:`b = \beta` into the
        equation of the wall.

        TESTS::

            sage: from stability_conditions import *

            sage: var('a, b', domain = RR)
            (a, b)
            sage: tilt.TiltWall(empty=True).equation(a, b)
            -1
            sage: tilt.TiltWall(s=-1, vertical=True).equation(a, b)
            b + 1
            sage: tilt.TiltWall(s=-3, p=1).equation(a, b)
            a^2 + (b + 3)^2 - 1
        """
        if self.is_empty:
            return -1
        elif self.is_vertical:
            return b - self.s
        else:
            return a**2 + (b - self.s)**2 - self.p


def bar_beta(v):
    r"""
    Returns the number

    .. math:: \overline{\beta}(v) :=
              \begin{cases}
              \frac{v_1 - \sqrt{\Delta(v)}}{v_0} &\text{, if } v_0 \neq 0 \\
              \frac{v_2}{v_1} &\text{, if } v_0 = 0.
              \end{cases}

    TESTS::

        sage: from stability_conditions import *
        sage: var('r, c, d', domain=RR)
        (r, c, d)

        sage: v = Element([r, c, d])
        sage: bool(tilt.bar_beta(v) == (c - sqrt(c^2 - 2*d*r))/r)
        True

        sage: v = Element([0, c, d])
        sage: bool(tilt.bar_beta(v) == d/c)
        True

        sage: v = Element([1, 0, -2])
        sage: tilt.bar_beta(v)
        -2
        sage: tilt.bar_beta(-v)
        2

        sage: v = Element([3, -2, 0, 2/3])
        sage: tilt.bar_beta(v)
        -4/3
        sage: tilt.bar_beta(-v)
        0

        sage: v = Element([0, 1, -1/2])
        sage: tilt.bar_beta(v)
        -1/2
        sage: tilt.bar_beta(-v)
        -1/2
    """
    if v[0] == 0:
        return v[2] / v[1]
    else:
        return (v[1] - sqrt(delta(v))) / v[0]


def beta_minus(v):
    r"""
    Returns the smallest solutions to the equation
    :math:`\nu_{0, \beta}(v) = 0`.


    TESTS::

        sage: from stability_conditions import *
        sage: var('r, c, d', domain=RR)
        (r, c, d)

        sage: v = Element([r, c, d])
        sage: bool(tilt.beta_minus(v) == (c - sqrt(c^2 - 2*d*r))/r)
        Traceback (most recent call last):
        ...
        ValueError: The rank needs to have a sign or be zero.

        sage: assume(r > 0)
        sage: bool(tilt.beta_minus(v) == (c - sqrt(c^2 - 2*d*r))/r)
        True

        sage: forget()
        sage: assume(r < 0)
        sage: bool(tilt.beta_minus(v) == (c + sqrt(c^2 - 2*d*r))/r)
        True

        sage: v = Element([0, c, d])
        sage: bool(tilt.beta_minus(v) == d/c)
        True

        sage: v = Element([1, 0, -2])
        sage: tilt.beta_minus(v)
        -2

        sage: v = Element([3, -2, 0, 2/3])
        sage: tilt.beta_minus(-v)
        -4/3

        sage: v = Element([0, 1, -1/2])
        sage: tilt.beta_minus(v)
        -1/2

        sage: forget()
    """
    if v[0] == 0:
        return v[2] / v[1]
    elif v[0] > 0:
        return bar_beta(v)
    elif v[0] < 0:
        return bar_beta(-v)
    else:
        raise ValueError("The rank needs to have a sign or be zero.")


def beta_plus(v):
    r"""
    Returns the largest solutions to the equation
    :math:`\nu_{0, \beta}(v) = 0`.

    TESTS::

        sage: from stability_conditions import *
        sage: var('r, c, d', domain=RR)
        (r, c, d)

        sage: v = Element([r, c, d])
        sage: bool(tilt.beta_plus(v) == (c - sqrt(c^2 - 2*d*r))/r)
        Traceback (most recent call last):
        ...
        ValueError: The rank needs to have a sign or be zero.

        sage: assume(r > 0)
        sage: bool(tilt.beta_plus(v) == (c + sqrt(c^2 - 2*d*r))/r)
        True

        sage: forget()
        sage: assume(r < 0)
        sage: bool(tilt.beta_plus(v) == (c - sqrt(c^2 - 2*d*r))/r)
        True

        sage: v = Element([0, c, d])
        sage: bool(tilt.beta_plus(v) == d/c)
        True

        sage: v = Element([1, 0, -2])
        sage: tilt.beta_plus(v)
        2

        sage: v = Element([3, -2, 0, 2/3])
        sage: tilt.beta_plus(-v)
        0

        sage: v = Element([0, 1, -1/2])
        sage: tilt.beta_plus(v)
        -1/2

        sage: forget()
    """
    if v[0] == 0:
        return v[2] / v[1]
    elif v[0] > 0:
        return bar_beta(-v)
    elif v[0] < 0:
        return bar_beta(v)
    else:
        raise ValueError("The rank needs to have a sign or be zero.")


def hyperbola(v, a, b):
    r"""
    Computes the real part of the central charge of tilt stability.

    TESTS::

        sage: from stability_conditions import *

        sage: v = Element((1, -1, 1/2, 5, 5, 5))
        sage: tilt.hyperbola(v, 1, -2)
        0
        sage: tilt.hyperbola(v, 0, 0)
        -1/2

        sage: var('a, b, r, c, d', domain = RR)
        (a, b, r, c, d)
        sage: p = tilt.hyperbola(Element((r, c, d)), a, b)
        sage: bool(p == 1/2*a^2*r - 1/2*b^2*r + b*c - d)
        True
    """
    return real_part(z(v, a, b))


def nu(v, a, b):
    r"""
    Computes the slope :math:`\nu_{a, b}(v)`.

    TESTS::

        sage: from stability_conditions import *

        sage: v = Element([1, 0, 0])
        sage: tilt.nu(v, 1, 0)
        +Infinity
        sage: tilt.nu(v, 1, -1)
        0
        sage: tilt.nu(v, 2, 1)
        3/2

        sage: v = Element([3, -5, 1/2, 1/6, 1/24])
        sage: tilt.nu(v, 5, 6)
        -47/23

        sage: var('a, b, r, c, d', domain = RR)
        (a, b, r, c, d)
        sage: p = tilt.nu(Element([r, c, d]), a, b)
        sage: bool(p == 1/2*(a^2*r - b^2*r + 2*b*c - 2*d)/(b*r - c))
        True
    """
    central_charge = z(v, a, b)
    if imag_part(central_charge) == 0:
        return infinity
    return -real_part(central_charge) / imag_part(central_charge)


def q(v, w=None, a=0, b=0):
    r"""
    Computes the quadratic form :math:`Q_{a, b}(v, w)`.

    INPUT:

    - ``v`` -- Element of the numerical Chow ring
    - ``w`` -- Element of the numerical Chow ring or None. If `w` is None,
               then it assumes `w = v`. (default: `None`)
    - ``a`` -- positive real number (default: `0`)
    - ``b`` -- arbitrary real number (default: `0`)

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
        return q(v, v, a, b)
    v_tw = ch(v, b)
    w_tw = ch(w, b)
    return (a ** 2 * delta(v, w) + 4 * v_tw[2] * w_tw[2]
            - 3 * v_tw[1] * w_tw[3] - 3 * v_tw[3] * w_tw[1])


def q_wall(v):
    r"""
    Computes the numerical wall defined by :math:`Q_{a, b}(v) = 0`.

    TESTS::

        sage: from stability_conditions import *

        sage: tilt.q_wall(Element([1, 0, -4, 8]))
        TiltWall(s=-3, p=1)

        sage: tilt.q_wall(Element([3, -2, 0, 2/3]))
        TiltWall(s=-3/2, p=1/4)

        sage: tilt.q_wall(Element([1, -10, 48, -431/3]))
        TiltWall(s=-49/4, p=17/16)

        sage: tilt.q_wall(Element([2, 10, 20, 50/3]))
        TiltWall(empty=True)
    """
    return wall(v, [v[1], 2 * v[2], 3 * v[3]])


def wall(v, w):
    r"""
    Returns wall between `v` and `w`.

    TESTS::

        sage: from stability_conditions import *

        sage: v = Element([1, 0, -4, 8])
        sage: w = Element([0, -8, 24])
        sage: tilt.wall(v, w)
        TiltWall(s=-3, p=1)
        sage: tilt.wall(w, v)
        TiltWall(s=-3, p=1)

        sage: v = Element([1, -1, 1/2])
        sage: w = Element([1, -1, -1/2])
        sage: tilt.wall(v, w)
        TiltWall(s=-1, vertical=True)

        sage: v = Element([0, 0, 1])
        sage: w = Element([2, 2, 1])
        sage: tilt.wall(v, w)
        TiltWall(s=1, vertical=True)
        sage: tilt.wall(w, v)
        TiltWall(s=1, vertical=True)

        sage: v = Element([1, 0, -1])
        sage: w = Element([1, -1, -1/2])
        sage: tilt.wall(v, w)
        TiltWall(empty=True)
    """
    if v[1] * w[0] - v[0] * w[1] == 0:
        if v[0] == 0:
            if w[0] == 0:
                return TiltWall(empty=True)
            else:
                return TiltWall(s=mu(w), vertical=True)
        else:
            return TiltWall(s=mu(v), vertical=True)

    s = (v[2] * w[0] - v[0] * w[2]) / (v[1] * w[0] - v[0] * w[1])
    p = ((v[2] ** 2 * w[0] ** 2 + v[0] ** 2 * w[2] ** 2
          - 2 * v[0] * v[2] * w[0] * w[2] - 2 * v[1] * v[2] * w[0] * w[1]
          - 2 * v[0] * v[1] * w[1] * w[2] + 2 * v[0] * v[2] * w[1] ** 2
          + 2 * v[1] ** 2 * w[0] * w[2]) / (v[1] * w[0] - v[0] * w[1]) ** 2)

    return TiltWall(s, p)


def walls_left(v, var, s=None, b=None, bounds=None):
    r"""
    Computes walls left of the vertical wall in tilt stability.

    This function assumes v[0] != 0. It computes those walls for `v`
    that are larger than or equal to the numerical wall with center
    :math:`\beta = s`. If `s` is None, it computes those walls that intersect
    the vertical ray :math:`\beta = b`. If both `b` and `s` are None and
    :math:`\beta^{-}(v)` is rational, it computes all walls to the left of the
    vertical wall.

    If `bounds` is None, then we will simply use the Bogomolov inequality to
    bound the second Chern characters of destabilizing subobjects. For more
    sophisticated bounds, one can use a different custom `Bounds` class. See
    the documentation for the `bounds` module for more details.

    INPUT:

    - ``v`` -- Element in the numerical Chow ring.
    - ``var`` -- Variety.
    - ``s`` -- rational number or None (default: `None`).
    - ``b`` -- rational number or None (default: `None`).
    - ``bounds`` -- Instance of bounds class that specifies Chern character
      bounds.

    If both `s` and `b` are different from None, then `b` is ignored.

    OUTPUT: Returns a dictionary where keys are TiltWall objects and the
    items are sets of those potentially destabilizing semistable
    subobjects or quotients which have positive rank and whose slope is
    smaller than the slope of v.

    ALGORITHM:

    WARNING:

    If ``bounds`` is None, then a new instance of it is created. If you use a
    lot of wall computations, it might be wise to create an instance manually
    once and keep passing it as a parameter.

    TESTS::

        sage: from stability_conditions import *

        sage: output = {tilt.TiltWall(s=-5/2, p=9/4): {Element((1, -1, 1/2))}}
        sage: tilt.walls_left(Element((1, 0, -2)), p(2)) == output
        True

        sage: w = Element((1, -1, 1/2))
        sage: output = {tilt.TiltWall(s=-5/2, p=9/4): {w}}
        sage: tilt.walls_left(Element((1, 0, -2)), p(2), b=-1) == output
        True

        sage: tilt.walls_left(Element((3, -1, -1/2, -1/6)), p(3))
        {}

        sage: v = Element((3, -2, 0, 2/3))
        sage: w = Element((1, -1, 1/2))
        sage: output =  {tilt.TiltWall(s=-3/2, p=1/4): {w, 2*w, 3*w, 4*w}}
        sage: tilt.walls_left(v, p(3), s=-4/3) == output
        True

        sage: var = Variety(3, (1, 1, 1/2, 1/6), (1, 0, 0, 0))
        sage: v = Element((2, 0, -3/2))
        sage: w = Element((2, -1, 0))
        sage: onlyWall = tilt.TiltWall(s=-3/2, p=3/4)
        sage: tilt.walls_left(v, var, s=-3/2) == {onlyWall: {w}}
        True

        sage: tilt.walls_left(Element((1, 0, -3)), p(2))
        Traceback (most recent call last):
        ...
        ValueError: There are infinitely many such walls.

        sage: tilt.walls_left(Element((0, 1, 1/2)), p(2))
        Traceback (most recent call last):
        ...
        ValueError: Use the appropriate function for rank zero objects.

        sage: tilt.walls_left(Element((1, 0, -3)), p(2), s=1)
        Traceback (most recent call last):
        ...
        ValueError: There is no numerical wall with this center on this
        side of the vertical wall.
    """
    if v[0] == 0:
        raise ValueError("Use the appropriate function for rank zero objects.")

    if v[0] < 0:
        v = -Element(v[0:3])
    else:
        v = Element(v[0:3])

    if bounds is None:
        bounds = Bounds(var)

    # First we determine s if necessary and check the input.
    if s is None:
        if b is None:
            s = beta_minus(v)
            if s not in QQ:
                raise ValueError("There are infinitely many such walls.")
        elif b >= mu(v):
            # Clearly, no walls to the left of the vertical wall can
            # intersect `\beta = b` in this case.
            return {}
        else:
            s = wall(v, var.o(b)).s
    elif s > beta_minus(v):
        raise ValueError("There is no numerical wall with this center on "
                         "this side of the vertical wall.")

    # Find the maximal rank for a destabilizing subobject or quotient and
    # the minimal radius squared of our walls.
    if s == beta_minus(v):
        p = 0
        n = Integer(s.denominator())
        m = Integer(1/var.gens[2])
        r_max = (n * ch(v, s)[1] - 1) ** 2 * m / gcd(2 * n ** 2, m)
    else:
        p = (ch(v, s)[1] ** 2 - delta(v)) / v[0] ** 2
        r_max = floor(v[0] / 2 + sqrt(v[0] ** 2 / 4 + delta(v) / 4 / p))

    walls = {}
    min_wall = TiltWall(s=s, p=p)
    m = previous_farey(mu(v), r_max)
    while m >= s and (m - s) ** 2 >= p:
        r = m.denominator()
        c = m.numerator()
        while r <= r_max:
            d = bounds.ch_max(r, c)
            w = Element([r, c, d])
            u = v - w
            current_wall = wall(v, w)
            while min_wall <= current_wall and current_wall.p > 0:
                if ((u[0] > 0 and current_wall.s <= mu(u)) or
                        u[0] == 0 or
                        (u[0] < 0 and current_wall.s >= mu(u))):
                    if bounds.satisfies_bounds(*u):
                        if current_wall in walls.keys():
                            walls[current_wall].add(w)
                        else:
                            walls[current_wall] = {w}
                d -= var.gens[2]
                w = Element([r, c, d])
                u = v - w
                current_wall = wall(v, w)
            r += m.denominator()
            c += m.numerator()
        m = previous_farey(m, r_max)
    return walls


def walls_right(v, var, s=None, b=None, bounds=None):
    r"""
    Computes walls right of the vertical wall in tilt stability.

    This function assumes v[0] != 0. It computes those walls for `-v`
    that are larger than or equal to the numerical wall with center
    :math:`\beta = s`. If `s` is None, it computes those walls that intersect
    the vertical ray :math:`\beta = b`. If both `b` and `s` are None and
    :math:`\beta^{-}(v)` is rational, it computes all walls to the left of the
    vertical wall.

    If `bounds` is None, then we will simply use the Bogomolov inequality to
    bound the second Chern characters of destabilizing subobjects. For more
    sophisticated bounds, one can use a different custom `Bounds` class. See
    the documentation for the `bounds` module for more details.

    INPUT:

    - ``v`` -- Element in the numerical Chow ring.
    - ``var`` -- Variety.
    - ``s`` -- rational number or None (default: `None`).
    - ``b`` -- rational number or None (default: `None`).
    - ``bounds`` -- Instance of bounds class that specifies Chern character
      bounds.

    If both `s` and `b` are different from None, then `b` is ignored.

    OUTPUT: Returns a dictionary where keys are TiltWall objects and the
    items are sets of those potentially destabilizing semistable
    subobjects or quotients which have negative rank and whose slope is
    smaller than the slope of -v.

    ALGORITHM:

    TESTS::

        sage: from stability_conditions import *

        sage: output = {tilt.TiltWall(s=5/2, p=9/4): {-Element((1, 1, 1/2))}}
        sage: tilt.walls_right(Element((1, 0, -2)), p(2)) == output
        True

        sage: output = {tilt.TiltWall(s=5/2, p=9/4): {-Element((1, 1, 1/2))}}
        sage: tilt.walls_right(Element((1, 0, -2)), p(2), b=1) == output
        True

        sage: w = Element((-1, 0, 0))
        sage: only_wall = tilt.TiltWall(s=1/2, p=1/4)
        sage: output = {only_wall: {w, 2*w, 3*w, 4*w}}
        sage: tilt.walls_right(Element((3, -1, -1/2, -1/6)), p(3)) == output
        True

        sage: tilt.walls_right(Element((3, -2, 0, 2/3)), p(3), s=0)
        {}

        sage: var = Variety(3, (1, 1, 1/2, 1/6), (1, 0, 0, 0))
        sage: v = Element((2, 0, -3/2))
        sage: w = -Element((2, 1, 0))
        sage: only_wall = tilt.TiltWall(s=3/2, p=3/4)
        sage: tilt.walls_right(v, var, s=3/2) == {only_wall: {w}}
        True

        sage: tilt.walls_right(Element((2, 0, -3)), p(2))
        Traceback (most recent call last):
        ...
        ValueError: There are infinitely many such walls.

        sage: tilt.walls_right(Element((0, 10, 3)), p(2))
        Traceback (most recent call last):
        ...
        ValueError: Use the appropriate function for rank zero objects.

        sage: tilt.walls_right(Element((2, 0, -3)), p(2), s=-1)
        Traceback (most recent call last):
        ...
        ValueError: There is no numerical wall with this center on this
        side of the vertical wall.

        sage: var = Variety(3, (1, 1, 1/3, 1/3), (1, 1, 2/3, 1/3))
        sage: wall = tilt.TiltWall(s=1/2, p=1/4)
        sage: w = Element((-1, 0, 0))
        sage: output = {wall: {w, 2*w, 3*w, 4*w}}
        sage: tilt.walls_right(Element((3, -1, -1/2)), var) == output
        True
    """
    if bounds is None:
        bounds = Bounds(var)

    if s is None:
        if b is None:
            walls_dual = walls_left(v.dual(), var, bounds=bounds)
        else:
            walls_dual = walls_left(v.dual(), var, None, -b, bounds)
    elif b is None:
        walls_dual = walls_left(v.dual(), var, -s, bounds=bounds)
    else:
        walls_dual = walls_left(v.dual(), var, -s, -b, bounds)

    walls = {}
    for key, val in walls_dual.items():
        new_wall = TiltWall(-key.s, key.p)
        walls[new_wall] = set()
        for w in val:
            walls[new_wall].add(-w.dual())

    return walls


def walls_torsion(v, var, p=None, b=None, bounds=None):
    r"""
    Compute walls for v in tilt stability if v[0] == 0.

    This function computes those walls for `v` that are larger than or
    equal to the numerical wall with radius squared `p`. If `p` is None,
    it computes those walls that intersect the vertical ray :math:`\beta = b`.
    If both `b` and `p` are None it computes all walls to the left of
    the vertical wall.

    If `bounds` is None, then we will simply use the Bogomolov inequality to
    bound the second Chern characters of destabilizing subobjects. For more
    sophisticated bounds, one can use a different custom `Bounds` class. See
    the documentation for the `bounds` module for more details.

    INPUT:

    - ``v`` -- Element in the numerical Chow ring.
    - ``var`` -- Variety.
    - ``p`` -- rational number or None (default: `None`).
    - ``b`` -- rational number or None (default: `None`).
    - ``bounds`` -- Instance of bounds class that specifies Chern character
      bounds.

    If both `p` and `b` are different from None, then `b` is ignored. If
    `p < 0`, then the function assumes p is 0.

    OUTPUT: Returns a dictionary where keys are TiltWall objects and the
    items are sets of those potentially destabilizing semistable
    subobjects or quotients which have positive rank.

    ALGORITHM:

    TESTS::

        sage: from stability_conditions import *

        sage: v = Element((0, 1, -1/2))
        sage: w = Element((1, 0, 0))
        sage: only_wall = tilt.TiltWall(s=-1/2, p=1/4)
        sage: tilt.walls_torsion(v, p(2)) == {only_wall: {w}}
        True

        sage: v = Element((0, 1, -1/2))
        sage: tilt.walls_torsion(v, p(2), p=1)
        {}

        sage: v = Element((0, 1, -1/2))
        sage: w = Element((1, 0, 0))
        sage: only_wall = tilt.TiltWall(s=-1/2, p=1/4)
        sage: tilt.walls_torsion(v, p(2), b=0) == {only_wall: {w}}
        True

        sage: v = Element((0, 4, 8))
        sage: wall1 = tilt.TiltWall(s=2, p=4)
        sage: wall2 = tilt.TiltWall(s=2, p=2)
        sage: wall3 = tilt.TiltWall(s=2, p=1)
        sage: output = {wall1: set(), wall2: set(), wall3: set()}
        sage: output[wall1].add(Element((1, 4, 8)))
        sage: output[wall2].add(Element((1, 4, 7)))
        sage: output[wall3].add(Element((1, 5, 17/2)))
        sage: output[wall3].add(Element((1, 3, 9/2)))
        sage: output[wall3].add(Element((2, 6, 9)))
        sage: tilt.walls_torsion(v, p(2)) == output
        True

        sage: var = Variety(3, (1, 1, 1/2, 1/6), (1, 0, 0, 0))
        sage: v = Element((0, 2, 1/2))
        sage: wall1 = tilt.TiltWall(s=1/4, p=9/16)
        sage: wall2 = tilt.TiltWall(s=1/4, p=1/16)
        sage: output = {wall1: set(), wall2: set()}
        sage: output[wall1].add(Element((1, 1, 1/2)))
        sage: output[wall2].add(Element((1, 2, 1/2)))
        sage: output[wall2].add(Element((2, 2, 1/2)))
        sage: output[wall2].add(Element((3, 2, 1/2)))
        sage: output[wall2].add(Element((4, 2, 1/2)))
        sage: tilt.walls_torsion(v, var, p=1/16) == output
        True

        sage: tilt.walls_torsion(Element((1, 0, 0)), p(2))
        Traceback (most recent call last):
        ...
        ValueError: Use the appropriate function for non-torsion objects.
    """
    if v[0] != 0:
        raise ValueError("Use the appropriate function for non-torsion "
                         "objects.")

    if v[1] < 0:
        raise ValueError("There are no tilt-semistable objects with this "
                         "Chern character.")
    elif v[1] == 0:
        return {}

    if bounds is None:
        bounds = Bounds(var)

    if len(v) > 3:
        return walls_torsion(Element(v[0:3]), var, p, b, bounds)

    # Determine p if necessary.
    if p is None:
        if b is None:
            p = 0
        else:
            p = wall(v, var.o(b)).p
    elif p < 0:
        p = 0

    s = v[2]/v[1]

    # Find the maximal rank for a destabilizing subobject or quotient.
    if p == 0:
        n = Integer(s.denominator())
        m = Integer(1 / var.gens[2])
        r_max = (n * v[1] - 1) ** 2 * m / gcd(2 * n ** 2, m)
    else:
        r_max = floor(sqrt(delta(v) / 4 / p))

    if r_max == 0:
        return {}

    # Find maximal slope for a destabilizing subobject or quotient.
    m = v[1] + s
    while m.denominator() > r_max:
        m = previous_farey(m, m.denominator())

    walls = {}
    min_wall = TiltWall(s=s, p=p)

    while m >= s and (m - s) ** 2 >= p:
        r = m.denominator()
        c = m.numerator()
        while r <= r_max:
            d = bounds.ch_max(r, c)
            w = Element((r, c, d))
            u = v - w
            current_wall = wall(v, w)
            while min_wall <= current_wall and current_wall.p > 0:
                if s >= mu(u) and bounds.satisfies_bounds(*u):
                    if current_wall in walls.keys():
                        walls[current_wall].add(w)
                    else:
                        walls[current_wall] = {w}
                d -= var.gens[2]
                w = Element([r, c, d])
                u = v - w
                current_wall = wall(v, w)
            r += m.denominator()
            c += m.numerator()
        m = previous_farey(m, r_max)
    return walls


def walls_vertical(v, var, bounds=None):
    r"""
    Computes vertical walls for `v` in tilt stability.

    This function assumes v[0] != 0.

    OUTPUT: Returns a set of potentially destabilizing semistable
    subobjects or quotients `w`. It gives exactly half of these objects
    as follows. If `w` has rank zero, then instead `-v - w` is in the
    set. Assume that both `w` and `-v - w` have negative rank. If `-w`
    would not numerically destabilize a sheaf with class `-v` in
    Gieseker stability, then `w` is in the set. If it does, then
    `v + w` would numerically destabilize a sheaf with class `-v` in
    Gieseker stability and `-v - w` is added to the set.

    If `bounds` is None, then we will simply use the Bogomolov inequality to
    bound the second Chern characters of destabilizing subobjects. For more
    sophisticated bounds, one can use a different custom `Bounds` class. See
    the documentation for the `bounds` module for more details.

    INPUT:

    - ``v`` -- Element in the numerical Chow ring.
    - ``var`` -- Variety.
    - ``bounds`` -- Instance of bounds class that specifies Chern character
      bounds.

    ALGORITHM:

    TESTS::

        sage: from stability_conditions import *

        sage: output = set()
        sage: output.add(-Element((2, 0, -1)))
        sage: output.add(-Element((2, 0, -2)))
        sage: output.add(-Element((2, 0, 0)))
        sage: output.add(-Element((1, 0, -1)))
        sage: output.add(-Element((1, 0, 0)))
        sage: tilt.walls_vertical(Element((2, 0, -2)), p(2)) == output
        True

        sage: output = set()
        sage: output.add(-Element((4, -1, -1/2)))
        sage: output.add(-Element((4, -1, -3/2)))
        sage: tilt.walls_vertical(Element((4, -1, -3/2)), p(2)) == output
        True

        sage: output = set()
        sage: output.add(-Element((3, -4, 2)))
        sage: output.add(-Element((6, -8, 4)))
        sage: output.add(-Element((6, -8, 5)))
        sage: output.add(-Element((9, -12, 6)))
        sage: output.add(-Element((9, -12, 7)))
        sage: output.add(-Element((9, -12, 8)))
        sage: output.add(-Element((12, -16, 8)))
        sage: output.add(-Element((12, -16, 9)))
        sage: output.add(-Element((12, -16, 10)))
        sage: tilt.walls_vertical(Element((12, -16, 8, 1)), p(2)) == output
        True

        sage: output = set()
        sage: output.add(-Element((1, 0, 0)))
        sage: output.add(-Element((1, 0, -1/3)))
        sage: output.add(-Element((2, 0, 0)))
        sage: output.add(-Element((2, 0, -1/3)))
        sage: output.add(-Element((2, 0, -2/3)))
        sage: output.add(-Element((2, 0, -1)))
        sage: var = Variety(3, (1, 1, 1/3, 1/3), (1, 1, 2/3, 1/3))
        sage: tilt.walls_vertical(Element((2, 0, -1)), var) == output
        True

        sage: tilt.walls_vertical(Element((0, 10, 3)), p(2))
        Traceback (most recent call last):
        ...
        ValueError: Use the appropriate function for rank zero objects.
    """
    if v[0] == 0:
        raise ValueError("Use the appropriate function for rank zero objects.")
    elif v[0] < 0:
        v = -Element(v[0:3])
    else:
        v = Element(v[0:3])

    if bounds is None:
        bounds = Bounds(var)

    walls = set()
    for n in range(1, gcd(v[0], v[1]) + 1):
        r = mu(v).denominator() * n
        c = mu(v).numerator() * n
        d = bounds.ch_max(r, c)
        w = Element((r, c, d))
        u = v - w
        if r == v[0]:
            while d >= v[2] - v[1] ** 2 / 2 / v[0] + c ** 2 / 2 / r:
                walls.add(-w)
                d -= var.gens[2]
                w = Element((r, c, d))
        else:
            while bounds.satisfies_bounds(*u):
                if delta(w) / r ** 2 <= delta(v) / v[0] ** 2:
                    walls.add(-w)
                d -= var.gens[2]
                w = Element((r, c, d))
                u = v - w

    # for n in range(1, gcd(v[0], v[1]) + 1):
    #     s = mu(v).denominator() * n
    #     x = mu(v).numerator() * n
    #     y = v[2] - v[1] ** 2 / 2 / v[0] + x ** 2 / 2 / s
    #     y = var.floor(s, x, y)
    #     y_max = bounds.ch_max(s, x)
    #     while y <= y_max:
    #         w = Element((s, x, y))
    #         if delta(w) / s ** 2 <= delta(v) / v[0] ** 2:
    #             walls.add(-w)
    #         y += var.gens[2]

    return walls


def z(v, a, b):
    r"""
    Computes the central charge :math:`Z_{a, b}(v)` for tilt stability.

    TESTS::

        sage: from stability_conditions import *

        sage: v = Element([1, 0, 0])
        sage: tilt.z(v, 1, 0)
        1/2
        sage: tilt.z(v, 1, -1)
        I
        sage: tilt.z(v, 2, 1)
        -I + 3/2

        sage: v = Element([3, -5, 1/2, 1/6, 1/24])
        sage: tilt.z(v, 5, 6)
        -23*I - 47

        sage: var('a, b, r, c, d', domain = RR)
        (a, b, r, c, d)
        sage: p = tilt.z(Element([r, c, d]), a, b)
        sage: bool(p == 1/2*a^2*r - 1/2*b^2*r + b*c - I*b*r + I*c - d)
        True
    """
    twisted_ch = ch(v, b)
    return -twisted_ch[2] + a ** 2 / 2 * twisted_ch[0] + i * twisted_ch[1]
