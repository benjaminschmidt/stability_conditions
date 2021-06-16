r"""
Module for things specifically regarding :math:`\mathbb{P}^3`.

EXAMPLES::

<Lots and lots of examples>

.. TODO::
    - Add examples.
    - Write anything at all.
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
from sage.all import ceil, Expression, infinity, Integer, Rational, RR
# noinspection PyUnresolvedReferences
from sage.all import SageObject, SR, var

import json

from .bounds import Bounds
from .variety import Element, o, p
from .. import slope, tilt


class BoundsP3(Bounds):
    r"""
    Manages explicit Chern character bounds saved in a file.
    """
    def __init__(self, data_json=None, compute=False):
        super().__init__(p(3))
        self.compute = compute
        self.d = var('d', domain=RR)
        if data_json is None:
            # Max ch2 values.
            self.ch2 = {(Integer(1), Integer(0)): Rational(0)}

            # Max ch3 values for special cases
            self.ch3_special = {}

            # Before which ch2 the ch3 special values are required.
            self.ch2_general = {(Integer(1), Integer(0)): Rational(0)}

            # General ch3 bounds below the special cases
            self.ch3 = {(Integer(1), Integer(0)): SR('d^2/2 - d/2')}

            # List of exceptional objects
            self.exceptional = [Element((Integer(1), Integer(0),
                                         Rational(0), Rational(0)))]

            # Up to which rank all bounds have been computed
            self.rank = Integer(1)
        else:
            data = json.loads(data_json)
            ch2_json = data[0]
            ch3_special_json = data[1]
            ch2_general_json = data[2]
            ch3_json = data[3]
            exceptional_json = data[4]
            self.rank = data[5]

            self.ch2 = {}
            for k_json, v in ch2_json.items():
                k = json.loads(k_json)
                k = (Integer(k[0]), Integer(k[1]))
                self.ch2[k] = Rational(v)

            self.ch3_special = {}
            for k_json, v in ch3_special_json.items():
                k = json.loads(k_json)
                k = (Integer(k[0]), Integer(k[1]), Rational(k[2]))
                self.ch3_special[k] = Rational(v)

            self.ch2_general = {}
            for k_json, v in ch2_general_json.items():
                k = json.loads(k_json)
                k = (Integer(k[0]), Integer(k[1]))
                self.ch2_general[k] = Rational(v)

            self.ch3 = {}
            for k_json, v in ch3_json.items():
                k = json.loads(k_json)
                k = (Integer(k[0]), Integer(k[1]))
                self.ch3[k] = SR(v)

            self.exceptional = []
            for v in exceptional_json:
                self.exceptional.append(Element((Integer(v[0]),
                                                 Integer(v[1]),
                                                 Rational(v[2]),
                                                 Rational(v[3]))))

    def json(self):
        """
        Convertes the class into JSON data.
        """
        ch2_json = {}
        for k, v in self.ch2.items():
            k_json = json.dumps((str(k[0]), str(k[1])))
            ch2_json[k_json] = str(v)

        ch3_special_json = {}
        for k, v in self.ch3_special.items():
            k_json = json.dumps((str(k[0]), str(k[1]), str(k[2])))
            ch3_special_json[k_json] = str(v)

        ch2_general_json = {}
        for k, v in self.ch2_general.items():
            k_json = json.dumps((str(k[0]), str(k[1])))
            ch2_general_json[k_json] = str(v)

        ch3_json = {}
        for k, v in self.ch3.items():
            k_json = json.dumps((str(k[0]), str(k[1])))
            ch3_json[k_json] = v.json()

        exceptional_json = []
        for v in self.exceptional:
            exceptional_json.append((str(v[0]), str(v[1]),
                                     str(v[2]), str(v[3])))

        return json.dumps((ch2_json, ch3_special_json, ch2_general_json,
                           ch3_json, exceptional_json, self.rank))

    def ch_max(self, *args):
        if len(args) <= 1 or len(args) >= 4:
            return super().ch_max(*args)
        elif len(args) == 2:
            return self.ch2_max(args[0], args[1])
        elif len(args) == 3:
            return self.ch3_max(args[0], args[1], args[2])
        else:
            raise ValueError()

    def ch_min(self, *args):
        if len(args) <= 1 or len(args) >= 4:
            return super().ch_min(*args)
        elif len(args) == 2:
            return self.ch2_min(args[0], args[1])
        elif len(args) == 3:
            return -infinity
        else:
            raise ValueError()

    def ch2_max(self, r, c):
        if r == 0:
            return infinity
        elif r < 0:
            return infinity

        # Now we know r > 0.
        n = ceil(c / r)
        if (r, c - n * r) in self.ch2.keys():
            d_new = self.ch2[(r, c - n * r)]
            v = Element((r, c - n * r, d_new))
            return (v * o(n, 2))[2]
        else:
            if self.compute:
                raise NotImplementedError('Computing unknown ch_2 bounds ' +
                                          'is not yet implemented.')
            else:
                return self.bogomolov_max(r, c)

    def ch2_min(self, r, c):
        return -self._ch2_max(-r, c)

    def ch3_max(self, r, c, d=var('d', domain=RR)):
        if r == 0:
            if type(d) == Expression:
                if c != 1:
                    raise ValueError("This bound is not of polynomial form.")
                else:
                    return Integer(1)/Integer(24) + d ** 2 / 2
            else:
                f = Integer(-(d + c ** 2 / 2)) % c
                epsilon = f / 2 * (c - f - 1 + f / c)
                return c ** 3 / 24 + d ** 2 / 2 / c - epsilon
        elif r < 0:
            return self._ch3_max(-r, c, -d)

        n = ceil(c / r)
        if n != 0:
            # Tensor v with O(-n).
            v = Element((r, c, d))*o(-n, 2)
            e_new = self.ch3_max(v[0], v[1], v[2])
            # Tensor back.
            w = Element((v[0], v[1], v[2], e_new))*o(n, 3)
            return w[3]

        if (r, c, d) in self.ch3_special.keys():
            # Never happens if d is a symbolic expression.
            return self.ch3_special[(r, c, d)]

        if (r, c) in self.ch2_general.keys():
            if d > self.ch2_max(r, c):
                # Never happens if d is a symbolic expression.
                raise ValueError("There is no such object.")
            else:  # This means we are in the general range.
                return self.ch3[(r, c)].subs(self.d == d)

        # The remaining case means either (r, c, d) is invalid or
        # has not been computed yet.
        if self.compute:
            raise NotImplementedError('Computing unknown ch_3 bounds is ' +
                                      'not yet implemented.')
        else:
            raise IndexError('This bound is unknown.')

    def satisfies_bounds(self, *args):
        if len(args) <= 3:
            return super().satisfies_bounds(*args)
        if len(args) == 4:
            if args[3] >= self.ch_min(args[0], args[1], args[2]):
                if args[3] <= self.ch_max(args[0], args[1], args[2]):
                    return True
                else:
                    return False
            else:
                return False
        else:
            raise ValueError()

    def _find_ch2_max(self, r, c):
        r"""Computes the largest possible ch_2 for given rank r and c_1 = c
        and saves it.

        WARNING:

        Assumes r > 0, c \in [1 - r, 0]. If self.compute is False, bad things
        will happen from invoking this function. Don't do it!"""
        d = self.bogomolov_max(r, c)
        while self.ch3_max(r, c, d) < -self.ch3_max(r, -c, d):
            d -= 1
        self.ch2[(r, c)] = d

    def _find_ch3_max_special(self, r, c, d):
        r"""Computes the largest possible ch_3 for given rank r, c_1 = c,
        ch_2 = d and saves it.

        WARNING: r > 0, c \in [1 - r, 0], \Delta > 0. If self.compute is
        False, bad things will happen from invoking this function.
        Don't do it!"""
        if c == 0:
            if d == 0:
                self.ch3_special[(r, c, d)] = 0
                return

        # We start by

        v = Element((r, c, d))
        e_var = var('e', domain=RR)
        e = infinity
        for exceptional in self.exceptional:
            for i in range(3):
                exc = exceptional*self.var.o(Integer(i))

                if slope.mu(exc) < slope.mu(v):
                    continue

                if (slope.mu(exc) == slope.mu(v) and
                        tilt.delta(exc) / exc[0] ** 2
                        >= tilt.delta(v) / v[0] ** 2):
                    continue

                min_wall = tilt.wall(v, exc*self.var.o(-4))

                if min_wall.is_empty:
                    continue

                if min_wall.s > slope.mu(v):
                    continue

                if min_wall.s < slope.mu(exc) - 4:
                    continue

                v_complete = Element((r, c, d, e_var))
                e_max_euler = exc.chi(v_complete).solve(e)[0].rhs()
                e_max_euler = self.var.floor(r, c, d, e_max_euler)

                if e <= e_max_euler:
                    # In this case, we will not be able to get a better bound,
                    # no matter what the bounds from the walls say.
                    continue

                walls = tilt.walls_left(v, self.var, min_wall.s, bounds=self)

                e_max_walls = -infinity
                for wall in walls:
                    for w in walls[wall]:
                        u = v - w
                        e_max_walls = max(e_max_walls,
                                          (self.ch3_max(w[0], w[1], w[2])
                                           + self.ch3_max(u[0], u[1], u[2])))

                if e_max_walls < e:
                    e = max(e_max_euler, e_max_walls)
        self.ch3Special[(r, c, d)] = e

    #     if d == Polynomial((0, 1), 'd'):
    #         # Subobject at wall: (s, x, y, z)
    #         subSlope = previousFarey(Fraction(c, r), r)
    #         s = subSlope.denominator
    #         x = subSlope.numerator
    #
    #         # Problem: multiples of s, x might give better bounds.
    #         # Thus, we maximize k*s <= r, i.e., k <= r/s.
    #         k = r//s
    #         s *= k
    #         x *= k
    #
    #         # Subobject at wall: (s, x, y, z)
    #         y = self.ch2Max(s, x)
    #         z = self.ch3Max(s, x, y)
    #
    #         # e - z <= ch3Max(Quotient)
    #         self.ch3[(r, c)] = z + self.ch3Max(r - s, c - x, d)(d - y)
    #         return
    #
    #     elif c == 0:
    #         if d == 0:
    #             self.ch3Special[(r, c, d)] = 0
    #             return
    #
    #         n = 1
    #         while 2*r > n^2 + n:
    #             n += 1
    #
    #         if d > -n:
    #             self.ch3Special[(r, c, d)] = -2*d - r
    #             return
    #
    #         elif d <= -n:
    #             self.ch3Special[(r, c, d)] = (Fraction(d**2, 2) +
    #                                           (n - Fraction(3, 2))*d +
    #                                           Fraction(n*(n + 1), 2) - r)
    #             return
    #
    #     # elif (c, d) == (0, 0):
    #     #     self.ch3Special[(r, c, d)] = 0
    #     #     return
    #
    #     v = Element((r, c, d))
    #     walls = WallStructure(v)
    #     eVar = Polynomial((0, 1), 'e')
    #     e = inf
    #     for i in range(3):
    #         # exc = O(i)
    #         minWall = v.wall(O(i-4))
    #         if minWall.s < i - 4: # or minWall.s > i:
    #             continue
    #         walls.compute(minWall.s, self)
    #
    #         # We are under the hypothesis that the largest wall always
    #         # gives the weakest bound!
    #         if len(walls.walls) != 0:
    #             w, wall = walls.walls[0]
    #             u = v - w
    #             ePotential = (self.ch3Max(w[0], w[1], w[2])
    #                  + self.ch3Max(u[0], u[1], u[2]))
    #         else:
    #             vComplete = Element((r, c, d, eVar))
    #             ePotential = O(i).chi(vComplete).roots()[0]
    #             ePotential = ch3Floor(c, d, ePotential)
    #         if ePotential < e:
    #             e = ePotential
    #
    #     # If no bound has been found so far, assume general bound is correct.
    #     if e == inf:
    #         # Subobject at wall: (s, x, y, z)
    #         subSlope = previousFarey(Fraction(c, r), r)
    #         s = subSlope.denominator
    #         x = subSlope.numerator
    #
    #         # Problem: multiples of s, x might give better bounds.
    #         # Thus, we maximize k*s <= r, i.e., k <= r/s.
    #         k = r//s
    #         s *= k
    #         x *= k
    #
    #         # Subobject at wall: (s, x, y, z)
    #         y = self.ch2Max(s, x)
    #         z = self.ch3Max(s, x, y)
    #
    #         # e - z <= ch3Max(Quotient)
    #         e = z + self.ch3Max(r - s, c - x, d - y)
    #
    #     self.ch3Special[(r, c, d)] = e

    # def _findch3Max(self, r, c, d = Polynomial((0, 1), 'd')):
    #     """Computes the largest possible ch_3 for given rank r, c_1 = c,
    #     ch_2 = d and saves it."""
    #     if d == Polynomial((0, 1), 'd'):
    #         # Subobject at wall: (s, x, y, z)
    #         subSlope = previousFarey(Fraction(c, r), r)
    #         s = subSlope.denominator
    #         x = subSlope.numerator
    #
    #         # Problem: multiples of s, x might give better bounds.
    #         # Thus, we maximize k*s <= r, i.e., k <= r/s.
    #         k = r//s
    #         s *= k
    #         x *= k
    #
    #         # Subobject at wall: (s, x, y, z)
    #         y = self.ch2Max(s, x)
    #         z = self.ch3Max(s, x, y)
    #
    #         # e - z <= ch3Max(Quotient)
    #         self.ch3[(r, c)] = z + self.ch3Max(r - s, c - x, d)(d - y)
    #         return
    #
    #     # elif (c, d) == (0, 0):
    #     #     self.ch3Special[(r, c, d)] = 0
    #     #     return
    #

    def _find_exceptional(self, r):
        """Finds all exceptional bundles of rank r. Assumes that all bounds
        have been determined for objects of rank <= r."""
        if r % 2 == 0:
            return

        for c in range(1 - r, 1):
            d = (2 * c ** 2 - r**2 + 1)/r
            v = Element((r, c, d))
            if self.var.valid(v):
                if d <= self.ch2_max(r, c):
                    v = Element((r, c, d, self.ch3_max(r, c, d)))
                    self.exceptional.append(v)
