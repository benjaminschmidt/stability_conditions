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


class BoundsP3(Bounds):
    r"""
    Manages explicit Chern character bounds saved in a file.
    """
    def __init__(self, data_json=None):
        super().__init__(p(3))
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

        if type(d) == Expression:
            if (r, c) in self.ch2_general.keys():
                return self.ch3[(r, c)]
            else:
                raise IndexError('This bound is unknown.')
        else:
            if (r, c) in self.ch2_general.keys():
                if d <= self.ch2_general[(r, c)]:
                    return self.ch3[(r, c)]
            if d > self.ch2_max(r, c):
                raise ValueError("There is no such object.")
            elif (r, c, d) in self.ch3_special.keys():
                return self.ch3_special[(r, c, d)]
            else:
                # The remaining case means either (r, c, d) is invalid or
                # has not been computed yet.
                raise IndexError('This bound is unknown.')
