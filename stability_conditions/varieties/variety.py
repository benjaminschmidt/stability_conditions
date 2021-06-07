r"""
Module for the numerical Chow ring of some special varieties.

This module is only dealing with smooth projective varieties :math:`X` whose
numerical Chow ring :math:`N^*(X)` has rank one in each dimension. For all
purposes of this module one should only work with varieties for which the
Bogomolov inequality for semistable sheaves holds. This is always true in
characteristic :math:`0`. Let :math:`H` be the class of the ample generator of
the Neron-Severi group :math:`N^1(X)`.

EXAMPLES::

<Lots and lots of examples>

.. TODO::
    - Add examples.
    - Document attributes of Variety.
    - Implement _latex for both classes.
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
from sage.functions.other import ceil, factorial, floor
# noinspection PyUnresolvedReferences
from sage.rings.all import Integer
# noinspection PyUnresolvedReferences
from sage.structure.all import SageObject


class Element(SageObject):
    r"""
    An element in the numerical Chow ring of a variety.

    INPUT:

    - ``vec`` -- list of rational numbers. vec[i] is the coefficient in
      front of :math:`H^i`.
    """
    def __init__(self, vec):
        self.vec = tuple(vec)

    def __add__(self, other):
        r"""
        Returns this object plus other.

        TESTS::

            sage: from stability_conditions import *

            sage: Element([1, 0, -4, 9]) + Element([1, -1/2, 1/2, -1])
            (2, -1/2, -7/2, 8)

            sage: Element([]) + Element([])
            ()

            sage: Element([1, 2]) + Element([1, 2, 3])
            Traceback (most recent call last):
            ...
            ValueError: Two Elements need to have the same length for adding.
        """
        if len(self) != len(other):
            raise ValueError("Two Elements need to have the same length for "
                             "adding.")

        return Element([self[j] + other[j] for j in range(len(self))])

    def __eq__(self, other):
        r"""
        Two elements are equal if they have identical coefficients.

        TESTS::

            sage: from stability_conditions import *

            sage: Element([1, 2, 3]) == Element([1, 2, 3])
            True

            sage: Element([1, 2, 3]) == Element([1, 2, 4])
            False

            sage: Element([5, 2, 3]) == Element([5, 2, 3, 4])
            False

            sage: Element([]) == Element([])
            True

            sage: Element([2, 2, 3]) == 1
            False
        """
        try:
            return self.vec == other.vec
        except AttributeError:
            return False

    def __getitem__(self, j):
        r"""
        Returns the coefficient in front of :math:`H^j`.

        TESTS::

            sage: from stability_conditions import *

            sage: v = Element([3/2, 2, 1])
            sage: v[0]
            3/2
            sage: v[2]
            1
        """
        return self.vec[j]

    def __hash__(self):
        r"""
        Returns a hash value for this object.

        TESTS::

            sage: from stability_conditions import *

            sage: v = Element((1, 2, 3))
            sage: w = Element((1, 2, 3))
            sage: v.__hash__() == w.__hash__()
            True

            sage: v = Element((1, 2, 3))
            sage: w = Element((1, 2, 4))
            sage: v.__hash__() == w.__hash__()
            False

            sage: v = Element((5, 2, 3))
            sage: w = Element((5, 2, 3, 4))
            sage: v.__hash__() == w.__hash__()
            False

            sage: v = Element(())
            sage: w = Element(())
            sage: v.__hash__() == w.__hash__()
            True

            sage: Element([2, 2, 3]).__hash__() == hash(1)
            False
        """
        return hash(self.vec)

    def __len__(self):
        r"""
        Returns the length of this object.

        TESTS::

            sage: from stability_conditions import *

            sage: len(Element([8, 3/2, 4, 5]))
            4

            sage: len(Element([]))
            0
        """
        return len(self.vec)

    def __ne__(self, other):
        r"""
        Returns the opposite of self == other.

        TESTS::

            sage: from stability_conditions import *

            sage: Element([1, 2, 3]) != Element([1, 2, 3])
            False

            sage: Element([1, 2, 3]) != Element([1, 2, 4])
            True

            sage: Element([5, 2, 3]) != Element([5, 2, 3, 4])
            True

            sage: Element([]) != Element([])
            False

            sage: Element([2, 2, 3]) != 1
            True
        """
        return not self == other

    def __mul__(self, other):
        r"""
        Returns this object multiplied with other.

        TESTS::

            sage: from stability_conditions import *
            sage: Element([1, 0, -1]) * Element([1, -1, 1/2])
            (1, -1, -1/2)

            sage: Element([]) * Element([])
            ()

            sage: Element([0, 1, -3/2, 7/6]) * 2
            (0, 2, -3, 7/3)

            sage: var('x')
            x
            sage: Element((1, 0, 0))*x
            (x, 0, 0)
        """
        if not isinstance(other, Element):
            return Element([other * v for v in self])

        vec = [0]*len(self)

        for k in range(len(self)):
            for j in range(k + 1):
                vec[k] += self[j]*other[k - j]

        return Element(vec)

    def __neg__(self):
        r"""
        Returns the additive inverse of this object.

        TESTS::

            sage: from stability_conditions import *

            sage: -Element([2, -1/2, 3, 5/6])
            (-2, 1/2, -3, -5/6)

            sage: -Element([])
            ()
        """
        return Element([-v for v in self])

    def __pos__(self):
        r"""
        Returns a new copy of this object.

        TESTS::

            sage: from stability_conditions import *

            sage: +Element([1, 1/2, 1/3, 1/4, 1/5])
            (1, 1/2, 1/3, 1/4, 1/5)

            sage: +Element([])
            ()
        """
        return Element(self.vec)

    def __pow__(self, k):
        r"""
        Returns this object to the power k.

        INPUT:

        - ``k`` - positive integer.

        TESTS::

            sage: from stability_conditions import *
            sage: Element([1, -1, 1/2])**6
            (1, -6, 18)

            sage: Element([2, 2, 1, 1, 1])**7
            (128, 896, 3136, 7616, 15008)

            sage: Element([])**100
            ()

            sage: Element([1, -1, 1/2, -1/6])**(-1)
            (1, 1, 1/2, 1/6)

            sage: Element([2, -1, -1/2])**(-1)
            (1/2, 1/4, 1/4)

            sage: Element([4, -1, 3, 6, 7])**(-4)
            (1/256, 1/256, -19/2048, -151/4096, -2957/65536)

            sage: Element([0, 1])**(-1)
            Traceback (most recent call last):
            ...
            ValueError: Element is not invertible.
        """
        if k == 0:
            return Element([1] + [0]*(len(self) - 1))
        elif k == 1:
            return Element(self.vec)
        elif k > 1:
            return self*self**(k-1)
        elif k == -1:
            if self[0] == 0:
                raise ValueError("Element is not invertible.")
            vec = [1/self[0]]
            for m in range(1, len(self)):
                vec.append(-vec[0]*sum([self[i]*vec[m - i]
                                        for i in range(1, m + 1)]))
            return Element(vec)
        elif k < -1:
            return (self**(-1))**(-k)

    def __rmul__(self, other):
        r"""
        Returns other multiplied with this object.

        TESTS::

            sage: from stability_conditions import *
            sage: 2 * Element([0, 1, -3/2, 7/6])
            (0, 2, -3, 7/3)

            sage: 2 * Element([])
            ()

            sage: var('x')
            x
            sage: x*Element((1, 0, 0))
            (x, 0, 0)
        """
        return self.__mul__(other)

    def __sub__(self, other):
        r"""
        Returns this object minus other.

        TESTS::

            sage: from stability_conditions import *

            sage: Element([1, 0, -4, 9]) - Element([1, -1/2, 1/2, -1])
            (0, 1/2, -9/2, 10)

            sage: Element([]) - Element([])
            ()

            sage: Element([1, 2, 3, 4, 5]) - Element([0, 1, 2, 3])
            Traceback (most recent call last):
            ...
            ValueError: Two Elements need to have the same length for adding.
        """
        return self + (-other)

    def __truediv__(self, k):
        r"""
        Divides this object by a rational number `k`.

        TESTS::

            sage: from stability_conditions import *
            sage: Element([10, 2, -1, -5])/5
            (2, 2/5, -1/5, -1)

            sage: Element([])/100
            ()
        """
        return self*(1/k)

    def _repr_(self):
        r"""
        String representation of this object.

        TESTS::

            sage: from stability_conditions import *

            sage: Element([1, 2, 3])
            (1, 2, 3)
        """
        return self.vec.__str__()

    def c(self, k=None):
        r"""
        Computes Chern classes from Chern characters.

        More precisely, if this object is :math:`\operatorname{ch}(E)` for
        some object :math:`E`, then this method computes :math:`c(E)`.

        INPUT:

        - ``k`` -- integer between 0 and len(self) or None (default: `None`).

        OUTPUT:

        A rational number describing the k-th Chern class or the total
        Chern class, if k is None. Note that by convention :math:`c_0 = 1`.

        ALGORITHM:

        We are using the Newton Identities that say

        .. math:: k c_k = \sum_{l = 1}^k (-1)^{l - 1} l! c_{k - l}
           \operatorname{ch}_l.

        (https://en.wikipedia.org/wiki/Newton%27s_identities)

        TESTS::

            sage: from stability_conditions import *

            sage: v = Element([1, -1, 1/2, -1/6])
            sage: v.c(1)
            -1
            sage: v.c(0)
            1
            sage: v.c(2)
            0
            sage: v.c(3)
            0
            sage: v.c()
            (1, -1, 0, 0)

            sage: v = Element([2, -1, -1/2])
            sage: v.c(2)
            1
            sage: v.c()
            (1, -1, 1)
        """
        if k is None:
            k_max = len(self) - 1
        else:
            k_max = k

        c = [1]

        for j in range(1, k_max + 1):
            c.append(sum([(-1)**(m - 1)*c[j - m]*factorial(m)*self[m]
                          for m in range(1, j + 1)])/j)

        if k is None:
            return Element(c)
        else:
            return c[-1]

    def ch(self, k=None, rk=1):
        r"""
        Computes Chern characters from Chern classes.

        More precisely, if this object is :math:`c(E)` for some object
        :math:`E`, then this method computes :math:`\operatorname{ch}(E)`.
        Note that the Chern character encodes the rank, while the total Chern
        class does not. Therefore, you can set the rank `rk`. On default we
        assume rank one.

        INPUT:

        - ``k`` -- integer between 0 and len(self) or None (default: `None`).
        - ``rk`` -- rank of the object (default: `1`).

        OUTPUT:

        A rational number describing the k-th Chern character or the total
        Chern character, if k == None.

        ALGORITHM:

        We are using the Newton Identities that say

        .. math:: k! \operatorname{ch}_k = (-1)^{k - 1} k c_k -
           \sum_{l = 1}^{k - 1} (-1)^{k + l} l! c_{k - l} \operatorname{ch}_l.

        (https://en.wikipedia.org/wiki/Newton%27s_identities)

        TESTS::

            sage: from stability_conditions import *

            sage: v = Element([1, -1, 0, 0])
            sage: v.ch(0)
            1
            sage: v.ch(0, rk=2)
            2
            sage: v.ch(1)
            -1
            sage: v.ch(2)
            1/2
            sage: v.ch(3)
            -1/6
            sage: v.ch()
            (1, -1, 1/2, -1/6)

            sage: v = Element([1, -1, 1])
            sage: v.ch(0)
            1
            sage: v.ch(0, rk=2)
            2
            sage: v.ch(1)
            -1
            sage: v.ch(2)
            -1/2
            sage: v.ch(rk=2)
            (2, -1, -1/2)
        """
        if k is None:
            k_max = len(self) - 1
        else:
            k_max = k

        chern = [rk]

        for j in range(1, k_max + 1):
            s = sum([(-1)**(j + m)*factorial(m)*self[j - m]*chern[m]
                     for m in range(1, j)])/Integer(factorial(j))
            chern.append((-1) ** (j - 1) * j * self[j] / factorial(j) - s)

        if k is None:
            return Element(chern)
        else:
            return chern[-1]

    def dual(self):
        r"""
        If this object is :math:`\operatorname{ch}(E)`, then this returns
        :math:`\operatorname{ch}(E^{\vee})`.

        TESTS::

            sage: from stability_conditions import *
            sage: Element([5, -54, 6/10, 345, 34/5, 67/2]).dual()
            (5, 54, 3/5, -345, 34/5, -67/2)

            sage: Element([]).dual()
            ()
        """
        return Element([(-1)**j*self[j] for j in range(len(self))])

    # def _latex(self):
    #     pass


class Variety(SageObject):
    r"""
    A smooth projective variety.

    INPUT:

    - ``dim_var`` -- positive integer describing the dimension of the variety.
    - ``generators`` -- list of rational numbers. :math:`gens[i]H^i`
      generates :math:`N^i(X)` for `i` between `0` and dim. By assumption
      :math:`gens[dim] H^{dim}` is the class of a point. Therefore,
      1/gens[dim] has to be an integer equal to :math:`H^{dim}`.
    - ``td`` -- list of rational numbers describing the Todd class of
      :math:`X` given as :math:`td[i]H^i = \operatorname{td}_i(X)`.

    TESTS::

        sage: from stability_conditions import *

        sage: Variety(-1, [], [])
        Traceback (most recent call last):
        ...
        ValueError: Varieties must have positive dimension.

        sage: Variety(2, [1, 1], [1, 3/2, 1])
        Traceback (most recent call last):
        ...
        ValueError: Incorrect number of generators for the numerical Chow
        group.

        sage: Variety(3, [2, 1, 1, 1], [1, 2, 11/6, 1])
        Traceback (most recent call last):
        ...
        ValueError: By definition `1` has to generate `N^0(X)`.

        sage: Variety(3, [1, 2, 1, 1], [1, 2, 11/6, 1])
        Traceback (most recent call last):
        ...
        ValueError: By definition `H` has to generate `N^1(X)`.

        sage: Variety(3, [1, 1, 1, 2/3], [1, 2, 11/6, 1])
        Traceback (most recent call last):
        ...
        ValueError: By definition 1/gens[dim] = H^dim is an integer.

        sage: Variety(2, [1, 1, 1], [2, 2, 11/6])
        Traceback (most recent call last):
        ...
        ValueError: `\td_0(X)` is always `1`.

        sage: Variety(2, [1, 1, 1], [1, 3/2])
        Traceback (most recent call last):
        ...
        ValueError: Todd class does not have the correct length.
    """

    def __init__(self, dim_var, generators, td):
        if dim_var <= 0:
            raise ValueError("Varieties must have positive dimension.")
        self.dim = dim_var

        if len(generators) != dim_var + 1:
            raise ValueError("Incorrect number of generators for the "
                             "numerical Chow group.")
        if generators[0] != 1:
            raise ValueError(r"By definition `1` has to generate `N^0(X)`.")
        if generators[1] != 1:
            raise ValueError(r"By definition `H` has to generate `N^1(X)`.")
        if generators[dim_var].numerator() != 1:
            raise ValueError(r"By definition 1/gens[dim] = H^dim "
                             "is an integer.")
        self.gens = tuple(generators)

        if td[0] != 1:
            raise ValueError(r"`\td_0(X)` is always `1`.")
        if len(td) != dim_var + 1:
            raise ValueError("Todd class does not have the correct length.")
        self.td = Element(td)

    def _repr_(self):
        r"""
        Returns a string that shows the constructor for this object.

        TESTS::

            sage: from stability_conditions import *

            sage: Variety(2, [1, 1, 1], [1, 3/2, 1])
            Variety(2, (1, 1, 1), (1, 3/2, 1))
        """
        return "Variety(%d, %s, %s)" % (self.dim, self.gens, self.td)

    def o(self, m):
        r"""
        Returns the Chern character of the line bundle :math:`\mathcal{O}(mH)`.

        TESTS::

            sage: from stability_conditions import *

            sage: X = Variety(2, [1, 1, 1], [1, 3/2, 1])
            sage: X.o(-5)
            (1, -5, 25/2)

            sage: X = Variety(3, [1, 1, 1, 1], [1, 2, 11/6, 1])
            sage: X.o(3)
            (1, 3, 9/2, 9/2)
        """
        return o(m, self.dim)

    def chi(self, v, w=None):
        r"""
        Returns the Euler characteristic :math:`\chi(v, w)` or :math:`\chi(v)`.

        INPUT:

        - ``v`` -- Element of the numerical Chow ring of this variety.
        - ``w`` -- Element of the numerical Chow ring of this variety or None.

        OUTPUT:

        If w == None, this returns the Euler characteristic
        :math:`\chi(v, w)`. Otherwise, it returns :math:`\chi(v)`.

        TESTS::

            sage: from stability_conditions import *

            sage: X = Variety(3, [1, 1, 1, 1], [1, 2, 11/6, 1])
            sage: v = Element([1, -1, 1/2, -1/6])
            sage: w = Element([1, -5, 25/2, -125/6])
            sage: X.chi(v, w)
            -1

            sage: X = Variety(3, [1, 1, 1/3, 1/3], [1, 1, 2/3, 1/3])
            sage: v = Element([1, 1, 1/2, 1/6])
            sage: X.chi(v)
            5

            sage: X = Variety(3, [1, 1, 1/3, 1/3], [1, 1, 2/3, 1/3])
            sage: v = Element([3, -1, -1/2, 1/6])
            sage: X.chi(v, X.o(0))
            3

            sage: X = Variety(2, [1, 1, 1], [1, 3/2, 1])
            sage: v = Element([1, -1, 1/2])
            sage: w = Element([2, -1, -1/2])
            sage: X.chi(v, w)
            3
        """
        # 1/self.gens[self.dim] = H^{self.dim}
        # Rest of the formula is Hirzebruch-Riemann-Roch
        if w is None:
            return (v*self.td)[self.dim]/self.gens[self.dim]
        else:
            return (v.dual()*w*self.td)[self.dim]/self.gens[self.dim]

    def valid(self, v):
        r"""Checks whether the Element v could occur as a Chern character.

        This method returns True if the rank, the Euler characteristic
        :math:`\chi(v)`, and the Chern classes corresponding to
        :math:`v = \operatorname{ch}(E)` are all integral.

        TESTS::

            sage: from stability_conditions import *

            sage: X = Variety(2, [1, 1, 1], [1, 3/2, 1])
            sage: v = Element([1, 1, 1/2])
            sage: X.valid(v)
            True
            sage: v = Element([1, 1, 1])
            sage: X.valid(v)
            False
            sage: v = Element([1, 0, 1/6])
            sage: X.valid(v)
            False

            sage: X = Variety(3, [1, 1, 1, 1], [1, 2, 11/6, 1])
            sage: v = Element([1, 1, 1/2, 1/6])
            sage: X.valid(v)
            True
            sage: v = Element([2, -1, -1/2, 1/3])
            sage: X.valid(v)
            False
            sage: v = Element([2, -1, -1/2, -1/6])
            sage: X.valid(v)
            True
            sage: v = Element([3/2, -1, -1/2, 1/3])
            sage: X.valid(v)
            False

            sage: X = Variety(3, [1, 1, 1/2, 1/6], [1, 0, 0, 0])
            sage: v = Element([1, -1, 1/2, -1/3])
            sage: X.valid(v)
            True
            sage: v = Element([1, -1, 1/3, 0])
            sage: X.valid(v)
            False
        """
        if v[0].denominator() != 1:
            return False

        if self.chi(v).denominator() != 1:
            return False

        c = v.c()
        for j in range(len(v)):
            if (c[j]/self.gens[j]).denominator() != 1:
                return False

        return True

    def floor(self, *args):
        r"""Rounds down Chern characters.

        Let n = len(args) - 1 and assume that for some object :math:`E` we
        have :math:`\operatorname{ch}_i(E) = args[i]` for all `i` in between
        `0` and `n - 1`. This method then returns the largest rational number
        :math:`x < args[n]` such that :math:`\operatorname{ch}_n(E) = x` is
        a possible value.

        Basically, this function simply rounds down based on the
        fact that Chern classes are integral.

        INPUT:

        - ``args`` -- list of rational numbers with length smaller than or
          equal to the dimension of the variety plus one.

        WARNING:

        It is not checked whether the values in args are valid inputs.

        TESTS::

            sage: from stability_conditions import *

            sage: X = Variety(2, [1, 1, 1], [1, 3/2, 1])
            sage: X.floor(1, 2, 3, 4)
            Traceback (most recent call last):
            ...
            IndexError: Too many arguments given the dimension of the variety.

            sage: X = Variety(2, [1, 1, 1], [1, 3/2, 1])
            sage: X.floor()
            Traceback (most recent call last):
            ...
            IndexError: There has to be at least one argument.

            sage: X = Variety(2, [1, 1, 1], [1, 3/2, 1])
            sage: X.floor(2, -1, 0)
            -1/2

            sage: X = Variety(3, [1, 1, 1/2, 1/6], [1, 0, 0, 0])
            sage: X.floor(1, 0, 0, 1/6)
            1/6
            sage: X.floor(1, 0, 0, 1/7)
            1/12
            sage: X.floor(1, 0, 0, 1/13)
            0
        """
        if len(args) > self.dim + 1:
            raise IndexError("Too many arguments given the dimension of the "
                             "variety.")

        if len(args) == 0:
            raise IndexError("There has to be at least one argument.")

        v = Element(args)
        j = len(args) - 1  # i-th Chern character has to be rounded

        if j == 0:
            return floor(args[0])
        else:
            c = list(v.c().vec)
            if j % 2 == 0:
                c[j] = ceil(c[j]/self.gens[j])*self.gens[j]
            else:  # if j % 2 == 1
                c[j] = floor(c[j]/self.gens[j])*self.gens[j]
            return Element(c).ch(j)

    # def _latex(self):
    #     pass


def ch(v, b):
    r"""
    Computes the twisted Chern character :math:`\operatorname{ch}^{b}(v)`.

    The twisted Chern character is defined as
    :math:`\operatorname{ch}(v) e^{-bH}`.

    TESTS::

        sage: from stability_conditions import *

        sage: variety.ch(Element([1, 0, 0]), 1)
        (1, -1, 1/2)

        sage: variety.ch(Element([3, -1, -1/2, -1/6]), 2)
        (3, -7, 15/2, -31/6)

        sage: variety.ch(Element([]), 101)
        ()

        sage: var('b', domain = RR)
        b
        sage: variety.ch(Element([1, 0, 0, 0, 0, 0]), b)
        (1, -b, 1/2*b^2, -1/6*b^3, 1/24*b^4, -1/120*b^5)
    """
    return v * o(-b, len(v))


def o(m, n):
    r"""
    Returns the Chern character of the line bundle :math:`\mathcal{O}(mH)` on
    a variety of dimension `n`.

    INPUT:

    - ``m`` -- integer
    - ``n`` -- positive integer describing the dimension of the variety

    TESTS::

        sage: from stability_conditions.varieties import o

        sage: o(-5, 2)
        (1, -5, 25/2)

        sage: o(3, 3)
        (1, 3, 9/2, 9/2)
    """
    return Element([m**j/factorial(j) for j in range(n + 1)])


def p(n):
    r"""
    Returns projective space of dimension `n`.

    ALGORITHM:

    We use that the Todd class of :math:`\mathbb{P}^n` is given by
    :math:`\left(\frac{H}{1 - e^{-H}}\right)^{n + 1}`, where :math:`H` is
    the class of a hyperplane.

    TESTS::

            sage: from stability_conditions import *

            sage: p(1)
            Variety(1, (1, 1), (1, 1))

            sage: p(2)
            Variety(2, (1, 1, 1), (1, 3/2, 1))

            sage: p(3)
            Variety(3, (1, 1, 1, 1), (1, 2, 11/6, 1))

            sage: p(4)
            Variety(4, (1, 1, 1, 1, 1), (1, 5/2, 35/12, 25/12, 1))

            sage: p(27)
            Variety(27, (1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), (1, 14, 581/6, 441,
            1070657/720, 158319/40, 37400899/4320, 767921/48,
            2201402821/86400, 56805833/1600, 749183964913/17107200,
            339766721/7040, 3577239914712089/74724249600,
            127498710046229/2965248000, 2019789295599781/57480192000,
            10370744131087/395366400, 62157476672898653/3464487936000,
            34433122825921/3055104000,
            742857089765534754163/114042281633280000,
            1882907719410379703/543058483968000,
            1015770864598729418381/597364332364800000,
            1213894246448392783/1580328921600000,
            2558719956867520777483/8014638125894400000,
            86109807874052383/706758212160000,
            103452650387006862847/2431106898187968000,
            816088653136373/60295309974900, 312536252003/80313433200, 1))

    """
    structure_sheaf = Element([Integer(1)] + [Integer(0)]*(n + 1))
    o_minus_one = Element([Integer(-1)**j/factorial(j) for j in range(n + 2)])
    v = Element((structure_sheaf - o_minus_one)[1:])
    td = v ** (-n - 1)
    return Variety(n, [Integer(1)]*(n + 1), td.vec)
