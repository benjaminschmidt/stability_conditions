r"""
Module for dealing with Chern character bounds.

EXAMPLES::

<Lots and lots of examples>

.. TODO::
    - Add examples.
    - Write BoundsFromFile class
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
from sage.all import infinity, SageObject


class Bounds(SageObject):
    r"""
    Class for managing Chern character bounds.

    This basic class is just using the Bogomolov inequality. For more
    advanced behavior one can derive a class from this one and overwrite
    required funtions. Examples of this can be found in the varieties
    subpackage or in the class BoundsFromFile.

    INPUT:

        - ``var`` -- Variety of dimension at least two.

    TESTS::

            sage: from stability_conditions import *
    """

    def __init__(self, var):
        self.var = var

    def bogomolov_max(self, r, c):
        r"""
        Returns the maximal ch_2 satisfying Bogomolov's inequality or infinity.

        INPUT:

        - ``r`` -- rank
        - ``c`` -- first Chern class

        WARNING:

        This does not check whether the inputs are valid, e.g., fractional
        ranks are probably a bad idea to input.

        TESTS::

            sage: from stability_conditions import *

            sage: X = p(3)
            sage: bounds = varieties.bounds.Bounds(X)
            sage: bounds.bogomolov_max(0, 1)
            +Infinity
            sage: bounds.bogomolov_max(-10, 100)
            +Infinity
            sage: bounds.bogomolov_max(1, 0)
            0
            sage: bounds.bogomolov_max(3, 7)
            15/2

            sage: X = Variety(3, [1, 1, 1/2, 1/6], [1, 0, 0, 0])
            sage: bounds = varieties.bounds.Bounds(X)
            sage: bounds.bogomolov_max(2, -1)
            0
        """
        if r <= 0:
            return infinity
        elif r > 0:
            return self.var.floor(r, c, c**2/(2*r))
        else:
            raise ValueError("r = " + str(r) + " is not a valid input.")

    def bogomolov_min(self, r, c):
        r"""
        Returns the minimal ch_2 satisfying Bogomolov's inequality or
        -infinity.

        INPUT:

        - ``r`` -- rank
        - ``c`` -- first Chern class

        WARNING:

        This does not check whether the inputs are valid, e.g., fractional
        ranks are probably a bad idea to input.

        TESTS::

            sage: from stability_conditions import *

            sage: X = p(3)
            sage: bounds = varieties.bounds.Bounds(X)
            sage: bounds.bogomolov_min(0, 1)
            -Infinity
            sage: bounds.bogomolov_min(-3, 7)
            -15/2
            sage: bounds.bogomolov_min(-1, 0)
            0
            sage: bounds.bogomolov_min(3, 7)
            -Infinity

            sage: X = Variety(3, [1, 1, 1/2, 1/6], [1, 0, 0, 0])
            sage: bounds = varieties.bounds.Bounds(X)
            sage: bounds.bogomolov_min(-2, 1)
            0
        """
        if r >= 0:
            return -infinity
        elif r < 0:
            return self.var.ceil(r, c, c**2/(2*r))
        else:
            raise ValueError("r = " + str(r) + " is not a valid input.")

    def ch_max(self, *args):
        r"""
        Returns the maximal i-th Chern character, where i = len(args).

        INPUT:

        - ``args`` -- list of the first Chern characters up to i - 1.

        WARNING:

        This basic implementation only works up to the second Chern character.
        It simply uses the Bogomolov inequality. For more advanced
        functionality check classes derived from this one. It also does not
        check whether the inputs are valid, e.g., fractional ranks are
        probably a bad idea to input.

        TESTS::

            sage: from stability_conditions import *

            sage: X = p(3)
            sage: bounds = varieties.bounds.Bounds(X)
            sage: bounds.ch_max(0, 1)
            +Infinity
            sage: bounds.ch_max(-10, 100)
            +Infinity
            sage: bounds.ch_max(1, 0)
            0
            sage: bounds.ch_max(3, 7)
            15/2
            sage: bounds.ch_max()
            +Infinity
            sage: bounds.ch_max(1)
            +Infinity

            sage: X = Variety(3, [1, 1, 1/2, 1/6], [1, 0, 0, 0])
            sage: bounds = varieties.bounds.Bounds(X)
            sage: bounds.ch_max(2, -1)
            0
            sage: bounds.ch_max(1, 2, 3)
            Traceback (most recent call last):
            ...
            NotImplementedError: Third Chern character or higher.
        """
        if len(args) <= 1:
            return infinity
        elif len(args) == 2:
            return self.bogomolov_max(args[0], args[1])
        else:
            raise NotImplementedError('Third Chern character or higher.')

    def ch_min(self, *args):
        r"""
        Returns the minimal i-th Chern character, where i = len(args).

        INPUT:

        - ``args`` -- list of the first Chern characters up to i - 1.

        WARNING:

        This basic implementation only works up to the second Chern character.
        It simply uses the Bogomolov inequality. For more advanced
        functionality check classes derived from this one. It also does not
        check whether the inputs are valid, e.g., fractional ranks are
        probably a bad idea to input.

        TESTS::

            sage: from stability_conditions import *

            sage: X = p(3)
            sage: bounds = varieties.bounds.Bounds(X)
            sage: bounds.ch_min(0, 1)
            -Infinity
            sage: bounds.ch_min(-3, 7)
            -15/2
            sage: bounds.ch_min(-1, 0)
            0
            sage: bounds.ch_min(3, 7)
            -Infinity
            sage: bounds.ch_min()
            -Infinity
            sage: bounds.ch_min(1)
            -Infinity

            sage: X = Variety(3, [1, 1, 1/2, 1/6], [1, 0, 0, 0])
            sage: bounds = varieties.bounds.Bounds(X)
            sage: bounds.ch_min(-2, 1)
            0
            sage: bounds.ch_min(1, 2, 3)
            Traceback (most recent call last):
            ...
            NotImplementedError: Third Chern character or higher.
        """
        if len(args) <= 1:
            return -infinity
        elif len(args) == 2:
            return self.bogomolov_min(args[0], args[1])
        else:
            raise NotImplementedError('Third Chern character or higher.')

    def satisfies_bounds(self, *args):
        r"""
        For a list of Chern characters checks whether the last one satisfies
        the bounds of this class for the previous values.

        INPUT:
        - ``args`` -- list of Chern characters

        WARNING:

        This basic implementation only works up to the second Chern character.
        It simply uses the Bogomolov inequality. For more advanced
        functionality check classes derived from this one. It also does not
        check whether the inputs are valid, e.g., fractional ranks are
        probably a bad idea to input.

        TESTS::

            sage: from stability_conditions import *

            sage: X = p(2)
            sage: bounds = varieties.bounds.Bounds(X)
            sage: bounds.satisfies_bounds(1)
            True
            sage: bounds.satisfies_bounds(1)
            True
            sage: bounds.satisfies_bounds(1, 2)
            True
            sage: bounds.satisfies_bounds(1, 0, 1)
            False
            sage: bounds.satisfies_bounds(1, 0, 0)
            True
            sage: bounds.satisfies_bounds(1, 0, -1)
            True
            sage: bounds.satisfies_bounds(3, 7, 13/2)
            True
            sage: bounds.satisfies_bounds(-3, 7, -9/2)
            True
            sage: bounds.satisfies_bounds(10, 10, 6)
            False

            sage: X = Variety(3, [1, 1, 1/2, 1/6], [1, 0, 0, 0])
            sage: bounds = varieties.bounds.Bounds(X)
            sage: bounds.satisfies_bounds(-2, 1, 1/2)
            True
            sage: bounds.satisfies_bounds(1, 2, 3, 4)
            Traceback (most recent call last):
            ...
            NotImplementedError: Third Chern character or higher.
        """
        if len(args) <= 2:
            return True
        if len(args) == 3:
            if args[2] >= self.ch_min(args[0], args[1]):
                if args[2] <= self.ch_max(args[0], args[1]):
                    return True
                else:
                    return False
            else:
                return False

        else:
            raise NotImplementedError('Third Chern character or higher.')
