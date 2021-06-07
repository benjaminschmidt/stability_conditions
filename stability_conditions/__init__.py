r"""
Sage library for computations with (very weak) Bridgeland stability
conditions.


.. TODO::

    - Add examples.
    - Add tests for checking imports.


EXAMPLES::

<Lots and lots of examples>


TESTS::

    Check whether imports are working.

    sage: from stability_conditions import *

    sage: type(library)
    <class 'module'>

    sage: type(slope)
    <class 'module'>

    sage: type(tilt)
    <class 'module'>

    sage: type(tilt_two)
    <class 'module'>
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
# from .all import *

from . import library
from . import slope
from . import tilt
from . import tilt_two
from .varieties import variety
from .varieties.variety import Element, p, Variety
from . import varieties
