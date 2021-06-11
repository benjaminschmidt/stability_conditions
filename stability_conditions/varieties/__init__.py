r"""
Module for the numerical Chow ring of some special varieties. Additionally,
contains submodules for with additional features for specific varieties.

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

from . import p3
from . import bounds
