"""
Make sure we show a decent error if a to old version of sage is being used

    ====================================================================

    This file is part of Isogeny Primes.

    Copyright (C) 2022 Barinder S. Banwait and Maarten Derickx

    Isogeny Primes is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

    The authors can be reached at: barinder.s.banwait@gmail.com and
    maarten@mderickx.nl.

    ====================================================================

"""

try:
    from sagemath.check_version import check_version

    check_version(">=9.4")
except ModuleNotFoundError as err:
    raise RuntimeError(
        """
        Not all requirements of Isogeny Primes are installed. Please do
            sage -pip install -r requirements.txt
        Before continuing.
        """
    ) from ModuleNotFoundError

from .monkey_patch_sage import monkey_patch

monkey_patch()
