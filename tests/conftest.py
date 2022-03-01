"""conftest.py

Test configs.

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

from pytest import fixture

import logging
import cProfile
import os

from sage_code.config import PROFILING_DIR

logger = logging.getLogger(__name__)


@fixture(scope="module", autouse=True)
def profile_module(request):
    module = request.module.__name__
    scope = os.environ.get("PROFILE_SCOPE", "").lower()
    allowed_scopes = ["module", "function", ""]
    assert (
        scope in allowed_scopes
    ), f"PROFILE_SCOPE should be module or function not {scope}"

    if not scope == "module":
        yield module
        return

    with cProfile.Profile() as pr:
        yield module

    test_dir = PROFILING_DIR.joinpath(module)
    try:
        os.mkdir(test_dir)
    except FileExistsError:
        pass

    pr.dump_stats(test_dir.joinpath(f"module_tests.pstats"))


@fixture(autouse=True)
def profile_function(request):
    function = request.function
    test_names = [key for key in request.keywords if key.startswith(function.__name__)]
    test_name = test_names[0] if test_names else function.__name__

    scope = os.environ.get("PROFILE_SCOPE", "").lower()

    if not scope == "function":
        yield test_name
        return

    with cProfile.Profile() as pr:
        yield test_name

    test_dir = PROFILING_DIR.joinpath(function.__module__)
    try:
        os.mkdir(test_dir)
    except FileExistsError:
        pass

    pr.dump_stats(test_dir.joinpath(f"{test_name}.pstats"))
