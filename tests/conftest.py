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

    pr.dump_stats(test_dir.joinpath(f"{module}.pstats"))


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
