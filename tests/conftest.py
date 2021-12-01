from pytest import fixture

import logging
import cProfile
import os

from sage_code.config import PROFILING_DIR

logger = logging.getLogger(__name__)


@fixture(autouse=True)
def profile(request):
    function = request.function
    test_names = [key for key in request.keywords if key.startswith(function.__name__)]
    test_name = test_names[0] if test_names else function.__name__

    if not os.environ.get("PROFILE"):
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
