"""
Make sure we show a decent error if a to old version of sage is being used
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
