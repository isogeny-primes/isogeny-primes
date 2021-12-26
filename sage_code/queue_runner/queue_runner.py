import importlib
import logging
import time
from typing import Optional

from .redis_queue import RedisFifoQueue

logger = logging.getLogger(__name__)


class QueueRunner:
    _backends = {"redis": RedisFifoQueue}

    def __init__(self, backend: str = "redis"):
        self._backend = self._backends[backend]()

    def fill_queue(self, queue_name, file_name):
        self._backend.push_file(queue_name, file_name)

    def clear_queue(self, queue_name):
        self._backend.clear(queue_name)

    def len_queue(self, queue_name):
        return self._backend.len(queue_name)

    def run_function(
        self,
        queue_name,
        module_name,
        function_name,
        out_queue_name=None,
        err_queue_name=None,
        default_arguments: Optional[dict] = None,
    ):
        module = importlib.import_module(module_name)
        function = getattr(module, function_name)
        if not default_arguments:
            default_arguments = {}
        if not out_queue_name:
            out_queue_name = f"{queue_name}_out"
        if not err_queue_name:
            err_queue_name = f"{queue_name}_err"
        while True:
            try:
                arguments = self._backend.pop(queue_name)
            except IndexError:
                break
            summary = {
                "arguments": arguments,
                "default_arguments": default_arguments,
            }
            arguments = {**default_arguments, **arguments}
            t = time.process_time()
            try:
                result = function(**arguments)
            except Exception as err:
                t = time.process_time() - t
                logger.warning(repr(err))
                summary["error"] = repr(err)
                summary["time"] = t
                self._backend.push(err_queue_name, summary)
                continue

            t = time.process_time() - t
            summary["result"] = result
            summary["time"] = t
            self._backend.push(out_queue_name, summary)

    def dump(self, queue_name, file_name):
        self._backend.dump(queue_name, file_name)
