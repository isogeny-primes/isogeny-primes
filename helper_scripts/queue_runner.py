# main.py
import sys
from sage_code.queue_runner.queue_runner import QueueRunner

if __name__ == "__main__":
    command = sys.argv[1]
    runner = QueueRunner()
    if command == "fill":
        runner.fill_queue(*sys.argv[2:4])
    elif command == "clear":
        runner.clear_queue(sys.argv[2])
    elif command == "len":
        print(runner.len_queue(sys.argv[2]))
    elif command == "run":
        runner.run_function(sys.argv[2], "sage_code.lmfdb_isogeny_primes", "isogeny_primes")
    elif command == "dump":
        runner.dump(*sys.argv[2:4])
