import time
import logging
import os
import subprocess
import threading

dim = "dim8/classic"  # Edit to test for other dimensions
dirs = os.listdir(f"../ea-equivalence/resources/TT_library/{dim}")


def affine(name):
    filename = f"../ea-equivalence/resources/TT_library/{dim}/{name}"
    dest_name = f"./results/{dim}/{name[:-3]}.txt"
    subprocess.run(["./a.out", filename, dest_name])


if __name__ == '__main__':
    time_format = "%(asctime)s: %(message)s"
    logging.basicConfig(format=time_format, level=logging.INFO, datefmt="%H:%M:%S")
    start = time.perf_counter()
    threads = list()
    for file in dirs:
        x = threading.Thread(target=affine, args=(file,))
        threads.append(x)
        x.start()
    for thread in threads:
        thread.join()
    end = time.perf_counter()
    print(end - start)
