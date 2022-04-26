import time
import logging
import os
import subprocess
import threading

dim = "dim8/classic"  # Edit to test for other dimensions
dirs = os.listdir(f"../TT_library/{dim}")  # The path to the directory with all the files to run


def affine(name):
    filename = f"../TT_library/{dim}/{name}"
    dest_name = f"./results/{dim}/{name[:-3]}.txt"  # The path to where you want to store the results
    subprocess.run(["./a.out", filename, dest_name])


if __name__ == '__main__':
    time_format = "%(asctime)s: %(message)s"
    logging.basicConfig(format=time_format, level=logging.INFO, datefmt="%H:%M:%S")
    start = time.perf_counter()
    # Use threading
    threads = list()
    for file in dirs:
        x = threading.Thread(target=affine, args=(file,))
        threads.append(x)
        x.start()
    for thread in threads:
        thread.join()

    end = time.perf_counter()
    print(end - start)
