import logging
import subprocess
import time

# dim = "dim8"
# filenames = ["classic_1.tt", "classic_2.tt", "classic_3.tt", "classic_4.tt"]
# root = f"../TT_library/{dim}/classic/"
dim = "dim6"
filenames = ["q_6_1.tt", "q_6_2.tt", "q_6_3.tt", "q_6_4.tt"]
root = f"../TT_library/{dim}/"
count = 10


def affine(name):
    destination = f"results/computationally_testing/{dim}/{name}"
    start = time.perf_counter()
    for i in range(count):
        subprocess.run(["./affine", "-w", f"{destination[:-3]}_{i}.txt", f"{root}{name}"])
    total_time = time.perf_counter() - start
    average = total_time / count
    print(average)
    wf = open(f"results/computationally_testing/{dim}_{name[:-3]}.txt", "w+")
    wf.write(f"{average}")
    wf.close()


if __name__ == '__main__':
    time_format = "%(asctime)s: %(message)s"
    logging.basicConfig(format=time_format, level=logging.INFO, datefmt="%H:%M:%S")
    for f in filenames:
        affine(f)
