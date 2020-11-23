import numpy as np
import matplotlib.pyplot as plt

filename = "/home/luke/CLionProjects/MCGRP_debug/instance/debug_test/log/2020-Nov-25 10:04:36/TWe1_t.dat.prob"

if __name__ == '__main__':
    arr = np.loadtxt(filename)

    plt.imshow(arr)
    plt.colorbar()
    plt.show()

