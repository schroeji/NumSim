import numpy as np
import matplotlib.pyplot as plt

def plot_res_level():
    name = "build/mg_residuals.txt"
    f = open(name, 'r')
    res  = []
    levels = 4
    i = 0
    l = []
    x = []
    for line in f:
        vals = line.split(",")
        l.append(int(vals[0]))
        res.append(float(vals[1]))
        if(int(vals[2]) > 0):
            x.append(-int(vals[0]) + 2*levels + 1)
        else:
            x.append(int(vals[0]))
        i += 1
        if i >= 2*levels + 2:
           break
    print(l)
    print(res)
    plt.vlines(levels + 0.5, 0, max(res), linestyles="dashed", color="r")
    plt.plot(x, res)
    plt.xticks(x, l)
    plt.xlabel("Level")
    plt.ylabel("Residuum")
    plt.show()

def konvergenz_analyse():
    prefix = "build/counts_3_"
    sizes = [16,32,64,128]
    avgs = np.zeros_like(sizes)
    for i,size in enumerate(sizes):
        name = prefix + str(size)
        f = open(name, 'r')
        line_count = 0
        for line in f:
            line_count += 1
            avgs[i] += int(line)
        avgs[i] /= line_count
    print avgs
    plt.plot(sizes, avgs)
    plt.xticks(sizes)
    plt.xlabel("# Zellen pro Dimension")
    plt.ylabel("# Iterationen fr res 10-3")
    plt.title("MG Konvergenzverhalten")
    plt.show()
# plot_res_level()
konvergenz_analyse()