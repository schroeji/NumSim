from os import rename, listdir
import numpy as np

def expected():
    runs = 2000
    re = 0
    exp_u_120_5 = np.zeros(12501)
    var_u_120_5 = np.zeros(12501)
    exp_u_5_120 = np.zeros(12501)
    var_u_5_120 = np.zeros(12501)
    exp_u_64_64 = np.zeros(12501)
    var_u_64_64 = np.zeros(12501)
    for i in range(runs):
        name = "build/uvalues/run_" + str(i)
        f = open(name, 'r')
        re += float(f.readline())
        j = 0
        for line in f:
            u = line.split(",")
            exp_u_120_5[j] += float(u[0])
            exp_u_64_64[j] += float(u[1])
            exp_u_5_120[j] += float(u[2])
            j += 1
    exp_u_120_5 /= runs
    exp_u_5_120 /= runs
    exp_u_64_64 /= runs
    re /= runs
    for i in range(runs):
        name = "build/uvalues/run_" + str(i)
        f = open(name, 'r')
        f.readline()
        j = 0
        for line in f:
            u = line.split(",")
            var_u_120_5[j] += (float(u[0]) - exp_u_120_5)**2
            var_u_64_64[j] += (float(u[1]) - exp_u_64_64)**2
            var_u_5_120[j] += (float(u[2]) - exp_u_5_120)**2
            j += 1
    var_u_120_5 /= runs - 1
    var_u_5_120 /= runs - 1
    var_u_64_64 /= runs - 1

expected()
