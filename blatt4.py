from os import rename, listdir
import numpy as np
import matplotlib.pyplot as plt

def expected(runs):
    # left = -0.0552
    # right = -0.0325
    # step = (right - left) / 500
    # x = np.arange(left, right, step)
    # intervals_120_5 = np.zeros_like(x)
    # intervals_5_120 = np.zeros_like(x)
    # intervals_64_64 = np.zeros_like(x)

    res = np.zeros(runs)
    re = 0
    exp_u_120_5 = np.zeros(12501)
    exp_u_5_120 = np.zeros(12501)
    exp_u_64_64 = np.zeros(12501)

    var_u_5_120 = np.zeros(12501)
    var_u_120_5 = np.zeros(12501)
    var_u_64_64 = np.zeros(12501)

    stdabw_u_120_5 = np.zeros(12501)
    stdabw_u_5_120 = np.zeros(12501)
    stdabw_u_64_64 = np.zeros(12501)

    u_120_5 = np.zeros(runs)
    u_5_120 = np.zeros(runs)
    u_64_64 = np.zeros(runs)
    for i in range(runs):
        name = "build/uvalues/run_" + str(i)
        f = open(name, 'r')
        res[i] = float(f.readline())
        re += res[i]

        j = 0
        for line in f:
            u = line.split(",")
            exp_u_120_5[j] += float(u[0])
            exp_u_64_64[j] += float(u[1])
            exp_u_5_120[j] += float(u[2])
            j += 1
            u_120_5[i] = u[0]
            u_64_64[i] = u[1]
            u_5_120[i] = u[2]
        # intervals_120_5[int((u_120_5[i] - left)/step)] += 1
        # intervals_5_120[int((u_5_120[i] - left)/step)] += 1
        # intervals_64_64[int((u_64_64[i] - left)/step)] += 1

    # print(min(u_120_5))
    # print(max(u_120_5))

    # print(min(u_5_120))
    # print(max(u_5_120))

    # print(min(u_64_64))
    # print(max(u_64_64))

    exp_u_120_5 /= runs
    exp_u_5_120 /= runs
    exp_u_64_64 /= runs
    re /= runs
    # plt.title("Verteilung u 120x5 {}er Intervalle".format(step))
    # plt.plot(x, intervals_120_5)
    # plt.show()

    # plt.title("Verteilung u 5x120 {}er Intervalle".format(step))
    # plt.plot(x, intervals_5_120)
    # plt.show()

    # plt.title("Verteilung u 64x64 {}er Intervalle".format(step))
    # plt.plot(x, intervals_64_64)
    # plt.show()

    for i in range(runs):
        name = "build/uvalues/run_" + str(i)
        f = open(name, 'r')
        f.readline()
        j = 0
        for line in f:
            u = line.split(",")
            var_u_120_5[j] += (float(u[0]) - exp_u_120_5[j])**2
            var_u_64_64[j] += (float(u[1]) - exp_u_64_64[j])**2
            var_u_5_120[j] += (float(u[2]) - exp_u_5_120[j])**2
            j += 1

    var_u_120_5 /= runs - 1
    var_u_5_120 /= runs - 1
    var_u_64_64 /= runs - 1


    print("runs: {} exp_u_64_64: {}, {}".format(runs, exp_u_64_64[12500], exp_u_64_64[11999]) )
    return exp_u_5_120[-1], exp_u_64_64[-1], exp_u_120_5[-1], var_u_5_120[-1], var_u_64_64[-1], var_u_120_5[-1]


    # stdabw_u_120_5 = np.sqrt(var_u_120_5)
    # stdabw_u_5_120 = np.sqrt(var_u_5_120)
    # stdabw_u_64_64 = np.sqrt(var_u_64_64)

    # plt.xlabel(r"t")
    # plt.ylabel(r"u")
    # plt.title("Erwartungswert und Standardabweichung bei 120x5")
    # x = np.linspace(0,50,12501)
    # plt.plot(x, exp_u_120_5)
    # plt.plot(x, exp_u_120_5 - stdabw_u_120_5/2, "k--")
    # plt.plot(x, exp_u_120_5 + stdabw_u_120_5/2, "k--")
    # plt.show()

    # plt.xlabel(r"t" )
    # plt.ylabel(r"u" )
    # plt.title("Erwartungswert und Standardabweichung bei 5x120")
    # plt.plot(x, exp_u_5_120)
    # plt.plot(x, exp_u_5_120 - stdabw_u_5_120/2, "k--")
    # plt.plot(x, exp_u_5_120 + stdabw_u_5_120/2, "k--")
    # plt.show()

    # plt.xlabel(r"t" )
    # plt.ylabel(r"u" )
    # plt.title("Erwartungswert und Standardabweichung bei 64x64")
    # plt.plot(x, exp_u_64_64)
    # plt.plot(x, exp_u_64_64 + stdabw_u_64_64/2, "k--")
    # plt.plot(x, exp_u_64_64 - stdabw_u_64_64/2, "k--")
    # plt.show()

# ------------------------- Konvergenzplots

e_5_120 = []
e_64_64 = []
e_120_5 = []

v_5_120 = []
v_64_64 = []
v_120_5 = []
x = [500, 1000, 2000]

for runs in x:
    res = expected(runs)
    e_5_120.append(res[0])
    e_64_64.append(res[1])
    e_120_5.append(res[2])

    v_5_120.append(res[3])
    v_64_64.append(res[4])
    v_120_5.append(res[5])

print("e0=" + str(abs(e_64_64[0] -  e_64_64[2])))
print("e1=" + str(abs(e_64_64[1] -  e_64_64[2])))
print(np.log2( (abs(e_5_120[0] -  e_5_120[2])) / (abs(e_5_120[1] - e_5_120[2])) ))
print(np.log2( (abs(e_120_5[0] -  e_120_5[2])) / (abs(e_120_5[1] - e_120_5[2])) ))
print(np.log2( (abs(e_64_64[0] -  e_64_64[2])) / (abs(e_64_64[1] - e_64_64[2])) ))


print(np.log2( (abs(v_5_120[0] -  v_5_120[2])) / (abs(v_5_120[1] - v_5_120[2])) ))
print(np.log2( (abs(v_120_5[0] -  v_120_5[2])) / (abs(v_120_5[1] - v_120_5[2])) ))
print(np.log2( (abs(v_64_64[0] -  v_64_64[2])) / (abs(v_64_64[1] - v_64_64[2])) ))

# plt.xlabel(r"runs")
# plt.ylabel(r"E[u]")
# plt.title("Erwartungswertkonvergenz bei 5x120 t=50 mit Montecarlo")
# plt.plot(x, e_5_120, label="E[u] 5x120")
# plt.show()

# plt.xlabel(r"runs")
# plt.ylabel(r"E[u]")
# plt.title("Erwartungswertkonvergenz bei 120x5 t=50 mit Montecarlo")
# plt.plot(x, e_120_5, label="E[u] 120x5")
# plt.show()

# plt.xlabel(r"runs")
# plt.ylabel(r"E[u]")
# plt.title("Erwartungswertkonvergenz bei 64x64 t=50 mit Montecarlo")
# plt.plot(x, e_64_64, label="E[u] 64x64")
# plt.show()

# plt.xlabel(r"runs")
# plt.ylabel(r"V[u]")
# plt.title("Varianzkonvergenz bei 5x120 t=50 mit Montecarlo")
# plt.plot(x, v_5_120, label="V[u] 5x120")
# plt.show()

# plt.xlabel(r"runs")
# plt.ylabel(r"V[u]")
# plt.title("Varianzkonvergenz bei 120x5 t=50 mit Montecarlo")
# plt.plot(x, v_120_5, label="V[u] 120x5")
# plt.show()

# plt.xlabel(r"runs")
# plt.ylabel(r"V[u]")
# plt.title("Varianzkonvergenz bei 64x64 t=50 mit Montecarlo")
# plt.plot(x, v_64_64, label="V[u] 64x64")
# plt.show()
# expected(2000)
