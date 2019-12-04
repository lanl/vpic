import numpy as np
import sys


# Begin functions
def mu(x):
    ss = np.linspace(0, x, 1000)
    dss = ss[1]-ss[0]
    mux = (2/np.sqrt(np.pi))*np.sum(np.exp(-ss)*np.sqrt(ss))*dss
    return mux


def load_slice(file):
    fstr = file
    fd = open(fstr, "rb")
    arr = np.fromfile(fd, dtype=np.float64)
    fd.close()
    return arr
# end functions


datafile = "./JX_data"

dtcoll = 5e-3  # nu_0*dt_jx from info file
E = 0.5*1  # 0.5*me*vdre^2/Te

nuei = mu(E)*np.power(E, -1.5) / np.sqrt(2)

jx = load_slice(datafile)
jx = jx/jx[0]

numsteps = 80  # hard coded to match deck
t = np.linspace(0, numsteps*dtcoll, numsteps)
theory = jx[0]*np.exp(-t*nuei)

show_plot = 0
if show_plot:
    import matplotlib.pyplot as plt
    fig, (ax1) = plt.subplots(nrows=1)
    im = ax1.plot(t, jx, 'b')
    im2 = ax1.plot(t, theory, 'r')
    # plt.axis([0, 0.1, 0, 0.04])
    plt.show()

# Check data is "close" to theory, tracking max error
max_rel = 0.0
max_ab = 0.0
for i in range(numsteps):
    a = jx[i]
    b = theory[i]
    rel = (1.0 - (a/b)) * 100.0
    ab = a-b

    if rel > max_rel:
        max_rel = rel

    if ab > max_ab:
        max_ab = ab
    # print("{} vs {} diff rel {} abs {}".format(a, b, rel, ab) )
print("max errors -- relative {}% absolute {}".format(max_rel, max_ab))

rel_tol = 1.0  # 1%
if max_rel > rel_tol:
    # Error out
    print("=> Fail")
    sys.exit(1)
else:
    print("=> Pass")
