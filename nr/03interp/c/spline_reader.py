import numpy as np
import matplotlib.pyplot as plt

def splineplot(fignum, bctype, xa, ya, x, y_int):
    """
    Matplotlib script to plot dataset and interpolation.
    """
    plt.figure(figsize=(5,5))
    plt.title(r'Figure %d. $Dataset \hspace{0.5}\#%d$'%(fignum, fignum),
              fontsize='xx-large', fontweight='bold')
    plt.xlabel(r'$x$', fontsize='large')
    plt.ylabel(r'$y$', fontsize='large')
    plt.axvline(color='k')
    plt.axhline(color='k')
    plt.plot(xa, ya, 'o', label='$y_i(x_i)$')
    plt.plot(x, y_int, label='Cubic spline\n(%s)'%(bctype))
    plt.legend(loc='best', fontsize='x-large')
    plt.grid()
    plt.tight_layout()
    plt.show()

def read(filename):
    f = open(filename, 'r')
    data = []
    x = np.array([])
    y = np.array([])
    for line in f:
        data.append(line)
    for i in data:
        x = np.append(x, (i.split(('\t\t'))[0]))
        y = np.append(y, (i.split(('\t\t'))[1]))
    x = np.delete(x, 0)
    y = np.delete(y, 0)
    return x.astype(np.float), y.astype(np.float)

xa, ya = read("dataset.txt")
x, y_int = read("interp.txt")
splineplot(1, "clamped", xa, ya, x, y_int)
