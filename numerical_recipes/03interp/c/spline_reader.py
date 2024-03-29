import numpy as np
import matplotlib.pyplot as plt

def splineplot(fignum, dset, bctype, xa, ya, x, y_int):
    """
    Matplotlib script to plot dataset and interpolation.

    Parameters
    ----------
    fignum : integer
        Figure number to be displayed on plot.

    dset : integer
        Dataset number to be displayed on plot.

    bctype : string
        Boundary condition type to be displayed on plot.

    xa : array_like
        Monotonic array of independent values.
    ya : 
        Array of dependent values.
    x : 
        Monotonic array of independent values.
    y_int : 
        Array of interpolated values.
    """
    plt.figure(figsize=(5,5))
    plt.title(r'Figure %d. $Dataset \hspace{0.5}\#%d$'%(fignum, dset),
              fontsize='xx-large', fontweight='bold')
    plt.xlabel(r'$x$', fontsize='large')
    plt.ylabel(r'$y$', fontsize='large')
    plt.axvline(color='k')
    plt.axhline(color='k')
    plt.plot(xa, ya, 'o', label='$y_i(x_i)$')
    plt.plot(x, y_int, label='Cubic spline (%s)'%(bctype))
    plt.legend(loc='best', fontsize='x-large')
    plt.grid()
    plt.tight_layout()
    plt.show()

def read(filename):
    f = open(filename, 'r')
    data = []
    x = np.aray([])
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
splineplot(1, "", xa, ya, x, y_int)
