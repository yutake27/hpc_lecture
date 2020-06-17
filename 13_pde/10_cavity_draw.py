import numpy as np
import pandas as pd
from matplotlib import pyplot, cm
from mpl_toolkits.mplot3d import Axes3D


def read_txt(path, nx):
    df = pd.read_csv(path, sep=' ', header=None, usecols=lambda x: x != nx)
    return np.array(df)


def draw():
    nx = 41
    ny = 41

    x = np.linspace(0, 2, nx)
    y = np.linspace(0, 2, ny)
    X, Y = np.meshgrid(x, y)

    u = read_txt('u.txt', nx)
    v = read_txt('v.txt', nx)
    p = read_txt('p.txt', nx)


    fig = pyplot.figure(figsize=(11,7), dpi=100)
    # plotting the pressure field as a contour
    pyplot.contourf(X, Y, p, alpha=0.5, cmap=cm.viridis)  
    pyplot.colorbar()
    # plotting the pressure field outlines
    pyplot.contour(X, Y, p, cmap=cm.viridis)  
    # plotting velocity field
    pyplot.quiver(X[::2, ::2], Y[::2, ::2], u[::2, ::2], v[::2, ::2]) 
    pyplot.xlabel('X')
    pyplot.ylabel('Y')
    pyplot.savefig('10_cavity.png')
    # pyplot.show()lsl


if __name__ == '__main__':
    draw()