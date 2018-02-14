# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys


def solve(phi, step):
    N = phi.shape[0] - 1
    omega = 2 / (1 + np.pi / N)
    phi_new = phi
    count = 0
    while count < step:
        for i in range(1, N):
            for j in range(1, N):
                phi_new[i, j] = phi[i, j] +\
                    (omega / 4) *\
                    (phi[i+1, j] +
                     phi_new[i-1, j] +
                     phi[i, j+1] +
                     phi_new[i, j-1] -
                     4*phi[i, j])
        phi = phi_new
        count = count + 1
    return phi


def main():
    N = 25
    Top = 0
    Right = 0
    Bottom = 0
    Left = 0
    num = 100

    for a in sys.argv:
        if a[0:2] == 'N=':
            N = int(a[2:])
        if a[0:2] == 'T=':
            Top = float(a[2:])
        if a[0:2] == 'R=':
            Right = float(a[2:])
        if a[0:2] == 'B=':
            Bottom = float(a[2:])
        if a[0:2] == 'L=':
            Left = float(a[2:])
        if a[0:4] == 'num=':
            num = int(a[4:])
    phi0 = np.zeros([N+1, N+1], np.float32)
    for i in range(N+1):
        phi0[0, i] = Top

    for i in range(N+1):
        phi0[N, i] = Bottom
    for i in range(N+1):
        phi0[i, N] = Right
    for i in range(N+1):
        phi0[i, 0] = Left
    solution = solve(phi0, num)
    # print(solution)
    # print(solution.shape)
    fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3D')
    ax = Axes3D(fig)
    axis = np.linspace(0, N+1, N+1)
    X, Y = np.meshgrid(axis, axis)
    if abs(max(solution.flatten())) >= abs(min(solution.flatten())):
        vmax = abs(max(solution.flatten()))
        vmin = 0-abs(max(solution.flatten()))
    else:
        vmax = abs(min(solution.flatten()))
        vmin = 0-abs(min(solution.flatten()))
    ax.plot_surface(X, Y, solution, cmap='seismic', vmax=vmax, vmin=vmin)
    plt.show()


main()
# vim: cc=42
