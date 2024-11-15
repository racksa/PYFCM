import numpy as np
import fcm
from utils import *
import matplotlib.pyplot as plt
import random as random

# System parameters
dt = 1
frames = 1

# FCM parameters
N = 10
rh = 1.
alpha = 1.
beta = 10.
eta = 12.
nx = 128
ny = 128
nz = 128
boxsize = 100
repeat = 1
checkerror = False
rotation = False

# Global arrays
Y = np.zeros((N, 3))
F = np.zeros((N, 3))
T = np.zeros((N, 3))
V = np.zeros((N, 3))
W = np.zeros((N, 3))
# Test case
for i in range(N):
    Y = np.random.uniform(0, boxsize, size=(N, 3))
    F = np.random.uniform(-1, 1, size=(N, 3))
    F[:,2] = 0

    Y[:,0] = np.linspace(0.1, boxsize, N)
    Y[:,1] = np.linspace(0.1, boxsize/2, N)
    F[:,0] = np.linspace(0.1, 1, N)
    F[:,1] = np.linspace(0.1, 4, N)

# Initialise pars
pars = Pars(N, rh, alpha, beta, eta, nx, ny, nz, boxsize, repeat, checkerror, rotation)

# Initialise FCM solver
fcm_solver = fcm.FCM(pars)
fcm_solver.print_pars()

for t in range(frames):
    fcm_solver.hydrodynamic_solve(Y,F,T,V,W)


# Plot for info
fig = plt.figure()
ax = fig.add_subplot()
ax.vlines(np.linspace(0, fcm_solver.Lx, fcm_solver.Mx), 0, fcm_solver.Ly, colors='lightgray', linestyles='--', linewidth=0.7)
ax.hlines(np.linspace(0, fcm_solver.Ly, fcm_solver.My), 0, fcm_solver.Lx, colors='lightgray', linestyles='--', linewidth=0.7)
print(V)
ax.quiver(Y[:,0], Y[:,1], F[:,0], F[:,1], color='b')
ax.quiver(Y[:,0], Y[:,1], V[:,0], V[:,1], color='r')
ax.set_aspect('equal')
ax.set_xlim(0, boxsize)
ax.set_ylim(0, boxsize)
plt.show()