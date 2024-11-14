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
eta = 1
nx = 128
ny = 128
nz = 128
boxsize = 100
repeat = 1
checkerror = False

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

# Initialise pars
pars = Pars(N, rh, alpha, beta, eta, nx, ny, nz, boxsize, repeat, checkerror)

# Initialise FCM solver
fcm_solver = fcm.FCM(pars)
fcm_solver.print_pars()

for t in range(frames):
    fcm_solver.hydrodynamic_solve(Y,F,T,V,W)


fig = plt.figure()
ax = fig.add_subplot()
ax.quiver(Y[:,0], Y[:,1], F[:,0], F[:,1], color='b')
ax.quiver(Y[:,0], Y[:,1], V[:,0], V[:,1], color='r')
ax.set_aspect('equal')
ax.set_xlim(0, boxsize)
ax.set_ylim(0, boxsize)
plt.show()