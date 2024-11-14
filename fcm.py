import numpy as np
from utils import *

class FCM:
    
    def __init__(self, pars):
        self.pars = pars
        for field in pars.__dataclass_fields__:
            # Use setattr to assign new attributes in MyClass with names from the dataclass
            setattr(self, field, getattr(pars, field))
        
        self.grid_size = self.nx*self.ny*self.nz
        self.fft_grid_size = (self.nx/2+1)*self.ny*self.nz
        self.dx = self.boxsize/self.nx
        self.ngd = round(self.alpha*self.beta)

        self.init_grid()

    def print_pars(self):
        for arr, val in self.pars:
            print(f"{arr:<12} =  {val}")

    def init_grid(self):
        self.hx = np.zeros(self.grid_size)
        self.hy = np.zeros(self.grid_size)
        self.hz = np.zeros(self.grid_size)

        self.fx = np.zeros(self.fft_grid_size)
        self.fy = np.zeros(self.fft_grid_size)
        self.fz = np.zeros(self.fft_grid_size)

    def reset_grid(self):
        self.hx.fill(0)
        self.hy.fill(0)
        self.hz.fill(0)

        self.fx.fill(0)
        self.fy.fill(0)
        self.fz.fill(0)

    def box_particle(self):
        pass

    def spatial_hasing(self):
        pass
    

    def hydrodynamic_solve(self, Y, F, T, V, W):
        self.Y = Y
        self.F = F
        self.T = T
        self.V = V
        self.W = W



        self.stokes_drag()
        return 0

    def stokes_drag(self):
        factor = 1./(6*np.pi*self.rh)
        for i, f in enumerate(self.F):
            self.V[i] = factor*self.F[i]

