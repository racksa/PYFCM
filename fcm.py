import numpy as np
from utils import *
from data import *

class FCM:
    
    def __init__(self, pars):
        self.pars = pars

        self.init_config()

        self.init_fcm_vars()

        self.init_arrays()

        self.init_grid()

    def init_config(self):
        for field in self.pars.__dataclass_fields__:
            # Use setattr to assign new attributes in MyClass with names from the dataclass
            setattr(self, field, getattr(self.pars, field))

        self.grid_size = int(self.nx*self.ny*self.nz)
        self.fft_grid_size = int((self.nx/2+1)*self.ny*self.nz)
        self.dx = self.boxsize/self.nx
        self.ngd = round(self.alpha*self.beta)

    def init_fcm_vars(self):
        self.sigma = self.rh/np.pi**.5
        self.sigmadip = self.rh/(6.*np.pi)**(1./3)

        self.Sigma_grid = self.sigma

        self.Rc_fac = self.eta*self.alpha
        self.Rc = self.Rc_fac*self.dx
        self.Rcsq = self.Rc**2

        self.Lx = self.boxsize
        self.Ly = self.boxsize/self.nx*self.ny
        self.Lz = self.boxsize/self.nx*self.ny
        self.Lmin = min(min(self.Lx, self.Ly), self.Lz)
        self.Lmax = max(max(self.Lx, self.Ly), self.Lz)
        self.Mx = int(max(self.Lx/self.Rc, 3))
        self.My = int(max(self.Ly/self.Rc, 3))
        self.Mz = int(max(self.Lz/self.Rc, 3))

        self.ncell = self.Mx*self.My*self.Mz
        self.mapsize = 13*self.ncell

        self.volume_frac = (self.N*4/3*np.pi*self.rh**3) / (self.Lx*self.Ly*self.Lz)

    def init_arrays(self):
        self.particle_cellhash = np.zeros(self.N, dtype=int)
        self.particle_index = np.zeros(self.N, dtype=int)
        self.sortback_index = np.zeros(self.N, dtype=int)
        self.key_buf = np.zeros(self.N, dtype=int)
        self.index_buf = np.zeros(self.N, dtype=int)

        self.cell_start = np.zeros(self.ncell, dtype=int)
        self.cell_end = np.zeros(self.ncell, dtype=int)

        self.map = np.zeros(self.mapsize, dtype=int)

        bulkmap(self.map, self.Mx, self.My, self.Mz)

        

    def print_pars(self):
        for arr, val in self.pars:
            print(f"{arr:<12} =  {val}")
        
        print(f"{'dx':<12} =  {self.dx}")
        print(f"{'Boxsize':<12} =  {self.boxsize}")
        print(f"{'Grid points':<12} =  ({self.nx}, {self.ny}, {self.nz})")
        print(f"{'Cell number':<12} =  ({self.Mx}, {self.My}, {self.Mz})")
        print(f"{'Rc/a':<12} =  {self.Rc/self.rh}")

    def init_grid(self):
        self.h =np.zeros((self.grid_size, 3))

    def reset_grid(self):
        self.h.fill(0)

    def box_particle(self):
        box(self.Y, self.N, self.Lx, self.Ly, self.Lz )

    def spatial_hasing(self):
        create_hash(self.particle_cellhash, self.Y, self.N, self.Mx, self.My, self.Mz, self.Lx, self.Ly, self.Lz)
        self.particle_index = np.arange(self.N)

    def sort_particles(self):
        sorted_index = np.argsort(self.particle_cellhash)
        self.particle_index = self.particle_index[sorted_index]
        self.particle_cellhash = self.particle_cellhash[sorted_index]
        self.Y[:] = self.Y[sorted_index]
        self.F[:] = self.F[sorted_index]
        self.T[:] = self.T[sorted_index]
        self.V[:] = self.V[sorted_index]
        self.W[:] = self.W[sorted_index]

        self.cell_start.fill(0)
        self.cell_end.fill(0)

        create_cell_list(self.particle_cellhash, self.cell_start, self.cell_end, self.N)

    def spread(self):
        spread_kernel(self.h,\
                  self.Y, self.T, self.F,\
                    self.N, self.ngd,\
                        self.sigma, self.sigmadip, self.Sigma_grid,\
                            self.dx, self.nx, self.ny, self.nz,\
                                self.particle_index, self.cell_start, self.cell_end,
                                self.rotation)

    def fft_solve(self):
        pass

    def gather(self):
        pass

    def correction(self):
        pass

    def sortback(self):
        self.sortback_index = np.arange(self.N)
        sorted_index = np.argsort(self.particle_index)
        self.particle_index = self.particle_index[sorted_index]
        self.sortback_index = self.sortback_index[sorted_index]

        self.Y[:] = self.Y[sorted_index]
        self.F[:] = self.F[sorted_index]
        self.T[:] = self.T[sorted_index]
        self.V[:] = self.V[sorted_index]
        self.W[:] = self.W[sorted_index]
    

    def hydrodynamic_solve(self, Y, F, T, V, W):
        self.Y = Y
        self.F = F
        self.T = T
        self.V = V
        self.W = W

        self.reset_grid()

        self.box_particle()

        self.spatial_hasing()

        self.sort_particles()

        self.spread()

        # self.fft_solve()

        # self.gather()

        # self.correction()

        self.stokes_drag()

        self.sortback()

        return 0

    def stokes_drag(self):
        factor = 1./(6*np.pi*self.rh)
        for i, f in enumerate(self.F):
            self.V[i] = factor*self.F[i]

