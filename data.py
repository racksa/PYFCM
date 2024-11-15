import numpy as np

def icell(Mx, My, Mz, x, y, z):
    xi = int((x+Mx) % Mx)
    yi = int((y+My) % My)
    zi = int((z+Mz) % Mz)

    return xi + (yi + zi*My)*Mx

def bulkmap(map, Mx, My, Mz):
    for iz in range(Mz):
        for iy in range(My):
            for ix in range(Mx):
                tempmap = icell(Mx, My, Mz, ix, iy, iz)
                imap = tempmap*13
                map[imap] = icell(Mx, My, Mz, ix+1, iy, iz)
                map[imap+1] = icell(Mx, My, Mz, ix+1, iy+1, iz)
                map[imap+2] = icell(Mx, My, Mz, ix, iy+1, iz)
                map[imap+3] = icell(Mx, My, Mz, ix-1, iy+1, iz)
                map[imap+4] = icell(Mx, My, Mz, ix+1, iy, iz-1)
                map[imap+5] = icell(Mx, My, Mz, ix+1, iy+1, iz-1)
                map[imap+6] = icell(Mx, My, Mz, ix, iy+1, iz-1)
                map[imap+7] = icell(Mx, My, Mz, ix-1, iy+1, iz-1)
                map[imap+8] = icell(Mx, My, Mz, ix+1, iy, iz+1)
                map[imap+9] = icell(Mx, My, Mz, ix+1, iy+1, iz+1)
                map[imap+10] = icell(Mx, My, Mz, ix, iy+1, iz+1)
                map[imap+11] = icell(Mx, My, Mz, ix-1, iy+1, iz+1)
                map[imap+12] = icell(Mx, My, Mz, ix, iy, iz+1)

def create_cell_list(particle_cellhash, cell_start, cell_end, N):
    c1, c2 = 0, 0
    for n in range(N):
        c2 = particle_cellhash[n]
        if (n>0):
            c1 = particle_cellhash[n-1]
        if (c1!=c2 or n==0):
            cell_start[c2] = n
            if(n>0):
                cell_end[c1] = n
        if (n==n-1):
            cell_end[c2] = N

def box(Y, N, Lx, Ly, Lz):
    for n in range(N):
        images(Y[n][0], Lx)
        images(Y[n][1], Ly)
        images(Y[n][2], Lz)

def images(x, boxsize):
    x -= np.floor(x/boxsize)*boxsize

    if (x==boxsize):
        x = 0
    
def create_hash(hash_table, Y, N, Mx, My, Mz, Lx, Ly, Lz):
    for n in range(N):
        xc = int(Y[n][0]/Lx*Mx)
        yc = int(Y[n][1]/Ly*My)
        zc = int(Y[n][2]/Lz*Mz)

        hash_table[n] = xc + (yc + zc*My)*Mx

def spread_kernel(f,\
                  Y, T, F,\
                    N, ngd,\
                        sigma, sigmadip, Sigma,\
                            dx, nx, ny, nz,\
                                particle_index, start, end,
                                rotation):
    pass

def gather_kernel(ux, uy, uz,\
                  Y, YTEMP, WTEMP,\
                    N, nga,\
                        sigma, sigmadip, Sigma,\
                            dx, nx, ny, nz,\
                                particle_index, start, end,\
                                    rotation):
    pass

def flow_solve(fk_x, fk_y, fk_z, uk_x, uk_y, uk_z, nx, ny, nz, Lx, Ly, Lz):
    pass
