import numpy as np

with open('trim_xyz', 'r') as tr:
    x_grid = tr.readline().split()
    y_grid = tr.readline().split()
    z_grid = tr.readline().split()

x_size = int(x_grid[1])-int(x_grid[0])+1
y_size = int(y_grid[1])-int(y_grid[0])+1
z_size = int(z_grid[1])-int(z_grid[0])+1

with open('matsubara_params', 'r') as mp:
    mp.readline()
    mp.readline()
    T = float(mp.readline().split()[0])

with open('delta_xyz') as dl:
    box = dl.readline().split()
    box_grid = dl.readline().split()

ax = float(box[0])/(int(box_grid[0]))
ay = float(box[1])/(int(box_grid[1]))
az = float(box[2])/(int(box_grid[2]))

pre = (ax*ay*az*T)/(8*np.pi*51.4220674763)

def row_3(x, y, z, nx, ny, nz, ax, ay, az):
    x -= 1
    y -= 1
    z -= 1
    ith_row = np.zeros(nx*ny*nz)
    ith_row[np.mod(x+1,nx)*ny*nz + y*nz + z] = 1/ax**2
    ith_row[np.mod(x-1,nx)*ny*nz + y*nz + z] = 1/ax**2
    ith_row[x*ny*nz + np.mod(y+1,ny)*nz + z] = 1/ay**2
    ith_row[x*ny*nz + np.mod(y-1,ny)*nz + z] = 1/ay**2
    ith_row[x*ny*nz + y*nz + np.mod(z+1,nz)] = 1/az**2
    ith_row[x*ny*nz + y*nz + np.mod(z-1,nz)] = 1/az**2
    ith_row[x*ny*nz + y*nz + z] = -(2/ax**2 + 2/ay**2 + 2/az**2)
    return ith_row

# Set lattice size
N = x_size*y_size*z_size
lap = np.zeros((N, N))
i = 0
for x in range(1, x_size+1):
    for y in range(1, y_size+1):
        for z in range(1, z_size+1):
            lap[i, :] = pre*row_3(x, y, z, x_size, y_size, z_size,ax, ay, az)
            i += 1

# Extract sparse entries
rows, cols = np.nonzero(lap)
values = lap[rows, cols]
sparse_index = len(values)

# Fortran expects 1-based indexing
rows_1based = rows + 1
cols_1based = cols + 1

# Write binary file in the format Fortran expects
with open("laplacian.bin", "wb") as f:
    np.array([sparse_index], dtype=np.int64).tofile(f)
    values.astype(np.float64).tofile(f)
    rows_1based.astype(np.int64).tofile(f)
    cols_1based.astype(np.int64).tofile(f)

