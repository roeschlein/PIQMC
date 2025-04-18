import numpy as np

with open('active_window', 'r') as AW:
    lines = AW.readlines()
    HO_index = int(lines[0])
    AWS = int(lines[1])
    AWE = int(lines[2])

print(f'Lowest KS orbital index in active window read from active_window:{AWS}')
print(f'Highest KS orbital index in active window read from active_window:{AWE}')



active_window_start = AWS - 1 #KS index - 1
active_window_end = AWE - 1 #KS index - 1
active_window_size = active_window_end-active_window_start+1

with open('./wfc1.dat', 'rb') as f:
    # Moves the cursor 4 bytes to the right
    f.seek(4)

    ik = np.fromfile(f, dtype='int32', count=1)[0]
    xk = np.fromfile(f, dtype='float64', count=3)
    ispin = np.fromfile(f, dtype='int32', count=1)[0]
    gamma_only = bool(np.fromfile(f, dtype='int32', count=1)[0])
    scalef = np.fromfile(f, dtype='float64', count=1)[0]

    # Move the cursor 8 byte to the right
    f.seek(8, 1)

    ngw = np.fromfile(f, dtype='int32', count=1)[0]
    igwx = np.fromfile(f, dtype='int32', count=1)[0]
    npol = np.fromfile(f, dtype='int32', count=1)[0]
    nbnd = np.fromfile(f, dtype='int32', count=1)[0]

    # Move the cursor 8 byte to the right
    f.seek(8, 1)

    b1 = np.fromfile(f, dtype='float64', count=3)
    b2 = np.fromfile(f, dtype='float64', count=3)
    b3 = np.fromfile(f, dtype='float64', count=3)

    f.seek(8,1)
    
    mill = np.fromfile(f, dtype='int32', count=3*igwx)
    mill = mill.reshape( (igwx, 3) ) 

    evc = np.zeros( (active_window_size, npol*igwx), dtype="complex128")

    f.seek(8,1)
    for i in range(active_window_start):
        np.fromfile(f, dtype='complex128', count=npol*igwx)
        f.seek(8,1)
    for i in range(active_window_size):
        evc[i,:] = np.fromfile(f, dtype='complex128', count=npol*igwx)
        f.seek(8, 1)



for i in range(len(evc)):
    with open('./evc'+str(i+1+active_window_start),'w') as evc_write:
        for j in range(len(evc[i])):
            evc_write.write(str(np.real(evc[i][j]))+'   '+str(np.imag(evc[i][j]))+'\n')

with open('./gkvectors','w') as gkvec:
    for i in range(len(evc[1])):
        gkvec.write(str(b1[0]*mill[i][0])+'   '+str(b2[1]*mill[i][1])+'   '+str(b3[2]*mill[i][2])+'\n')

