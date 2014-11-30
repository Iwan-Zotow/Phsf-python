#!/usr/bin/python

# requires python 3.x. tested with 3.4.2

import math
import struct

import matplotlib.pyplot as plt

from H1Du import H1Du

def read_record_short(phsf):
    """
    Read MODE0 PhSF particle record
    """
    LATCH = struct.unpack('i', phsf.read(4))[0]
    E = struct.unpack('f', phsf.read(4))[0]
    X = struct.unpack('f', phsf.read(4))[0]
    Y = struct.unpack('f', phsf.read(4))[0]
    U = struct.unpack('f', phsf.read(4))[0]
    V = struct.unpack('f', phsf.read(4))[0]
    W = math.sqrt(1.0 - (U*U + V*V))
    WT = struct.unpack('f', phsf.read(4))[0]
    if (WT < 0.0): # weight sign is in fact Z directional cosine sign
        WT = -WT
        W  = -W
        
    return (LATCH, E, X, Y, U, V, W, WT)
    
def read_record_long(phsf):
    """
    Read MODE2 PhSF particle record
    """
    LATCH = struct.unpack('i', phsf.read(4))[0]
    E = struct.unpack('f', phsf.read(4))[0]
    X = struct.unpack('f', phsf.read(4))[0]
    Y = struct.unpack('f', phsf.read(4))[0]
    U = struct.unpack('f', phsf.read(4))[0]
    V = struct.unpack('f', phsf.read(4))[0]
    W = math.sqrt(1.0 - (U*U + V*V))
    WT = struct.unpack('f', phsf.read(4))[0]
    ZLAST = struct.unpack('f', phsf.read(4))[0]
    if (WT < 0.0): # weight sign is in fact Z directional cosine sign
        WT = -WT
        W  = -W

    return (LATCH, E, X, Y, U, V, W, WT, ZLAST)

with open("C25.egsphsp1", "rb") as phsf:
    mode = phsf.read(5)
    print(mode) # shall be MODE0 or MODE2

    if mode == b"MODE0":
        print("SHORT BEAMNRC file found")
    elif mode == b"MODE2":
        print("LONG BEAMNRC file found")
    else:
        print("Unknown phase space file format")
        exit()
    
    NPPHSP = struct.unpack('i', phsf.read(4))[0]
    print(NPPHSP) # total nof of particle records

    NPHOTPHSP = struct.unpack('i', phsf.read(4))[0]
    print(NPHOTPHSP) # total nof of photon records

    EKMAX = struct.unpack('f', phsf.read(4))[0]
    print(EKMAX) # max kinetic energy, MeV
    
    EKMIN = struct.unpack('f', phsf.read(4))[0]
    print(EKMIN) # min electron kinetic energy, MeV, shall be ECUT-0.511
    
    NINCP = struct.unpack('i', phsf.read(4))[0]
    print(NINCP) # nof original incident particles    
    print("================ End of PhSF header ======================")

    if mode == b"MODE2":
        dummy = phsf.read(7) # skip fortran garbage, for direct access of record with size 32bytes
    else: # MODE0 file
        dummy = phsf.read(3) # skip fortran garbage, for direct access of record with size 28bytes

    he = H1Du(50, 0.01, 1.33)
    hx = H1Du(31, -15.5, 15.5)
    hy = H1Du(31, -15.5, 15.5)

    for i in range (0, NPPHSP):
        (LATCH, E, X, Y, U, V, W, WT, ZLAST) = read_record_long(phsf)
        if LATCH == 8388608: #8388608=2^23, this is photon, see PIRS-509, page #96
            if E<0:
                E = -E

            he.fill(E, WT)
            hx.fill(X, WT)
            hy.fill(Y, WT)

    print(he.nof_events(), he.integral())

    X = []
    Y = []

    step = he.step()
    for i in range (-1, he.size()+1):
        x = he.lo() + (float(i) + 0.5)*step
        d = he[i] # data from bin with index i
        y = d[0]  # first part of bin is collected weights
        X.append(x)
        Y.append(y)

    width = 0.8*step
    p1 = plt.bar(X, Y, width, color='r')

    plt.xlabel('Energy(MeV)')
    plt.ylabel('N of photons')
    plt.title('Energy distribution')
    
    plt.grid(True);
    plt.tick_params(axis='x', direction='out')
    plt.tick_params(axis='y', direction='out')

    plt.show()
