#!/usr/bin/env python

# requires python 3.x. tested with 3.4.2

import math
import struct

import matplotlib.pyplot as plt

from H1Dn import H1Dn

def make_scale(nof_bins_hi2me):
    """
    """
    lo = 0.01
    me = 1.170001
    hi = 1.330001
    
    # bins in hi-to-me region
    step = (hi - me)/float(nof_bins_hi2me)
    
    scale = []

    # filling lo-to-me region with the same bin size
    prev = me
    while True:
        scale.append(prev)
        prev -= step
        if prev < lo:
            break
            
    if scale[-1] != lo:
        scale.append(lo)
        
    # make it ascending
    scale = sorted(scale)
    
    # finally fill me-to-hi
    for k in range(1,nof_bins_hi2me):
        scale.append(me + float(k)*step)
        
    scale.append(hi)
    
    return scale


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
    WT = struct.unpack('f', phsf.read(4))[0]
    ZLAST = struct.unpack('f', phsf.read(4))[0]
    
    W = 1.0 - (U*U + V*V)
    if W < 0.0:
        W = 0.0
    W = math.sqrt(W)

    if (WT < 0.0): # weight sign is in fact Z directional cosine sign
        WT = -WT
        W  = -W

    return (LATCH, E, X, Y, U, V, W, WT, ZLAST)
    
    
def load_events(filename, nof_events = -1):
    """
    load all photon events
    """
    
    with open(filename, "rb") as phsf:
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
    
        NINCP = struct.unpack('f', phsf.read(4))[0]
        print(NINCP) # nof original incident particles    
        print("================ End of PhSF header ======================")
        
        if nof_events < 0:
            nof_events = NPPHSP
        if nof_events > NPPHSP:
            nof_events = NPPHSP

        if mode == b"MODE2":
            dummy = phsf.read(7) # skip fortran garbage, for direct access of record with size 32bytes
        else: # MODE0 file
            dummy = phsf.read(3) # skip fortran garbage, for direct access of record with size 28bytes

        bit29 = 0b0100000000000000000000000000000
        bit30 = 0b1000000000000000000000000000000
        
        nof_photons   = 0
        nof_electrons = 0
        nof_positrons = 0
        
        events = []
        for k in range (0, nof_events):
            (LATCH, E, X, Y, U, V, W, WT, ZLAST) = read_record_long(phsf)
            
            IQ = 0
            if (LATCH & bit30) != 0:
                IQ = -1
                nof_electrons += 1
            elif (LATCH & bit29) != 0:
                IQ = 1
                nof_positrons += 1
            
            if IQ == 0:
                nof_photons += 1
                if E < 0.0:
                    E = -E
                    
                event = [WT, E, X, Y, ZLAST, U, V, W]
                events.append(event)
                
        return (events, nof_photons, nof_electrons, nof_positrons)
        
def write_all_events(filename , events):
    """
    """
    with open(filename, "wt") as f:
        for e in events:
            f.write('{0:15.4e} {1:15.4e} {2:15.4e} {3:15.4e} {4:15.4e} {5:15.4e} {6:15.4e} {7:15.4e}\n'.format(e[0], e[1], e[2], e[3], e[4], e[5], e[6], e[7]))
            
    return None

scale = make_scale(5)
#print(len(scale))
#print(scale)

(events, nof_photons, nof_electrons, nof_positrons) = load_events("../C25x.egsphsp1")

print(len(events), nof_photons, nof_electrons, nof_positrons)

#write_all_events("QQQ", events)

he = H1Dn(scale)

for e in events:
    WT = e[0]
    E  = e[1]
    he.fill(E, WT)

print(he.nof_events(), he.integral())

print(he.underflow())   
print(he.overflow())

X = []
Y = []
W = []

scale = he.x()
n     = len(scale)
norm  = 1.0/he.integral()

sum = 0.0

for k in range (-1, he.size()+1):
    x = 0.0
    w = (he.lo() - x)
    if k == he.size():
        w = (scale[-1]-scale[-2])
        x = he.hi()
    elif k >= 0:
        w = (scale[k+1] - scale[k])
        x = scale[k]
        
    d = he[k]     # data from bin with index k
    y = d[0] / w  # first part of bin is collected weights
    y = y * norm
    X.append(x)
    Y.append(y)
    W.append(w)
    sum += y*w

print(sum)

p1 = plt.bar(X, Y, W, color='r')

plt.xlabel('Energy(MeV)')
plt.ylabel('PDF of the photons')
plt.title('Energy distribution')
    
plt.grid(True);
plt.tick_params(axis='x', direction='out')
plt.tick_params(axis='y', direction='out')

plt.show()
