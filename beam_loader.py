#!/usr/bin/env python3

import math
import struct


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

                event = (WT, E, X, Y, ZLAST, U, V, W)
                events.append(event)

        return (events, nof_photons, nof_electrons, nof_positrons)
