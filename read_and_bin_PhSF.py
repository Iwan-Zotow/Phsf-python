#!/usr/bin/env python3

import sys

import matplotlib.pyplot as plt

import beam_loader
import text_loader

from H1Dn import H1Dn

def make_scale(nof_bins_hi2me):
    """
    FFiven number of bins, make energy scale
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
    for k in range(1, nof_bins_hi2me):
        scale.append(me + float(k)*step)

    scale.append(hi)

    return scale


def write_all_events(filename , events):
    """
    Given filename, write all events
    """
    with open(filename, "wt") as f:
        for e in events:
            f.write('{0:15.4e} {1:15.4e} {2:15.4e} {3:15.4e} {4:15.4e} {5:15.4e} {6:15.4e} {7:15.4e}\n'.format(e[0], e[1], e[2], e[3], e[4], e[5], e[6], e[7]))

    return None


def main(pshf_name):

    scale = make_scale(5)
    #print(scale)

    (events, nof_photons, nof_electrons, nof_positrons) = text_loader.load_events(pshf_name)
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

    for k in range(-1, he.size()+1):
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

if __name__ == "__main__":

    #fname = "/home/beamuser/C25.egsphsp1"
    fname = "../PHSF"

    main(fname)

    sys.exit(0)