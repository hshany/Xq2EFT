#!/usr/bin/env python2

import numpy as np
from time import time
import heapq
from matplotlib import pyplot as plt

from eft_calculator import EFT_calculator, Water
import tools


def load_coordinates(name):
    lines = open('random/'+name).readlines()[-7:-1]
    coors = [[float(item) for item in line.split()[2:5]] for line in lines]
    return np.array(coors)

def test_random_set():
    ener = []
    force = []
    torque = []
    t0 = time()
    for i in range(2, 2000):
        # load atomic coor 
        name = 'test%04d.inp' % i
        coors = load_coordinates(name)
        # convert atomic coor to r, phi, theta... 
        X0, q0 = calculator.mol.atomic2Xq(coors[:3])
        X1, q1 = calculator.mol.atomic2Xq(coors[3:])
        # evaluate with calculator
        eft = calculator.eval(X0, q0, X1, q1)
        ener.append(eft[0])
        force.append(eft[1:4])
        torque.append(eft[4:7])
    t1 = time()
    print 'took %.1f s to evaluate the random set' % (t1 - t0)
    return ener, force, torque

if __name__ == '__main__':
    order = 2
    calculator = EFT_calculator(order)
    t0 = time()
    calculator.setup('grid.dat')
    t1 = time()
    print 'took %.1f s to fill the grid' % (t1 - t0)
    ener, force, torque = test_random_set()

