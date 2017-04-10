#!/usr/bin/env python

import numpy as np
import itertools

from grid import Grid
import tools


# A class that carries the logic of evaluating the energy, force and torque 
# of a pair of rigid molecules. The coordinates of each molecule are given
# in the form of Xcom and q, with Xcom being the Cartesian coordinates of the 
# center of mass, q being the quaternion representation of its orientation 
# wrt a reference pose. The class evaluates the EFTs for a pair of such 
# coordinates by 
#   1. Apply translational and rotational operations to the pair to align the 
#      COM of the first molecule with the origin and its orientation to the 
#      reference pose. 
#   2. Convert the modified Xcom and q of the second molecule into spherical 
#      coordinates.
#   3. Use the resulted six-dimensional coordinate to query a six-dimensional 
#      grid that stores precomputed EFTs.
#   4. Unapply rotation in step 1 to obtain correctly oriented forces and torques
class EFT_calculator:
    def __init__(self, order=2):
        self.mol = Water()
        self.grid = Grid()
        self.order = order # order of the interpolant, 1 for linear

    # Setup the grid structure. If provided with a data file, load it
    def setup(self, filename=None):
        if not filename:
            self.grid.setup()
        else:
            self.grid.load(filename) 

    # Evaluate the Xcom and q for a pair of mols by querying the grid
    def eval(self, coors0, coors1):
        Xcom0 = self.mol.getCOM(coors0)
        Xcom1 = self.mol.getCOM(coors1)
        R0 = self.mol.getR(coors0)
        q0 = tools.R2q(R0)
        R1 = self.mol.getR(coors1)
        q1 = tools.R2q(R1)
        # move COM of mol0 to origin
        X = Xcom1 - Xcom0
        # reorient to align mol0 with refCoor
        X = np.dot(X, R0)
        q = tools.qdiv(q1, q0)
        # Use mirror symmetry of mol0 to move mol1 such that its COM has positive y and z values
        reflections = []
        qsub = q[1:]
        for i in self.mol.refl_axes:
            if X[i] < 0:
                X[i] = -X[i]
                # the following operation on q is equivalent to changing R to MRM 
                # i.e., the probe mol is reflected twice, once in the reference frame,
                # once in the molecular frame.
                qsub[i] = -qsub[i]
                qsub[:] = -qsub
                reflections.append(i)
        # Use mirror symmetry of mol1 to orient it such that it has positive q[0] and q[1] values
        if q[0] < 0:
            q = -q
        if q[1] < 0:
            q[0], q[1], q[2], q[3] = -q[1], q[0], q[3], -q[2]
        # convert X, q to polar coordinates
        r, phi, theta = tools.xyz2spherical(X)
        ophi1, ophi2, otheta = tools.q2spherical(q)
        coor = [r, phi, theta, ophi1, ophi2, otheta]
        # use the grid to obtain results
        eft = self.grid.interpolate(coor, self.order)
        ener = eft[0]
        force = eft[1:4]
        torque = eft[4:7]
        # Reverse the operations for mol0 mirror symmetry back
        for i in reflections:
            force[i] = -force[i]
            torque[i] = -torque[i]
            torque[:] = -torque
        # Reverse the reorientation applied to align mol0 with refCoor
        force[:] = np.dot(force, R0.T)
        torque[:] = np.dot(torque, R0.T)
        return eft


# A class that holds information related to the atomic structure of a water
# molecule. It also includes several methods that carries out operations 
# related to the atomic coordinates.
class Water:
    def __init__(self):
        self.mass = np.array([15.99900, 1.00800, 1.00800])
        refCoor = np.array([ [-0.06556939,   0.00000000,    0.00000000],
                             [0.52035943,    0.76114632,    0.00000000],
                             [0.52035943,   -0.76114632,    0.00000000] ])
        # The following code ensures that refCoor has COM at origin and orientation
        # aligned with the getR() method
        refCoor = refCoor - self.getCOM(refCoor)
        R = self.getR(refCoor) 
        refCoor = np.dot(refCoor, R)
        self.refCoor = refCoor
        # refl_axes is a list of indices of axes. Reflection along each of these
        # axes corresponds to a mirror symmetry of the molecule
        self.refl_axes = [1, 2]

    # Calculate the rotation matrix R that relates coors to self.refCoor
    # vec_in_reference_frame = R \dot vec_in_body_frame 
    # R(coors) \dot refCoor = coors - COM(coors)
    # This function defines the orientation of self.refCoor
    # Need to be consistent with self.refl_axes
    def getR(self, coors):
        coors = np.copy(coors)
        offset = coors[0]
        coors -= offset
        xvec = coors[1] + coors[2]
        zvec = np.cross(coors[1], coors[2])
        yvec = np.cross(zvec, xvec)
        xvec /= np.linalg.norm(xvec)
        yvec /= np.linalg.norm(yvec)
        zvec /= np.linalg.norm(zvec)
        R = np.array([xvec, yvec, zvec]).T 
        return R

    # Calculate the center of mass
    def getCOM(self, coors):
        return np.dot(self.mass, coors) / self.mass.sum()


