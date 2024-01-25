#!/usr/bin/python
import pyPPSTM.basUtils as Bu
import os
import numpy as np
import math
import matplotlib as mpl
# mpl.use('Agg') # Force matplotlib to not use any Xwindows backend. ## !!! important for working on clusters !!!!
import matplotlib.pyplot as plt
from pylab import genfromtxt
import matplotlib.pylab as pl
import matplotlib.colors as colors
import matplotlib.cm as cmx
from pyPPSTM import ReadSTM as RS
from collections import Counter
# ===========================================================================================================
def getSislGeom(geom='input_plot.xyz', lvs = None):
    Ratin, tmp1, tmp2 = Bu.loadAtoms(geom);
    del tmp1, tmp2;
    print(" Number of atoms: ", len(Ratin[0]))
    global num_at_
    num_at_ = len(Ratin[0])
    print(Ratin)
    return Ratin

eigEnz = []
def getSislEnz(fEnz='enz.dat'):
    with open(fEnz, 'r') as input_file:
        # Read the number of atoms
        # numOrbz = int(input_file.readline())
        # Read the atomic coordinates and negate the x-coordinate
        for line in input_file:
            parts = line.split()
            eigEnz.append(parts)
            eigEnz2 = np.reshape(eigEnz, -1)
    return eigEnz2

def getSislCoefs(fn='coefsAsMatrix.npy'):
    coefs = np.load(fn)
    return coefs

#    with open(fn, 'r') as inputF:
#        # Parse the coordinates and atomic elems for each atom
#        data = []
#        for line in data:
#        fields = line.split()
#        elems.append(fields[0])
#        coordinates.append([float(fields[1]), float(fields[2]), float(fields[3])])


def readSislAll(geom='input_plot.xyz', fEnz= 'enz.dat', fn='coefs.dat.npy', lvs = None):
    print(" reading the geometry from the xyz file...")
    Ratin = getSislGeom(geom)
    Ratin = RS.for_PBC(Ratin, lvs)
    print(" geometry loaded successfully!")
    
    eigEnz = getSislEnz(fEnz)
    print(" Energy eigenvalues loaded successfully!")

    print(" loading the LCAO coefficients...")
    coef = getSislCoefs(fn)
    print("coef's type = ", type(coef))
    print("coef's shape = ", coef.shape)

    from pyPPSTM.ReadSTM import n_max_, n_min_, Ynum_
    print("n_max_", n_max_)
    print("n_min_", n_min_)
    print("num_at_", num_at_)
    print("Ynum_", Ynum_)
    n_bands = len(eigEnz)

    coeff = np.reshape(coef, (n_bands, num_at_ * Ynum_))
    print("coeff's shape = ", coeff.shape)

    print("eigen-energies read")
    #print("DEBUG: type of EIG  = ", type(eig))
    return eigEnz.copy(), coeff.copy(), Ratin.copy()
'''

   # print ("PRE EIG = ", pre_eig)
   # eig = np.loadtxt(name, skiprows=1, usecols=(0,))
    num_at_ = int(pre_eig[0])
    n_bands = 1
    NoOrb = n_bands
    eig = []
    eig.append(pre_eig[2])
    print(type(eig))

    nc = Ratin[0].count('C')
    nh = Ratin[0].count('H')
    nn = Ratin[0].count('N')
    nb = Ratin[0].count('B')

    if szp == True:
        noCoefz = 5*nb + 5*nn + 5*nc + 1*nh
    elif dzp == True:
        noCoefz = 13*nb + 13*nn + 13*nc + 5*nh
    else:
        print("the basis set is not defined!")

    print ("number of coefz per energy", noCoefz)

    print("n orbs / bands: ", n_bands)
    print("eigen en read = ", eig)

    j = 0
    pointerz = []
    for i in np.arange(len(Ratin[0])):
        pointerz.append(j)
        if Ratin[0][i] == 'B':
            j += 4
        if Ratin[0][i] == 'C':
            j += 4
        elif Ratin[0][i] == 'H':
            j += 1
        elif Ratin[0][i] == 'N':  # This should be debugged, we need to make a list of species, maybe it exists somewhere.
            j += 4
        else:
            print("Species more than B, C, N and H detected, quitting now...")
        print("j = ", j)

    print("Element pointers = ", pointerz)
    # the next line removes the possible duplicates, coz of random lattice parameters or comments in the xyz file.
    temp = Counter(pointerz)
    res = [*temp]
    print("after removal of possible duplicates = ", res)

    coefs = np.zeros((n_bands, num_at_, Ynum_))
    print("Ynum_ = ", Ynum_)

    cc = np.loadtxt('lumoLCAO.dat', skiprows=1) #, usecols=(0,))

    for i in np.arange(n_bands):
        for j in np.arange(len(Ratin[0])):
            k = res[j]
            # print("species: ", spec[j])
            # print("j = ", j)
            # =print("k = ", k)
            if Ratin[0][j] == 'B' or Ratin[0][j] == 'C' or Ratin[0][j] == 'N':
                coefs[i, j, 0] = cc[k + 0] # S-orbital supposedly
                coefs[i, j, 3] = cc[k + 3]  # Px-orbital supposedly
                coefs[i, j, 1] = cc[k + 1]  # Py-orbital supposedly
                coefs[i, j, 2] = cc[k + 2] # Pz-orbital supposedly
            elif Ratin[0][j] == 'H':
                coefs[i, j, 0] = cc[k + 0]  # S-orbital supposedly
            else:
                print("Species other than carbon, nitrogen and hydrogen detected. Can't handle it, sorry.")


    if (num_at_ > 1):
        coef[:, :, 2] = np.loadtxt(name, skiprows=1, usecols=(0,))
       # coef[:, :, 2] = np.loadtxt(name, skiprows=1,usecols=tuple(range(1, num_at_*2+1, 2)) )
# this line depends on the format of the file produced by Diego
        print(coef)

   # coeff = coef.flatten()
   # coeff.reshape((n_bands, num_at_ * Ynum_))
   # print("coeff's shape", coeff.shape)
'''
   # print("DEBUG: te shape = ", te.shape)
 #   eig = np.ndarray(1)
    #print("eig = ", eig)
#    eig = np.repeat(1, num_at_) 
  #  print("eig = ", eig)
  #  print(type(eig))
  #  print("len(eig) = ", len(eig))
  #  assert (len(eig)==n_bands), "number of bands wrongly specified"
