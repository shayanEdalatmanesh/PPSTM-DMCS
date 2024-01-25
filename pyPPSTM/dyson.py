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
def getDysonGeom(geom='input_plot.xyz', lvs = None):
    Ratin, tmp1, tmp2 = Bu.loadAtoms(geom);
    del tmp1, tmp2;
    print(Ratin)
    return Ratin

def readDysonAll(name= 'dyson.dat', geom='input_plot.xyz', lvs = None):

    Ratin = getDysonGeom(geom)
    Ratin = RS.for_PBC(Ratin, lvs)

    print(" loading the LCAO coefficients")
    filein = open(name)
    pre_eig = filein.readline().split()
    filein.close()
   # print ("PRE EIG = ", pre_eig)
#    eig = np.loadtxt(name, skiprows=1, usecols=(0,))
    num_at_ = int(pre_eig[0])
    n_bands = 1
    NoOrb = n_bands
    eig = []
    eig.append(pre_eig[2])
    print(type(eig))

#    if eig == int(0):
#        eigs = [i for i in range(1)]
#        print(eigs)

    print("n orbs: ", n_bands)
    print("eigen en read = ", eig)
    from pyPPSTM.ReadSTM import n_max_, n_min_, Ynum_
    coef = np.zeros((n_bands,num_at_,Ynum_))
    print("coef's shape", coef.shape)
    if (num_at_ > 1):
        coef[:, :, 2] = np.loadtxt(name, skiprows=1, usecols=(0,))
       # coef[:, :, 2] = np.loadtxt(name, skiprows=1,usecols=tuple(range(1, num_at_*2+1, 2)) )
# this line depends on the format of the file produced by Diego
        print(coef)

    #coeff = coef.flatten()
   # coeff.reshape((n_bands, num_at_ * Ynum_))
   # print("coeff's shape", coeff.shape)

    te = np.reshape(coef, (n_bands, num_at_ * Ynum_))
   # print("DEBUG: te shape = ", te.shape)

 #   eig = np.ndarray(1)
    #print("eig = ", eig)
#    eig = np.repeat(1, num_at_) 
  #  print("eig = ", eig)
  #  print(type(eig))
  #  print("len(eig) = ", len(eig))
  #  assert (len(eig)==n_bands), "number of bands wrongly specified"

    print("eigen-energies read")
    #print("DEBUG: type of EIG  = ", type(eig))
    return eig.copy(), te.copy(), Ratin.copy()

