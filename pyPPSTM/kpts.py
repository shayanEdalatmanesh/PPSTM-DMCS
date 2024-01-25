#!/env/python
# Use this script to plot multiple ".dat" files in one figure.
import os
from . import ReadSTM as RS
from matplotlib import pyplot as plt
import numpy as np
from pylab import genfromtxt
#==================================================================
# T
#==================================================================
# 1. set the k (kpoint number x from the list, starting from 0)
# 2. readfireball the k-file no x
# 3. read the k number x from the *.kpts file
# 3.1. load the geometry
# 4. do the dot product k.r in a loop over the atoms -> list 
# 5. calculate the cos(k.r) and sin(k.r), one pair for each atom
# 6. for n atoms we have 2*n coeffs (n real and n imaginaries)
# 7. 
#==================================================================
#lat vec needs to be used, k = n pi / a? Should it be multiplied or pi / a = 1? latVec=lvs
def kGetterFire(nk=1, kName= 'temp.kpts'):
    kfil = np.loadtxt(kName, skiprows=1, usecols=(0,1,2))
    print ("The selected k-point is: ", kfil[nk-1]) 
    kCoord = kfil[nk-1]
    print ("#DEBUG#", kCoord)
    return kCoord


def atomiCoords(geo='input_plot.xyz'):
    geo00 = np.loadtxt(geo, skiprows=2, usecols=(1,2,3))
    print ("geometry file read: ", geo00)
    return geo00


def dotter(geo00, kCoord):
    dotProds=[]
    for i in np.arange(len(geo00)):
        prod = np.dot(geo00[i], kCoord)
        dotProds.append(prod)
    return dotProds


def phaser(dotProds):
    cos0 = []
    sin0 = []
    for j in np.arange(len(dotProds)):
        c = np.cos(dotProds[j])
        cos0.append(c)
        s = np.sin(dotProds[j])
        sin0.append(s)
    return cos0, sin0


def kCoeFireball(name = 'phi_' , geom='input_plot.xyz', fermi=None, orbs = 'sp', pbc=(1,1), cut_min=-15.0, cut_max=5.0, cut_at=-1, lvs = None, lower_atoms=[], lower_coefs=[]):
    #RS.initial_check(orbs=orbs, pbc=pbc, imaginary=imaginary, cut_min=cut_min, cut_max=cut_max, cut_at=cut_at, lower_atoms=lower_atoms, lower_coefs=lower_coefs)
    # obtaining the geometry :
    Ratin, num_at_ = RS.get_FIREBALL_geom(geom=geom, lvs=lvs)
    # getting eigen-energies
    filein = open(name+'s.dat' )
    pre_eig = filein.readline().split()
    filein.close()
    #print "DEBUG: num_at_  and pre_eig[0] " , num_at_ ,pre_eig[0], type(num_at_), type(pre_eig[0][0]), "int(pre_eig[0])", int(pre_eig[0]) , "(num_at_<=int(pre_eig[0][0]))", (num_at_<=int(pre_eig[0]))
    assert (num_at_<=int(pre_eig[0])),"coefficients for lower amount of atoms, that atoms in geometry file";
    n_bands= int(pre_eig[1]);
    eig = np.loadtxt(name+'s.dat',skiprows=1, usecols=(0,))
    assert (len(eig)==n_bands), "number of bands wrongly specified"
    eig = RS.to_fermi(eig, fermi, orig_fermi=float(pre_eig[2]))
    del pre_eig;
    eig2 = RS.cut_eigenenergies(eig)
    print("Eigenenergies in the selected range are: ", eig2)
    print("Eigen shape", np.shape(eig2))

    from pyPPSTM.ReadSTM import n_max_, n_min_, Ynum_
    print("n_max_", n_max_)
    print("n_min_", n_min_)
    print("num_at_", num_at_)
    print("Ynum_", Ynum_)

    print(" loading the LCAO coefficients")
    coefRe = np.zeros((n_bands,num_at_,Ynum_))
    coefIm = np.zeros((n_bands,num_at_,Ynum_))
    print(" loading the LCAO coefficients - imaginary active!")
    if (num_at_ > 1):
        coefRe[:,:,0] = np.loadtxt(name+'s.dat',skiprows=1,usecols=tuple(range(1, num_at_*2+1, 2)) )
        coefRe[:,:,1] = np.loadtxt(name+'py.dat',skiprows=1,usecols=tuple(range(1, num_at_*2+1, 2)) )
        coefRe[:,:,2] = np.loadtxt(name+'pz.dat',skiprows=1,usecols=tuple(range(1, num_at_*2+1, 2)) )
        coefRe[:,:,3] = np.loadtxt(name+'px.dat',skiprows=1,usecols=tuple(range(1, num_at_*2+1, 2)) )
        # imaginary coefs
        coefIm[:,:,0] = np.loadtxt(name+'s.dat',skiprows=1,usecols=tuple(range(2, num_at_*2+2, 2)) )
        coefIm[:,:,1] = np.loadtxt(name+'py.dat',skiprows=1,usecols=tuple(range(2, num_at_*2+2, 2)) )
        coefIm[:,:,2] = np.loadtxt(name+'pz.dat',skiprows=1,usecols=tuple(range(2, num_at_*2+2, 2)) )
        coefIm[:,:,3] = np.loadtxt(name+'px.dat',skiprows=1,usecols=tuple(range(2, num_at_*2+2, 2)) )
        if (orbs =='spd'):
            coefRe[:,:,4] = np.loadtxt(name+'dxy.dat',skiprows=1,usecols=tuple(range(1, num_at_*2+1, 2)) )
            coefRe[:,:,5] = np.loadtxt(name+'dyz.dat',skiprows=1,usecols=tuple(range(1, num_at_*2+1, 2)) )
            coefRe[:,:,6] = np.loadtxt(name+'dz2.dat',skiprows=1,usecols=tuple(range(1, num_at_*2+1, 2)) )
            coefRe[:,:,7] = np.loadtxt(name+'dxz.dat',skiprows=1,usecols=tuple(range(1, num_at_*2+1, 2)) )
            coefRe[:,:,8] = np.loadtxt(name+'dx2y2.dat',skiprows=1,usecols=tuple(range(1, num_at_*2+1, 2)) )
        # imaginary coefs
            coefIm[:,:,4] = np.loadtxt(name+'dxy.dat',skiprows=1,usecols=tuple(range(2, num_at_*2+2, 2)) )
            coefIm[:,:,5] = np.loadtxt(name+'dyz.dat',skiprows=1,usecols=tuple(range(2, num_at_*2+2, 2)) )
            coefIm[:,:,6] = np.loadtxt(name+'dz2.dat',skiprows=1,usecols=tuple(range(2, num_at_*2+2, 2)) )
            coefIm[:,:,7] = np.loadtxt(name+'dxz.dat',skiprows=1,usecols=tuple(range(2, num_at_*2+2, 2)) )
            coefIm[:,:,8] = np.loadtxt(name+'dx2y2.dat',skiprows=1,usecols=tuple(range(2, num_at_*2+2, 2)) )
    else:
        coefRe[:,0,0] = np.loadtxt(name+'s.dat',skiprows=1,usecols=tuple(range(1, num_at_*2+1, 2)) )
        coefRe[:,0,1] = np.loadtxt(name+'py.dat',skiprows=1,usecols=tuple(range(1, num_at_*2+1, 2)) )
        coefRe[:,0,2] = np.loadtxt(name+'pz.dat',skiprows=1,usecols=tuple(range(1, num_at_*2+1, 2)) )
        coefRe[:,0,3] = np.loadtxt(name+'px.dat',skiprows=1,usecols=tuple(range(1, num_at_*2+1, 2)) )
        #
        coefIm[:,0,0] = np.loadtxt(name+'s.dat',skiprows=1,usecols=tuple(range(2, num_at_*2+2, 2)) )
        coefIm[:,0,1] = np.loadtxt(name+'py.dat',skiprows=1,usecols=tuple(range(2, num_at_*2+2, 2)) )
        coefIm[:,0,2] = np.loadtxt(name+'pz.dat',skiprows=1,usecols=tuple(range(2, num_at_*2+2, 2)) )
        coefIm[:,0,3] = np.loadtxt(name+'px.dat',skiprows=1,usecols=tuple(range(2, num_at_*2+2, 2)) )

        if (orbs =='spd'):
            coefRe[:,0,4] = np.loadtxt(name+'dxy.dat',skiprows=1,usecols=tuple(range(1, num_at_*2+1, 2)) )
            coefRe[:,0,5] = np.loadtxt(name+'dyz.dat',skiprows=1,usecols=tuple(range(1, num_at_*2+1, 2)) )
            coefRe[:,0,6] = np.loadtxt(name+'dz2.dat',skiprows=1,usecols=tuple(range(1, num_at_*2+1, 2)) )
            coefRe[:,0,7] = np.loadtxt(name+'dxz.dat',skiprows=1,usecols=tuple(range(1, num_at_*2+1, 2)) )
            coefRe[:,0,8] = np.loadtxt(name+'dx2y2.dat',skiprows=1,usecols=tuple(range(1, num_at_*2+1, 2)) )
            #
            coefIm[:,0,4] = np.loadtxt(name+'dxy.dat',skiprows=1,usecols=tuple(range(2, num_at_*2+2, 2)) )
            coefIm[:,0,5] = np.loadtxt(name+'dyz.dat',skiprows=1,usecols=tuple(range(2, num_at_*2+2, 2)) )
            coefIm[:,0,6] = np.loadtxt(name+'dz2.dat',skiprows=1,usecols=tuple(range(2, num_at_*2+2, 2)) )
            coefIm[:,0,7] = np.loadtxt(name+'dxz.dat',skiprows=1,usecols=tuple(range(2, num_at_*2+2, 2)) )
            coefIm[:,0,8] = np.loadtxt(name+'dx2y2.dat',skiprows=1,usecols=tuple(range(2, num_at_*2+2, 2)) )
    
    print("coefRe =", coefRe)
    print("shape re =", coefRe.shape)
    print("coefIm =", coefIm)
    print("shape im =", coefIm.shape)

    return eig2.copy(), coefRe.copy(), coefIm.copy(), Ratin.copy();


def kHandlerFire(nk=1, kName= 'temp.kpts', imaginary=True, geo ='input_plot.xyz', name = 'phi_' , geom='input_plot.xyz', fermi=None, orbs = 'sp', pbc=(1,1), cut_min=-15.0, cut_max=5.0, cut_at=-1, lvs = None, lower_atoms=[], lower_coefs=[]):

    assert imaginary == True, f"imaginary = True was expected, got {imaginary}"
    print("Attention: you've chosen to read and include the imaginary parts of the coefficients in the calculations. If you'd like to do only Gamma-point calculations, set imaginary = False.")
    print("The selected k-point is k-p no: ", nk)
    print("reading the *.kpts file...")
    
    #nk should be used to automatically select the phik file with that number and open it. For now we're gonna do it manually in the input file.
    # eigEn, coefs, Ratin = RS.read_FIREBALL_all(name = files_path + 'phik_0001_', geom=files_path+geometry_file, lvs = cell, fermi=fermi, orbs = sample_orbs, pbc=pbc, cut_min=cut_min, cut_max=cut_max,
    
    # Reading the *.kpts file, finding the k-pt we're interested in.
    kCoord = kGetterFire(nk, kName)
    
    # Get the geometries (gotta make it work with the readfireball function later)
    geo00 = atomiCoords(geo)    
    
    # calculating exp(ik.r) - the phase that a particle picks up through the band
    print("Calculating the ik.r products (Bloch theorem) ...")

    dotProds = dotter(geo00, kCoord)
    #print ("Dot products = ", dotProds)
    print ("#DEBUG# Number of the prods = ", len(dotProds))
    #print ("dotProds is a list")
    #print ("Calculating the phases...")
    cos0, sin0 = phaser(dotProds)
    #print("cos = ", cos0)
    #print("sin = ", sin0)
    #print(len(cos0))
    #print(len(sin0))
    cos1 = np.array(cos0)
    sin1 = np.array(sin0)
    #print(np.shape(sin1))
    eig4, coefRe4, coefIm4, Ratin4 = kCoeFireball(name, geom, fermi, orbs, pbc, cut_min, cut_max, cut_at, lvs , lower_atoms, lower_coefs)

    #print("#DEBUG#",cos1[0]*coefRe4[0,0])
    cosRe = coefRe4*cos1[:,None]
    cosIm = coefIm4*cos1[:,None]
    sinRe = coefRe4*sin1[:,None]
    sinIm = coefIm4*sin1[:,None]

    print(cosRe.shape, cosIm.shape, sinRe.shape, sinIm.shape)
    #print("cosRe ", cosRe)
    # removing states (molecular orbitals) that are not wanted
    from pyPPSTM.ReadSTM import n_max_, n_min_
    cosRe2 = cosRe[n_min_:n_max_,:,:]
    cosIm2 = cosIm[n_min_:n_max_,:,:]
    sinRe2 = sinRe[n_min_:n_max_,:,:]
    sinIm2 = sinIm[n_min_:n_max_,:,:]

    print("including the phase in the coefs done...")
    print("#DEBUG# shape of the cosRe matrix: ", np.shape(cosRe2))
    
    # lowering over atoms and applying PBC
    cosRe3 = RS.handle_coef(cosRe2)
    cosIm3 = RS.handle_coef(cosIm2)
    sinRe3 = RS.handle_coef(sinRe2)
    sinIm3 = RS.handle_coef(sinIm2)

    print("All coefficients read.")
    print("EIGz shape: ", eig4.shape)

    return eig4, cosRe3, cosIm3, sinRe3, sinIm3, Ratin4

    # read the fireball coefs and multiply them

#=============================================================================================================================
def kGetterAIMS(kName='KS_eigenvectors.band_1.kpt_1.out'):
    k=[]
    with open(kName) as fh:
        for line in fh:
            if line.startswith("# Complex"):
                print(line)
                line = line.strip().split()
                print(line)
                print(line[15], line[16], line[17])
                k.append(line[15])
                k.append(line[16])
                k.append(line[17])
                #k = np.float(k)
                k = [float(x) for x in k]
    np.reshape(k,(1,3))
    return k


def kHandlerAIMS(kName = 'KS_eigenvectors.band_1.kpt_1.out', imaginary=True, geo ='input_plot.xyz', geom='geometry.in', fermi=None, orbs = 'sp', pbc=(0,0), cut_min=-15.0, cut_max=5.0, cut_at=-1, lower_atoms=[], lower_coefs=[]):

    kCoord = kGetterAIMS(kName)    
    print(kCoord)

    geo00 = atomiCoords(geo)  
    dotProds = dotter(geo00, kCoord)
    print("dotProds", dotProds)
    print(len(dotProds))
    
    cos0, sin0 = phaser(dotProds)
    cos1 = np.array(cos0)
    sin1 = np.array(sin0)
    #print("cos1= ", cos1)

    eig4, coeffs, Ratin4, origLvec = RS.read_AIMS_all(kName, geom, fermi, orbs, pbc, imaginary, cut_min, cut_max, cut_at, lower_atoms, lower_coefs)
    #np.delete(coeffs)
    #Ratin, at_num, origLvec = RS.get_AIMS_geom(geom=geom)

    coefRe4, coefIm4 = RS.read_AIMS_coefs(kName, at_num, imaginary)

    cosRe = coefRe4*cos1[:,None]
    cosIm = coefIm4*cos1[:,None]
    sinRe = coefRe4*sin1[:,None]
    sinIm = coefIm4*sin1[:,None]
    print(cosRe.shape, cosIm.shape, sinRe.shape, sinIm.shape)

    from pyPPSTM.ReadSTM import n_max_, n_min_
    cosRe2 = cosRe[n_min_:n_max_,:,:]
    cosIm2 = cosIm[n_min_:n_max_,:,:]
    sinRe2 = sinRe[n_min_:n_max_,:,:]
    sinIm2 = sinIm[n_min_:n_max_,:,:]

    print("including the phase in the coefs done...")
    print("#DEBUG# shape of the cosRe matrix: ", np.shape(cosRe2))
    
    # lowering over atoms and applying PBC
    cosRe3 = RS.handle_coef(cosRe2)
    cosIm3 = RS.handle_coef(cosIm2)
    sinRe3 = RS.handle_coef(sinRe2)
    sinIm3 = RS.handle_coef(sinIm2)
    print = ("DEBUG shape of coefs", np.shape(cosRe3))

    print("All coefficients read.")
    print("EIGz shape: ", np.shape(eig4))
    print (np.shape(Ratin4))

    return eig4, cosRe3, cosIm3, sinRe3, sinIm3, Ratin4, origLvec


