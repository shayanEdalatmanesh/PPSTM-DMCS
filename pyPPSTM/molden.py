#!/usr/bin/python
# To use on tarkil you should load python using: module load <your favorite python3 version> and uncomment the Agg line below (coming soon)

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
from . import ReadSTM as RS
from collections import Counter

# from RS import num_at_
# ===========================================================================================================
eHartree = float(27.2114)


# ===========================================================================================================

def moldenAtomz(name='molden.input'):
    linez = []
    with open(name) as infile:
        copy = False
        for line in infile:
            if line.strip() == "[Atoms] Angs":
                print("Coordinates in Angs, continuing to read the file...")
                copy = True
                continue
            elif line.strip() == "[Atoms] AU":
                print("Coordinates in Atomic Units instead of Angs, be careful!")
                copy = True
            elif line.strip() == "[GTO]" or line.strip() == "[5D7F]":
                copy = False
                continue
            elif copy:
                print(line)
                linez.append(line.split())
    print("Number of atoms found: ", len(linez))

    elem = []
    xs = []
    ys = []
    zs = []
    q = []
    for line in linez:
        elem.append(line[0])
        xs.append(line[3])
        ys.append(line[4])
        zs.append(line[5])
        q.append(0.0)

    print("Writing the atomic structure as an xyz file...")
    with open(name + '_.xyz', "w") as outfile:
        outfile.write(str(len(linez)) + "\n")
        outfile.write("\n")
        for i in np.arange(len(linez)):
            outfile.write(elem[i] + "   " + str(xs[i]) + "   " + str(ys[i]) + "   " + str(zs[i]) + "\n")
        outfile.close()
    print("Atomic structure written to file. You could check it now.")

    nDim = []
    lvec = []
    # print("")
    print(" Number of atoms: ", len(xs))
    global num_at_
    num_at_ = len(xs)

    return [elem, xs, ys, zs, q], nDim, lvec


# ===========================================================================================================

def getMOLDENgeom(name='molden.input'):
    '''
    Prepares geometry from the MOLDEN file format, no need to put some input_plot.xyz file
    '''
    print(" # ============ define atoms ============")
    atoms, nDim, tmp = moldenAtomz(name)
    # print "DEBUG: atoms", atoms
    del nDim, tmp  # for the future, if we had PBC or sth.
    atoms = RS.cut_atoms(atoms)

    all = Counter(atoms[0])
    print(all)

    # dict((i, atoms[0].count(i)) for i in atoms[0])
    # print(i)

    Ratin = RS.for_PBC(atoms, lvs=None)
    print("atomic geometry read")
    print("Ratin", Ratin)
    return Ratin, atoms[0]

def stripSpaces(astring):
    temp = astring.split(" ")
    return [element for element in temp if len(element) != 0]

# ===========================================================================================================
# Old procedure, didn't work for a case so I had to change it...
#
#def getMOLDENcoefs(fname, spin="closedShell"):
    #coefs = []
    #goody = []
    #with open(fname) as infile2:
    #    print("DEBUG: file opened for coefs")
    #    get = False
    #    for li in infile2:
    #        if li.startswith("Sym" or " Sym"):
                # occup.append(li.split()) # to store the occupation data. Is it accurate tho?
    #            get = False
    #            continue
    #        elif li.startswith("Ene" or " Ene"):  # Or ' Ene'
    #            get = False
    #            continue
    #        elif li.startswith("Occup" or " Occup" or " Occup=" or " Occup= "):
    #            get = True
    #            continue
    #        if get:
    #             print(li)
    #             coefs.append(li.split())
    #             print(coefs)
#    ===========================================================================================================
def getMOLDENcoefs(fname, spin="closedShell", noCoefz=1, impN=1):
    coLine = 0
    occupat = []
    coStart = [] 
    with open(fname) as infile2:
        for lin in infile2:
            coLine+=1
            if "Occup" in lin:
                occupat.append(lin.split())
                coStart.append(coLine)
                #print ("DEBUG", lin)
                
    #print ("coef first lines read")
    #print ("how many? ", len(coStart))
    #print ("occupat", occupat)
    print("Occupations data collected ")

    occupat2 = []
    for i in np.arange(len(occupat)):
        occupat2.append(occupat[i][1])

    print ("number of coefz per energy :", noCoefz)
    coStart2 = []
    for i in np.arange(len(coStart)):
        for j in np.arange(noCoefz):
            gline = coStart[i]+j
            coStart2.append(gline)

#    print ("coStart2", coStart2)
    coefs = []
    goody = []
    with open(fname) as infile3:
        for l, line in enumerate(infile3):
            if l in coStart2:
                coefs.append(line.strip())
            #for i in np.arange(len(coStart)):
                #for j in np.arange(noCoefz):
                #    if line [coStart[i]:coStart[i+j]]
                #print (x)
                #coefs.append(x)

#        for ls in infile3:
#            for c in np.arange(len(coStart)):
#                if ls == coStart[c]:
#                    coefs.append(ls.split())
#                    print coefs

            #coefs.append(coStart)
            #coefs.append(ls.split())
            #coStart[]
#
    #print(coefs)

    # assert len(coefs) == len(occup), "Reading the coefficients failed... "

    print("LEN COEFS = ", len(coefs))
    #print("coefs : ", coefs)
    print(np.shape(coefs))

    if spin == "closedShell":
        lenny = int(np.sqrt(len(coefs)))
    elif spin == "alpha" or spin == "beta":
        lenny = int(np.sqrt(len(coefs)/2))

    print("lenny", lenny)

    for i in np.arange(len(coefs)):
        print("#DEBUG# strip = ", stripSpaces(coefs[i]))
        g = stripSpaces(coefs[i])
        print (g[1]) 
        goody.append(g[1])  # removing the counter from the original file


    print("Len energies = ", len(goody))
    print(goody)
    print(np.shape(goody))
#    print(type(len(goody)))
    up = []
    dn = []
    if spin == "alpha":
        for i in np.arange(int(len(goody)/2)):
            up.append(goody[i])
    elif spin == "beta":
        for i in np.arange(int(len(goody)/2), int(len(goody))):
            dn.append(goody[i])
    else:
        print("No spin detected, closed-shell calculations then...")

    print("len up & dn = ", len(dn))

    if spin == "alpha":
        coco = np.reshape(up, (impN, noCoefz/2))
    elif spin == "beta":
        coco = np.reshape(dn, (impN, noCoefz/2))
    elif spin == "closedShell":
#        coco = np.reshape(goody, (lenny, lenny))
        coco = np.reshape(goody, (impN, noCoefz)) # only for the special case of Libor's butterfly calx.
    else:
        raise Exception

    #print("line by line", coefs)
    # print(coefs[13])
    # print(coefs[13][0])
    # print(coefs[13][0]+coefs[13][1]+coefs[13][2])
    return coco  # ?! pbc and stuff


# ===========================================================================================================

def getEigEnMolden(name='molden.input', spin="closedShell"):
    eigEn = []
    with open(name) as infile:
        for en in infile:
            if "Ene" in en:
#            if en.startswith("Ene" or " Ene"):
                eigEn.append(en.split())

    print("Eigen Energies read.")
    #print("#DEBUG# eigen = ", eigEn)
    eigEn2 = []
    for i in np.arange(len(eigEn)):
        eigEn2.append(eigEn[i][1])
        # eigEns need to be converted from Hartree to eV:
    print("eigEnz in Hartree = ", len(eigEn2))
    eigEn2 = [float(i) for i in eigEn2]
    # eigEn2 = np.sort(eigEn2, axis=None)
    eig = [x * eHartree for x in eigEn2]
    print("eigEnz in eV = ", eig)
    global n_bands
    n_bands = len(eig)
    print("number of molecular orbitals = ", n_bands)
    eigEn3 = []

    if spin == "alpha":
        for i in np.arange(int(n_bands/ 2)):
            eigEn3.append(eig[i])
    elif spin == "beta":
        for i in np.arange(int((n_bands / 2)), n_bands):
            eigEn3.append(eig[i])
    elif spin == "closedShell":
        eigEn3 = eig
    else:
        raise Exception
    print("final eigens selected: ", eigEn3)

    return eigEn3


# ===========================================================================================================
def getSpin(name='molden.input'):
    occup = []
    spin = []
    with open(name) as infile:
        for en in infile:
            if en.startswith(" Occup" or "Occup"):
                occup.append(en.split())
#            if en.startswith(" Spin"):
#                spin.append(en.split())

    print("Occupancy data read.")

    print("occup = ", occup)
    print("occup size = ", len(occup))

    print("spin = ", spin)
    print("spin size = ", len(spin))
 #   return occup, spin
    return occup

# ===========================================================================================================

def readMoldenAll(name='molden.input', fermi=None, orbs='sp', pbc=None, imaginary=False, cut_min=-50.0, cut_max=5.0,
                  cut_at=-1, lvs=None, lower_atoms=[], lower_coefs=[], spin="closedShell"):
    """
    spin= None or "alpha" or "beta" or "both" for spin-unrestricted calculations
    """
    RS.initial_check(orbs=orbs, pbc=pbc, imaginary=imaginary, cut_min=cut_min, cut_max=cut_max, cut_at=cut_at,
                     lower_atoms=lower_atoms, lower_coefs=lower_coefs)

    #occup, spin = getSpin(name)

    Ratin, spec = getMOLDENgeom(name)
    print("elems", spec)
    nc = spec.count('C')
    nh = spec.count('H')
    nn = spec.count('N')
    noCoefz = 14*nn + 14*nc + 5*nh
    print ("number of coefz per energy", noCoefz)
    el = enumerate(spec)
    print("el = ", el)
    j = 0
    pointerz = []
    for i in np.arange(len(spec)):
        pointerz.append(j)
        if spec[i] == 'C':
            j += 14
        elif spec[i] == 'H':
            j += 5
        elif spec[i] == 'N':  # This should be debugged, we need to make a list of species, maybe it exists somewhere.
            j += 14
        else:
            print("Species more than C, N and H detected, quitting now...")
        print("j = ", j)

    print("Element pointers = ", pointerz)
    Ratin = Ratin.astype(np.float64)
    eig = getEigEnMolden(name, spin)
    eig2 = RS.cut_eigenenergies(eig)
    print("Eigenenergies in the selected range are: ", eig2)
    print("Eigen shape", np.shape(eig2))
    impN = len(eig2)

    from pyPPSTM.ReadSTM import n_max_, n_min_, Ynum_
    print("n_max_", n_max_)
    print("n_min_", n_min_)
    print("num_at_", num_at_)
    print("Ynum_", Ynum_)

    cc = getMOLDENcoefs(name, spin, noCoefz=noCoefz, impN=impN)
    cc = cc.astype(float)
    n_bands = len(cc)
    print("cc has the shape = ", np.shape(cc))
    #print(cc)
    print(n_bands)
    coefs = np.zeros((n_bands, num_at_, Ynum_))
    print("Ynum_ = ", Ynum_)

    for i in np.arange(n_bands):
        for j in np.arange(len(spec)):
            k = pointerz[j]
            # print("species: ", spec[j])
            # print("j = ", j)
            # =print("k = ", k)
            if spec[j] == 'C' or spec[j] == 'N':
                coefs[i, j, 0] = cc[i][k + 0] + cc[i][k + 1] + cc[i][k + 2]  # S-orbital supposedly
                coefs[i, j, 3] = cc[i][k + 3] + cc[i][k + 6]  # Px-orbital supposedly
                coefs[i, j, 1] = cc[i][k + 4] + cc[i][k + 7]  # Py-orbital supposedly
                coefs[i, j, 2] = cc[i][k + 5] + cc[i][k + 8]  # Pz-orbital supposedly
            elif spec[j] == 'H':
                coefs[i, j, 0] = cc[i][k + 0] + cc[i][k + 1]  # S-orbital supposedly
                coefs[i, j, 3] = cc[i][k + 2]  # Px-orbital supposedly
                coefs[i, j, 1] = cc[i][k + 3]  # Py-orbital supposedly
                coefs[i, j, 2] = cc[i][k + 4]  # Pz-orbital supposedly
            else:
                print("Species other than carbon, nitrogen and hydrogen detected. Can't handle it, sorry.")
    # print("coefs = ", coefs)
    #print(np.shape(coefs))

    # coefs = np.zeros((n_max_-n_min_,num_at_,Ynum_))

    #  print(coefs)
    coef = coefs[n_min_:n_max_, :, :]
#    print(len(coefs))
#    print(coefs)
    coeffs = RS.handle_coef(coef)
    print("All coefficients read")
    # print("coefz2 = \n",coeffs)
    coeffs = np.array(coeffs)
    eig2 = np.array(eig2)

    return eig2.copy(), coeffs.copy(), Ratin.copy();
# ===========================================================================================================
# The end.
