#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
from . import elements
import math
#import matplotlib.pyplot as plt


# default variables:

default_atom_size     =  0.10

# procedures for loading geometry from different files:

def loadAtoms( name , sl=False):
    f = open(name,"r")
    n=0;
    l = f.readline()
    #print "--",l,"--"
    try:
        n=int(l)
    except:
                raise ValueError("First line of a xyz file should contain the "
                "number of atoms. Aborting...")
    if (n>0):
        n=int(l)
        e=[];x=[];y=[]; z=[]; q=[]
        i = 0;
        for line in f:
            if sl :
                print(" forced skipped line : ", line)
                sl = False
            else:
                words=line.split()
                nw = len( words)
                ie = None
            try:
                e.append( words[0] )
                x.append( float(words[1]) )
                y.append( float(words[2]) )
                z.append( float(words[3]) )
                if ( nw >=5 ):
                    q.append( float(words[4]) )
                else:
                    q.append( 0.0 )
            except:
                print(" skipped line : ", line)
    f.close()
    nDim = []
    lvec = [] 
    return [ e,x,y,z,q ], nDim, lvec

def loadGeometryIN( fname ):
    print("importin atoms from FHI-AIMS input")
    f = open(fname)
    e=[];x=[];y=[]; z=[]; q=[]
    lvec = [] 
    for i in np.arange(10000):
        ws = f.readline().split()
        print(ws)
        if (len(ws)>0):
            if (ws[0]=='atom'):
                e.append(ws[4]); x.append(float(ws[1])); y.append(float(ws[2])); z.append(float(ws[3])); q.append(0);
            elif (ws[0]=='lattice_vector'):
                lvec.append([float(ws[1]),float(ws[2]),float(ws[3])])
            elif (ws[0]=='trust_radius'):
                break
    f.close()
    print ("lvec", lvec)
    print ("e,x,y,z", e,x,y,z)
    nDim = []
    origLvec = lvec
    #print("DEBUG : readGeometry.in, lvec is = ", lvec)
    print("DEBUG : readGeometry.in, e is = ", e)
    return [ e,x,y,z,q ], nDim, lvec, origLvec

# other procedures for operating with geometries:

def multCell( xyz, cel, m=(2,2,1) ):
    n = len(xyz[0])
    mtot = m[0]*m[1]*m[2]*n
    es = [None] * mtot
    xs = [None] * mtot
    ys = [None] * mtot
    zs = [None] * mtot
    j  = 0
    for ia in range(m[0]):
        for ib in range(m[1]):
            for ic in range(m[2]):
                dx = ia*cel[0][0] + ib*cel[1][0] + ic*cel[2][0]
                dy = ia*cel[0][1] + ib*cel[1][1] + ic*cel[2][1]
                dz = ia*cel[0][2] + ib*cel[1][2] + ic*cel[2][2]
                for i in range(n):
                    es[j]=xyz[0][i]
                    xs[j]=xyz[1][i] + dx
                    ys[j]=xyz[2][i] + dy
                    zs[j]=xyz[3][i] + dz
                    j+=1
    return [es,xs,ys,zs]

# =========== Utils for plotting atoms =========================

XSF_HEAD_0='''ATOMS
'''

XSF_HEAD_2='''
BEGIN_BLOCK_DATAGRID_3D                        
   some_datagrid      
   BEGIN_DATAGRID_3D_whatever 
'''

#def At2XSF(atoms):
#    XSF_HEAD_1=XSF_HEAD_0
#    for i in range(len(atoms[0])):
#        XSF_HEAD_1 = XSF_HEAD_1+str(atoms[0][i])+" "+str(atoms[1][i])+" "+str(atoms[2][i])+" "+str(atoms[3][i])+"\n "
#    XSF_HEAD=XSF_HEAD_1 + XSF_HEAD_2
    #print "DEBUG: XSF_HEAD:"
    #print XSF_HEAD
#    return XSF_HEAD ;



def writeMatrix( fout, mat ):
    for v in mat:
        for num in v: fout.write(' %f ' %num )
        fout.write('\n')

#def writeAtoms( f, elems, xyzs ):
#    for i in range(len(elems)):
#        xyzsi = xyzs[i]
#        f.write( str(elems[i] ) );
#        f.write( " %10.10f %10.10f %10.10f\n" %(str(xyzsi[1][i]), str(xyzsi[2][i]), str(xyzsi[3][i]) )

'''
def saveGeomXSF( fname,elems,xyzs, primvec ):
    if convvec is None:
        primvec = convvec
    with open(fname,'w') as f:
        f.write('CRYSTAL\n')
        f.write('PRIMVEC\n')
        writeMatrix(f, primvec)
        f.write( 'CONVVEC\n' )
        writeMatrix( f, convvec )
        f.write( 'PRIMCOORD\n' )
        f.write( '%i %i\n' %(len(elems),1) )
        if bTransposed:
            writeAtomsTransposed( f, elems, xyzs )
        else:
            writeAtoms( f, elems, xyzs )
        f.write( '\n' )
'''
def At2XSF(fname,atoms,primvec,origLvec,ELEMENTS=elements.ELEMENTS):
    print("origLvec", origLvec)
    #convvec = primvec2
    print("atoms = ",atoms)
    spex=[]
    for e in atoms[0]:
        #print(speciesDict[e])
        ele = speciesDict[e]
        spex.append(ele)

    with open(fname,'w') as f:
        f.write('CRYSTAL\n')
        f.write('PRIMVEC\n')
        writeMatrix(f, origLvec)
        f.write('CONVVEC\n')
        writeMatrix(f, origLvec)
        f.write('PRIMCOORD\n')
        f.write( '%i %i\n' %(len(atoms[0]),1) )

        for i in range(len(spex)):
            f.write(str(spex[i]));
           # f.write(" %10.10f %10.10f %10.10f\n" % (str(atoms[1][i]), str(atoms[2][i]), str(atoms[3][i])))
            f.write(" %10.10f %10.10f %10.10f\n" % (atoms[1][i], atoms[2][i], atoms[3][i]))

        f.write("   BEGIN_BLOCK_DATAGRID_3D\n")
        f.write("   some_datagrid\n")
        f.write("BEGIN_DATAGRID_3D_whateve\n")

def specDict( ELEMENTS ):
    dic = { }
    for elem in  ELEMENTS:
        dic[ elem[1] ] = elem[0]
    return dic

speciesDict = specDict( ELEMENTS=elements.ELEMENTS )
# print("DEBUG: ELEMENT_DICT", speciesDict)

def flattenLst(input):
    new_list = []
    for i in input:
        for j in i:
            new_list.append(j)
    return new_list
'''
'''
def origVec(lvs):
    if dft_code != 'aims' or dft_code != 'AIMS' or dft_code != 'FHI-AIMS' or lvs != None:
        print ("CELL is = ", lvs)
        origLvec = lvs
      #  Rx =  cel[0][0] +  cel[1][0] +  cel[2][0]
      #  Ry =  cel[0][1] +  cel[1][1] +  cel[2][1]
      #  Rz =  cel[0][2] +  cel[1][2] +  cel[2][2]

    else :
        tmp1, tmp2, origLvec = get_AIMS_geom(geom=geom)
    return origLvec


def findBonds(xs, ys, zs, rmax):
   rr=float(rmax)
   bonds = []
   n = len( xs )
   for i in np.arange(n):
      for j in np.arange(i):
         dx=xs[j]-xs[i]
         dy=ys[j]-ys[i]
         dz=zs[j]-zs[i]
         r=math.sqrt(dx*dx+dy*dy+dz*dz)
         if (r<rr) :
            bonds.append( (i,j) )
#            plt.arrow(xs[i], ys[i], xs[j]-xs[i], ys[j]-ys[i], head_width=0.1, head_length=0.1,  fc='r', ec='C1', lw= 1.0,ls='solid',zorder=2 )
   return bonds
#   return plt.arrow(xs[i], ys[i], xs[j]-xs[i], ys[j]-ys[i], head_width=0.1, head_length=0.1,  fc='r', ec='C1', lw= 1.0,ls='solid',zorder=2 )

def plotBonds(xs, ys, bb):
   for (q1,q2) in bb:
#     i=q1[0]; j=q2[1];
       plt.arrow(xs[q1], ys[q1], xs[q2]-xs[q1], ys[q2]-ys[q1], head_width=0.0, head_length=0.0,  fc='none', ec='w', lw= 2.0,ls='dotted',zorder=2, alpha= 0.9)
     #  print (q1,q2)


def pltcolor(lst):
    cols=[]
    for l in lst:
        if l=='C':
            cols.append('g')
        elif l=='H':
            cols.append('w')
        else:
            cols.append('b')
    return cols