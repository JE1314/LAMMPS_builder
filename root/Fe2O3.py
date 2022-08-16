#!/usr/bin/env python

#07-11-2019

#Authors:Sebastian ECHEVERRI RESTREPO,
#	 	sebastian.echeverri.restrepo@skf.com, sebastianecheverrir@gmail.com
#	 James EWEN
#		j.ewen14@imperial.ac.uk, jimmyewen@gmail.com

#################################################################################3

#  This file contains a function to generate the smooth Fe2O3 surfaces.

#################################################################################3

from ase.spacegroup import crystal
import ase.io
import os
from ase.visualize import view
from ase.build import surface
from random import randint
from random import shuffle
from ase.neighborlist import *
import numpy as np
import copy
import random
import math

from ase.lattice.hexagonal import *
from ase.lattice.compounds import *
from ase import Atoms, Atom


def Fe2O3(FractalLevels,RMSin,H,boxLenghtX,boxLenghtY,boxLenghtZ,aFe,Separation):

  #converting the box lenght to angstroms
  aFe2O3 = 5.029
  bFe2O3 = 5.029
  cFe2O3 = 13.730

  boxLenghtXAngs=boxLenghtX*aFe2O3
  boxLenghtYAngs=boxLenghtY*bFe2O3*math.cos(30*math.pi/180)
  boxLenghtZAngs=boxLenghtZ*cFe2O3


  #####################################################################
  #generating a bulk Fe crystal
  print('#######################################')
  print('Generating the bulk Fe2O3 regions')

  atomsBulk = HEX_Fe2O3(symbol = ('Fe', 'O'),
                  latticeconstant={'a':5.029,'b':5.029, 'c':13.730,
                  'alpha':90,
                  'beta':90,
                  'gamma':120},
                  size=(boxLenghtX,boxLenghtY,boxLenghtZ))

  atomsBulkRough = copy.deepcopy(atomsBulk)

  atomsBulk2Rough = copy.deepcopy(atomsBulk)


  #Assembling the system
  print('#######################################')
  print('Assembling the Bulk1 and the Bulk2 Fe2O3 regions')

  atomsBulk2Rough.center(vacuum=0, axis=2)
  atomsBulkRough.center(vacuum=0, axis=2)

  atomsBulk2Rough.translate([0,0,0])
  atomsBulkRough.translate([0,0,boxLenghtZAngs+Separation])

  atomsWEA=atomsBulk2Rough+atomsBulkRough

  #############################################################
  #############################################################
  #############################################################

  print('#######################################')
  print('Writing file Fe2O3.lt for moltemplate')

  #Printing the .lt file for moltemplate
  f = open('Fe2O3.lt', 'w')

  f.write("FESurface inherits LOPLSAA {\n")
  f.write("write(\"Data Atoms\") {\n")

  for k in range(0, Atoms.get_number_of_atoms(atomsWEA)):
    if atomsWEA[k].symbol == 'Fe':
      f.write("$atom:FEX"+str(k)+" $mol:... @atom:10001 0.00 "+str(atomsWEA[k].x)+" "+str(atomsWEA[k].y)+" "+str(atomsWEA[k].z)+"\n")
    elif atomsWEA[k].symbol == 'O':
       f.write("$atom:OX"+str(k)+" $mol:... @atom:10002 0.00 "+str(atomsWEA[k].x)+" "+str(atomsWEA[k].y)+" "+str(atomsWEA[k].z)+"\n")
  f.write("} } \n")
  f.close()
  #    $atom:FE1 $mol:... @atom:10000  0.00   0.000  0.000             0.000
