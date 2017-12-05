#!/usr/bin/env python

#20-11-2017

#Authors:Sebastian ECHEVERRI RESTREPO
#	 	sebastian.echeverri.restrepo@skf.com, sebastianecheverrir@gmail.com
#	 James EWEN
#		j.ewen14@imperial.ac.uk, jimmyewen@gmail.com

#################################################################################3

#  This file contains all the information needed to define the dimensions
#    and structure of the system. This is the only file that should be
#    modified by the (regular) user

#################################################################################3

import os
import sys
sys.path.append("root")
from lopls import lopls
import math

#################################################################################3

#These are the variables that can be modified by the user. They relate to
#the structure of the system

#####General Inputs

#Size of the simulation box in the x and y direction. These
#	values have to be given in terms of the lattice parameter

boxLenghtX = 8
boxLenghtY = 9

#This value determines the Z size of the section of the
#	box that  will contain the alkane chains and OFMs.
#	The units are Angstrom
zhi = 100	#  Also, take into account the space for the alkane chains!


#####Inputs related to the Surfaces

#This flag determines if the the system will contain Surfaces.
#       if Surfaces = 0, do not generate surfaces
#	if Surfaces = 1, generate a rough Fe surface
#       if Surfaces = 2, generate a flat Fe2O3 surface
#
Surfaces = 2 

#Number of fractal levels used for the generation of the surfaces
#	The inputs are integers
#	Valid if Surfaces = 1
FractalLevels = 4

#Hurst exponent for the generation of the fractal surfaces
#       Valid if Surfaces = 1
H = 0.8

#Value of the desired RMS roughness of the surfaces.
#	The units are Angstrom
#       Valid if Surfaces = 1
RMSin = 8.0

#Thicknes of each of the Fe slabs that form the surfaces. This
#	is the thickness before the roughness is applied to the surface
#	This value does not have units. It corresponds to the number
#	of Fe (or Fe2O3) lattice constants
#	Exm: boxLenghtZ = 20 -> 20*aFe = 57.27320 Angstrom
boxLenghtZ = 1

#####Inputs related to the OFM's

#This flag determines if the the system will contain OFMs.
#	if OFM = 1, the generation of OFMs is activated
OFM = 1

#Defines the type of OFM that will be used in the simulation
#	the available options are: SA, SAm, GMS, OA, OAm, GMO
#	Exm: OFMtype = 'GMS'
OFMtype = 'GMO'

#OFM chains to be placed along the x and y directions.
#	The inputs are integers
#	Exm: 	OFMn_x = 5
#		OFMn_y = 7
OFMn_x = 2
OFMn_y = 2

#####Inputs related to the Alkanes

#This flag determines if the the system will contain Alkanes.
#	if Alkane = 1, the generation of Alkanes is activated
Alkane = 0

#Defines the number of monomers on each alkane chain.
#	The inputs are integers
#	Exm: 	nAlkane = 16
nAlkane = 8

#Number of alkane chains to be placed along each
#	direction of the simulation box
#	The inputs are integers
#	Exm: 	Alkanen_x = 5
#		Alkanen_y = 7
#		Alkanen_z = 3
Alkanen_x = 2
Alkanen_y = 3
Alkanen_z = 3

#####Inputs related to benzylbenzoate

#This flag determines if the the system will contain Alkanes.
#       if BZBZ = 1, the generation of benzylbenzoate is activated

BZBZ = 1

#Number of benzylbenzoate chains to be placed along each
#       direction of the simulation box
#       The inputs are integers
#       Exm:    BZBZn_x = 5
#               BZBZn_y = 7
#               BZBZn_z = 3

BZBZn_x = 3
BZBZn_y = 3
BZBZn_z = 3


#################################################################################3
#################################################################################3
#################################################################################3
#fixed input variables (Do not modify anything from here on)

#Lattice parameter of Fe in Angstrom
aFe = 2.86366

aFe2O3 = 5.029
bFe2O3 = 5.029
cFe2O3 = 13.730


xlo = 0                 #do not change!
ylo = 0                 #do not change!
zlo = 0 #Note that the size of the OFM needs to be taken into account (e.g. SA=22.8065)


if Surfaces == 1 or Surfaces == 0:
  xhi = boxLenghtX*aFe + xlo
  yhi = boxLenghtY*aFe + ylo


if Surfaces == 2 :
  xhi = boxLenghtX*aFe2O3  + xlo
  yhi = boxLenghtY*bFe2O3*math.cos(30*math.pi/180)  + ylo
  #zhi = boxLenghtZ*cFe2O3 + zlo



Separation = zhi-zlo+4
potential = 'lopls'


#if potential == 'lopls':
lopls(xlo,xhi,ylo,yhi,zlo,zhi,OFMn_x,OFMn_y,nAlkane, Alkanen_x,\
                Alkanen_y, Alkanen_z, Alkane, BZBZ, BZBZn_x, BZBZn_y, BZBZn_z, OFM,OFMtype, Surfaces,\
                FractalLevels,RMSin,H,boxLenghtX,boxLenghtY,boxLenghtZ,aFe,Separation)

#execfile('root/lopls.py')

print '#######################################'
print 'The folder ' + str(potential) + ' has been created'
print '#######################################'

