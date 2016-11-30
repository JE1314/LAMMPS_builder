#!/usr/bin/env python

#29-11-2016

#Authors:Sebastian ECHEVERRI RESTREPO,   
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

#################################################################################3

#These are the variables that can be modified by the user. They relate to
#the structure of the system

#####General Inputs

#Lattice parameter of Fe in Angstrom
aFe = 2.86366

#Size of the simulation box in the x and y direction. These
#	values have to be given in terms of the lattice parameter
#	of Fe. For example
#	xhi = 15*aFe
#	yhi = 15*aFe
xhi = 17*aFe #40*aFe
yhi = 17*aFe #40*aFe

#This value determines the Z size of the section of the 
#	box that  will contain the alkane chains and OFMs.
#	The units are Angstrom
zhi = 200	#  Also, take into account the space for the alkane chains! 

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
Alkane = 1

#Defines the number of monomers on each alkane chain.
#	The inputs are integers
#	Exm: 	nAlkane = 16
nAlkane= 16

#Number of alkane chains to be placed along each 
#	direction of the simulation box
#	The inputs are integers
#	Exm: 	Alkanen_x = 5
#		Alkanen_y = 7
#		Alkanen_z = 3
Alkanen_x = 2        
Alkanen_y = 4
Alkanen_z = 4

#####Inputs related to the Surfaces

#This flag determines if the the system will contain Surfaces.
#	if Surfaces = 1, the generation of Surfaces is activated
Surfaces = 1


#Number of fractal levels used for the generation of the surfaces
#	The inputs are integers
FractalLevels = 4

#Hurst exponent for the generation of the fractal surfaces
H = 0.8

#Value of the desired RMS roughness of the surfaces.
#	The units are Angstrom
RMSin = 8.0 		


#Thicknes of each of the Fe slabs that form the surfaces. This 
#	is the thickness before the roughness is applied to the surface
#	This value does not have units. It corresponds to the number 
#	of Fe lattice constants
#	Exm: boxLenghtZ = 20 -> 20*aFe = 57.27320 Angstrom
boxLenghtZ = 20

#################################################################################3
#################################################################################3
#################################################################################3

#fixed input variables (Do not modify anything from here on)
xlo = 0			#do not change!
ylo = 0			#do not change!
zlo = 0	#Note that the size of the OFM needs to be taken into account (e.g. SA=22.8065)
boxLenghtX = int((xhi-xlo)/aFe)
boxLenghtY = int((yhi-ylo)/aFe)
Separation = zhi-zlo+4
potential = 'lopls' 	


#if potential == 'lopls':
lopls(xlo,xhi,ylo,yhi,zlo,zhi,OFMn_x,OFMn_y,nAlkane, Alkanen_x,\
		Alkanen_y, Alkanen_z, Alkane, OFM, OFMtype, Surfaces, \
		FractalLevels,RMSin,H,boxLenghtX,boxLenghtY,boxLenghtZ,aFe,Separation)


#execfile('root/lopls.py')

print '#######################################'
print 'The folder ' + str(potential) + ' has been created'
print '#######################################'

