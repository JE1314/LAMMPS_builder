#!/usr/bin/env python

#15-09-2023

#Authors:Sebastian ECHEVERRI RESTREPO,   
#               sebastian.echeverri.restrepo@skf.com, sebastianecheverrir@gmail.com

#################################################################################3

#  This file contains a function to generate the rough Fe surfaces. The 
#    algorithm used is based on Tribol Lett 2011, 44, 279-85 and 
#    http://doi.org/10.1007/978-3-642-84574-1_34

#################################################################################3

from ase.spacegroup import crystal
import ase.io
import os
from ase import Atoms
from ase.visualize import view
from ase.build import surface
from ase import Atom
from random import randint
from random import shuffle
from ase.neighborlist import *
import numpy as np
import copy
import random 
import math
from ase.lattice.cubic import BodyCenteredCubic


def RoughCrack(FractalLevels,RMSin,H,boxLenghtX,boxLenghtY,boxLenghtZ,aFe,Separation,Orientation):

  #converting the box lenght to angstroms
  boxLenghtXAngs=boxLenghtX*aFe
  boxLenghtYAngs=boxLenghtY*aFe
  boxLenghtZAngs=boxLenghtZ*aFe

  #mean for the random number generator
  mu=0.0

  scaleFactorK=1.5*RMSin


  #####################################################################
  #Fractal Section For Bulk Fe
  #####################################################################
  print('#######################################')
  print('Generating the Fractal surface for bulk Fe')

  class Point:
    def __init__(self,x,y,z):
      self.x = x
      self.y = y
      self.z = z


  # Create an array
  PointsArray = [[0.0 for x in range(0,int(math.sqrt(4**FractalLevels)+1))] 
                    for y in range(0,int(math.sqrt(4**FractalLevels)+1))]
  xSize=len(PointsArray)
  ySize=len(PointsArray[0])

  #generatig the corner points
  zCorner=0.0 #random.gauss(mu,scaleFactorK*2**(-(1)*H))
  PointsArray[0][0]                     =Point(0.0,0.0,zCorner)
  PointsArray[xSize-1][0]               =Point(boxLenghtXAngs,0.0,zCorner)
  PointsArray[0][ySize-1]               =Point(0.0,boxLenghtYAngs,zCorner)
  PointsArray[xSize-1][ySize-1] =Point(boxLenghtXAngs,boxLenghtYAngs,zCorner)


  #Number of squares per side
  NSideSquares=2**(FractalLevels)

  RMS=0.0

  while (RMS > 1.05*RMSin or RMS < 0.95*RMSin):
    #Fractal cycle 
    for level in range(0, FractalLevels):
      deltaSideSquares=2**(FractalLevels-level)
      
      #calculate the central point
      for i in range(0,NSideSquares,deltaSideSquares):
        for j in range(0,NSideSquares,deltaSideSquares):
          
          PointsArray[i+int(deltaSideSquares/2)][j+int(deltaSideSquares/2)]=Point((PointsArray[i][j].x+PointsArray[i+deltaSideSquares][j].x)/2,
                                                                        (PointsArray[i][j].y+PointsArray[i][j+deltaSideSquares].y)/2, 
                                                                        (PointsArray[i][j].z+PointsArray[i][j+deltaSideSquares].z+PointsArray[i+deltaSideSquares][j].z+PointsArray[i+deltaSideSquares][j+deltaSideSquares].z)/4+
                                                                        random.gauss(mu,scaleFactorK*2**(-(level+1)*H)))
      #calculate the sides
      for i in range(0,NSideSquares,deltaSideSquares):
        for j in range(0,NSideSquares,deltaSideSquares):
          #point 0
          if i==0:
            PointsArray[i][j+int(deltaSideSquares/2)]=Point((i)*boxLenghtXAngs/NSideSquares, (j+int(deltaSideSquares/2))*boxLenghtYAngs/NSideSquares, 
                                                      (PointsArray[i-int(deltaSideSquares/2)-1][j+int(deltaSideSquares/2)].z+PointsArray[i][j].z+PointsArray[i+int(deltaSideSquares/2)][j+int(deltaSideSquares/2)].z+PointsArray[i][j+deltaSideSquares].z)/4+
                                                        random.gauss(mu,scaleFactorK*2**(-(level+1)*H)))
          else:
            PointsArray[i][j+int(deltaSideSquares/2)]=Point((i)*boxLenghtXAngs/NSideSquares, (j+int(deltaSideSquares/2))*boxLenghtYAngs/NSideSquares, 
                                                      (PointsArray[i-int(deltaSideSquares/2)][j+int(deltaSideSquares/2)].z+PointsArray[i][j].z+PointsArray[i+int(deltaSideSquares/2)][j+int(deltaSideSquares/2)].z+PointsArray[i][j+deltaSideSquares].z)/4+
                                                      random.gauss(mu,scaleFactorK*2**(-(level+1)*H)))
          #point 1
          if j==0:
            PointsArray[i+int(deltaSideSquares/2)][j]=Point((i+int(deltaSideSquares/2))*boxLenghtXAngs/NSideSquares, (j)*boxLenghtYAngs/NSideSquares, 
                                                      (PointsArray[i][j].z+PointsArray[i+int(deltaSideSquares/2)][j-int(deltaSideSquares/2)-1].z+PointsArray[i+deltaSideSquares][j].z+PointsArray[i+int(deltaSideSquares/2)][j+int(deltaSideSquares/2)].z)/4+
                                                      random.gauss(mu,scaleFactorK*2**(-(level+1)*H)))
          else:
            PointsArray[i+int(deltaSideSquares/2)][j]=Point((i+int(deltaSideSquares/2))*boxLenghtXAngs/NSideSquares, (j)*boxLenghtYAngs/NSideSquares, 
                                                      (PointsArray[i][j].z+PointsArray[i+int(deltaSideSquares/2)][j-int(deltaSideSquares/2)].z+PointsArray[i+deltaSideSquares][j].z+PointsArray[i+int(deltaSideSquares/2)][j+int(deltaSideSquares/2)].z)/4+
                                                      random.gauss(mu,scaleFactorK*2**(-(level+1)*H)))
          #point 2
          if i+deltaSideSquares== xSize-1:
            PointsArray[i+deltaSideSquares][j+int(deltaSideSquares/2)]=Point((i+deltaSideSquares)*boxLenghtXAngs/NSideSquares, (j+int(deltaSideSquares/2))*boxLenghtYAngs/NSideSquares, 
                                                                        PointsArray[0][j+int(deltaSideSquares/2)].z)
          else:
            PointsArray[i+deltaSideSquares][j+int(deltaSideSquares/2)]=Point((i+deltaSideSquares)*boxLenghtXAngs/NSideSquares, (j+int(deltaSideSquares/2))*boxLenghtYAngs/NSideSquares, 
                                                                        (PointsArray[i+int(deltaSideSquares/2)][j+int(deltaSideSquares/2)].z+PointsArray[i+deltaSideSquares][j].z+PointsArray[i+3*int(deltaSideSquares/2)][j+int(deltaSideSquares/2)].z+PointsArray[i+deltaSideSquares][j+deltaSideSquares].z)/4+
                                                                        random.gauss(mu,scaleFactorK*2**(-(level+1)*H)))
          #point 3
          if j+deltaSideSquares==ySize-1:
            PointsArray[i+int(deltaSideSquares/2)][j+deltaSideSquares]=Point((i+int(deltaSideSquares/2))*boxLenghtXAngs/NSideSquares, (j+deltaSideSquares)*boxLenghtYAngs/NSideSquares,
                                                                      PointsArray[i+int(deltaSideSquares/2)][0].z)#copy.deepcopy(PointsArray[i+int(deltaSideSquares/2)][j].z))
          else:
            PointsArray[i+int(deltaSideSquares/2)][j+deltaSideSquares]=Point((i+int(deltaSideSquares/2))*boxLenghtXAngs/NSideSquares, (j+deltaSideSquares)*boxLenghtYAngs/NSideSquares,
                                                                      (PointsArray[i][j+deltaSideSquares].z+PointsArray[i+int(deltaSideSquares/2)][j+int(deltaSideSquares/2)].z+PointsArray[i+deltaSideSquares][j+deltaSideSquares].z+PointsArray[i+int(deltaSideSquares/2)][j+3*int(deltaSideSquares/2)].z)/4+
                                                                      random.gauss(mu,scaleFactorK*2**(-(level+1)*H)))


    #forcing the min value to be at z=0
    #finding the min Z value

    minZ=PointsArray[0][0].z
    sumsq=0
    npts=0

    for i in range(0,xSize):
      for j in range(0,ySize):
        sumsq=(PointsArray[i][j].z-mu)**2+sumsq
        if minZ>PointsArray[i][j].z:
          minZ=PointsArray[i][j].z

    RMS=math.sqrt(sumsq/(xSize*ySize))

  print("RMS = ", RMS, "Angstrom")




  for i in range(0,xSize):
    for j in range(0,ySize):
      PointsArray[i][j].z=PointsArray[i][j].z-minZ

  #Making squares from the points
  class Square:
    def __init__(self, x0,y0,x1,y1,x2,y2,x3,y3,z):
      self.x0 = x0
      self.y0 = y0
      self.x1 = x1
      self.y1 = y1
      self.x2 = x2
      self.y2 = y2
      self.x3 = x3
      self.y3 = y3
      self.z  = z

  SquaresArray = [[0.0 for x in range(0,xSize-1)] 
                    for y in range(0,ySize-1)]

  for i in range(0,xSize-1):
    for j in range(0,ySize-1):
      SquaresArray[i][j]=Square(PointsArray[i][j].x, PointsArray[i][j].y,
                                PointsArray[i+1][j].x, PointsArray[i+1][j].y,
                                PointsArray[i+1][j+1].x, PointsArray[i+1][j+1].y,
                                PointsArray[i][j+1].x, PointsArray[i][j+1].y,
                                (PointsArray[i][j].z+PointsArray[i+1][j].z+PointsArray[i+1][j+1].z+PointsArray[i][j+1].z)/4)
      



  #Plotting the squares and the points  
  #fig = plt.figure()
  #ax = fig.add_subplot(111, projection='3d')

  x=[]
  y=[]
  z=[]
  for i in range(0,xSize):
    for j in range(0,ySize):
      x.append(PointsArray[i][j].x)
      y.append(PointsArray[i][j].y)
      z.append(PointsArray[i][j].z)
      #x.append(PointsArray[i][j].x+boxLenghtXAngs)
      #y.append(PointsArray[i][j].y)
      #z.append(PointsArray[i][j].z)
      #x.append(PointsArray[i][j].x)
      #y.append(PointsArray[i][j].y+boxLenghtYAngs)
      #z.append(PointsArray[i][j].z)
      #x.append(PointsArray[i][j].x+boxLenghtXAngs)
      #y.append(PointsArray[i][j].y+boxLenghtYAngs)
      #z.append(PointsArray[i][j].z)
      #print '%5.3f,%5.3f,%5.3f \t\t'%(PointsArray[i][j].x,PointsArray[i][j].y,PointsArray[i][j].z ),
    #print ""

  #for i in range(0,xSize-1):
  # for j in range(0,ySize-1):   
  #   ax.plot([SquaresArray[i][j].x0,SquaresArray[i][j].x1,SquaresArray[i][j].x2,SquaresArray[i][j].x3, SquaresArray[i][j].x0 ],
  #    [SquaresArray[i][j].y0,SquaresArray[i][j].y1,SquaresArray[i][j].y2,SquaresArray[i][j].y3, SquaresArray[i][j].y0], 0)

  #ax.scatter(x,y,z,c=z)
  #ax.plot([1,1], [2,2], 0)
  #ax.plot([SquaresArray[0][0].x0,SquaresArray[0][0].x1,SquaresArray[0][0].x2,SquaresArray[0][0].x3, SquaresArray[0][0].x0 ],
          #[SquaresArray[0][0].y0,SquaresArray[0][0].y1,SquaresArray[0][0].y2,SquaresArray[0][0].y3,SquaresArray[0][0].y0], 0)
  #ax.plot_trisurf(x, y, z, cmap=plt.cm.Spectral)
  #ax.set_xlabel('X axis')
  #ax.set_ylabel('Y axis')

  #plt.show()


  #####################################################################
  #generating a bulk Fe crystal 
  print('#######################################')
  print('Generating the bulk Fe region')

#  atomsBulk = crystal(spacegroup=229,
#                  symbols='Fe',
#                  basis=[0,0,0],
#                  cellpar=[aFe,aFe,aFe,90.0,90.0,90.0],
#                  size=(boxLenghtX,boxLenghtY,boxLenghtZ))
  atomsBulk = BodyCenteredCubic(directions=Orientation,
                                size=(boxLenghtX,boxLenghtY,boxLenghtZ),
                                symbol='Fe',
                                latticeconstant=aFe)

  #ase.io.write("FeBulk.cfg", atomsBulk, "cfg")
  #os.system("atomsk FeBulk.cfg lmp >& /dev/null")
  #os.system("mv FeBulk.lmp data.FeBulk")



  #####################################################################
  #Making the crystal surface rough

  print('#######################################')
  print('Applying the roughness to the bulk Fe region')

  atomsBulkRough = copy.deepcopy(atomsBulk)

  for k in reversed(range(0, Atoms.get_global_number_of_atoms(atomsBulkRough))):
    deleted=False
    for i in range(0,xSize-1):
      for j in range(0,ySize-1):
        if atomsBulkRough[k].x>=SquaresArray[i][j].x0 and atomsBulkRough[k].x<SquaresArray[i][j].x1 and atomsBulkRough[k].y>=SquaresArray[i][j].y0 and atomsBulkRough[k].y<SquaresArray[i][j].y2 and atomsBulkRough[k].z<SquaresArray[i][j].z:
          del atomsBulkRough[k]
          deleted=True
          break
      if deleted==True:
        break


  #ase.io.write("RoughFeBulk.cfg", atomsBulkRough, "cfg")

  #os.system("atomsk RoughFeBulk.cfg lmp >& /dev/null")
  #os.system("mv RoughFeBulk.lmp data.RoughFeBulk")




  #####################################################################
  #Fractal Section For Bulk2 Fe
  #####################################################################
  print('#######################################')
  print('Generating the Fractal surface for Bulk2 Fe')



  # Create an array
  PointsArray = [[0.0 for x in range(0,int(math.sqrt(4**FractalLevels)+1))] 
                    for y in range(0,int(math.sqrt(4**FractalLevels)+1))]
  xSize=len(PointsArray)
  ySize=len(PointsArray[0])

  #generatig the corner points
  zCorner=0.0#random.gauss(mu,scaleFactorK*2**(-(1)*H))
  PointsArray[0][0]                     =Point(0.0,0.0,zCorner)
  PointsArray[xSize-1][0]               =Point(boxLenghtXAngs,0.0,zCorner)
  PointsArray[0][ySize-1]               =Point(0.0,boxLenghtYAngs,zCorner)
  PointsArray[xSize-1][ySize-1] =Point(boxLenghtXAngs,boxLenghtYAngs,zCorner)


  #Number of squares per side
  NSideSquares=2**(FractalLevels)

  RMS=0.0

  while (RMS > 1.05*RMSin or RMS < 0.95*RMSin):

    #Fractal cycle 
    for level in range(0, FractalLevels):
      deltaSideSquares=2**(FractalLevels-level)
      
      #calculate the central point
      for i in range(0,NSideSquares,deltaSideSquares):
        for j in range(0,NSideSquares,deltaSideSquares):
          
          PointsArray[i+int(deltaSideSquares/2)][j+int(deltaSideSquares/2)]=Point((PointsArray[i][j].x+PointsArray[i+deltaSideSquares][j].x)/2,
                                                                        (PointsArray[i][j].y+PointsArray[i][j+deltaSideSquares].y)/2, 
                                                                        (PointsArray[i][j].z+PointsArray[i][j+deltaSideSquares].z+PointsArray[i+deltaSideSquares][j].z+PointsArray[i+deltaSideSquares][j+deltaSideSquares].z)/4+
                                                                        random.gauss(mu,scaleFactorK*2**(-(level+1)*H)))
      #calculate the sides
      for i in range(0,NSideSquares,deltaSideSquares):
        for j in range(0,NSideSquares,deltaSideSquares):
          #point 0
          if i==0:
            PointsArray[i][j+int(deltaSideSquares/2)]=Point((i)*boxLenghtXAngs/NSideSquares, (j+int(deltaSideSquares/2))*boxLenghtYAngs/NSideSquares, 
                                                      (PointsArray[i-int(deltaSideSquares/2)-1][j+int(deltaSideSquares/2)].z+PointsArray[i][j].z+PointsArray[i+int(deltaSideSquares/2)][j+int(deltaSideSquares/2)].z+PointsArray[i][j+deltaSideSquares].z)/4+
                                                        random.gauss(mu,scaleFactorK*2**(-(level+1)*H)))
          else:
            PointsArray[i][j+int(deltaSideSquares/2)]=Point((i)*boxLenghtXAngs/NSideSquares, (j+int(deltaSideSquares/2))*boxLenghtYAngs/NSideSquares, 
                                                      (PointsArray[i-int(deltaSideSquares/2)][j+int(deltaSideSquares/2)].z+PointsArray[i][j].z+PointsArray[i+int(deltaSideSquares/2)][j+int(deltaSideSquares/2)].z+PointsArray[i][j+deltaSideSquares].z)/4+
                                                      random.gauss(mu,scaleFactorK*2**(-(level+1)*H)))
          #point 1
          if j==0:
            PointsArray[i+int(deltaSideSquares/2)][j]=Point((i+int(deltaSideSquares/2))*boxLenghtXAngs/NSideSquares, (j)*boxLenghtYAngs/NSideSquares, 
                                                      (PointsArray[i][j].z+PointsArray[i+int(deltaSideSquares/2)][j-int(deltaSideSquares/2)-1].z+PointsArray[i+deltaSideSquares][j].z+PointsArray[i+int(deltaSideSquares/2)][j+int(deltaSideSquares/2)].z)/4+
                                                      random.gauss(mu,scaleFactorK*2**(-(level+1)*H)))
          else:
            PointsArray[i+int(deltaSideSquares/2)][j]=Point((i+int(deltaSideSquares/2))*boxLenghtXAngs/NSideSquares, (j)*boxLenghtYAngs/NSideSquares, 
                                                      (PointsArray[i][j].z+PointsArray[i+int(deltaSideSquares/2)][j-int(deltaSideSquares/2)].z+PointsArray[i+deltaSideSquares][j].z+PointsArray[i+int(deltaSideSquares/2)][j+int(deltaSideSquares/2)].z)/4+
                                                      random.gauss(mu,scaleFactorK*2**(-(level+1)*H)))
          #point 2
          if i+deltaSideSquares== xSize-1:
            PointsArray[i+deltaSideSquares][j+int(deltaSideSquares/2)]=Point((i+deltaSideSquares)*boxLenghtXAngs/NSideSquares, (j+int(deltaSideSquares/2))*boxLenghtYAngs/NSideSquares, 
                                                                        PointsArray[0][j+int(deltaSideSquares/2)].z)
          else:
            PointsArray[i+deltaSideSquares][j+int(deltaSideSquares/2)]=Point((i+deltaSideSquares)*boxLenghtXAngs/NSideSquares, (j+int(deltaSideSquares/2))*boxLenghtYAngs/NSideSquares, 
                                                                        (PointsArray[i+int(deltaSideSquares/2)][j+int(deltaSideSquares/2)].z+PointsArray[i+deltaSideSquares][j].z+PointsArray[i+3*int(deltaSideSquares/2)][j+int(deltaSideSquares/2)].z+PointsArray[i+deltaSideSquares][j+deltaSideSquares].z)/4+
                                                                        random.gauss(mu,scaleFactorK*2**(-(level+1)*H)))
          #point 3
          if j+deltaSideSquares==ySize-1:
            PointsArray[i+int(deltaSideSquares/2)][j+deltaSideSquares]=Point((i+int(deltaSideSquares/2))*boxLenghtXAngs/NSideSquares, (j+deltaSideSquares)*boxLenghtYAngs/NSideSquares,
                                                                      PointsArray[i+int(deltaSideSquares/2)][0].z)#copy.deepcopy(PointsArray[i+int(deltaSideSquares/2)][j].z))
          else:
            PointsArray[i+int(deltaSideSquares/2)][j+deltaSideSquares]=Point((i+int(deltaSideSquares/2))*boxLenghtXAngs/NSideSquares, (j+deltaSideSquares)*boxLenghtYAngs/NSideSquares,
                                                                      (PointsArray[i][j+deltaSideSquares].z+PointsArray[i+int(deltaSideSquares/2)][j+int(deltaSideSquares/2)].z+PointsArray[i+deltaSideSquares][j+deltaSideSquares].z+PointsArray[i+int(deltaSideSquares/2)][j+3*int(deltaSideSquares/2)].z)/4+
                                                                      random.gauss(mu,scaleFactorK*2**(-(level+1)*H)))


    #forcing the max value to be at z=boxLenghtZAngs
    #finding the min Z value

    maxZ=PointsArray[0][0].z
    sumsq=0.0

    for i in range(0,xSize):
      for j in range(0,ySize):
        sumsq=(PointsArray[i][j].z-mu)**2+sumsq

        if maxZ<PointsArray[i][j].z:
          maxZ=PointsArray[i][j].z
        
      RMS=math.sqrt(sumsq/(xSize*ySize))
    
    
  print("RMS = ", RMS, "Angstrom")
        
  for i in range(0,xSize):
    for j in range(0,ySize):
      PointsArray[i][j].z=PointsArray[i][j].z+boxLenghtZAngs-maxZ







  #Making squares from the points

  SquaresArray = [[0.0 for x in range(0,xSize-1)] 
                    for y in range(0,ySize-1)]

  for i in range(0,xSize-1):
    for j in range(0,ySize-1):
      SquaresArray[i][j]=Square(PointsArray[i][j].x, PointsArray[i][j].y,
                                PointsArray[i+1][j].x, PointsArray[i+1][j].y,
                                PointsArray[i+1][j+1].x, PointsArray[i+1][j+1].y,
                                PointsArray[i][j+1].x, PointsArray[i][j+1].y,
                                (PointsArray[i][j].z+PointsArray[i+1][j].z+PointsArray[i+1][j+1].z+PointsArray[i][j+1].z)/4)


  #Plotting the squares and the points  
  #fig = plt.figure()
  #ax = fig.add_subplot(111, projection='3d')

  x=[]
  y=[]
  z=[]
  for i in range(0,xSize):
    for j in range(0,ySize):
      x.append(PointsArray[i][j].x)
      y.append(PointsArray[i][j].y)
      z.append(PointsArray[i][j].z)


  #####################################################################
  #Making the crystal surface rough

  print('#######################################')
  print('Applying the roughness to the Bulk2 Fe region')

  atomsBulk2Rough = copy.deepcopy(atomsBulk)

  for k in reversed(range(0, Atoms.get_global_number_of_atoms(atomsBulk2Rough))):
    deleted=False
    for i in range(0,xSize-1):
      for j in range(0,ySize-1):
        if atomsBulk2Rough[k].x>=SquaresArray[i][j].x0 and atomsBulk2Rough[k].x<SquaresArray[i][j].x1 and atomsBulk2Rough[k].y>=SquaresArray[i][j].y0 and atomsBulk2Rough[k].y<SquaresArray[i][j].y2 and atomsBulk2Rough[k].z>SquaresArray[i][j].z:
          del atomsBulk2Rough[k]
          deleted=True
          break
      if deleted==True:
        break




  #Assembling the interface
  print('#######################################')
  print('Assembling the Bulk1 and the Bulk2 Fe regions')

#  atomsBulk2Rough.center(vacuum=0, axis=2)
#  atomsBulkRough.center(vacuum=0, axis=2)

  atomsBulk2Rough.translate([0,0,0])
  atomsBulkRough.translate([0,0,boxLenghtZAngs+aFe*Separation])

#  atomsWEA=atomsBulk2Rough+atomsBulkRough


  #############################################################
  #############################################################
  #############################################################


  #Creating a bulk slab to close the crack
  print('#######################################')
  print('Generating the bulk Fe region')

#  atomsBulkClose = crystal(spacegroup=229,
#                      symbols='Fe',
#                      basis=[0,0,0],
#                      cellpar=[aFe,aFe,aFe,90.0,90.0,90.0],
#                      size=(boxLenghtX,boxLenghtY,2*boxLenghtZ+Separation))

  atomsBulkClose = BodyCenteredCubic(directions=Orientation,
                                    size=(boxLenghtX,boxLenghtY,2*boxLenghtZ+Separation), 
                                    symbol='Fe', 
                                    latticeconstant=aFe)


  atomsBulkClose.translate([boxLenghtXAngs,0,0])

  atomsWEA=atomsBulk2Rough+atomsBulkRough+atomsBulkClose


  atomsWEA.center(vacuum=15, axis=0)
  atomsWEA.center(vacuum=15, axis=2)

  atomsWEA.write("data.RoughCrack", format="lammps-data")


if __name__ == "__main__":
        
    FractalLevels = 4
    RMSin         = 8.10
    H             = 0.8
    boxLenghtX    = 20
    boxLenghtY    = 10
    boxLenghtZ    = 20
    aFe           = 2.86366
    Separation    = 2
    Orientation   = [[1,0,0], [0,1,0], [0,0,1]] 

    RoughCrack(FractalLevels,RMSin,H,boxLenghtX,boxLenghtY,boxLenghtZ,aFe,Separation,Orientation)

