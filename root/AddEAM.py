#!/usr/bin/env python

#29-11-2016

#Authors:Sebastian ECHEVERRI RESTREPO,
#	 	sebastian.echeverri.restrepo@skf.com, sebastianecheverrir@gmail.com
#	 James EWEN
#		j.ewen14@imperial.ac.uk, jimmyewen@gmail.com

#################################################################################3

#  This file adds the information needed to define the EAM interaction 
#     for the Fe atoms in the lammps input files

#################################################################################3

import math


def AddEAM():


  class PairCoeff:
    def __init__(self,Type1,Type2,Epsilon,Sigma):
      self.Type1 = Type1
      self.Type2 = Type2
      self.Epsilon = Epsilon
      self.Sigma = Sigma

##################################################33
  #Reading the file lopls.in.settings 
  f=open('lopls.in.settings',"r")
  lines = f.readlines()

  PairCoeffs = []
  OtherCoeffs = []
  i = 0
  for line in lines:
    if line.find("pair_coeff") != -1:
      line=line.split()
      PairCoeffs.append(PairCoeff(int(line[1]), int(line[2]), float(line[4]), float(line[5])))
      #print PairCoeffs[i].Type1, PairCoeffs[i].Type2, PairCoeffs[i].Epsilon, PairCoeffs[i].Sigma
      i += 1
    else :
      OtherCoeffs.append(line)
  f.close()

  #using gemetric combination rules to calculate the interaction between dissimilar atoms
  Ntypes = len(PairCoeffs)
  for i in range(0, Ntypes):
    for j in range(i, Ntypes):
      if i != j:
	type1 = i + 1
	type2 = j + 1
	Epsilon = math.sqrt(PairCoeffs[i].Epsilon*PairCoeffs[j].Epsilon)
	Sigma = math.sqrt(PairCoeffs[i].Sigma*PairCoeffs[j].Sigma)
	PairCoeffs.append(PairCoeff(type1,type2,Epsilon,Sigma))


  #Printing new coefficients
  f = open('lopls.in.settings','w')

  f.write("pair_coeff   * * eam/fs Fe_mm.eam.fs  ")
  for i in range(0, Ntypes-1):
    f.write("NULL ")
  f.write("Fe \n")


  for i in range(0, len(PairCoeffs)):
    if ((PairCoeffs[i].Type1 !=Ntypes) or (PairCoeffs[i].Type2 !=Ntypes)) :
      f.write("pair_coeff "+str(PairCoeffs[i].Type1) + " " +str(PairCoeffs[i].Type2) + " lj/cut/coul/long " + str(PairCoeffs[i].Epsilon) + " " + str(PairCoeffs[i].Sigma) + "\n")


  for i in range(0, len(OtherCoeffs)):
    f.write(OtherCoeffs[i])
  f.close()


############################################
  #adding eam/fs to the line pair_style hybrid lj/cut/coul/long 10.0 10.0

  f=open('lopls.in.init',"r")
  lines = f.readlines()
  LinesOut = []

  for line in lines:
     if line.find("pair_style") != -1:
       LinesOut.append("pair_style hybrid lj/cut/coul/long 10.0 10.0 eam/fs \n")
     else :
       LinesOut.append(line)
  f.close()


  f=open('lopls.in.init',"w")
  for i in range(0, len(LinesOut)):
    f.write(LinesOut[i])



