#!/usr/bin/env python

#20-11-2017

#Authors:Sebastian ECHEVERRI RESTREPO,
#	 	sebastian.echeverri.restrepo@skf.com, sebastianecheverrir@gmail.com
#	 James EWEN
#		j.ewen14@imperial.ac.uk, jimmyewen@gmail.com

#################################################################################3
#Appart from generating the Fe and the O atoms for the surfaces, the bonds between
#	them also have to be defined.
#	It seems that it is easier to let lammps do this at the begining of a
#	simulation.

#This script modifies the input files for lammps in order to define all the
#	variables that er needed for the geneartion of the bonds

#################################################################################3

import math
import fileinput


def AddFe2O3(name):
  #name = 'lopls'

#Reading the number of bond types
  f = open(name+'.data','r')
  lines = f.readlines()
  for line in lines:
    if line.find("bond types") != -1:
      line=line.split()
      nBondTypes = int(line[0])
  f.close()


#modifying (rewriting) the file in.lopls

  f = open('in.'+name,'wr+')

  f.write("#-------------- Initialization Section --------------------	\n")
  f.write("include       "+name+".in.CreateBonds             		\n")
  f.write("include       "+name+".in.init				\n")
  f.write("include       "+name+".in.settings                           \n")
  f.write("include       "+name+".in.charges	                	\n\n")

  f.write("dump           dump1 all atom 1000 lopls.dump			\n")
  f.write("thermo_style   custom step lx ly lz  density temp press etotal\n")
  f.write("thermo         1						\n")
  f.write("write_data     loplsInitial.data				\n\n")

  f.write("#--------------- Run Equilibriation --------------------   	\n")
  f.write("min_style       cg						\n")
  f.write("minimize        0.0 0.0 100000 100000			\n\n")

  f.close()

#generating a file to make the dummy start that lammps requires to generate bonds
#	see http://lammps.sandia.gov/threads/msg54748.html
  f = open(name+'.in.CreateBonds','wr+')

  f.write("#dummy start needed to create bonds in lammps		\n")
  f.write("units          metal 					\n")
  f.write("atom_style 	  full 						\n")
  f.write("read_data     "+name+".data extra/bond/types 6		\n")
  f.write("pair_style     lj/cut 10.0					\n")
  f.write("pair_coeff     * * 1.0 1.0					\n")
  f.write("bond_style harmonic						\n")
  f.write("bond_coeff     * 1.0 1.0					\n\n")

  f.write("# Creating groups						\n")
  f.write("group          fe               type      30			\n")
  f.write("group          ox               type      31			\n\n")

  f.write("#Ceating bonds\n")
  f.write("create_bonds   fe ox "+str(nBondTypes+1)+" 1.900 2.000				\n")
  f.write("create_bonds   fe ox "+str(nBondTypes+2)+" 2.000 2.500				\n")
  f.write("create_bonds   ox ox "+str(nBondTypes+3)+" 2.800 2.900				\n")
  f.write("create_bonds   ox ox "+str(nBondTypes+4)+" 2.700 2.799				\n")
  f.write("create_bonds   ox ox "+str(nBondTypes+5)+" 2.600 2.699				\n")
  f.write("create_bonds   fe fe "+str(nBondTypes+6)+" 2.900 3.000				\n")

  f.close()


#modifying (rewriting) the file lopls.in.init

  f = open(name+".in.init",'wr+')

  f.write("bond_style      hybrid harmonic				\n")
  f.write("angle_style     hybrid harmonic				\n")
  f.write("dihedral_style  hybrid opls multi/harmonic			\n")
  f.write("improper_style  hybrid harmonic				\n")
  f.write("pair_style      hybrid lj/cut/coul/long 10.0 10.0		\n")
  f.write("pair_modify     mix geometric				\n")
  f.write("special_bonds   lj/coul 0.0 0.0 0.5				\n")
  f.write("kspace_style    pppm 0.00001					\n")

  f.close()

#modifying (rewriting) the file lopls.in.settings
  replaced = False
  for line in fileinput.input(name+".in.settings", inplace=1):
    if line.startswith('    angle_coeff') and replaced  is False:
      print "    bond_coeff "+str(nBondTypes+1)+" harmonic 5.63738 1.945000"
      print "    bond_coeff "+str(nBondTypes+2)+" harmonic 5.63738 2.116000"
      print "    bond_coeff "+str(nBondTypes+3)+" harmonic 5.63738 2.888000"
      print "    bond_coeff "+str(nBondTypes+4)+" harmonic 5.63738 2.775000"
      print "    bond_coeff "+str(nBondTypes+5)+" harmonic 5.63738 2.669000"
      print "    bond_coeff "+str(nBondTypes+6)+" harmonic 5.63738 2.971000"
      replaced = True
    print line,


#modifying (rewriting) the file lopls.data
#	Adding the line "15  extra bond per atom" to allow more bonds per atom

  for line in fileinput.input(name+".data", inplace=1):
    print line,
    if line.endswith('improper types\n'):
      print "\t15  extra bond per atom"


