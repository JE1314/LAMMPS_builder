##Version: 1.0.0, 29-11-2016

##Authors - Sebastián ECHEVERRI RESTREPO: sebastian.echeverri.restrepo@skf.com, sebastianecheverrir@gmail.com; James EWEN: j.ewen14@imperial.ac.uk, jimmyewen@gmail.com

## GENERAL INFO 

 #This software generates a system to be run in LAMMPS(22 Jul 2016-ICMS)
 #containing:
 
 #-Two bcc Fe surfaces with nanoscale RMS roughness
 #-Monolayer films made of organic friction modifiers (OFMs)
 #-An oil region composed of n-alkane chains
 
#########################################################################
################## CITING THIS SOFTWARE #################################
#########################################################################
#If you use this software, please cite the following article:

#J. P. Ewen, S. Echeverri Restrepo, N. Morgan, D. Dini, Nonequilibrium 
#Molecular Dynamics Simulations of Stearic Acid Adsorbed on Iron Surfaces
#with Nanoscale Roughness, Tribology International (2016), 
#http://dx.doi.org/10.1016/j.triboint.2016.11.039

#########################################################################
################## SOFTWARE NEEDED ######################################
#########################################################################
#What is needed:
#-Bash shell
#-Python v2.17.12 
#-moltemplate version v1.34 2015-11-18 : http://www.moltemplate.org/
#-ASE v3.11.0 : https://wiki.fysik.dtu.dk/ase/

#########################################################################
################## HOW TO RUN ###########################################
#########################################################################
#In a bash shell:

#-Go to the directory where the file run.py and the folder root are located
#-Modify the file run.py as needed
#-Type:
  #python run.py
#-The folder "lopls" is generated containing the files needed to run 
  #simulation with lammps. To run the simulation, type something like (note
  #that the name of the executable of lammps might be different):
  #lmp_g++ < in.lopls

#########################################################################
################## FILES NEEDED #########################################
#########################################################################

#run.py
  #This file contains all the information needed to define the dimensions 
    #and structure of the system. This is the only file that needs to be 
    #modified by the (regular) user
    
#root/AddEAM.py
  #This file adds the information needed to define the EAM interaction 
     #for the Fe atoms in the lammps input files

#root/lopls.py
  #This file generates all the input files needed by moltemplate 
    #(.lt extension) and calls it to generate the input files needed 
    #by lammps
    
#root/loplsMETAL.lt
  #This file contains the information related to the LOPLS force field 
    #that is used for the behaviour of the polymer chains and the OFMs

#root/Rough.py
  #This file contains a function to generate the rough Fe surfaces. The 
    #algorithm used is based on Tribol Lett 2011;44:279–85 and 
    #http://doi.org/10.1007/978-3-642-84574-1_34
  
#Fe_mm.eam.fs
  #This file is provided with the installation of lammps. It contains the 
    #EAM parameters that define the Fe interactions. It needs to be in 
    #the PATH accessible by lammps

#########################################################################
###################### SCHEMATIC ########################################
#########################################################################

#The generated system has the following schematic structure: 
 
#     ---------------------------------------------
#     ---------------------------------------------	 SURFACE TOP
#     ---------------------------------------------
#      O     O     O     O     O     O     O     O
#      |     |     |     |     |     |     |     |    OFM TOP
#      |     |     |     |     |     |     |     |
#      |     |     |     |     |     |     |     |    
#        _______________________________________     _________
#       | c-c-c-c-c-c-c    c-c-c-c-c-c-c        |  /\
#       |                                       |  |
#       | c-c-c-c-c-c-c    c-c-c-c-c-c-c        |  |
#       |                                       |  |
#       | c-c-c-c-c-c-c    c-c-c-c-c-c-c        |  |
#       |                                 OIL   |  |  Z SEPARATION
#       | c-c-c-c-c-c-c    c-c-c-c-c-c-c        |  |
#       |                                       |  |
#       | c-c-c-c-c-c-c    c-c-c-c-c-c-c        |  |
#       |                                       |  |
#       | c-c-c-c-c-c-c    c-c-c-c-c-c-c        |  |
#        _______________________________________  \/ _________
#      |     |     |     |     |     |     |     |
#      |     |     |     |     |     |     |     |
#      |     |     |     |     |     |     |     |    OFM BOTTOM
#      O     O     O     O     O     O     O     O
#     ---------------------------------------------
#     ---------------------------------------------   SURFACE BOTTOM
#     ---------------------------------------------
#
#	       z
#         |
#         |
#         |_ _ _ _ x
#        /
#       y   

#########################################################################
###################### INPUT VARIABLES ##################################
#########################################################################
#The only file that should be modified by a regular user is:
#run.py

#These are the variables that can be modified by the user. They relate to
#the structure of the system

#####General Inputs

#aFe		Lattice parameter of Fe in Angstrom

#xhi	 	Size of the simulation box in the x and y direction. These
#yhi		  values have to be given in terms of the lattice parameter
		  #of Fe. For example
		    #xhi = 15*aFe
		    #yhi = 15*aFe

#zhi 	 	This value determines the Z size of the section of the 
		  #box that  will contain the alkane chains and OFMs.
		  #The units are Angstrom

#####Inputs related to the OFM's

#OFM 		This flag determines if the the system will contain OFMs.
		  #if OFM = 1, the generation of OFMs is activated

#OFMtype 	Defines the type of OFM that will be used in the simulation
		  #the available options are: SA, SAm, GMS, OA, OAm, GMO
      #for full details and chemical structures, see:
      #Ewen et al. Langmiur (2016) DOI:10.1021/acs.langmuir.6b00586
		  #Exm: OFMtype = 'GMS'
      
#OFMn_x	 	These two variables determine the number of independent 
#OFMn_y		  OFM chains to be placed along the x and y directions. 
		  #The inputs are integers
		  #Exm: 	OFMn_x = 5
			#	      OFMn_y = 7

#####Inputs related to the Alkanes

#Alkane 		This flag determines if the the system will contain Alkanes.
		  #if Alkane = 1, the generation of Alkanes is activated
      #Note only linear n-alkanes currently available

#nAlkane		Defines the number of CH2/CH3 monomers on each alkane chain.
		  #The inputs are integers
		  #Exm: 	nAlkane = 16

#Alkanen_x 	Number of alkane chains to be placed along each 
#Alkanen_y	  direction of the simulation box
#Alkanen_z	  The inputs are integers
		  #Exm: 	Alkanen_x = 5
			#	      Alkanen_y = 7
			#	      Alkanen_z = 3

#####Inputs related to the Surfaces

#Surfaces = 1	This flag determines if the the system will contain Surfaces.
		  #if Surfaces = 1, the generation of Surfaces is activated

#FractalLevels 	Number of fractal levels used for the generation of the surfaces
		  #The inputs are integers
      #Exm: 	FractalLevels = 4
      
#H 		Hurst exponent for the generation of the fractal surfaces
      #Spijker et al. Tribology Letters (2011) DOI: 10.1007/s11249-011-9846-y
      #Exm: 	H = 0.8

#RMSin 		Value of the desired RMS roughness of the surfaces.
		  #The units are Angstrom
      #Exm: 	RMSin = 8

#boxLenghtZ	Thicknes of each of the Fe slabs that form the surfaces. This 
		  #is the thickness before the roughness is applied to the surface
		  #This value does not have units. It corresponds to the number 
		  #of Fe lattice constants
		  #Exm: boxLenghtZ = 20 -> 20*aFe = 57.27320 Angstrom
