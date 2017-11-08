##Version: 08-11-2017

##Authors - Sebastián ECHEVERRI RESTREPO: sebastian.echeverri.restrepo@skf.com, sebastianecheverrir@gmail.com; James EWEN: j.ewen@imperial.ac.uk, jimmyewen@gmail.com

## GENERAL INFO

This software is suitable as a starting point to performing confined nonequilibrium molecular dynamics (NEMD) simulations of OFM films adsorbed to iron surfaces, separated by a layer of n-alkane molecules.

This software generates a LAMMPS datafile and basic input file for systems containing*:
 
 - Two bcc a-Fe or a-Fe2O3 slabs with/without nanoscale RMS roughness
 - Two organic friction modifier (OFM) monolayers above/below bottom/top slabs
 - An optional central linear alkane region
 
*Note that any of these components can be excluded by using the appropriate flags

## SCHEMATIC

The generated system has the following schematic structure: 
 
     ---------------------------------------------
     ---------------------------------------------	 SURFACE TOP
     ---------------------------------------------
     \_/```\_/```\_______/```\___/``\___/\___/\__/
      O     O     O     O     O     O     O     O
      |     |     |     |     |     |     |     |    OFM TOP
      |     |     |     |     |     |     |     |
      |     |     |     |     |     |     |     |    
        _______________________________________     _________
       | c-c-c-c-c-c-c    c-c-c-c-c-c-c        |  /\
       |                                       |  |
       | c-c-c-c-c-c-c    c-c-c-c-c-c-c        |  |
       |                                       |  |
       | c-c-c-c-c-c-c    c-c-c-c-c-c-c        |  |
       |                                ALKANE |  |  Z SEPARATION
       | c-c-c-c-c-c-c    c-c-c-c-c-c-c        |  |
       |                                       |  |
       | c-c-c-c-c-c-c    c-c-c-c-c-c-c        |  |
       |                                       |  |
       | c-c-c-c-c-c-c    c-c-c-c-c-c-c        |  |
        _______________________________________  \/ _________
      |     |     |     |     |     |     |    |
      |     |     |     |     |     |     |    |
      |     |     |     |     |     |     |    |    OFM BOTTOM
      O     O     O     O     O     O     O    O
     / \___/`\__/``\__/```\__/`\___/```````\__/`\_
     ---------------------------------------------
     ---------------------------------------------   SURFACE BOTTOM
     ---------------------------------------------

         z
         |
         |
         |_ _ _ _ x
        /
       y   

 
## CITING THIS SOFTWARE
If you use this software, please cite the following article: J. P. Ewen, S. Echeverri Restrepo, N. Morgan, D. Dini, Nonequilibrium Molecular Dynamics Simulations of Stearic Acid Adsorbed on Iron Surfaces with Nanoscale Roughness, Tribology International (2016), http://dx.doi.org/10.1016/j.triboint.2016.11.039

## SOFTWARE REQUIREMENTS
Required (all GNU):
 - Bash shell
 - Python v2.17.12 
 - moltemplate version v1.34 2015-11-18 : http://www.moltemplate.org/
 - ASE v3.11.0 : https://wiki.fysik.dtu.dk/ase/
 
Optional (all GNU):
 - LAMMPS 22 Jul 2016-ICMS : http://lammps.sandia.gov/ (to run MD simulation of resultant system)
 - VMD v1.9.1 : http://www.ks.uiuc.edu/Research/vmd/ (to visualise resultant datafile)
 
Note that versions given are those for which software has been tested, the code should still work with most older/newer versions.

## HOW TO RUN
In a bash shell:
 - Go to the directory where the file run.py and the folder root are located
 - Modify the file run.py as needed (see below)
 - To run the script, type:   $ python run.py
 - The folder "lopls" is generated, which contains the files needed to run a basic LAMMPS simulation
 - Check the resultant LAMMPS datafile with VMD -> Extensions, Tk Console, type: $ topo readlammpsdata lopls.data full
 - To test the LAMMPS simulation, type*:  $ mpirun -np X lmp_mpi < in.lopls       , where X is the number of processors
 - Extend/modify the in.lopls file in order to perform more complex MD simulations
   
*note that the name of the LAMMPS executable might be different

## FILES NEEDED
The below files are required to generate the LAMMPS datafile and input file, they are all provided in the distribution:

 - run.py

This file contains all the information needed to define the dimensions and structure of the system. This is the only file that needs to be modified by the (regular) user
    
 - root/AddEAM.py

This file adds the information needed to define the EAM interaction for the Fe atoms in the LAMMPS input files

 - root/lopls.py

This file generates all the input files needed by moltemplate (.lt extension) and calls it to generate the input files needed by LAMMPS
    
 - root/loplsMETAL.lt

This file contains the information related to the L-OPLS force field that is used to simulate the n-alkane chains and the OFMs. Ref: Siu et al. Journal of Chemical Theory and Computation (2012) http://pubs.acs.org/doi/abs/10.1021/ct200908r

 - root/Rough.py

This file contains a function to generate the rough Fe surfaces. For details of the RMD algorithm used, see Refs: Spijker et al. Tribology Letters (2011) http://link.springer.com/article/10.1007/s11249-011-9846-y and Voss in Fundamental Algorithms for Computer Graphics (1991) http://doi.org/10.1007/978-3-642-84574-1_34
  
 - Fe_mm.eam.fs
 
This file is provided with the installation of lammps. It contains the EAM parameters that define the Fe interactions. It needs to be in the PATH accessible by LAMMPS. Ref: Mendelev et al. Philosophical Magazine (2003) http://www.tandfonline.com/doi/abs/10.1080/14786430310001613264

## INPUT VARIABLES

The only file which needs to be modified by a regular user is:

 - run.py

## General Inputs

These are the variables that can be modified by the user, they relate to the structure of the system

 - aFe		

Lattice parameter of Fe in Angstrom

 - xhi
 - yhi

Size of the simulation box in the x and y direction. These values need to be given in terms of the lattice parameter of Fe (aFe).
 
For example: xhi = 15, yhi = 15

 - zhi

This value determines the Z size of the section of the box that will contain both the alkane chains and OFMs. The units are Angstrom.

For example: zhi = 60

## Inputs related to the OFMs

 - OFM 	

This flag determines if the the system will contain OFMs. If OFM = 1, the generation of OFMs is activated

 - OFMtype 	

Defines the type of OFM that will be used in the simulation the available options are: SA, SAm, GMS, OA, OAm, GMO. For full details and chemical structures, see: Ewen et al. Langmiur (2016) http://pubs.acs.org/doi/full/10.1021/acs.langmuir.6b00586
 
For example: OFMtype = GMS
      
 - OFMn_x
 - OFMn_y

These two variables determine the number of independent OFM chains to be placed along the x and y directions. The inputs are integers.

For example: OFMn_x = 5, OFMn_y = 7

## Inputs related to the Alkanes

 - Alkane 		
 
This flag determines if the the system will contain Alkanes. If Alkane = 1, the generation of Alkanes is activated. Note only linear n-alkanes are currently available.

 - nAlkane

Defines the number of CH2/CH3 monomers on each alkane chain. The inputs are integers.

For example: 	nAlkane = 16
 - Alkanen_x
 - Alkanen_y
 - Alkanen_z

Number of alkane chains to be placed along each direction of the simulation box. The inputs are integers

For example: 	Alkanen_x = 5, Alkanen_y = 7, Alkanen_z = 3

## Inputs related to the Surfaces

 - Surfaces

This flag determines if the the system will contain Surfaces. If Surfaces = 1, the generation of Surfaces is activated

 - FractalLevels 	

Number of fractal levels used for the generation of the surfaces. The inputs are integers.

For example: 	FractalLevels = 4
      
 - H

Hurst exponent for the generation of the fractal surfaces. Ref: Spijker et al. Tribology Letters (2011) http://link.springer.com/article/10.1007/s11249-011-9846-y

For example: 	H = 0.8

 - RMSin 

Value of the desired RMS roughness of the surfaces. The units are Angstrom.

For example:   RMSin = 8

 - boxLenghtZ

Thicknes of each of the Fe slabs that form the surfaces. This is the thickness before the roughness is applied to the surface This value corresponds to the number of Fe lattice constants (aFe)

For example:  boxLenghtZ = 20

## OUTPUT FILES

After running run.py, you should get the following files in the lopls folder, which can be used to run a LAMMPS simulation.

 - lopls/in.lopls
 
Basic LAMMPS input files for minimisation of lopls.data. We recommend changing this file to perform more complex simulations e.g. NEMD

 - lopls/lopls.data
 
LAMMPS datafile for slab-OFM-alkane-OFM-slab system.

 - lopls/lopls.in.charges
 
LAMMPS partial charges for L-OPLS atom types

 - lopls/lopls.in.init
 
LAMMPS units, atom_style, bond_style, angle_style, dihedral_style, pair_style, special_bonds, kspace_style etc.

 - lopls/lopls.in.settings
 
LAMMPS pair_coeff, bond_coeff, angle_coeff, dihedral_coeff etc.

## Notes

The script could be modified to simulate any confined system in LAMMPS using a range of force-fields. Please contact the authors if you would like help making these modifications.
