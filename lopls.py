#!/usr/bin/env python

#29-11-2016

#Authors:Sebastian ECHEVERRI RESTREPO,   
#	 	sebastian.echeverri.restrepo@skf.com, sebastianecheverrir@gmail.com
#	 James EWEN
#		j.ewen14@imperial.ac.uk, jimmyewen@gmail.com

#################################################################################3

#  This file generates all the input files needed by moltemplate 
#    (.lt extension) and calls it to generate the input files needed 
#    by lammps

#################################################################################3

import os
import sys

sys.path.append("root")

from Rough import Rough
from AddEAM import AddEAM



def lopls(xlo,xhi,ylo,yhi,zlo,zhi,OFMn_x,OFMn_y,nAlkane, Alkanen_x,\
		Alkanen_y, Alkanen_z, Alkane, OFM,OFMtype, Surfaces,\
		FractalLevels,RMSin,H,boxLenghtX,boxLenghtY,boxLenghtZ,aFe,Separation):

  f = open('lopls.lt','wr+')

  #############################################################
  f.write('import "root/loplsaaMETAL.lt"')
  f.write("\n")

  if Surfaces == 1:
  
    f.write('import "WEA.lt"')
    
  ## This part builds the basic CH, CH2, CH3, COOH, CONH2, RCOOR
  
  f.write("\n")
  f.write("\n")
  f.write('CH inherits LOPLSAA {')
  f.write("\n")
  f.write('  write("Data Atoms") {')
  f.write("\n")
  f.write('    $atom:C  $mol:... @atom:87  0.00   0.000  0.000             0.000')
  f.write("\n")
  f.write('    $atom:H $mol:... @atom:89  0.00    0.0 1.0 0.000')
  f.write("\n")
  f.write('  }')
  f.write("\n")
  f.write("\n")
  f.write('  write(\'Data Bond List\') {')
  f.write("\n")
  f.write('    $bond:CH $atom:C $atom:H')
  f.write("\n")
  f.write('}')
  f.write("\n")
  f.write('  } #CH')



  f.write("\n")
  f.write("\n")
  f.write('CH2 inherits LOPLSAA {')
  f.write("\n")
  f.write('  write("Data Atoms") {')
  f.write("\n")
  f.write('    $atom:C  $mol:... @atom:81  0.00   0.000  0.000             0.000')
  f.write("\n")
  f.write('    $atom:H1 $mol:... @atom:8500  0.00   0.892430762954  0.63104384422426  0.000')
  f.write("\n")
  f.write('    $atom:H2 $mol:... @atom:8500  0.00  	-0.892430762954  0.63104384422426 -0.000')
  f.write("\n")
  f.write('  }')
  f.write("\n")
  f.write("\n")
  f.write('  write(\'Data Bond List\') {')
  f.write("\n")
  f.write('    $bond:CH1 $atom:C $atom:H1')
  f.write("\n")
  f.write('    $bond:CH2 $atom:C $atom:H2')
  f.write("\n")
  f.write('}')
  f.write("\n")
  f.write('  } #CH2')


  f.write("\n")
  f.write("\n")
  f.write('CH3 inherits LOPLSAA {')
  f.write("\n")
  f.write('  write("Data Atoms") {')
  f.write("\n")
  f.write('    $atom:C  $mol:... @atom:80  0.00   0.000000    0.000000          0.000000')
  f.write("\n")
  f.write('    $atom:H1 $mol:... @atom:85  0.00   0.8924307629540046   0.6310438442242609 0.000000')
  f.write("\n")
  f.write('    $atom:H2 $mol:... @atom:85  0.00   -0.8924307629540046   0.6310438442242609 -0.000000')
  f.write("\n")
  f.write('    $atom:H3 $mol:... @atom:85  0.00  -0.000000 -0.6310438442242609 -0.8924307629540046')
  f.write("\n")
  f.write('  }')
  f.write("\n")
  f.write("\n")
  f.write('  write(\'Data Bond List\') {')
  f.write("\n")
  f.write('    $bond:CH1 $atom:C $atom:H1')
  f.write("\n")
  f.write('    $bond:CH2 $atom:C $atom:H2')
  f.write("\n")
  f.write('    $bond:CH3 $atom:C $atom:H3')
  f.write("\n")
  f.write('  }')
  f.write("\n")
  f.write("\n")
  f.write('} # CH3')


  f.write("\n")
  f.write("\n")
  f.write('COOH inherits LOPLSAA {')
  f.write("\n")
  f.write('  write("Data Atoms") {')
  f.write("\n")
  f.write('    $atom:C  $mol:... @atom:209  0.00   0.000000    0.000000          0.000000')
  f.write("\n")
  f.write('    $atom:O $mol:... @atom:210  0.00   0.00000   1.000000 0.000000')
  f.write("\n")
  f.write('    $atom:OH $mol:... @atom:211  0.00   -0.000000   0.000000 -1.000000')
  f.write("\n")
  f.write('    $atom:HO $mol:... @atom:212  0.00  -0.000000 0.400000 -2.00000')
  f.write("\n")
  f.write('  }')
  f.write("\n")
  f.write("\n")
  f.write('  write(\'Data Bond List\') {')
  f.write("\n")
  f.write('    $bond:CO $atom:C $atom:O')
  f.write("\n")
  f.write('    $bond:COH $atom:C $atom:OH')
  f.write("\n")
  f.write('    $bond:OHHO $atom:OH $atom:HO')
  f.write("\n")
  f.write('  }')
  f.write("\n")
  f.write("\n")
  f.write('} # COOH')



  f.write("\n")
  f.write("\n")
  f.write('CONH2 inherits LOPLSAA {')
  f.write("\n")
  f.write('  write("Data Atoms") {')
  f.write("\n")
  f.write('    $atom:C  $mol:... @atom:177  0.00   0.000000    0.000000          0.000000')
  f.write("\n")
  f.write('    $atom:O $mol:... @atom:178  0.00   0.00000   1.000000 0.000000')
  f.write("\n")
  f.write('    $atom:N $mol:... @atom:179  0.00   -0.000000   0.000000 -1.000000')
  f.write("\n")
  f.write('    $atom:H1 $mol:... @atom:182  0.00  -0.000000 -1.00000 -1.00000')
  f.write("\n")
  f.write('    $atom:H2 $mol:... @atom:182  0.00  -0.000000 0.000000 -2.00000')
  f.write("\n")
  f.write('  }')
  f.write("\n")
  f.write("\n")
  f.write('  write(\'Data Bond List\') {')
  f.write("\n")
  f.write('    $bond:CO $atom:C $atom:O')
  f.write("\n")
  f.write('    $bond:CN $atom:C $atom:N')
  f.write("\n")
  f.write('    $bond:NH1 $atom:N $atom:H1')
  f.write("\n")
  f.write('    $bond:NH2 $atom:N $atom:H2')
  f.write("\n")
  f.write('  }')
  f.write("\n")
  f.write("\n")
  f.write('} # CONH2')


  f.write("\n")
  f.write("\n")
  f.write('RCOOR inherits LOPLSAA {')
  f.write("\n")
  f.write('  write("Data Atoms") {')
  f.write("\n")
  f.write('    $atom:C  $mol:... @atom:406  0.00   0.000000    0.000000          0.000000')
  f.write("\n")
  f.write('    $atom:O $mol:... @atom:407  0.00   0.00000   1.000000 0.000000')
  f.write("\n")
  f.write('    $atom:OS $mol:... @atom:408  0.00   -0.000000   0.000000 -1.000000')
  f.write("\n")
  f.write('    $atom:CT $mol:... @atom:409  0.00  -0.000000 0.00000 -2.00000')
  f.write("\n")
  f.write('    $atom:H1 $mol:... @atom:8500  0.00   0.892430762954  0.63104384422426  -2.000')
  f.write("\n")
  f.write('    $atom:H2 $mol:... @atom:8500  0.00  	-0.892430762954  0.63104384422426 -2.000')
  f.write("\n")
  f.write('    $atom:CT2 $mol:... @atom:100  0.00  -0.000000 0.00000 -3.00000')
  f.write("\n")
  f.write('    $atom:HC $mol:... @atom:89  0.00   -0.892430762954  -0.63104384422426  -3.000')
  f.write("\n")
  f.write('    $atom:OH $mol:... @atom:96  0.00  	0.892430762954  -0.63104384422426 -3.000')
  f.write("\n")
  f.write('    $atom:HO $mol:... @atom:97  0.00  	1.6 -1.2 -3.000')
  f.write("\n")
  f.write('    $atom:CT3 $mol:... @atom:99  0.00  -0.000000 0.00000 -4.00000')
  f.write("\n")
  f.write('    $atom:H3 $mol:... @atom:8500  0.00   0.892430762954  0.63104384422426  -4.000')
  f.write("\n")
  f.write('    $atom:H4 $mol:... @atom:8500  0.00  	-0.892430762954  0.63104384422426 -4.000')
  f.write("\n")
  f.write('    $atom:OH2 $mol:... @atom:96  0.00  -0.000000 0.00000 -5.00000')
  f.write("\n")
  f.write('    $atom:HO2 $mol:... @atom:97  0.00  -0.000000 0.00000 -6.00000')
  f.write("\n")
  f.write('  }')
  f.write("\n")
  f.write("\n")
  f.write('  write(\'Data Bond List\') {')
  f.write("\n")
  f.write('    $bond:CO $atom:C $atom:O')
  f.write("\n")
  f.write('    $bond:COS $atom:C $atom:OS')
  f.write("\n")
  f.write('    $bond:OSCT $atom:OS $atom:CT')
  f.write("\n")
  f.write('    $bond:CTH1 $atom:CT $atom:H1')
  f.write("\n")
  f.write('    $bond:CTH2 $atom:CT $atom:H2')
  f.write("\n")
  f.write('    $bond:CTCT2 $atom:CT $atom:CT2')
  f.write("\n")
  f.write('    $bond:CT2HC $atom:CT2 $atom:HC')
  f.write("\n")
  f.write('    $bond:CT2OH $atom:CT2 $atom:OH')
  f.write("\n")
  f.write('    $bond:OHHO $atom:OH $atom:HO')
  f.write("\n")
  f.write('    $bond:CT2CT3 $atom:CT2 $atom:CT3')
  f.write("\n")
  f.write('    $bond:CT3H3 $atom:CT3 $atom:H3')
  f.write("\n")
  f.write('    $bond:CT3H4 $atom:CT3 $atom:H4')
  f.write("\n")
  f.write('    $bond:CT3OH2 $atom:CT3 $atom:OH2')
  f.write("\n")
  f.write('    $bond:OH2HO2 $atom:OH2 $atom:HO2')
  f.write("\n")
  f.write('  }')
  f.write("\n")
  f.write("\n")
  f.write('} # RCOOR')



  ######################################################################


  # This part ensembles the CH2 and CH3 into polymers, import is the data on line 95 which makes sure the final molecule is in the right direction

  #the OFM SA
  f.write("\n")
  f.write("\n")
  f.write('SA inherits LOPLSAA {')
  f.write("\n")
  f.write("\n")
  f.write('  create_var {$mol} ')
  f.write("\n")
  f.write('  OFMpolymers = new CH2 [')
  n_OFMpolymers = str(18)
  f.write(n_OFMpolymers)
  f.write('].rot(180,0,0,1).move(0,0,1.2533223)')
  f.write("\n")
  f.write('  delete OFMpolymers[0]')
  f.write("\n")
  f.write('  delete OFMpolymers[')
  OFMn_1= (17)
  f.write(str(OFMn_1))
  f.write("]")
  f.write("\n")
  f.write('  OFMpolymers[0] = new CH3')
  f.write("\n")
  f.write('  OFMpolymers[')
  f.write(str(OFMn_1))
  f.write('] = new COOH ')
  f.write("\n")
  f.write('  OFMpolymers[')
  f.write(str(OFMn_1))
  f.write("].rot(180,1,")
  f.write('0')
  f.write(',0).move(0,0,')
  OFMd_n_1 = (OFMn_1)*1.2533223
  f.write(str(OFMd_n_1))
  f.write(')')
  f.write("\n") 
  f.write("\n")
  f.write("  write('Data Bond List') {")
  f.write("\n")

  for a in range (0, (17)):
	  f.write("    $bond:b")	
	  f.write(str(a+1))
	  f.write("  $atom:OFMpolymers[")
	  f.write(str(a))
	  f.write("]/C $atom:OFMpolymers[")
	  f.write(str(a+1))
	  f.write("]/C")
	  f.write("\n")


  f.write("\n")
  f.write("  }")
  f.write("\n")
  f.write("} # SA")
  f.write("\n")


###############################################################
#the OFM SAm
  f.write("\n")
  f.write("\n")
  f.write('SAm inherits LOPLSAA {')
  f.write("\n")
  f.write("\n")
  f.write('  create_var {$mol} ')
  f.write("\n")
  f.write('  OFMpolymers = new CH2 [')
  n_OFMpolymers = str(18)
  f.write(n_OFMpolymers)
  f.write('].rot(180,0,0,1).move(0,0,1.2533223)')
  f.write("\n")
  f.write('  delete OFMpolymers[0]')
  f.write("\n")
  f.write('  delete OFMpolymers[')
  OFMn_1= (17)
  f.write(str(OFMn_1))
  f.write("]")
  f.write("\n")
  f.write('  OFMpolymers[0] = new CH3')
  f.write("\n")
  f.write('  OFMpolymers[')
  f.write(str(OFMn_1))
  f.write('] = new CONH2 ')
  f.write("\n")
  f.write('  OFMpolymers[')
  f.write(str(OFMn_1))
  f.write("].rot(180,1,")
  f.write('0')
  f.write(',0).move(0,0,')
  OFMd_n_1 = (OFMn_1)*1.2533223
  f.write(str(OFMd_n_1))
  f.write(')')
  f.write("\n") 
  f.write("\n")
  f.write("  write('Data Bond List') {")
  f.write("\n")

  for a in range (0, (17)):
	  f.write("    $bond:b")	
	  f.write(str(a+1))
	  f.write("  $atom:OFMpolymers[")
	  f.write(str(a))
	  f.write("]/C $atom:OFMpolymers[")
	  f.write(str(a+1))
	  f.write("]/C")
	  f.write("\n")


  f.write("\n")
  f.write("  }")
  f.write("\n")
  f.write("} # SAm")
  f.write("\n")
###############################################################
#the OFM GMS
  f.write("\n")
  f.write("\n")
  f.write('GMS inherits LOPLSAA {')
  f.write("\n")
  f.write("\n")
  f.write('  create_var {$mol} ')
  f.write("\n")
  f.write('  OFMpolymers = new CH2 [')
  n_OFMpolymers = str(18)
  f.write(n_OFMpolymers)
  f.write('].rot(180,0,0,1).move(0,0,1.2533223)')
  f.write("\n")
  f.write('  delete OFMpolymers[0]')
  f.write("\n")
  f.write('  delete OFMpolymers[')
  OFMn_1= (17)
  f.write(str(OFMn_1))
  f.write("]")
  f.write("\n")
  f.write('  OFMpolymers[0] = new CH3')
  f.write("\n")
  f.write('  OFMpolymers[')
  f.write(str(OFMn_1))
  f.write('] = new RCOOR ')
  f.write("\n")
  f.write('  OFMpolymers[')
  f.write(str(OFMn_1))
  f.write("].rot(180,1,")
  f.write('0')
  f.write(',0).move(0,0,')
  OFMd_n_1 = (OFMn_1)*1.2533223
  f.write(str(OFMd_n_1))
  f.write(')')
  f.write("\n") 
  f.write("\n")
  f.write("  write('Data Bond List') {")
  f.write("\n")

  for a in range (0, (17)):
	  f.write("    $bond:b")	
	  f.write(str(a+1))
	  f.write("  $atom:OFMpolymers[")
	  f.write(str(a))
	  f.write("]/C $atom:OFMpolymers[")
	  f.write(str(a+1))
	  f.write("]/C")
	  f.write("\n")


  f.write("\n")
  f.write("  }")
  f.write("\n")
  f.write("} # GMS")
  f.write("\n")

###############################################################
  
  #the OFM OA
  f.write("\n")
  f.write("\n")
  f.write('OA inherits LOPLSAA {')
  f.write("\n")
  f.write("\n")
  f.write('  create_var {$mol} ')
  f.write("\n")
  f.write('  OFMpolymers = new CH2 [')
  n_OFMpolymers = str(8)
  f.write(n_OFMpolymers)
  f.write('].rot(180,0,0,1).move(0,0,1.2533223)')
  f.write("\n")
  f.write('  delete OFMpolymers[0]')
  f.write("\n")
  f.write('  OFMpolymers[0] = new CH3')
  f.write("\n")
  f.write('  OFMpolymers[8] = new CH')
  f.write("\n")
  f.write('  OFMpolymers[8].rot(45,1,0,0).move(0,0,10.0)')
  f.write("\n")
  f.write('  OFMpolymers[9] = new CH')
  f.write("\n")
  f.write('  OFMpolymers[9].rot(45,1,0,0).move(0,0,10.0).move(0,-0.8,0.8)')
  f.write("\n")
  
  f.write('  OFMpolymers2 = new CH2 [7].rot(180,0,0,1).move(0,0,1.2533223)')
  f.write("\n")
  f.write('  OFMpolymers2[].rot(180,0,0,1).rot(45,1,0,0).move(0,0,10.0).move(0,-1.6,1.6)')
  f.write("\n")
  f.write('  OFMpolymers2[7] = new COOH')
  f.write("\n")
  f.write('  OFMpolymers2[7].rot(-180,0,0,1).rot(-125,1,0,0).move(0,0,10.0).move(0,-7.8,7.8)')
  f.write("\n")



  f.write("\n")
  f.write("  write('Data Bond List') {")
  f.write("\n")

  for a in range (0, (9)):
	  f.write("    $bond:b")	
	  f.write(str(a+1))
	  f.write("  $atom:OFMpolymers[")
	  f.write(str(a))
	  f.write("]/C $atom:OFMpolymers[")
	  f.write(str(a+1))
	  f.write("]/C")
	  f.write("\n")
  f.write(' $bond:bc  $atom:OFMpolymers[9]/C $atom:OFMpolymers2[0]/C ')
  f.write("\n")
	  
  for a in range (0, (7)):
	  f.write("    $bond:c")	
	  f.write(str(a+1))
	  f.write("  $atom:OFMpolymers2[")
	  f.write(str(a))
	  f.write("]/C $atom:OFMpolymers2[")
	  f.write(str(a+1))
	  f.write("]/C")
	  f.write("\n")


  f.write("\n")
  f.write("  }")
  f.write("\n")
  f.write("} # OA")
  f.write("\n")


###############################################################
   #the OFM OAm
  f.write("\n")
  f.write("\n")
  f.write('OAm inherits LOPLSAA {')
  f.write("\n")
  f.write("\n")
  f.write('  create_var {$mol} ')
  f.write("\n")
  f.write('  OFMpolymers = new CH2 [')
  n_OFMpolymers = str(8)
  f.write(n_OFMpolymers)
  f.write('].rot(180,0,0,1).move(0,0,1.2533223)')
  f.write("\n")
  f.write('  delete OFMpolymers[0]')
  f.write("\n")
  f.write('  OFMpolymers[0] = new CH3')
  f.write("\n")
  f.write('  OFMpolymers[8] = new CH')
  f.write("\n")
  f.write('  OFMpolymers[8].rot(45,1,0,0).move(0,0,10.0)')
  f.write("\n")
  f.write('  OFMpolymers[9] = new CH')
  f.write("\n")
  f.write('  OFMpolymers[9].rot(45,1,0,0).move(0,0,10.0).move(0,-0.8,0.8)')
  f.write("\n")
  
  f.write('  OFMpolymers2 = new CH2 [7].rot(180,0,0,1).move(0,0,1.2533223)')
  f.write("\n")
  f.write('  OFMpolymers2[].rot(180,0,0,1).rot(45,1,0,0).move(0,0,10.0).move(0,-1.6,1.6)')
  f.write("\n")
  f.write('  OFMpolymers2[7] = new CONH2')
  f.write("\n")
  f.write('  OFMpolymers2[7].rot(-180,0,0,1).rot(-125,1,0,0).move(0,0,10.0).move(0,-7.8,7.8)')
  f.write("\n")



  f.write("\n")
  f.write("  write('Data Bond List') {")
  f.write("\n")

  for a in range (0, (9)):
	  f.write("    $bond:b")	
	  f.write(str(a+1))
	  f.write("  $atom:OFMpolymers[")
	  f.write(str(a))
	  f.write("]/C $atom:OFMpolymers[")
	  f.write(str(a+1))
	  f.write("]/C")
	  f.write("\n")
  f.write(' $bond:bc  $atom:OFMpolymers[9]/C $atom:OFMpolymers2[0]/C ')
  f.write("\n")
	  
  for a in range (0, (7)):
	  f.write("    $bond:c")	
	  f.write(str(a+1))
	  f.write("  $atom:OFMpolymers2[")
	  f.write(str(a))
	  f.write("]/C $atom:OFMpolymers2[")
	  f.write(str(a+1))
	  f.write("]/C")
	  f.write("\n")


  f.write("\n")
  f.write("  }")
  f.write("\n")
  f.write("} # OAm")
  f.write("\n")


###############################################################
 
###############################################################
   #the OFM GMO
  f.write("\n")
  f.write("\n")
  f.write('GMO inherits LOPLSAA {')
  f.write("\n")
  f.write("\n")
  f.write('  create_var {$mol} ')
  f.write("\n")
  f.write('  OFMpolymers = new CH2 [')
  n_OFMpolymers = str(8)
  f.write(n_OFMpolymers)
  f.write('].rot(180,0,0,1).move(0,0,1.2533223)')
  f.write("\n")
  f.write('  delete OFMpolymers[0]')
  f.write("\n")
  f.write('  OFMpolymers[0] = new CH3')
  f.write("\n")
  f.write('  OFMpolymers[8] = new CH')
  f.write("\n")
  f.write('  OFMpolymers[8].rot(45,1,0,0).move(0,0,10.0)')
  f.write("\n")
  f.write('  OFMpolymers[9] = new CH')
  f.write("\n")
  f.write('  OFMpolymers[9].rot(45,1,0,0).move(0,0,10.0).move(0,-0.8,0.8)')
  f.write("\n")
  
  f.write('  OFMpolymers2 = new CH2 [7].rot(180,0,0,1).move(0,0,1.2533223)')
  f.write("\n")
  f.write('  OFMpolymers2[].rot(180,0,0,1).rot(45,1,0,0).move(0,0,10.0).move(0,-1.6,1.6)')
  f.write("\n")
  f.write('  OFMpolymers2[7] = new RCOOR')
  f.write("\n")
  f.write('  OFMpolymers2[7].rot(-180,0,0,1).rot(-125,1,0,0).move(0,0,10.0).move(0,-7.8,7.8)')
  f.write("\n")



  f.write("\n")
  f.write("  write('Data Bond List') {")
  f.write("\n")

  for a in range (0, (9)):
	  f.write("    $bond:b")	
	  f.write(str(a+1))
	  f.write("  $atom:OFMpolymers[")
	  f.write(str(a))
	  f.write("]/C $atom:OFMpolymers[")
	  f.write(str(a+1))
	  f.write("]/C")
	  f.write("\n")
  f.write(' $bond:bc  $atom:OFMpolymers[9]/C $atom:OFMpolymers2[0]/C ')
  f.write("\n")
	  
  for a in range (0, (7)):
	  f.write("    $bond:c")	
	  f.write(str(a+1))
	  f.write("  $atom:OFMpolymers2[")
	  f.write(str(a))
	  f.write("]/C $atom:OFMpolymers2[")
	  f.write(str(a+1))
	  f.write("]/C")
	  f.write("\n")


  f.write("\n")
  f.write("  }")
  f.write("\n")
  f.write("} # GMO")
  f.write("\n")


###############################################################
  
###############################################################


  ####
  #The AlKane
  f.write("\n")
  f.write("\n")
  f.write('Hexadecane inherits LOPLSAA {')
  f.write("\n")
  f.write("\n")
  f.write('  create_var {$mol} ')
  f.write("\n")
  f.write('  AlkanePolymer = new CH2 [')
  n_AlkanePolymer = str(nAlkane)
  f.write(n_AlkanePolymer)
  f.write('].rot(180,0,0,1).move(0,0,1.2533223)')
  f.write("\n")
  f.write('  delete AlkanePolymer[0]')
  f.write("\n")
  f.write('  delete AlkanePolymer[')
  Alkanen_1= (nAlkane-1)
  f.write(str(Alkanen_1))
  f.write("]")
  f.write("\n")
  f.write('  AlkanePolymer[0] = new CH3')
  f.write("\n")
  f.write('  AlkanePolymer[')
  f.write(str(Alkanen_1))
  f.write('] = new CH3 ')
  f.write("\n")
  f.write('  AlkanePolymer[')
  f.write(str(Alkanen_1))
  f.write("]")

  if nAlkane%2==0:
	  f.write('.rot(180,0,0,1).rot(180,0,1,0).')
  else:
	  f.write('.rot(180,0,1,0).')
  f.write('move(0,0,')
  Alkaned_n_1 = (Alkanen_1)*1.2533223

  f.write(str(Alkaned_n_1))

  f.write(')')
  f.write("\n") 
  f.write("\n")
  f.write("  write('Data Bond List') {")
  f.write("\n")

  for a in range (0, (nAlkane-1)):
	  f.write("    $bond:b")	
	  f.write(str(a+1))
	  f.write("  $atom:AlkanePolymer[")
	  f.write(str(a))
	  f.write("]/C $atom:AlkanePolymer[")
	  f.write(str(a+1))
	  f.write("]/C")
	  f.write("\n")


  f.write("\n")
  f.write("  }")
  f.write("\n")
  f.write("} # Hexadecane")
  f.write("\n")

  ######################################################################
  #rough Iron surfaces
  if Surfaces == 1:


    Rough(FractalLevels,RMSin,H,boxLenghtX,boxLenghtY,boxLenghtZ,aFe,Separation)




  ######################################################################


  # This calculates the amount of molecules in each direction and the total amount of molecules



  ####
  #The OFM

  #OFMn_x = int((xhi-xlo)/OFMs_x)
  #OFMn_y = int((yhi-ylo)/OFMs_y)
  OFMn_z = 1 #int((zhi-zlo)/OFMs_z)
  #N_total = OFMn_x*OFMn_y*OFMn_z

  # This determines how far apart all OFMpolymers will be placed
  OFMs_x = (xhi-xlo)/OFMn_x
  OFMs_y = (yhi-ylo)/OFMn_y
  OFMs_z = 0 #(1.5*n)+3


  ####
  #The AlKane
  # This determines how far apart all Alkanes polymers will be placed
  Alkanes_x = (xhi-xlo)/Alkanen_x #(1.2533223*(nAlkane-1))+5
  Alkanes_y = (yhi-ylo)/Alkanen_y
  Alkanes_z = ((zhi-23.3065-4)-(zlo+23.3065+4))/(Alkanen_z-1)


  ######################################################################
  #Placing the polymers in the box


  f.write("\n")
  f.write("# Periodic boundary conditions:")
  f.write("\n")
  f.write("write_once(\"Data Boundary\") {")
  f.write("\n")
  f.write(str(xlo)+"  "+str(xhi))
  f.write("  xlo xhi")
  f.write("\n")

  f.write(str(ylo)+"  "+str(yhi))
  f.write("  ylo yhi")
  f.write("\n")

  f.write(str(zlo-boxLenghtZ*aFe-20)+"  "+str(zhi+boxLenghtZ*aFe+20))
  f.write("  zlo zhi")
  f.write("\n")
  f.write("}")





  #####
  #The OFMs

  if OFM == 1:

    # Here the OFMpolymers are placed, using the numer of OFMpolymers in each direction and the set distance
    f.write("\n")
    f.write("\n")
    if OFMtype == 'SA':
      f.write("molecules = new SA [")
    elif OFMtype == 'SAm':
      f.write("molecules = new SAm [")
    elif OFMtype == 'GMS':
      f.write("molecules = new GMS [")
    elif OFMtype == 'OA':
      f.write("molecules = new OA [")
    elif OFMtype == 'OAm':
      f.write("molecules = new OAm [")
    elif OFMtype == 'GMO':
      f.write("molecules = new GMO [")
    f.write(str(OFMn_z))
    f.write("].move(0, 0,") 
    f.write(str(OFMs_z))
    f.write(")")
    f.write("\n")
    f.write("                           [")
    f.write(str(OFMn_y))
    f.write("].move(0, ")
    f.write(str(OFMs_y))
    f.write(", 0)")
    f.write("\n")
    f.write("                           [")
    f.write(str(OFMn_x))
    f.write("].move(")
    f.write(str(OFMs_x))
    f.write(", 0, 0)")
    f.write("\n")

    f.write("molecules[*][*][*].rot(180,1,0,0).move("+str(xlo)+","+str(ylo+(OFMn_y-1)*OFMs_y)+","+str(zlo+23.3065)+")")


    f.write("\n")
    f.write("\n")
    if OFMtype == 'SA':
      f.write("molecules2 = new SA [")
    elif OFMtype == 'SAm':
      f.write("molecules2 = new SAm [")
    elif OFMtype == 'GMS':
      f.write("molecules2 = new GMS [")
    elif OFMtype == 'OA':
      f.write("molecules2 = new OA [")
    elif OFMtype == 'OAm':
      f.write("molecules2 = new OAm [")
    elif OFMtype == 'GMO':
      f.write("molecules2 = new GMO [")
    f.write(str(OFMn_z))
    f.write("].move(0, 0,") 
    f.write(str(OFMs_z))
    f.write(")")
    f.write("\n")
    f.write("                           [")
    f.write(str(OFMn_y))
    f.write("].move(0, ")
    f.write(str(OFMs_y))
    f.write(", 0)")
    f.write("\n")
    f.write("                           [")
    f.write(str(OFMn_x))
    f.write("].move(")
    f.write(str(OFMs_x))
    f.write(", 0, 0)")
    f.write("\n")

    f.write("molecules2[*][*][*].move("+str(xlo)+","+str(ylo)+","+str(zlo+(zhi-zlo)-23.3065)+")")

  #####
  #The Alkanes

  if Alkane == 1:
    f.write("\n")
    f.write("\n")
    f.write("molecules3 = new Hexadecane.rot(90, 0, 1, 0) [")
    f.write(str(Alkanen_z))
    f.write("].move(0, 0,") 
    f.write(str(Alkanes_z))
    f.write(")")
    f.write("\n")
    f.write("                           [")
    f.write(str(Alkanen_y))
    f.write("].move(0, ")
    f.write(str(Alkanes_y))
    f.write(", 0)")
    f.write("\n")
    f.write("                           [")
    f.write(str(Alkanen_x))
    f.write("].move(")
    f.write(str(Alkanes_x))
    f.write(", 0, 0)")
    f.write("\n")

    f.write("molecules3[*][*][*].move("+str(xlo)+","+str(ylo)+","+str(zlo+23.3065+4)+")")



  ######
  #The Surfaces
  if Surfaces == 1:

    f.write("\n")
    f.write("\n")
    f.write("molecules4 = new FESurface.move("+str(xlo)+","+str(ylo)+","+str(zlo-boxLenghtZ*aFe-1)+")")
    f.write("\n")
    f.write("\n")


  #f.write("molecules4 = new FESurface.move("+str(xlo)+","+str(ylo)+","+str(zlo-boxLenghtZ*aFe-1)+")")


  f.close()

  # Creates the name for the all the imput files
  name = 'lopls' #+str(N_total)
  os.rename('lopls.lt',name+'.lt')

  #runs Moltemplate
  os.system('moltemplate.sh '+name+'.lt')


  # Builds the .in file
  ########################################################################
  f = open('in.'+name,'wr+')
  f.write("# ------------------------------- Initialization Section --------------------")
  f.write("\n")
  f.write("include         	"+name+".in.init")
  f.write("\n")
  f.write("read_data        	" + name +".data")
  f.write("\n")
  f.write("include         	"+name+".in.settings")
  f.write("\n")
  f.write("include         	"+name+".in.charges")
  f.write("\n")
  f.write("\n")
  f.write("dump            dump1 all atom 1 "+name+".dump")
  f.write("\n")
  f.write("thermo_style    custom step lx ly lz  density temp press etotal")
  f.write("\n")
  f.write("thermo          1")
  f.write("\n")
  f.write("\n")
  f.write("# ------------------Run Equilibriation ---------------------------")
  f.write("\n")
  f.write("\n")
  f.write("min_style       cg")
  f.write("\n")
  f.write("minimize        0.0 0.0 100000 100000")
  f.write("\n")
  f.write("\n")

  f.close()

  ############################################################################################################
  if Surfaces == 1:

    AddEAM()
    os.system('rm WEA.lt')

  # Moves all files to a seperate folder 
  os.system('rm -r lopls')
  os.system('mkdir lopls')
  os.system('rm -r output_ttree')
  os.system('rm '+name+'.in')
  os.system('mv '+name+'.in.init lopls')
  os.system('mv '+name+'.in.settings lopls')
  os.system('mv '+name+'.in.charges lopls')
  os.system('rm '+name+'.lt')
  os.system('mv '+name+'.data lopls')

  # moves all new input files to the folder
  os.system("mv in."+name+" lopls")

