#!/usr/bin/env python

#22-04-2025

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
from Fe2O3 import Fe2O3
from AddFe2O3 import AddFe2O3
from RoughFe2O3 import RoughFe2O3

def lopls(xlo,xhi,ylo,yhi,zlo,zhi,OFMn_x,OFMn_y,nAlkane, Alkanen_x,\
		Alkanen_y, Alkanen_z, Alkane, BZBZ, BZBZn_x, BZBZn_y, BZBZn_z,\
		Squalane, Squalanen_x, Squalanen_y, Squalanen_z, R123, R123n_x, R123n_y, R123n_z,\
    OFM ,OFMtype, Surfaces,\
		FractalLevels,RMSin,H,boxLenghtX,boxLenghtY,boxLenghtZ,aFe,Separation):

  f = open('lopls.lt','w+')

  #############################################################
  f.write('import "root/loplsaaMETAL.lt"')
  f.write("\n")

  if Surfaces == 1:

    f.write('import "WEA.lt"')

  if Surfaces == 2:

    f.write('import "Fe2O3.lt"')

  if Surfaces == 3:

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
  f.write('    $atom:HC $mol:... @atom:8500  0.00   -0.892430762954  -0.63104384422426  -3.000')
#  f.write('    $atom:HC $mol:... @atom:89  0.00   -0.892430762954  -0.63104384422426  -3.000')
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

  ###############################################################
  ####
  #Benzyl Benzoate
  f.write("BZBZ inherits LOPLSAA {\n")

  f.write("  # atomID      molID   atomType     charge   X       Y        Z\n")
  f.write("  write('Data Atoms') {\n")
  f.write("        $atom:C1        $mol:...        @atom:90        0.00    1.617   2.278    2.109  \n")
  f.write("        $atom:C2        $mol:...        @atom:90        0.00    2.453   2.282    0.996  \n")
  f.write("        $atom:C3        $mol:...        @atom:90        0.00    3.836   2.153    1.161  \n")
  f.write("        $atom:C4        $mol:...        @atom:90        0.00    4.380   2.026    2.443  \n")
  f.write("        $atom:C5        $mol:...        @atom:90        0.00    3.546   2.027    3.561  \n")
  f.write("        $atom:C6        $mol:...        @atom:842       0.00    2.158   2.153    3.398  \n")
  f.write("        $atom:C7        $mol:...        @atom:40600     0.00    1.215   2.158    4.553  \n")
  f.write("        $atom:C8        $mol:...        @atom:81        0.00    0.980   2.119    6.918  \n")
  f.write("        $atom:C9        $mol:...        @atom:90        0.00    2.283   3.365    8.682  \n")
  f.write("        $atom:C10       $mol:...        @atom:90        0.00    1.854   2.146    8.142  \n")
  f.write("        $atom:C11       $mol:...        @atom:90        0.00    2.261   0.951    8.751  \n")
  f.write("        $atom:C12       $mol:...        @atom:90        0.00    3.083   0.974    9.879  \n")
  f.write("        $atom:C13       $mol:...        @atom:90        0.00    3.507   2.195   10.410  \n")
  f.write("        $atom:C14       $mol:...        @atom:90        0.00    3.105   3.392    9.810  \n")
  f.write("        $atom:H1        $mol:...        @atom:91        0.00    0.541   2.374    2.000  \n")
  f.write("        $atom:H2        $mol:...        @atom:91        0.00    2.030   2.381    0.000  \n")
  f.write("        $atom:H3        $mol:...        @atom:91        0.00    4.490   2.152    0.293  \n")
  f.write("        $atom:H4        $mol:...        @atom:91        0.00    5.454   1.927    2.572  \n")
  f.write("        $atom:H5        $mol:...        @atom:91        0.00    3.965   1.933    4.557  \n")
  f.write("        $atom:H6        $mol:...        @atom:8500      0.00    0.349   1.226    6.888  \n")
  f.write("        $atom:H7        $mol:...        @atom:8500      0.00    0.332   2.997    6.858  \n")
  f.write("        $atom:H8        $mol:...        @atom:91        0.00    1.968   4.298    8.219  \n")
  f.write("        $atom:H9        $mol:...        @atom:91        0.00    1.928   0.000    8.341  \n")
  f.write("        $atom:H10       $mol:...        @atom:91        0.00    3.388   0.040   10.346  \n")
  f.write("        $atom:H11       $mol:...        @atom:91        0.00    4.143   2.214   11.292  \n")
  f.write("        $atom:H12       $mol:...        @atom:91        0.00    3.429   4.344   10.223  \n")
  f.write("        $atom:O1        $mol:...        @atom:40700     0.00    0.000   2.208    4.444  \n")
  f.write("        $atom:O2        $mol:...        @atom:40800     0.00    1.844   2.107    5.743  \n")
  f.write("        }\n")
  f.write(" write('Data Bond List') {\n")
  f.write(" # Aromatic Carbon Cycles\n")
  f.write("        $bond:CC1       $atom:C1        $atom:C2\n")
  f.write("        $bond:CC2       $atom:C2        $atom:C3\n")
  f.write("        $bond:CC3       $atom:C3        $atom:C4\n")
  f.write("        $bond:CC4       $atom:C4        $atom:C5\n")
  f.write("        $bond:CC5       $atom:C5        $atom:C6\n")
  f.write("        $bond:CC6       $atom:C6        $atom:C1\n")
  f.write("        $bond:CC7       $atom:C9        $atom:C10\n")
  f.write("        $bond:CC8       $atom:C10       $atom:C11\n")
  f.write("        $bond:CC9       $atom:C11       $atom:C12\n")
  f.write("        $bond:CC10      $atom:C12       $atom:C13\n")
  f.write("        $bond:CC11      $atom:C13       $atom:C14\n")
  f.write("        $bond:CC12      $atom:C14       $atom:C9\n")
  f.write("# Non Aromatic Carbons\n")
  f.write("        $bond:CC13      $atom:C7        $atom:C6\n")
  f.write("        $bond:CC14      $atom:C8        $atom:C10\n")
  f.write("# Hydrogen - Aromatic Carbon\n")
  f.write("        $bond:CH1       $atom:C1        $atom:H1\n")
  f.write("        $bond:CH2       $atom:C2        $atom:H2\n")
  f.write("        $bond:CH3       $atom:C3        $atom:H3\n")
  f.write("        $bond:CH4       $atom:C4        $atom:H4\n")
  f.write("        $bond:CH5       $atom:C5        $atom:H5\n")
  f.write("        $bond:CH6       $atom:C9        $atom:H8\n")
  f.write("        $bond:CH7       $atom:C11       $atom:H9\n")
  f.write("        $bond:CH8       $atom:C12       $atom:H10\n")
  f.write("        $bond:CH9       $atom:C13       $atom:H11\n")
  f.write("        $bond:CH10      $atom:C14       $atom:H12\n")
  f.write("# Hydrogen - Non aromatic Carbon\n")
  f.write("        $bond:CH11      $atom:C8        $atom:H6\n")
  f.write("        $bond:CH12      $atom:C8        $atom:H7\n")
  f.write("# Oxygen\n")
  f.write("        $bond:CO1   $atom:C7    $atom:O1\n")
  f.write("        $bond:CO2       $atom:C7        $atom:O2\n")
  f.write("        $bond:CO3       $atom:C8        $atom:O2\n")
  f.write("        }\n")
  f.write("}\n")


  ###############################################################
  ####
  #Squalane

  f.write("squalane inherits LOPLSAA {  \n")

  f.write("  # atomID      molID   atomType     charge   X       Y        Z\n")
  f.write("  write('Data Atoms') {\n")
  f.write("     $atom:H33       $mol:...        @atom:85        0.00    -13.888 -2.118  0.507\n")
  f.write("     $atom:C1        $mol:...        @atom:80        0.00    -13.557 -2.018  -0.555\n")
  f.write("     $atom:H31       $mol:...        @atom:85        0.00    -14.464 -1.954  -1.204\n")
  f.write("     $atom:H32       $mol:...        @atom:85        0.00    -12.986 -2.935  -0.835\n")
  f.write("     $atom:C2        $mol:...        @atom:8100      0.00    -12.698 -0.782  -0.729\n")
  f.write("     $atom:H34       $mol:...        @atom:8500      0.00    -12.351 -0.740  -1.801\n")
  f.write("     $atom:C3        $mol:...        @atom:80        0.00    -13.512 0.462   -0.441\n")
  f.write("     $atom:H35       $mol:...        @atom:85        0.00    -13.877 0.457   0.614\n")
  f.write("     $atom:H36       $mol:...        @atom:85        0.00    -12.896 1.380   -0.598\n")
  f.write("     $atom:H37       $mol:...        @atom:85        0.00    -14.397 0.510   -1.121\n")
  f.write("     $atom:C4        $mol:...        @atom:81        0.00    -11.472 -0.866  0.169\n")
  f.write("     $atom:H38       $mol:...        @atom:8500      0.00    -11.770 -0.645  1.228\n")
  f.write("     $atom:H39       $mol:...        @atom:8500      0.00    -11.073 -1.915  0.148\n")
  f.write("     $atom:C5        $mol:...        @atom:81        0.00    -10.377 0.088   -0.262\n")
  f.write("     $atom:H40       $mol:...        @atom:8500      0.00    -10.070 -0.144  -1.315\n")
  f.write("     $atom:H41       $mol:...        @atom:8500      0.00    -10.769 1.139   -0.254\n")
  f.write("     $atom:C6        $mol:...        @atom:81        0.00    -9.176  -0.018  0.654\n")
  f.write("     $atom:H42       $mol:...        @atom:8500      0.00    -9.473  0.252   1.701\n")
  f.write("     $atom:H43       $mol:...        @atom:8500      0.00    -8.829  -1.086  0.672\n")
  f.write("     $atom:C7        $mol:...        @atom:8100      0.00    -8.025  0.876   0.209\n")
  f.write("     $atom:H44       $mol:...        @atom:8500      0.00    -7.883  0.745   -0.902\n")
  f.write("     $atom:C8        $mol:...        @atom:80        0.00    -8.336  2.333   0.485\n")
  f.write("     $atom:H45       $mol:...        @atom:85        0.00    -8.488  2.503   1.578\n")
  f.write("     $atom:H46       $mol:...        @atom:85        0.00    -7.496  2.984   0.140\n")
  f.write("     $atom:H47       $mol:...        @atom:85        0.00    -9.263  2.644   -0.054\n")
  f.write("     $atom:C9        $mol:...        @atom:81        0.00    -6.741  0.452   0.910\n")
  f.write("     $atom:H48       $mol:...        @atom:8500      0.00    -6.777  0.778   1.982\n")
  f.write("     $atom:H49       $mol:...        @atom:8500      0.00    -6.671  -0.669  0.907\n")
  f.write("     $atom:C10       $mol:...        @atom:81        0.00    -5.505  1.020   0.245\n")
  f.write("     $atom:H50       $mol:...        @atom:8500      0.00    -5.489  0.723   -0.837\n")
  f.write("     $atom:H51       $mol:...        @atom:8500      0.00    -5.537  2.141   0.280\n")
  f.write("     $atom:C11       $mol:...        @atom:81        0.00    -4.247  0.523   0.926\n")
  f.write("     $atom:H52       $mol:...        @atom:8500      0.00    -4.258  0.825   2.006\n")
  f.write("     $atom:H53       $mol:...        @atom:8500      0.00    -4.238  -0.599  0.895\n")
  f.write("     $atom:C12       $mol:...        @atom:8100      0.00    -2.982  1.057   0.265\n")
  f.write("     $atom:H54       $mol:...        @atom:8500      0.00    -3.107  0.981   -0.853\n")
  f.write("     $atom:C13       $mol:...        @atom:80        0.00    -2.757  2.511   0.626\n")
  f.write("     $atom:H55       $mol:...        @atom:85        0.00    -2.590  2.622   1.725\n")
  f.write("     $atom:H56       $mol:...        @atom:85        0.00    -1.864  2.913   0.090\n")
  f.write("     $atom:H57       $mol:...        @atom:85        0.00    -3.644  3.126   0.340\n")
  f.write("     $atom:C14       $mol:...        @atom:81        0.00    -1.789  0.203   0.672\n")
  f.write("     $atom:H58       $mol:...        @atom:8500      0.00    -1.543  0.396   1.749\n")
  f.write("     $atom:H59       $mol:...        @atom:8500      0.00    -2.064  -0.882  0.584\n")
  f.write("     $atom:C15       $mol:...        @atom:81        0.00    -0.569  0.468   -0.184\n")
  f.write("     $atom:H60       $mol:...        @atom:8500      0.00    -0.830  0.332   -1.267\n")
  f.write("     $atom:H61       $mol:...        @atom:8500      0.00    -0.236  1.532   -0.055\n")
  f.write("     $atom:C16       $mol:...        @atom:81        0.00    0.568   -0.463  0.183\n")
  f.write("     $atom:H62       $mol:...        @atom:8500      0.00    0.830   -0.327  1.266\n")
  f.write("     $atom:H63       $mol:...        @atom:8500      0.00    0.236   -1.526  0.054\n")
  f.write("     $atom:C17       $mol:...        @atom:81        0.00    1.788   -0.197  -0.673\n")
  f.write("     $atom:H64       $mol:...        @atom:8500      0.00    1.542   -0.389  -1.751\n")
  f.write("     $atom:H65       $mol:...        @atom:8500      0.00    2.064   0.887   -0.584\n")
  f.write("     $atom:C18       $mol:...        @atom:8100      0.00    2.982   -1.053  -0.267\n")
  f.write("     $atom:H66       $mol:...        @atom:8500      0.00    3.105   -0.979  0.851\n")
  f.write("     $atom:C19       $mol:...        @atom:80        0.00    2.756   -2.505  -0.631\n")
  f.write("     $atom:H67       $mol:...        @atom:85        0.00    2.591   -2.615  -1.730\n")
  f.write("     $atom:H68       $mol:...        @atom:85        0.00    3.642   -3.122  -0.345\n")
  f.write("     $atom:H69       $mol:...        @atom:85        0.00    1.862   -2.908  -0.097\n")
  f.write("     $atom:C20       $mol:...        @atom:81        0.00    4.247   -0.518  -0.927\n")
  f.write("     $atom:H70       $mol:...        @atom:8500      0.00    4.257   -0.816  -2.008\n")
  f.write("     $atom:H71       $mol:...        @atom:8500      0.00    4.239   0.604   -0.892\n")
  f.write("     $atom:C21       $mol:...        @atom:81        0.00    5.505   -1.018  -0.248\n")
  f.write("     $atom:H72       $mol:...        @atom:8500      0.00    5.489   -0.725  0.836\n")
  f.write("     $atom:H73       $mol:...        @atom:8500      0.00    5.537   -2.138  -0.287\n")
  f.write("     $atom:C22       $mol:...        @atom:81        0.00    6.741   -0.447  -0.910\n")
  f.write("     $atom:H74       $mol:...        @atom:8500      0.00    6.777   -0.768  -1.984\n")
  f.write("     $atom:H75       $mol:...        @atom:8500      0.00    6.672   0.673   -0.901\n")
  f.write("     $atom:C23       $mol:...        @atom:8100      0.00    8.025   -0.877  -0.211\n")
  f.write("     $atom:H76       $mol:...        @atom:8500      0.00    7.883   -0.752  0.900\n")
  f.write("     $atom:C24       $mol:...        @atom:80        0.00    8.334   -2.332  -0.496\n")
  f.write("     $atom:H77       $mol:...        @atom:85        0.00    9.261   -2.647  0.041\n")
  f.write("     $atom:H78       $mol:...        @atom:85        0.00    7.494   -2.984  -0.156\n")
  f.write("     $atom:H79       $mol:...        @atom:85        0.00    8.487   -2.495  -1.590\n")
  f.write("     $atom:C25       $mol:...        @atom:81        0.00    9.176   0.020   -0.651\n")
  f.write("     $atom:H80       $mol:...        @atom:8500      0.00    9.472   -0.243  -1.700\n")
  f.write("     $atom:H81       $mol:...        @atom:8500      0.00    8.829   1.088   -0.662\n")
  f.write("     $atom:C26       $mol:...        @atom:81        0.00    10.378  -0.093  0.263\n")
  f.write("     $atom:H82       $mol:...        @atom:8500      0.00    10.071  0.131   1.319\n")
  f.write("     $atom:H83       $mol:...        @atom:8500      0.00    10.769  -1.144  0.248\n")
  f.write("     $atom:C27       $mol:...        @atom:81        0.00    11.473  0.863   -0.162\n")
  f.write("     $atom:H84       $mol:...        @atom:8500      0.00    11.770  0.650   -1.222\n")
  f.write("     $atom:H85       $mol:...        @atom:8500      0.00    11.074  1.912   -0.133\n")
  f.write("     $atom:C28       $mol:...        @atom:8100      0.00    12.699  0.773   0.735\n")
  f.write("     $atom:H86       $mol:...        @atom:8500      0.00    12.353  0.723   1.807\n")
  f.write("     $atom:C29       $mol:...        @atom:80        0.00    13.513  -0.470  0.438\n")
  f.write("     $atom:H87       $mol:...        @atom:85        0.00    14.398  -0.523  1.116\n")
  f.write("     $atom:H88       $mol:...        @atom:85        0.00    12.897  -1.388  0.587\n")
  f.write("     $atom:H89       $mol:...        @atom:85        0.00    13.877  -0.457  -0.618\n")
  f.write("     $atom:C30       $mol:...        @atom:80        0.00    13.559  2.009   0.570\n")
  f.write("     $atom:H90       $mol:...        @atom:85        0.00    12.988  2.925   0.857\n")
  f.write("     $atom:H91       $mol:...        @atom:85        0.00    14.466  1.940   1.217\n")
  f.write("     $atom:H92       $mol:...        @atom:85        0.00    13.889  2.117   -0.492\n")
  f.write("     }\n")

  f.write(" write('Data Bond List') {\n")
  f.write("     $bond:B1        $atom:H33       $atom:C1\n")
  f.write("     $bond:B2        $atom:C1        $atom:H31\n")
  f.write("     $bond:B3        $atom:C1        $atom:H32\n")
  f.write("     $bond:B4        $atom:C1        $atom:C2\n")
  f.write("     $bond:B5        $atom:C2        $atom:H34\n")
  f.write("     $bond:B6        $atom:C2        $atom:C3\n")
  f.write("     $bond:B7        $atom:C2        $atom:C4\n")
  f.write("     $bond:B8        $atom:C3        $atom:H35\n")
  f.write("     $bond:B9        $atom:C3        $atom:H36\n")
  f.write("     $bond:B10       $atom:C3        $atom:H37\n")
  f.write("     $bond:B11       $atom:C4        $atom:H38\n")
  f.write("     $bond:B12       $atom:C4        $atom:H39\n")
  f.write("     $bond:B13       $atom:C4        $atom:C5\n")
  f.write("     $bond:B14       $atom:C5        $atom:H40\n")
  f.write("     $bond:B15       $atom:C5        $atom:H41\n")
  f.write("     $bond:B16       $atom:C5        $atom:C6\n")
  f.write("     $bond:B17       $atom:C6        $atom:H42\n")
  f.write("     $bond:B18       $atom:C6        $atom:H43\n")
  f.write("     $bond:B19       $atom:C6        $atom:C7\n")
  f.write("     $bond:B20       $atom:C7        $atom:H44\n")
  f.write("     $bond:B21       $atom:C7        $atom:C8\n")
  f.write("     $bond:B22       $atom:C7        $atom:C9\n")
  f.write("     $bond:B23       $atom:C8        $atom:H45\n")
  f.write("     $bond:B24       $atom:C8        $atom:H46\n")
  f.write("     $bond:B25       $atom:C8        $atom:H47\n")
  f.write("     $bond:B26       $atom:C9        $atom:H48\n")
  f.write("     $bond:B27       $atom:C9        $atom:H49\n")
  f.write("     $bond:B28       $atom:C9        $atom:C10\n")
  f.write("     $bond:B29       $atom:C10       $atom:H50\n")
  f.write("     $bond:B30       $atom:C10       $atom:H51\n")
  f.write("     $bond:B31       $atom:C10       $atom:C11\n")
  f.write("     $bond:B32       $atom:C11       $atom:H52\n")
  f.write("     $bond:B33       $atom:C11       $atom:H53\n")
  f.write("     $bond:B34       $atom:C11       $atom:C12\n")
  f.write("     $bond:B35       $atom:C12       $atom:H54\n")
  f.write("     $bond:B36       $atom:C12       $atom:C13\n")
  f.write("     $bond:B37       $atom:C12       $atom:C14\n")
  f.write("     $bond:B38       $atom:C13       $atom:H55\n")
  f.write("     $bond:B39       $atom:C13       $atom:H56\n")
  f.write("     $bond:B40       $atom:C13       $atom:H57\n")
  f.write("     $bond:B41       $atom:C14       $atom:H58\n")
  f.write("     $bond:B42       $atom:C14       $atom:H59\n")
  f.write("     $bond:B43       $atom:C14       $atom:C15\n")
  f.write("     $bond:B44       $atom:C15       $atom:H60\n")
  f.write("     $bond:B45       $atom:C15       $atom:H61\n")
  f.write("     $bond:B46       $atom:C15       $atom:C16\n")
  f.write("     $bond:B47       $atom:C16       $atom:H62\n")
  f.write("     $bond:B48       $atom:C16       $atom:H63\n")
  f.write("     $bond:B49       $atom:C16       $atom:C17\n")
  f.write("     $bond:B50       $atom:C17       $atom:H64\n")
  f.write("     $bond:B51       $atom:C17       $atom:H65\n")
  f.write("     $bond:B52       $atom:C17       $atom:C18\n")
  f.write("     $bond:B53       $atom:C18       $atom:H66\n")
  f.write("     $bond:B54       $atom:C18       $atom:C19\n")
  f.write("     $bond:B55       $atom:C18       $atom:C20\n")
  f.write("     $bond:B56       $atom:C19       $atom:H67\n")
  f.write("     $bond:B57       $atom:C19       $atom:H68\n")
  f.write("     $bond:B58       $atom:C19       $atom:H69\n")
  f.write("     $bond:B59       $atom:C20       $atom:H70\n")
  f.write("     $bond:B60       $atom:C20       $atom:H71\n")
  f.write("     $bond:B61       $atom:C20       $atom:C21\n")
  f.write("     $bond:B62       $atom:C21       $atom:H72\n")
  f.write("     $bond:B63       $atom:C21       $atom:H73\n")
  f.write("     $bond:B64       $atom:C21       $atom:C22\n")
  f.write("     $bond:B65       $atom:C22       $atom:H74\n")
  f.write("     $bond:B66       $atom:C22       $atom:H75\n")
  f.write("     $bond:B67       $atom:C22       $atom:C23\n")
  f.write("     $bond:B68       $atom:C23       $atom:H76\n")
  f.write("     $bond:B69       $atom:C23       $atom:C24\n")
  f.write("     $bond:B70       $atom:C23       $atom:C25\n")
  f.write("     $bond:B71       $atom:C24       $atom:H77\n")
  f.write("     $bond:B72       $atom:C24       $atom:H78\n")
  f.write("     $bond:B73       $atom:C24       $atom:H79\n")
  f.write("     $bond:B74       $atom:C25       $atom:H80\n")
  f.write("     $bond:B75       $atom:C25       $atom:H81\n")
  f.write("     $bond:B76       $atom:C25       $atom:C26\n")
  f.write("     $bond:B77       $atom:C26       $atom:H82\n")
  f.write("     $bond:B78       $atom:C26       $atom:H83\n")
  f.write("     $bond:B79       $atom:C26       $atom:C27\n")
  f.write("     $bond:B80       $atom:C27       $atom:H84\n")
  f.write("     $bond:B81       $atom:C27       $atom:H85\n")
  f.write("     $bond:B82       $atom:C27       $atom:C28\n")
  f.write("     $bond:B83       $atom:C28       $atom:H86\n")
  f.write("     $bond:B84       $atom:C28       $atom:C29\n")
  f.write("     $bond:B85       $atom:C28       $atom:C30\n")
  f.write("     $bond:B86       $atom:C29       $atom:H87\n")
  f.write("     $bond:B87       $atom:C29       $atom:H88\n")
  f.write("     $bond:B88       $atom:C29       $atom:H89\n")
  f.write("     $bond:B89       $atom:C30       $atom:H90\n")
  f.write("     $bond:B90       $atom:C30       $atom:H91\n")
  f.write("     $bond:B91       $atom:C30       $atom:H92\n")
  f.write("     }\n")
  f.write("}\n")

  ###############################################################
  ####
  #R123
  f.write("R123 inherits LOPLSAA {\n")
  
  f.write("  # atomID      molID   atomType     charge   X       Y        Z\n")
  f.write("  write('Data Atoms') {\n")
  f.write("        $atom:C00        $mol:...        @atom:80000        0.00    1.000000    1.000000    0.000000  \n")
  f.write("        $atom:C01        $mol:...        @atom:80100        0.00   -0.536000    1.000000    0.000000  \n")
  f.write("        $atom:F02        $mol:...        @atom:80200        0.00   -1.065000    1.000000    1.247000  \n")
  f.write("        $atom:F03        $mol:...        @atom:80300        0.00   -1.003000    2.123000   -0.608000  \n")
  f.write("        $atom:F04        $mol:...        @atom:80400        0.00   -1.077000   -0.051000   -0.652000  \n")
  f.write("        $atom:Cl0        $mol:...        @atom:80500        0.00    1.711000    0.698000    1.615000  \n")
  f.write("        $atom:Cl1        $mol:...        @atom:80600        0.00    1.701000   -0.178000   -1.155000  \n")
  f.write("        $atom:H07        $mol:...        @atom:80700        0.00    1.361000    1.984000   -0.311000  \n")
  f.write("        }\n")

  f.write(" write('Data Bond List') {\n")
  f.write("        $bond:C01-C00       $atom:C01        $atom:C00\n")
  f.write("        $bond:F02-C01       $atom:F02        $atom:C01\n")
  f.write("        $bond:F03-C01       $atom:F03        $atom:C01\n")
  f.write("        $bond:F04-C01       $atom:F04        $atom:C01\n")
  f.write("        $bond:Cl0-C00       $atom:Cl0        $atom:C00\n")
  f.write("        $bond:Cl1-C00       $atom:Cl1        $atom:C00\n")
  f.write("        $bond:H07-C00       $atom:H07        $atom:C00\n")
  f.write("        }\n")
  
  f.write("}\n")


  ######################################################################
  #rough Iron surfaces
  if Surfaces == 1:


    Rough(FractalLevels,RMSin,H,boxLenghtX,boxLenghtY,boxLenghtZ,aFe,Separation)


  ######################################################################
  #Flat Fe2O3 surfaces

  if Surfaces == 2:


    Fe2O3(FractalLevels,RMSin,H,boxLenghtX,boxLenghtY,boxLenghtZ,aFe,Separation)


  if Surfaces == 3:

    RoughFe2O3(FractalLevels,RMSin,H,boxLenghtX,boxLenghtY,boxLenghtZ,aFe,Separation)

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
  if Alkanen_z == 1:
    Alkanes_z = 0.0
  else:
    Alkanes_z = ((zhi-23.3065-4)-(zlo+23.3065+4))/(Alkanen_z-1)


  ####
  #BZBZ
  # This determines how far apart all BZBZ molecules will be placed
  BZBZ_x = (xhi-xlo)/BZBZn_x #(1.2533223*(nAlkane-1))+5
  BZBZ_y = (yhi-ylo)/BZBZn_y
  if BZBZn_z == 1:
    BZBZ_z = 0.0
  else:
    BZBZ_z = ((zhi-23.3065-5)-(zlo+23.3065+5))/(BZBZn_z-1)

  ####
  #Squalane
  # This determines how far apart all Squalane molecules will be placed
  Squalane_x = (xhi-xlo)/Squalanen_x #(1.2533223*(nAlkane-1))+5
  Squalane_y = (yhi-ylo)/Squalanen_y
  if Squalanen_z == 1:
    Squalane_z = 0.0
  else:
    Squalane_z = ((zhi-23.3065-5)-(zlo+23.3065+5))/(Squalanen_z-1)


  ####
  #R123
  # This determines how far apart all R123 molecules will be placed
  R123_x = (xhi-xlo)/R123n_x #(1.2533223*(nAlkane-1))+5
  R123_y = (yhi-ylo)/R123n_y
  if Surfaces == 0:
    if R123n_z == 1:
      R123_z = 0.0
    else:
      R123_z = (zhi)/(R123n_z)
  else:

    if R123n_z == 1:
      R123_z = 0.0
    else:
      R123_z = ((zhi-23.3065-5)-(zlo+23.3065+5))/(R123n_z-1)


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

  #no surfaces (making a triclinic box)
  if Surfaces == 0 :
    f.write(str(zlo)+"  "+str(zhi))
    f.write("  zlo zhi")
    f.write("\n")
    f.write("0.0 0.0 0.0 xy xz yz")
    f.write("\n")
    f.write("}")
  #rough Iron surfaces
  if Surfaces == 1 :
    f.write(str(zlo-boxLenghtZ*aFe-20)+"  "+str(zhi+boxLenghtZ*aFe+20))
    f.write("  zlo zhi")
    f.write("\n")
    f.write("}")
  #Flat Fe2O3 surfaces
  if Surfaces == 2 :
    f.write(str(zlo-boxLenghtZ*13.730-20)+"  "+str(zhi+boxLenghtZ*13.730+20))
    f.write("  zlo zhi")
    f.write("\n")
    f.write("}")
  if Surfaces == 3 :
    f.write(str(zlo-boxLenghtZ*13.730-20)+"  "+str(zhi+boxLenghtZ*13.730+20))
#    f.write(str(zlo-boxLenghtZ*aFe-20)+"  "+str(zhi+boxLenghtZ*aFe+20))
    f.write("  zlo zhi")
    f.write("\n")
    f.write("}")



  #####
  #The OFMs

  if OFM == 1 and Surfaces!=3:

    # Here the OFMpolymers are placed, using the number of OFMpolymers in each direction and the set distance
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

  if Alkane == 1 and Surfaces!=3:
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
  #BZBZ
  if BZBZ == 1 and Surfaces!=3:
    f.write("\n")
    f.write("\n")
    f.write("molecules4 = new BZBZ.rot(90, 0, 1, 0) [")
    f.write(str(BZBZn_z))
    f.write("].move(0, 0,")
    f.write(str(BZBZ_z))
    f.write(")")
    f.write("\n")
    f.write("                           [")
    f.write(str(BZBZn_y))
    f.write("].move(0, ")
    f.write(str(BZBZ_y))
    f.write(", 0)")
    f.write("\n")
    f.write("                           [")
    f.write(str(BZBZn_x))
    f.write("].move(")
    f.write(str(BZBZ_x))
    f.write(", 0, 0)")
    f.write("\n")

    f.write("molecules4[*][*][*].move("+str(xlo)+","+str(ylo)+","+str(zlo+23.3065+8)+")")

  #####
  #Squalane

  if Squalane == 1 and Surfaces!=3:
    f.write("\n")
    f.write("\n")
    f.write("molecules6 = new squalane.move(15, 0, 0) [")
    f.write(str(Squalanen_z))
    f.write("].move(0, 0,")
    f.write(str(Squalane_z))
    f.write(")")
    f.write("\n")
    f.write("                           [")
    f.write(str(Squalanen_y))
    f.write("].move(0, ")
    f.write(str(Squalane_y))
    f.write(", 0)")
    f.write("\n")
    f.write("                           [")
    f.write(str(Squalanen_x))
    f.write("].move(")
    f.write(str(Squalane_x))
    f.write(", 0, 0)")
    f.write("\n")

    f.write("molecules6[*][*][*].move("+str(xlo)+","+str(ylo)+","+str(zlo+23.3065+4)+")")

  ######
  #R123
  if R123 == 1 and Surfaces!=3:
    f.write("\n")
    f.write("\n")
    f.write("molecules7 = new R123.rot(0, 0, 1, 0) [")
    f.write(str(R123n_z))
    f.write("].move(0, 0,")
    f.write(str(R123_z))
    f.write(")")
    f.write("\n")
    f.write("                           [")
    f.write(str(R123n_y))
    f.write("].move(0, ")
    f.write(str(R123_y))
    f.write(", 0)")
    f.write("\n")
    f.write("                           [")
    f.write(str(R123n_x))
    f.write("].move(")
    f.write(str(R123_x))
    f.write(", 0, 0)")
    f.write("\n")

    if Surfaces == 0:
      f.write("molecules7[*][*][*].move("+str(xlo)+","+str(ylo)+","+str(zlo)+")")
    else:
      f.write("molecules7[*][*][*].move("+str(xlo)+","+str(ylo)+","+str(zlo+23.3065+8)+")")



  ######
  #The Surfaces
  if Surfaces == 1:

    f.write("\n")
    f.write("\n")
    f.write("molecules5 = new FESurface.move("+str(xlo)+","+str(ylo)+","+str(zlo-boxLenghtZ*aFe-1)+")")
    f.write("\n")
    f.write("\n")

  if Surfaces == 2:

    f.write("\n")
    f.write("\n")
    f.write("molecules5 = new FESurface.move("+str(xlo)+","+str(ylo)+","+str(zlo-boxLenghtZ*13.730-1)+")")
    f.write("\n")
    f.write("\n")

  if Surfaces == 3:

    f.write("\n")
    f.write("\n")
    f.write("molecules5 = new FESurface.move("+str(xlo)+","+str(ylo)+","+str(zlo-boxLenghtZ*13.730-1)+")")
#    f.write("molecules5 = new FESurface.move("+str(xlo)+","+str(ylo)+","+str(zlo-boxLenghtZ*aFe-1)+")")
    f.write("\n")
    f.write("\n")

  #f.write("molecules4 = new FESurface.move("+str(xlo)+","+str(ylo)+","+str(zlo-boxLenghtZ*aFe-1)+")")


  f.close()

  # Creates the name for the all the imput files
  name = 'lopls' #+str(N_total)
  os.rename('lopls.lt',name+'.lt')

  #runs Moltemplate
  os.system('moltemplate.sh '+name+'.lt')
  #os.system('moltemplate_2016-12-18.sh '+name+'.lt')
  #os.system('moltemplate_2017-2-10.sh '+name+'.lt')

  # Builds the .in file
  ########################################################################
  f = open('in.'+name,'w+')
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
  f.write("dump            dump1 all atom 1000 "+name+".dump")
  f.write("\n")
  f.write("thermo_style    custom step lx ly lz  density temp press etotal")
  f.write("\n")
  f.write("thermo          1")
  f.write("\n")
  f.write("write_data		"+ name +"Initial.data")
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

  if Surfaces == 2:
    AddFe2O3(name)
    os.system('rm Fe2O3.lt')

  if Surfaces == 3:
    AddEAM()
    os.system('rm WEA.lt')

  # Moves all files to a seperate folder 
  os.system('rm -r lopls')
  os.system('mkdir lopls')
  os.system('rm -r output_ttree')
  os.system('rm '+name+'.in')
  os.system('mv '+name+'.in.init lopls')
  os.system('mv '+name+'.in.CreateBonds lopls')
  os.system('mv '+name+'.in.settings lopls')
  os.system('mv '+name+'.in.charges lopls')
  os.system('rm '+name+'.lt')
  os.system('mv '+name+'.data lopls')

  # moves all new input files to the folder
  os.system("mv in."+name+" lopls")
  
  #if Surfaces=3 make a rough surface of Fe2O3, ONLY the surface. The rest is trash
  if Surfaces == 3:
    os.system('rm lopls/in.lopls lopls/lopls.in.charges lopls/lopls.in.init lopls/lopls.in.settings')

  #copy the eam potential file in case that Surfaces == 1
  if Surfaces == 1:
    os.system('cp root/Fe_mm.eam.fs lopls')

