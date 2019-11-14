# python calc_scat_from_previous_calc.py  bamlet15_80D2O  bamlet15_100D2O 200 0.8 A -2.05 17202 B -4.83 10485 
# python calc_scat_from_previous_calc.py  bamlet15_XRAY   bamlet15_100D2O 200 0.8 A 3.007 17202 B -0.32 10485 

# python calc_scat_from_previous_calc.py [output scattering file name prefix] [input_scattering_file_prefix] [approx_num_points_in_output] [q_max] [chain-id-1] [contrast-1] [volume-1] [chain-id-2] [contrast-2] [volume-2] ... 

# This program reads in the scattering curves for each scattering component and cross-term between scattering components, 
# and reads in the contrasts and volumes, calculates the new scattering. 
# A previous program calculated the scattering curves, and it may have taken a while for that to be calculated. 
# This program allows the contrast to be varied and scattering to be recalculated with the varied contrast 
# without having to rerun the timeconsuming calculating for component scattering. 

# This program reads in : 
#	[input_scattering_file_prefix]_qI_AA.dat 
#	[input_scattering_file_prefix]_qI_BB.dat 
#	[input_scattering_file_prefix]_qI_AB.dat 
# etc. 
#	[input_scattering_file_prefix]_qI_CC.dat 
#	[input_scattering_file_prefix]_qI_AC.dat 
#	[input_scattering_file_prefix]_qI_BC.dat 
# etc. 
# according to how many chains are provided on the input line. 

# The previous python program started of as a python version of a Steen Hanson MCsimul program (written in FORTRAN) 
# that was download from http://www.matfys.kvl.dk/~steen/MCsimul.zip May 2012. 
# It was changed to read input in the PDB file format. 
# It was extended to process multiple densities according to Whitten and Trewhella 2008. 
# The Guinier processing in this python program has not yet been tested to ensure that it is still correct after the extensions made to the python program. 

# The previous program has read in a file in pdb format, containing points that represent scattering points. 
# For x-ray scattering, the points represent electrons. 
# The previous program program reads the co-ordinates of the points, 
# not the atom type or any other information in the pdb format file that might be used for display. 
# The input is not a pdb file describing atomic co-ordinates, 
# but can conveniently be used by pdb programs to view the scattering points as if they were atoms. 
# The pdb format does not have a field for contrast (or scattering length density). 
# This field will be read from columns 81-88. If the field is not filled in the input pdb file, 
# then it will be defaulted to the value of 1 for that scattering point. 

# The user enters the contrasts and volumes of the components on the command line. 
# Any contrasts or volumes not entered will be set to 1 (which could be reasonable for contrast and probably not reasonable for volume). 

# The previous program outputs the total P(r) and the P(r) for each component. 
# The previous program outputs the total expected scattering, and the scattering for each component, 
# at the relative levels that they contribute to total scattering and normalised. 

# This program does include contrast (scattering length density) and volume in the equations, 
# so that the output for each Q-point is the final scattering according to the following equation 
# that includes the form factors and structure factor(s) : 
# To get the final scattering for a two component system, the equation calculated by this program is, for each Q-point : 
# 
# 	Itotal = (p1*p1 * V1*V1 * P11) + (p2*p2 * V2*V2 * P22) + (2 * p1*V1*p2*V2 * P12) 
# 
# where : 
#	p1	= contrast of component 1 
#	V1	= volume of component 1 
#	P11	= form factor for component 1 
#	p2	= contrast of component 2 
#	V2	= volume of component 2 
#	P22	= form factor for component 2 
#	P12	= structure factor for the inter-component scattering 
# 
# That formula was taken from : 
# "MULCh: modules for the analysis of small-angle neutron contrast variation data from biomolecular assemblies" 
# Whitten AE, Cai S, Trewhella J (2008) J Appl Cryst 41, 222-226. 

# Other papers that have implemented the expected scattering from molecules are : 
# CRYSOL : Svergun D, Barberato C, Koch MHJ (1995) J Appl Cryst 28:768-773, and 
# FOXS : Schneidman-Duhovny D, Hammel M, Sali A (2010) Nucleic Acids Research 38 Web Server issue: W540-W544. 
# These programs don't seem to handle negative contrast 
# which is macromolecule scattering length density that is less than the solvent scattering length density that the macromolecule is in. 
# This python program does handle that, as shown by the fact that this program gives the same expected scattering result 
# as does the classic core shell form factor equation, when the core or shell have negative contrast. 
# However, this python program does not account for exclude volume 
# (volume of solvent that is not contributing to scattering due to the fact that the macromolecule is occupying that volume - 
# which affects the scattering at high Q of Q > 0.25 ang^-1). 
# CRYSOL and FOXS do account for excluded volume. 

# The output files that are always produced are : 
#	*_qI.dat		output scattering file 
#	*_Pr.dat		output pair distribution file 
#	*_info.txt		output information file 
# The chain field in the pdb formatted input file marks the components of the scattering particle. 
# For each component, a file is outputted containing the form factor for that component, 
# which is the scattering from pairs of scattering points within the component. 
# A file is also outputted for each component with each other component, 
# containing the cross scattering from pairs of scattering points such that one point is in one component 
# and the second point is in the other component. 
# All the scattering points within a given component should have the same contrast (or scattering length density) 
# as each other. However, the program doesn't check for this. 
# The output files for form factors and cross scattering structure factors are : 
#	*_qI_11.dat		scattering for form factor for component 1 
#	*_qI_22.dat		scattering for form factor for component 2 
#	*_qI_12.dat		scattering for structure factor for cross scattering between components 1 and 2 
#	*_qI_33.dat		scattering for form factor for component 3 
#	*_qI_13.dat		scattering for structure factor for cross scattering between components 1 and 3 
#	*_qI_23.dat		scattering for structure factor for cross scattering between components 2 and 3 
#	etc. 
#	*_Pr_11.dat		pair distribution function for form factor for component 1 
#	*_Pr_22.dat		pair distribution function for form factor for component 2 
#	*_Pr_12.dat		pair distribution function for structure factor for cross scattering between components 1 and 2 
#	*_Pr_33.dat		pair distribution function for form factor for component 3 
#	*_Pr_13.dat		pair distribution function for structure factor for cross scattering between components 1 and 3 
#	*_Pr_23.dat		pair distribution function for structure factor for cross scattering between components 2 and 3 
#	etc. 

# Example of input pdb format : 
# ATOM      1  N   GLU A   1      20.657  31.671 -10.811  1.00 53.57           N  
# ATOM      2  CA  GLU A   1      20.644  33.109 -10.368  1.00 55.43           C  
# ATOM      3  C   GLU A   1      21.812  33.850 -11.031  1.00 53.43           C  
# ATOM      4  O   GLU A   1      22.961  33.778 -10.594  1.00 55.32           O  
# ATOM      5  CB  GLU A   1      20.729  33.181  -8.839  1.00 56.27           C  
# ATOM      6  CG  GLU A   1      21.358  34.449  -8.267  1.00 68.39           C  
# ATOM      7  CD  GLU A   1      20.612  35.733  -8.632  1.00 73.96           C  
# ATOM      8  OE1 GLU A   1      20.371  36.576  -7.721  1.00 74.73           O  
# ATOM      9  OE2 GLU A   1      20.278  35.900  -9.829  1.00 76.18           O1- 
# HETATM  991  O   HOH A 124      30.563  28.980  -5.128  1.00 20.12           O  
# HETATM  992  O   HOH A 125      27.155  30.562  -5.694  1.00 25.32           O  

# PDB format from http://deposit.rcsb.org/adit/docs/pdb_atom_format.html 
# COLUMNS        DATA TYPE       CONTENTS                            
# -------------------------------------------------------------------------------- 
#  1 -  6        Record name     "ATOM  "                                            
#  7 - 11        Integer         Atom serial number.                   
# 13 - 16        Atom            Atom name.                            
# 17             Character       Alternate location indicator.         
# 18 - 20        Residue name    Residue name.                         
# 22             Character       Chain identifier.                     
# 23 - 26        Integer         Residue sequence number.              
# 27             AChar           Code for insertion of residues.       
# 31 - 38        Real(8.3)       Orthogonal coordinates for X in Angstroms.                       
# 39 - 46        Real(8.3)       Orthogonal coordinates for Y in Angstroms.                            
# 47 - 54        Real(8.3)       Orthogonal coordinates for Z in Angstroms.                            
# 55 - 60        Real(6.2)       Occupancy.                            
# 61 - 66        Real(6.2)       Temperature factor (Default = 0.0).                   
# 73 - 76        LString(4)      Segment identifier, left-justified.   
# 77 - 78        LString(2)      Element symbol, right-justified.      
# 79 - 80        LString(2)      Charge on the atom.       
# Example: 
#          1         2         3         4         5         6         7         8 
# 12345678901234567890123456789012345678901234567890123456789012345678901234567890 
# ATOM    145  N   VAL A  25      32.433  16.336  57.540  1.00 11.92      A1   N 
# ATOM    146  CA  VAL A  25      31.132  16.439  58.160  1.00 11.85      A1   C 
# ATOM    147  C   VAL A  25      30.447  15.105  58.363  1.00 12.34      A1   C 
# ATOM    148  O   VAL A  25      29.520  15.059  59.174  1.00 15.65      A1   O 
# ATOM    149  CB AVAL A  25      30.385  17.437  57.230  0.28 13.88      A1   C 
# ATOM    150  CB BVAL A  25      30.166  17.399  57.373  0.72 15.41      A1   C 
# ATOM    151  CG1AVAL A  25      28.870  17.401  57.336  0.28 12.64      A1   C 
# ATOM    152  CG1BVAL A  25      30.805  18.788  57.449  0.72 15.11      A1   C 
# ATOM    153  CG2AVAL A  25      30.835  18.826  57.661  0.28 13.58      A1   C 
# ATOM    154  CG2BVAL A  25      29.909  16.996  55.922  0.72 13.25      A1   C 


##################################################### 
# importing the required libraries 
import sys 
import os 
import math 
import re 
#import numpy 
import string 


###################################################### 
# global variables 
x_array = [] 
y_array = [] 
z_array = [] 
chain_array = [] 
id_for_unique_chain_array = [] 
contrast_for_unique_chain_array = [] 
volume_for_unique_chain_array = [] 
corr_array = [] 
form_or_struc_factor_cor_array = [] 
n_total = 0 
r_max = 0 
r_step = 0 
total_q_array = [] 
total_intensity_array = [] 
total_err_array = [] 
scattering_chains_not_in_chain_parameters = [] 


###################################################### 
def is_numeric(s): 
    try: 
        float(s) 
        return True 
    except ValueError: 
        return False 


###################################################### 
def init_arrays( array_length ): 

	global corr_array 
	for i in range( 0, (array_length+2) ): 
		corr_array.append( float(0) ) 

	return 


###################################################### 
def main(argv=None): 

	# read the input parameters from the command line 

	if argv is None: 
		argv = sys.argv 
	output_file_name_prefix = sys.argv[1] 
	input_file_name_prefix = sys.argv[2] 
	n_points = sys.argv[3] 
	n_points = int(n_points) 
	q_max = sys.argv[4] 
	for j in range( 5, len(sys.argv), 3 ): 
		id_for_unique_chain_array.append( sys.argv[j] ) 
		contrast_for_unique_chain_array.append( float(sys.argv[j+1]) ) 
		volume_for_unique_chain_array.append( float(sys.argv[j+2]) ) 

	# initialise the array for total scattering, getting the q value from the first scattering file 

	# this_input_file_name = input_file_name_prefix + '_qI' + '.dat' 
	this_input_file_name = input_file_name_prefix + '_qI_AA' + '.dat' 
	infile = open( this_input_file_name, "r" ) 
	inlines = infile.readlines() 
	infile.close() 

	for inline in inlines: 
		inline = inline.strip() 
		if (inline != ''): 
			infields = inline.split() 
			this_q = float(infields[0]) 
			this_I = float(infields[1]) 
			this_err = float(infields[2]) 

			total_q_array.append( this_q ) 
			total_intensity_array.append( float(0) ) 
			total_err_array.append( float(0) ) 

	# for each component in the complex, read in it's scattering form factor, 
	# output the component scattering (adjusted for volume and contrast), 
	# and add it to the total scattering. 

	for i in range( 0, len(id_for_unique_chain_array) ): 
		for j in range( i, len(id_for_unique_chain_array) ): 
			this_input_file_name = input_file_name_prefix + '_qI_' + id_for_unique_chain_array[i] + id_for_unique_chain_array[j] + '.dat' 
			this_output_file_name = output_file_name_prefix + '_qI_' + id_for_unique_chain_array[i] + id_for_unique_chain_array[j] + '_XpV.dat' 

			outfile = open( this_output_file_name, "w" ) 

			this_multiplier = contrast_for_unique_chain_array[i] * contrast_for_unique_chain_array[j] * volume_for_unique_chain_array[i] * volume_for_unique_chain_array[j] 
			if ( i != j ): 
				this_multiplier = this_multiplier * 2 

			infile = open( this_input_file_name, "r" ) 
			inlines = infile.readlines() 
			infile.close() 

			k = 0 
			for inline in inlines: 
				inline = inline.strip() 
				if (inline != ''): 
					infields = inline.split() 
					this_q = float(infields[0]) 
					this_I = float(infields[1]) 
					this_err = float(infields[2]) 

					this_I = this_I * this_multiplier 
					this_err = this_err * this_multiplier 

					outline = ' ' + str(this_q) + ' ' + str(this_I) + ' ' + str(this_err) + "\r\n" 
					outfile.write( outline ) 

					total_intensity_array[k] = total_intensity_array[k] + this_I 
					total_err_array[k] = total_err_array[k] + this_err 
					if (this_q != total_q_array[k]): 
						print "ERROR : There is a mismatch between q's" 
					k = k + 1 

			outfile.close() 

	# write out the total scattering now that it has been calculated by adding scattering for all the components 

	this_output_file_name = output_file_name_prefix + '_qI_XpV.dat' 
	outfile = open( this_output_file_name, "w" ) 
	for i in range( 0, len(total_q_array) ): 
		q = total_q_array[i] 
		intensity = total_intensity_array[i] 
		err_value = total_err_array[i] 
		outline = ' ' + str(q) + ' ' + str(intensity) + ' ' + str(err_value) + "\r\n" 
		outfile.write( outline ) 
	outfile.close() 


if __name__ == "__main__": 
    sys.exit(main()) 


