# python calc_scat_with_without_params.py bamlet15_100D2O bamlet15.pdb    200 0.8 A -3.18 17202 B -6.21 10485 

# python calc_scat_with_without_params.py [output scattering file name prefix] [input PDB file of macromolecule] [approx_num_points_in_output] [q_max] [chain-id-1] [contrast-1] [volume-1] [chain-id-2] [contrast-2] [volume-2] ... 

# This python program started of as a python version of a Steen Hanson MCsimul program (written in FORTRAN) 
# that was download from http://www.matfys.kvl.dk/~steen/MCsimul.zip May 2012. 
# It was changed to read input in the PDB file format. 
# It was extended to process multiple densities according to Whitten and Trewhella 2008. 
# The Guinier processing in this python program has not yet been tested to ensure that it is still correct after the extensions made to the python program. 

# This program reads in a file in pdb format, containing points that represent scattering points. 
# For x-ray scattering, the points represent electrons. 
# This program reads the co-ordinates of the points, 
# not the atom type or any other information in the pdb format file that might be used for display. 
# The input is not a pdb file describing atomic co-ordinates, 
# but can conveniently be used by pdb programs to view the scattering points as if they were atoms. 
# The pdb format does not have a field for contrast (or scattering length density). 
# This field will be read from columns 81-88. If the field is not filled in the input pdb file, 
# then it will be defaulted to the value of 1 for that scattering point. 

# The user enters the contrasts and volumes of the components on the command line. 
# Any contrasts or volumes not entered will be set to 1 (which could be reasonable for contrast and probably not reasonable for volume). 

# This program outputs the total P(r) and the P(r) for each component. 
# This program outputs the total expected scattering, and the scattering for each component, 
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
def read_file( this_file, output_file_name_prefix ): 

	global x_array, y_array, z_array, chain_array, id_for_unique_chain_array, scattering_chains_not_in_chain_parameters 

	infile = open( this_file, "r" ) 
	inlines = infile.readlines() 
	infile.close() 

	n_total = 0 

	for inline in inlines: 
		inline = inline.strip() 
		if (inline != ''): 
			if ( len(inline) > 6 ): 
				atom_type = inline[0:6] 
				if ((atom_type == 'ATOM  ') or (atom_type == 'HETATM')): 
					x = inline[30:38] 
					y = inline[38:46] 
					z = inline[46:54] 
					chain = inline[21:22] 
					x = x.strip() 
					y = y.strip() 
					z = z.strip() 
					x = float(x) 
					y = float(y) 
					z = float(z) 
					x_array.append( x ) 
					y_array.append( y ) 
					z_array.append( z ) 
					chain_array.append( chain ) 
					n_total = n_total + 1 
	x_array.append( 0 ) 
	y_array.append( 0 ) 
	z_array.append( 0 ) 
	chain_array.append( 0 ) 

	for i in range( 0, (len(chain_array)-1) ): 
		this_chain = chain_array[i] 
		found_it = 0 
		j = 0 
		while ( j < len(id_for_unique_chain_array) ): 
			this_unique_chain = id_for_unique_chain_array[j] 
			if ( this_unique_chain == this_chain ): 
				found_it = 1 
				j = len(id_for_unique_chain_array) 
			j = j + 1 
		if ( found_it == 0 ): 
			found_err_msg = 0 
			j = 0 
			while ( j < len(scattering_chains_not_in_chain_parameters) ): 
				this_err_chain = scattering_chains_not_in_chain_parameters[j] 
				if ( this_err_chain == this_chain ): 
					found_err_msg = 1 
					j = len(scattering_chains_not_in_chain_parameters) 
				j = j + 1 
			if ( found_err_msg == 0 ): 
				scattering_chains_not_in_chain_parameters.append( this_chain ) 
				err_msg = 'Input scattering points having chain ' + this_chain + ' will be ignored because no information was given for this chain in the command line input parameters.' + "\n" 
				print err_msg 

	return n_total 


###################################################### 
def calculate_parameters( output_info_file_name, n_total ): 

	global x_array, y_array, z_array, chain_array 
	global r_max, r_step 

	init_arrays( n_total+2 ) 
	x_max = -1000000000000000000000000000000 
	x_min = 1000000000000000000000000000000 
	y_min = x_min 
	y_max = x_max 
	z_min = x_min 
	z_max = x_max 
	for i in range( 0, n_total ): 
		x = x_array[i] 
		y = y_array[i] 
		z = z_array[i] 
		if (x < x_min): 
			x_min = x 
		if (x > x_max): 
			x_max = x 
		if (y < y_min): 
			y_min = y 
		if (y > y_max): 
			y_max = y 
		if (z < z_min): 
			z_min = z 
		if (z > z_max): 
			z_max = z 

	outfile = open( output_info_file_name, "w" ) 
	outline = 'No. of scatterers = ' + str(n_total) + "\r\n" 
	outfile.write( outline ) 
	outline = 'Dimensions = ' + str(x_min) + ' ' + str(x_max) + ' ' + str(y_min) + ' ' + str(y_max) + ' ' + str(z_min) + ' ' + str(z_max) + ' ' + "\r\n" 
	outfile.write( outline ) 

	r_max = math.sqrt( (x_max - x_min)**2 + (y_max - y_min)**2 + (z_max - z_min)**2 ) 
	r_step = r_max / float(n_total) 

	outline = 'Step length for dd-function = ' + str(r_step) + ' ' + "\r\n" 
	outfile.write( outline ) 
	outfile.close() 

	return 


###################################################### 
def guinier_radius( output_info_file_name, n_total ): 

	global x_array, y_array, z_array, chain_array 

	cmx = 0 
	cmy = 0 
	cmz = 0 
	b_total = 0 
	for i in range( 1, (n_total+1) ): 
		cmx = cmx + 1 # b_array[i] * x_array[i] 
		cmy = cmy + 1 # b_array[i] * y_array[i] 
		cmz = cmz + 1 # b_array[i] * z_array[i] 
		b_total = 1 + 1 # b_total + b_array[i] 
	cmx = cmx / b_total 
	cmy = cmy / b_total 
	cmz = cmz / b_total 

	outfile = open( output_info_file_name, "a" ) 
	outline = 'Center of (scattering-)mass = ' + str(cmx) + ' ' + str(cmy) + ' ' + str(cmz) + "\r\n" 
	outfile.write( outline ) 

	rg = 0 
	for i in range( 1, (n_total+1) ): 
		# rg = rg + b_array[i] * ( (x_array[i] - cmx)**2 + (y_array[i] - cmy)**2 + (z_array[i] - cmz)**2 ) 
		rg = rg + ( (x_array[i] - cmx)**2 + (y_array[i] - cmy)**2 + (z_array[i] - cmz)**2 ) 
	rg_div_b = rg / b_total 
	if ( rg_div_b >= 0 ) : 
		rg = math.sqrt( rg / b_total ) 
	else: 
		rg = 'NaN' 

	outline = 'Guinier-radius = ' + str(rg) + "\r\n" 
	outfile.write( outline ) 
	outfile.close() 

	return 


###################################################### 
def distance_distribution_function( output_Pr_file_name, output_info_file_name, output_file_name_prefix, n_total, r_step, q_max ): 

	global x_array, y_array, z_array, chain_array, id_for_unique_chain_array 
	global corr_array, form_or_struc_factor_cor_array 

	for i in range( 0, len(id_for_unique_chain_array) ): 
		temp_new_array = [] 
		for j in range( 0, len(id_for_unique_chain_array) ): 
			temp_new_array_2 = [] 
			for k in range( 0, (n_total+2) ): 
				temp_new_array_2.append(0) 
			temp_new_array.append( temp_new_array_2 ) 
		form_or_struc_factor_cor_array.append( temp_new_array ) 

	for i in range( 0, (n_total+2) ): 
		corr_array[i] = 0 

	#n_div = n_total % 2 
	#n_total = n_total - n_div 
	for j in range( 0, n_total ): 
		for k in range( 0, n_total ): 
			if ( j != k ): 
				dxx = x_array[j] - x_array[k] 
				dyy = y_array[j] - y_array[k] 
				dzz = z_array[j] - z_array[k] 
				d = math.sqrt( dxx*dxx + dyy*dyy + dzz*dzz ) 
				#print x_array[j], x_array[k], y_array[j], y_array[k], z_array[j], z_array[k] 
				l = int(round( d / r_step )) 
				corr_array[l] = corr_array[l] + 1 

				this_chain_1 = chain_array[j] 
				this_chain_2 = chain_array[k] 
				found_it_1 = 0 
				found_it_2 = 0 
				chain_1_idx = -1 
				chain_2_idx = -1 
				m = 0 
				while ( m < len(id_for_unique_chain_array) ): 
					this_unique_chain = id_for_unique_chain_array[m] 
					if ( this_unique_chain == this_chain_1 ): 
						found_it_1 = 1 
						chain_1_idx = m 
					if ( this_unique_chain == this_chain_2 ): 
						found_it_2 = 1 
						chain_2_idx = m 
					if ( (found_it_1 == 1) and (found_it_2 == 1) ): 
						m = len(id_for_unique_chain_array) 
					m = m + 1 
				if ( (chain_1_idx == -1) or (chain_2_idx == -1) ): 
					ignore_this_pair = 1 
				else: 
					form_or_struc_factor_cor_array[chain_1_idx][chain_2_idx][l] = form_or_struc_factor_cor_array[chain_1_idx][chain_2_idx][l] + 1 

	n_max = 0 
	for i in range( 1, (n_total+1) ): 
		if ( corr_array[i] != 0 ): 
			n_max = i 
	outfile = open( output_info_file_name, "a" ) 
	max_length_molecule = n_max * r_step 
	outline = 'Maximum length of molecule = ' + str(max_length_molecule) + "\r\n" 
	outfile.write( outline ) 
	outline = 'Qmax at Fourier transformation = ' + str(q_max) + "\r\n" 
	outfile.write( outline ) 
	outfile.close() 

	corr_array[0] = 2 * corr_array[0] 
	corr_array[n_max] = 2 * corr_array[n_max] 
	outfile = open( output_Pr_file_name, "w" ) 
	n_total = n_max 
	for j in range( 1, (n_max+1) ): 
		r = j * r_step 
		pr = corr_array[j] 
		if ( corr_array[j] >= 0): 
			pr_sqrt = math.sqrt( corr_array[j] ) 
		else: 
			pr_sqrt = 'NaN' 
		outline = ' ' + str(r) + ' ' + str(pr) + ' ' + str(pr_sqrt) + "\r\n" 
		outfile.write( outline ) 
	outfile.close() 

	for i in range( 0, len(id_for_unique_chain_array) ): 
		for j in range( i, len(id_for_unique_chain_array) ): 
			output_file_name = output_file_name_prefix + '_Pr_' + id_for_unique_chain_array[i] + id_for_unique_chain_array[j] + '.dat' 
			outfile = open( output_file_name, "a" ) 
			for k in range( 1, (n_max+1) ): 
				r = k * r_step 
				pr = form_or_struc_factor_cor_array[i][j][k] 
				if ( form_or_struc_factor_cor_array[i][j][k] >= 0): 
					pr_sqrt = math.sqrt( form_or_struc_factor_cor_array[i][j][k] ) 
				else: 
					pr_sqrt = 'NaN' 
				outline = ' ' + str(r) + ' ' + str(pr) + ' ' + str(pr_sqrt) + "\r\n" 
				outfile.write( outline ) 
			outfile.close() 

	return n_total 


###################################################### 
def fourier_transform( output_file_name, output_file_name2, q_step, r_step, n_total, n_points, corr_array, pV_multiplier ): 

	f_intensity = [] 
	err_array = [] 
	for i in range( 0, max((n_total+1), (n_points+1)) ): 
		f_intensity.append( float(0) ) 
		err_array.append( float(0) ) 

	for k in range( 1, (n_points+1) ): 
		q = k * q_step 
		#print 'q =', q 
		fadd = 0 

		eadd = 0 
		for i in range( 1, (n_total+1) ): 
			r = i * r_step 
			debeye = math.sin( q * r ) / (q * r ) 
			fadd = corr_array[i] * debeye 
			f_intensity[k] = f_intensity[k] + fadd 
			#print 'q =', q, '     r =', r, '     Pr =', corr_array[i], '     sin(qr) =', math.sin( q * r ), '     qr =', ( q * r ), '     sin(qr)/qr =', (math.sin( q * r ) / (q * r )), '     intensity to add =', fadd 
			eadd = debeye * corr_array[i] * debeye 
			err_array[k] = err_array[k] + eadd 
		f_intensity[k] = f_intensity[k] + corr_array[0] * 0.5 
		err_array[k] = err_array[k] + corr_array[0] * 0.5 
	for i in range( 1, (n_total+1) ): 
		f_intensity[0] = f_intensity[0] + corr_array[i] 
		err_array[0] = err_array[0] + corr_array[i] 
	f_intensity[0] = f_intensity[0] + 0.5 * corr_array[0] 
	err_array[0] = err_array[0] + 0.5 * corr_array[0] 
	for i in range( 1, (n_points+1) ): 
		err_array[i] = math.sqrt( abs(err_array[i]) ) 

	outfile2 = open( output_file_name2, "w" ) 
	Io = f_intensity[0] 
	for i in range( 0, (n_points+1) ): 
		q = i * q_step 
		f_intensity[i] = f_intensity[i] / Io 
		err_array[i] = err_array[i] / Io 
		intensity_value = f_intensity[i] 
		err_value = err_array[i] 
		outline2 = ' ' + str(q) + ' ' + str(intensity_value) + ' ' + str(err_value) + "\r\n" 
		outfile2.write( outline2 ) 
	outfile2.close() 

	outfile = open( output_file_name, "w" ) 

	for i in range( 0, (n_points+1) ): 
		q = i * q_step 
		f_intensity[i] = f_intensity[i] * pV_multiplier 
		err_array[i] = err_array[i] * pV_multiplier 
		intensity_value = f_intensity[i] 
		err_value = err_array[i] 
		outline = ' ' + str(q) + ' ' + str(intensity_value) + ' ' + str(err_value) + "\r\n" 
		outfile.write( outline ) 
		total_q_array[i] = q 
		total_intensity_array[i] = total_intensity_array[i] + intensity_value 
		total_err_array[i] = total_err_array[i] + err_value 
	outfile.close() 

	return 

###################################################### 
def main(argv=None): 

	# read the input parameters from the command line 

	if argv is None: 
		argv = sys.argv 
	output_file_name_prefix = sys.argv[1] 
	input_file_name = sys.argv[2] 
	n_points = sys.argv[3] 
	n_points = int(n_points) 
	q_max = sys.argv[4] 
	for j in range( 5, len(sys.argv), 3 ): 
		id_for_unique_chain_array.append( sys.argv[j] ) 
		contrast_for_unique_chain_array.append( float(sys.argv[j+1]) ) 
		volume_for_unique_chain_array.append( float(sys.argv[j+2]) ) 

	# initialise some of the output files 

	output_qI_file_name = output_file_name_prefix + '_qI_XpV.dat' 
	output_Pr_file_name = output_file_name_prefix + '_Pr.dat' 
	output_info_file_name = output_file_name_prefix + '_info.txt' 
	for i in range( 0, len(id_for_unique_chain_array) ): 
		for j in range( i, len(id_for_unique_chain_array) ): 
			output_file_name = output_file_name_prefix + '_Pr_' + id_for_unique_chain_array[i] + id_for_unique_chain_array[j] + '.dat' 
			outfile = open( output_file_name, "w" ) 
			outfile.close() 

	# read the input file that contains the coordinates of the scattering points 

	n_total = read_file( input_file_name, output_file_name_prefix ) 

	# output some of the information calculated about these scattering points 

	calculate_parameters( output_info_file_name, n_total ) 

	# calculate the guinier radius of the scattering complex 
	# this code hasn't been tested or verified since copying it from the original fortran scat program 

	guinier_radius( output_info_file_name, n_total ) 

	# calculate the pair distance distribution functions for this scattering complex of points 

	n_total = distance_distribution_function( output_Pr_file_name, output_info_file_name, output_file_name_prefix, n_total, r_step, q_max ) 

	# initialise the arrays for total scattering 

	q_step = float(q_max) / n_points 
	for i in range( 0, (n_points+1) ): 
		q = i * q_step 
		total_q_array.append( float(0) ) 
		total_intensity_array.append( float(0) ) 
		total_err_array.append( float(0) ) 

	# for each component in the complex, calculate its pair distance distribution function, 
	# then calculate it's scattering as calculated by the fourier transform of the distribution function. 
	# finally, add the component scattering to the total scattering. 

	for i in range( 0, len(id_for_unique_chain_array) ): 
		for j in range( i, len(id_for_unique_chain_array) ): 
			this_output_file_name = output_file_name_prefix + '_qI_' + id_for_unique_chain_array[i] + id_for_unique_chain_array[j] + '_XpV.dat' 
			this_output_file_name2 = output_file_name_prefix + '_qI_' + id_for_unique_chain_array[i] + id_for_unique_chain_array[j] + '.dat' 
			this_Pr = form_or_struc_factor_cor_array[i][j] 
			#print 'i =', i, ', chain =', id_for_unique_chain_array[i] 
			#print 'j =', j, ', chain =', id_for_unique_chain_array[j] 
			#print 'p1 =', contrast_for_unique_chain_array[i] 
			#print 'V1 =', volume_for_unique_chain_array[i] 
			#print 'p2 =', contrast_for_unique_chain_array[j] 
			#print 'V2 =', volume_for_unique_chain_array[j] 
			this_multiplier = contrast_for_unique_chain_array[i] * contrast_for_unique_chain_array[j] * volume_for_unique_chain_array[i] * volume_for_unique_chain_array[j] 
			if ( i != j ): 
				this_multiplier = this_multiplier * 2 
				#print '2 = 2' 
			#print 'multiplier =', this_multiplier 
			#print 
			#print 'before fourier_transform for form/structure factor,', i, j, this_multiplier, contrast_for_unique_chain_array[i], contrast_for_unique_chain_array[j], volume_for_unique_chain_array[i], volume_for_unique_chain_array[j] 
			fourier_transform( this_output_file_name, this_output_file_name2, q_step, r_step, n_total, n_points, this_Pr, this_multiplier ) 

	# write out the total scattering now that it has been calculated by adding scattering for all the components 

	outfile = open( output_qI_file_name, "w" ) 
	for i in range( 0, len(total_q_array) ): 
		q = total_q_array[i] 
		intensity = total_intensity_array[i] 
		err_value = total_err_array[i] 
		outline = ' ' + str(q) + ' ' + str(intensity) + ' ' + str(err_value) + "\r\n" 
		outfile.write( outline ) 
	outfile.close() 


if __name__ == "__main__": 
    sys.exit(main()) 

