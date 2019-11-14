/****************************************************************** 
 ***                            atom2.c                         *** 
 ******************************************************************
 *
 * Additional Info:
 *
 * Date: 9/4/06
 * Author: A.E. Whitten
 * Description:
 * 
 * Variety of atomic information and functions, for extracting 
 * atomic properties
 *
 * Version 1.0: File created
 *
 * Oct 2016 : Create compost2.c from compost.c
 *
 * Oct 2016 : Changes to atom.c to create atom2.c include :
 *            Include cesium and others.
 */

#include <stdio.h>
#include "atom2.h"

/* Atomic data */

const char ELEMENT[TABLESIZE][4] ={
	
	/* Labels for all elements, so that if need be atom labels can be searched to find the
	 * position of that element in the array. Deuterium has been added */
	
	"D\0",
	"H\0",                                                                                                                 "He\0",
	"Li\0","Be\0",                                                                      "B\0" ,"C\0" ,"N\0" ,"O\0" ,"F\0" ,"Ne\0",
	"Na\0","Mg\0",                                                                      "Al\0","Si\0","P\0" ,"S\0" ,"Cl\0","Ar\0",
	"K\0" ,"Ca\0","Sc\0","Ti\0","V\0" ,"Cr\0","Mn\0","Fe\0","Co\0","Ni\0","Cu\0","Zn\0","Ga\0","Ge\0","As\0","Se\0","Br\0","Kr\0",
	"Rb\0","Sr\0","Y\0" ,"Zr\0","Nb\0","Mo\0","Tc\0","Ru\0","Rh\0","Pd\0","Ag\0","Cd\0","In\0","Sn\0","Sb\0","Te\0","I\0" ,"Xe\0",
	"Cs\0","Ba\0",
	              "La\0","Ce\0","Pr\0","Nd\0","Pm\0","Sm\0","Eu\0","Gd\0","Tb\0","Dy\0","Ho\0","Er\0","Tm\0","Yb\0","Lu\0",
	                     "Hf\0","Ta\0","W\0", "Re\0","Os\0","Ir\0","Pt\0","Au\0","Hg\0","Tl\0","Pb\0","Bi\0","Po\0","At\0","Rn\0",
	"Fr\0","Ra\0",
	              "Ac\0","Th\0","Pa\0","U\0", "Np\0","Pu\0","Am\0","Cm\0","Bk\0","Cf\0","Es\0","Fm\0","Md\0","No\0","Lr\0",
	                     "Rf\0","Db\0","Sg\0","Bh\0","Hs\0","Mt\0","Ds\0","Rg\0","Cn\0","Uut\0","Uuq\0","UUp\0","Uuh\0","Uus\0","Uuo\0"

};

const double ATOMICMASS[TABLESIZE] ={   /* Atomic masses to 4 sig. figs */

	2.014,
	1.008,                                                                                                4.003,
	6.941,9.012,                                                            10.81,12.01,14.01,16.00,19.00,20.18,
	22.99,24.31,                                                            26.98,28.09,30.97,32.07,35.45,39.95,
	39.10,40.08,44.96,47.88,50.94,52.00,54.94,55.85,58.93,58.69,63.55,65.39,69.72,72.59,74.92,78.96,79.90,83.80,
	85.47,87.62,88.91,91.22,92.91,95.94,98.00,101.1,102.9,106.4,107.9,112.4,114.8,118.7,121.8,127.6,126.9,131.3,
	132.9054519,137.327,
		138.90547,140.116,140.90765,144.242,145,150.36,151.964,157.25,158.92535,162.5,164.93032,167.259,168.93421,173.054,174.9668,
			178.49,180.94788,183.84,186.207,190.23,192.217,195.084,196.966569,200.59,204.3833,207.2,208.9804,209,210,222,
	223,226,
		227,232.03806,231.03588,238.02891,237,244,243,247,247,251,252,257,258,259,262,
			267,268,271,272,270,276,281,280,285,284,289,288,293,294,294

};

const double NSL[TABLESIZE] ={  /* Bound, average, coherant scattering lengths */ 

	/* Values for H and D refer to specific isotopes, and those with a scattering length of zero have complex magnitudes */
	/* Values for D and values from H to Xe were in the original Mulch program */
	/* Some of the values from Cs to Uuo were taken from http://www.ncnr.nist.gov/resources/n-lengths/ 
	   but most were set to 0 and need to be filled in with the correct value if they are to be used. */

	 6.671,
	-3.741,                                                                                                                 3.260,
	-1.900, 7.790,                                                                       0.000, 6.646, 9.360, 5.803, 5.654, 4.566,
	 3.630, 5.375,                                                                       3.449, 4.149, 5.130, 2.847, 9.577, 1.909,
	 3.670, 4.700,12.290,-3.438,-0.382, 3.635,-3.730, 9.450, 2.490,10.300, 7.718, 5.680, 7.288, 8.185, 6.580, 7.970, 6.795, 7.810,
	 7.900, 7.020, 7.750, 7.160, 7.054, 6.715, 6.800, 7.030, 5.880, 5.910, 5.922, 0.000, 0.000, 6.225, 5.570, 5.800, 5.280, 4.920,
	 5.42,  5.07,
	               0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
	                      0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
	 0,     0,
	               0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
	                      0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0

};

const double ATOMICVOL[TABLESIZE] ={  /* approximate, atomic volumes - very rough */ 

	/* These values are approximate - Na, Mg, K, Ca, Cl, Br and I are calculated from ionic radii (Shannon radii)
	 * H, C, N and O are taken from Fraser, MacRae and Suzuki 
	 * P, S, Fe, Co, Ni, Cu and Zn are calculated from atomic radii (periodic.lanl.gov - originally from CRC-HCP) */

	/* van der Waals volume = 4/3 PI radius^3 		4/3 PI = 4.188786667 */

	/*
		http://abulafia.mt.ic.ac.uk/shannon/radius.php?Element=Li 	: coordination VI, ionic radius 0.76
		http://abulafia.mt.ic.ac.uk/shannon/radius.php?Element=Na 	: coordination VI, ionic radius 1.38
		http://abulafia.mt.ic.ac.uk/shannon/radius.php?Element=K 	: coordination VI, ionic radius 1.02
		http://abulafia.mt.ic.ac.uk/shannon/radius.php?Element=Cs 	: coordination VI, ionic radius 1.67

		All data presented here has been taken from "Revised Effective Ionic Radii and Systematic Studies of Interatomic Distances in Halides and Chalcogenides" 
		By R. D. Shannon. Central Research and Development Department, Experimental Station, E. I. Du Pont de Nemours and Company, Wilmington, Delaware 19898, U.S.A.
		Published in Acta Crystallographica. (1976). A32, Pages 751-767.
	*/

	  5.15,
	  5.15,                                                                                                                 0.000,
	  1.84, 0.000,                                                                       0.000, 16.44,  2.49,  9.13, 0.000, 0.000,
	  4.45,  1.56,                                                                       0.000, 0.000,  3.37, 26.09, 24.84, 0.000,
	 11.01,  4.19, 0.000, 0.000, 0.000, 0.000,  0.00,  7.99,  7.99,  8.18,  8.78,  9.85, 0.000, 0.000, 0.000, 0.000, 31.54, 0.000,
	 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 44.60, 0.000,
	 19.51, 0,
	               0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
	                      0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
	 0,     0,
	               0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
	                      0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0

};

/* *********** */

/* Functions involving the manipulation of the atomic data above */

int getAtomicNum(char *atom){   /* Returns the atomic number (technically the array postion), given the element label.
                                 * The exception is D, which lies in the 0 position of the array, even though Z = 1 */ 
	int i;
	
	for(i = 0; i < TABLESIZE; i ++) {
		//fprintf(stdout,"getAtomicNum : i = %d, ELEMENT[i] = %s<br>\n", i, ELEMENT[i]);
		if(strcasecmp(ELEMENT[i], atom) == 0) return(i);
	}
	
	return(-1);  /* If the function does not find a match */
	
}

double getNSL(int atomicNum){   /* Returns the neutron scattering length given the position of that element in the array */
	
	if(atomicNum < TABLESIZE) return(NSL[atomicNum]);
	
	return(0.0);
	
}

double getXSL(int atomicNum){   /* Returns the x-ray scattering length, given the atomic number */
	
	if(atomicNum < TABLESIZE){
	
		if(atomicNum == 0) return(THOMPSON_RADIUS);   /* Special case for deuterium */
		else return(THOMPSON_RADIUS*(double)atomicNum);
		
	}
	
	return(0.0);
	
}

double getAtomicMass(int atomicNum){   /* Returns the atomic mass, given the atomic number */

	if(atomicNum < TABLESIZE) return(ATOMICMASS[atomicNum]);
	
	return(0.0);
	
}

double getAtomicVol(int atomicNum){   /* Returns the atomic volume, given the atomic number */

	if(atomicNum < TABLESIZE) return(ATOMICVOL[atomicNum]);
	
	return(0.0);
	
}

char *getElementLabel(int atomicNum){   /* Returns the element label, given the atomic number */

	if(atomicNum < TABLESIZE) return((char *)ELEMENT[atomicNum]);
	
	return((char *)NULL);
	
}
