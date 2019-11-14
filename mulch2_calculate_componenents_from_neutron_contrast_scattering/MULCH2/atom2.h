/****************************************************************** 
 ***                            atom2.h                         *** 
 ******************************************************************
 *
 * Additional Info:
 *
 * Date: 9/4/06
 * Author: A.E. Whitten
 * Description:
 * 
 * Function prototypes and definitions from atom.c
 *
 * Version 1.0: File created
 *
 * Oct 2016 : Create compost2.c from compost.c
 *
 * Oct 2016 : Changes to atom.c to create atom2.c include :
 *            Include cesium and others.
 */

/* Definitions */

#define TABLESIZE 119
#define THOMPSON_RADIUS 2.82   /* x 10**-13 cm*/
#define NAVAGADRO 6.0221415E+23   /* Atoms per mole */

/* *********** */

 /* Data structures */
 
typedef enum{

	 D, H, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar, K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu,
	 Zn, Ga, Ge, As, Se, Br, Kr, Rb, Sr, Y, Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, Sn, Sb, Te, I, Xe,
	 Cs, Ba,
		La, Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu,
			Hf, Ta, W, Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi, Po, At, Rn,
	 Fr, Ra,
		Ac, Th, Pa, U, Np, Pu, Am, Cm, Bk, Cf, Es, Fm, Md, No, Lr,
			Rf, Db, Sg, Bh, Hs, Mt, Ds, Rg, Cn, Uut, Uuq, Uup, Uuh, Uus, Uuo

}enumAtoms;   /* Structure allowing specification of a specific element, without knowing its location in the array */

 /* *************** */
 
/* Function prototypes */

int getAtomicNum(char *atom);   /* returns the atomic number of a given atom */

double getNSL(int atomicNum);   /* returns the neutron scattering length of a given atom */

double getXSL(int atomicNum);   /* returns the X-ray scattering length of a given atom */

double getAtomicMass(int atomicNum);   /* returns the atomic mass of a given atom */

double getAtomicVol(int atomicNum);  /* returns the atomic volume of a given atom */

char *getElementLabel(int atomicNum);   /* returns the symbol of an element */ 

