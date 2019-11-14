/****************************************************************** 
 ***                          molecule2.h                       *** 
 ******************************************************************
 *
 * Additional Info:
 *
 * Date: 9/4/06 last modfied on 25/9/08
 * Author: A.E. Whitten
 * Description:
 * 
 * Function prototypes and definitions from molecule.c
 *
 * Version 1.0: File created
 *
 * Version 1.1: allocateMolecule function changed to a macro
 *
 * Oct 2016 : Create molecule2.h from molecule.h
 */
 
#include "atom2.h"

#define MOLECULE(mol, oldSize, newSize, name, func_name, i, j){								\
															\
	if(oldSize == 0){												\
	        													\
		if (( mol = (molecule *)malloc((newSize)*sizeof(molecule))) == NULL){					\
															\
			fprintf(stderr, "\nError: Can't allocate memory for *%s in function %s.\n", name, func_name);	\
			return(0);											\
															\
		}													\
	}														\
    															\
	else if(oldSize != newSize){   /* If more memory needs to be allocated */					\
															\
		 if (( mol = (molecule *)realloc(mol, (newSize)*sizeof(molecule))) == NULL){				\
															\
			fprintf(stderr, "\nError: Can't reallocate memory for *%s in function %s.\n", name, func_name);	\
			return(0);											\
															\
		}													\
	} 														\
															\
	mol->size = newSize;   /* Update the number of elements */							\
															\
	for(i = oldSize; i < newSize; i ++){  /* Initialise componenets of the structure */				\
															\
		mol[i].vol = 0.0;											\
		mol[i].exH = 0.0;											\
															\
		for(j = 0; j < TABLESIZE; j ++)										\
			mol[i].composition[j] = 0.0;									\
															\
		for(j = 0; j < 9999; j ++)										\
			mol[i].formula[j] = '\0';									\
															\
	}														\
															\
}															\
/* Allocates (or reallocates) memory for molecule data structures, and initialises the elements */

/* Data structures */
 
typedef enum{   /* Enumerated data structure matching for interpreting the formula of a molecule*/

	 PROTEIN, RNA, DNA, MOLECULE

}enumSequence;
 
typedef struct{   /* Data stucture that contains amino acid and protein data*/

	int size, status;
	char formula[9999];
	double n, exH, vol, composition[TABLESIZE];
	enumSequence type;
	
	/* size is the number of elements in a molecule array, and is set when allocate_molecule is called
	 * fragment is a flag defining which fragment the molecule belongs to
	 * formula indicates the formula of a molcule e.g. C10H10 for a molecule or ATGC for an amino or nucleic acid sequence
	 * exH is the number of exchangable protons associated with the molecule
	 * vol is the volume of the molecule
	 * composition contains information regarding the number of atoms
	 * type flags how the formula is to be interpreted */

}molecule;

 /* *************** */
 
/* Function prototypes */

int getComposition(molecule *mol);   /* converts a formula (e.g. Al3SO4) into a composition of atoms */

double getMolecularMass(molecule *mol);   /* given a composition, calculates the molecular mass */

int getMolecularVol(molecule *mol);   /* given a composition, calculates the molecular volume */

int addMoleculeComponents(molecule *molPlusMol, molecule *mol);   /* adds the characteristics of one molecule array to another */

double getNFScat(molecule *mol);   /* given a composition, calculates the neutron forward scattering */

double getXFScat(molecule *mol);   /* given a composition, calculates the X-ray forward scattering */


