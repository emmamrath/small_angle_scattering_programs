/****************************************************************** 
 ***                          molecule2.c                       *** 
 ******************************************************************
 *
 * Additional Info:
 *
 * Date: 9/4/06
 * Author: A.E. Whitten
 * Description:
 * 
 * Data structure, and functions for manipulating molecule data.
 *
 * Version 1.0: File created
 *
 * Oct 2016 : Create molecule2.c from molecule.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "molecule2.h"

/* Functions involving the manipulation of the molecule data */

int getComposition(molecule *mol){   /* Converts a formula into an atomic compostion */

	char element[3] = {'\0','\0','\0'}, n[3] = {'\0','\0','\0'}, *elementPos, *nPos;
	int i;
	
	/* This function requires some explanation, as it is not obvious what it is doing. There are three different
	 * tests carried out: Is the character - Upper case; Lower case; Number. An element always starts with a capital
	 * followed by a number, or a lower case character and then a number. The tests are carried out by checking
	 * their position on the ascii character table - which is why char is being compared to int in the if statements
	 * below.
	 * Upper case (A - Z): 65 - 90
	 * Lower case (a - z): 97 - 122
	 * Numbers (0 - 9): 48 - 57 */

	elementPos = element;   /* Set the pointers to the start of the element and n arrays */
	nPos = n;

	for(i = 0; mol->formula[i] != '\0'; i ++){
	
		/* If it is a character then it is part of the element label */
	
		if((mol->formula[i] >= 65 && mol->formula[i] <= 90) || (mol->formula[i] >= 97 && mol->formula[i] <= 122)){
			
			*elementPos = mol->formula[i];
			elementPos ++;   /* Increment the pointer */
			
		}
			
		/* Otherwise it is a number */
		
		else if(mol->formula[i] >= 48 && mol->formula[i] <= 57){
			
			*nPos = mol->formula[i];
			nPos ++;
			
		}
		
		/* If the next character is a capital or if the end of the string has been reached, information can be processed for that element */
		
		if((mol->formula[i+1] >= 65 && mol->formula[i+1] <= 90) || (mol->formula[i+1] == '\0')){
			
			if(n[0] == '\0') n[0] = '1'; /* If n is blank, then 1 is implied in the formula */
			mol->composition[getAtomicNum(element)] += strtod(n, (char **)NULL); /* Convert n to a double */
			//fprintf(stdout,"element = %s, getAtomicNum(element) = %d, mol->composition[getAtomicNum(element)] = %f<br>\n", element, getAtomicNum(element), mol->composition[getAtomicNum(element)]);

			element[0] = '\0'; element[1] = '\0'; element[2] = '\0'; n[0] = '\0'; n[1] = '\0'; n[2] = '\0';
			elementPos = element;   /* Reset all values, to start again */
			nPos = n;
		
		}
		
	}

	return(1);
	
}

double getMolecularMass(molecule *mol){   /* Returns the molecular mass of a molecule */
	
	int i;
	double mass = 0.0;
	
	for(i = 0; i < TABLESIZE; i ++)
		mass += mol->composition[i]*getAtomicMass(i);
		
	return(mass);
	
}

int getMolecularVol(molecule *mol){   /* Updates the vol component of the structure, and returns an integer flag */
	
	int i;
		
	for(i = 0; i < TABLESIZE; i ++)
		mol->vol += getAtomicVol(i)*mol->composition[i];
	
	return(1);
	
}

int addMoleculeComponents(molecule *molPlusMol, molecule *mol){   /* Routine for adding the components of 2 molecules together */
	
	/* molPlusMol is a pointer to the structure that is to be added to, mol is what is to be added to molPlusMol
	 * nMol can be considered as the concentration or number of mol units. It is esentially a multiplier for the
	 * data structure components */

	int i;

	molPlusMol->vol += mol->vol;
	molPlusMol->exH += mol->exH;

	for(i = 0; i < TABLESIZE; i ++)
		molPlusMol->composition[i] += mol->composition[i];

	return(1);
			
}

double getNFScat(molecule *mol){   /* Calculates the forward scattering of neutrons for a given molecule */
	
	int i;
	double nForwardScat = 0.0;
	 	
	for(i = 0; i < TABLESIZE; i ++)
		nForwardScat += getNSL(i)*mol->composition[i];
		
	return(nForwardScat);
	
}

double getXFScat(molecule *mol){   /* Calculates the forward scattering of X-rays for a given molecule */
	
	int i;
	double XForwardScat = 0.0;
	 		
	for(i = 0; i < TABLESIZE; i ++)
		XForwardScat += getXSL(i)*mol->composition[i]; 
	
	return(XForwardScat);
	
}

