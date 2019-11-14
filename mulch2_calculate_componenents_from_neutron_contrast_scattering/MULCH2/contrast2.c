/****************************************************************** 
 ***                          contrast2.c                       *** 
 ******************************************************************
 *
 * Additional Info:
 *
 * Date: 9/4/06
 * Author: A.E. Whitten
 * Description:
 * 
 * Various functions and routines specific to the Contrast tool.
 *
 * Version 1.0: File created
 *
 * Oct 2016 : Create contrast2.c from contrast.c
 *
 * Oct 2016 : Changes include file name and parameter changes
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mathFunctions2.h"
#include "molecule2.h"
 
 /* Definitions */

#define AMINOACIDS 20
#define NUCLEOTIDES 4

#define MOLESOFWATERPERLITRE (0.9982E+3/(1.008*2+16.00))   /* @ 20 degrees C */

#define CONTRAST(data, oldSize, newSize, name, func_name, i){								\
															\
	if(oldSize == 0){												\
	        													\
		if (( data = (contrast_t *)malloc((newSize)*sizeof(contrast_t))) == NULL){				\
															\
			fprintf(stderr, "\nError: Can't allocate memory for *%s in function %s.\n", name, func_name);	\
			return(0);											\
															\
		}													\
															\
	}														\
    															\
	else if(oldSize != newSize){											\
															\
		 if (( data = (contrast_t *)realloc(data, (newSize)*sizeof(contrast_t))) == NULL){			\
															\
			fprintf(stderr, "\nError: Can't reallocate memory for *%s in function %s.\n", name, func_name);	\
			return(0);											\
															\
		}													\
															\
	} 														\
															\
	data->size = newSize;												\
															\
	for(i = oldSize; i < newSize; i ++){										\
															\
		data[i].bX = 0.0;											\
		data[i].bN[0] = 0.0;											\
		data[i].bN[1] = 0.0;											\
		data[i].bN[2] = 0.0;											\
		data[i].XchH = 0.0;											\
		data[i].nH = 0.0;											\
		data[i].nD = 0.0;											\
		data[i].fAccessH = 0.0;											\
		data[i].mass[0] = 0.0;											\
		data[i].mass[1] = 0.0;											\
		data[i].mass[2] = 0.0;											\
		data[i].V = 0.0;											\
															\
	}														\
															\
}
/* Allocates (or reallocates) memory for contrast data structures, and initialises the elements */

/* *********** */

 /* Data structures */
 
typedef enum{

	 ALA, CYS, ASP, GLU, PHE, GLY, HIS, ILE, LYS, LEU, MET, ASN, PRO, GLN, ARG, SER, THR, VAL, TRP, TYR

}enumAminoAcid;   /*enumerated data structure matching a three letter amino acid code to a number*/

typedef enum{

	 DADENINE, DCYTOSINE, DGUANINE, DTHYMINE

}enumDNA;   /*enumerated data structure matching a nucleic acid code to a number*/

typedef enum{

	 RADENINE, RCYTOSINE, RGUANINE, RURACIL

}enumRNA;   /*enumerated data structure matching a nucleic acid code to a number*/

typedef struct{

	int size;

	/* size of the array, or number of fragments */

	double bX, bN[3], XchH, nH, nD, fAccessH, mass[3], V;

	/* bX is the sum of X-ray scattering lengths; bN is the sum of neutron scattering lengths:
	 * bN[0] is the normal value; bN[1] is the change due to the deuteration level; bN[2] is the
	 * contribution from exchangeables; XchH is the number of exchangeable H atoms; nH is the number of H atoms
	 * nD is the number of D (2H) atoms; fAccessH is the number of accessible exchangeable H atoms;
	 * mass is the mass: the 3 components are analogous to those for bN; V is the volume */

}contrast_t;   /* data structure for contrast information */

molecule *aminoAcidData;   /* these are defined externally as they are not variables as such */
molecule *RNAData;         /* they are initialised and then used throughout the program */
molecule *DNAData;

int initAminoAcidData(void){   /* initialise the amino acid data */
	
	int i, j;

	/* i and j are counters */

	/* Amino acid compositions are without a H2O in their composition, so they are consistent with
	 * those involved in peptide bonds */

	MOLECULE(aminoAcidData, 0, AMINOACIDS, "aminoAcidData", "initaminoAcidData", i, j);
	
	/* Volumes taken from Tsai, JMB, 290, 253-266 (1999) */

	aminoAcidData[ALA].formula[0] = 'A';
	aminoAcidData[ALA].exH = 1.0;
	aminoAcidData[ALA].composition[H] = 5.0;
	aminoAcidData[ALA].composition[C] = 3.0;
	aminoAcidData[ALA].composition[N] = 1.0;
	aminoAcidData[ALA].composition[O] = 1.0;
	aminoAcidData[ALA].composition[S] = 0.0;
	aminoAcidData[ALA].vol = 89.27;
	
	aminoAcidData[CYS].formula[0] = 'C';
	aminoAcidData[CYS].exH = 2.0;
	aminoAcidData[CYS].composition[H] = 5.0;
	aminoAcidData[CYS].composition[C] = 3.0;
	aminoAcidData[CYS].composition[N] = 1.0;
	aminoAcidData[CYS].composition[O] = 1.0;
	aminoAcidData[CYS].composition[S] = 1.0;
	aminoAcidData[CYS].vol = 112.84; /* CSS has a volume of 102.50 */
	
	aminoAcidData[ASP].formula[0] = 'D';   /* assumed to be acidic form pK=4.5 */
	aminoAcidData[ASP].exH = 1.0;
	aminoAcidData[ASP].composition[H] = 4.0;
	aminoAcidData[ASP].composition[C] = 4.0;
	aminoAcidData[ASP].composition[N] = 1.0;
	aminoAcidData[ASP].composition[O] = 3.0;
	aminoAcidData[ASP].composition[S] = 0.0;
	aminoAcidData[ASP].vol = 114.43;
	
	aminoAcidData[GLU].formula[0] = 'E';   /* assumed to be acidic form pK=4.5 */
	aminoAcidData[GLU].exH = 1.0;
	aminoAcidData[GLU].composition[H] = 6.0;
	aminoAcidData[GLU].composition[C] = 5.0;
	aminoAcidData[GLU].composition[N] = 1.0;
	aminoAcidData[GLU].composition[O] = 3.0;
	aminoAcidData[GLU].composition[S] = 0.0;
	aminoAcidData[GLU].vol = 138.81;
	
	aminoAcidData[PHE].formula[0] = 'F';
	aminoAcidData[PHE].exH = 1.0;
	aminoAcidData[PHE].composition[H] = 9.0;
	aminoAcidData[PHE].composition[C] = 9.0;
	aminoAcidData[PHE].composition[N] = 1.0;
	aminoAcidData[PHE].composition[O] = 1.0;
	aminoAcidData[PHE].composition[S] = 0.0;
	aminoAcidData[PHE].vol = 190.84;
	
	aminoAcidData[GLY].formula[0] = 'G';
	aminoAcidData[GLY].exH = 1.0;
	aminoAcidData[GLY].composition[H] = 3.0;
	aminoAcidData[GLY].composition[C] = 2.0;
	aminoAcidData[GLY].composition[N] = 1.0;
	aminoAcidData[GLY].composition[O] = 1.0;
	aminoAcidData[GLY].composition[S] = 0.0;
	aminoAcidData[GLY].vol = 63.76;
	
	aminoAcidData[HIS].formula[0] = 'H';   /* assumed to be intermediate pK=6.5 */
	aminoAcidData[HIS].exH = 1.5;
	aminoAcidData[HIS].composition[H] = 6.5;
	aminoAcidData[HIS].composition[C] = 6.0;
	aminoAcidData[HIS].composition[N] = 3.0;
	aminoAcidData[HIS].composition[O] = 1.0;
	aminoAcidData[HIS].composition[S] = 0.0;
	aminoAcidData[HIS].vol = 157.46;
	
	aminoAcidData[ILE].formula[0] = 'I';
	aminoAcidData[ILE].exH = 1.0;
	aminoAcidData[ILE].composition[H] = 11.0;
	aminoAcidData[ILE].composition[C] = 6.0;
	aminoAcidData[ILE].composition[N] = 1.0;
	aminoAcidData[ILE].composition[O] = 1.0;
	aminoAcidData[ILE].composition[S] = 0.0;
	aminoAcidData[ILE].vol = 163.01;
		
	aminoAcidData[LYS].formula[0] = 'K';  /* assumed to be in basic form pK=10 */
	aminoAcidData[LYS].exH = 4.0;  
	aminoAcidData[LYS].composition[H] = 13.0;
	aminoAcidData[LYS].composition[C] = 6.0;
	aminoAcidData[LYS].composition[N] = 2.0;
	aminoAcidData[LYS].composition[O] = 1.0;
	aminoAcidData[LYS].composition[S] = 0.0;
	aminoAcidData[LYS].vol = 165.08;
	
	aminoAcidData[LEU].formula[0] = 'L';
	aminoAcidData[LEU].exH = 1.0;
	aminoAcidData[LEU].composition[H] = 11.0;
	aminoAcidData[LEU].composition[C] = 6.0;
	aminoAcidData[LEU].composition[N] = 1.0;
	aminoAcidData[LEU].composition[O] = 1.0;
	aminoAcidData[LEU].composition[S] = 0.0;
	aminoAcidData[LEU].vol = 163.09;
	
	aminoAcidData[MET].formula[0] = 'M';
	aminoAcidData[MET].exH = 1.0;
	aminoAcidData[MET].composition[H] = 9.0;
	aminoAcidData[MET].composition[C] = 5.0;
	aminoAcidData[MET].composition[N] = 1.0;
	aminoAcidData[MET].composition[O] = 1.0;
	aminoAcidData[MET].composition[S] = 1.0;
	aminoAcidData[MET].vol = 165.82;
	
	aminoAcidData[ASN].formula[0] = 'N';
	aminoAcidData[ASN].exH = 3.0;
	aminoAcidData[ASN].composition[H] = 6.0;
	aminoAcidData[ASN].composition[C] = 4.0;
	aminoAcidData[ASN].composition[N] = 2.0;
	aminoAcidData[ASN].composition[O] = 2.0;
	aminoAcidData[ASN].composition[S] = 0.0;
	aminoAcidData[ASN].vol = 122.35;	

	aminoAcidData[PRO].formula[0] = 'P';
	aminoAcidData[PRO].exH = 0.0;
	aminoAcidData[PRO].composition[H] = 7.0;
	aminoAcidData[PRO].composition[C] = 5.0;
	aminoAcidData[PRO].composition[N] = 1.0;
	aminoAcidData[PRO].composition[O] = 1.0;
	aminoAcidData[PRO].composition[S] = 0.0;
	aminoAcidData[PRO].vol = 121.29;
	
	aminoAcidData[GLN].formula[0] = 'Q';
	aminoAcidData[GLN].exH = 3.0;
	aminoAcidData[GLN].composition[H] = 8.0;
	aminoAcidData[GLN].composition[C] = 5.0;
	aminoAcidData[GLN].composition[N] = 2.0;
	aminoAcidData[GLN].composition[O] = 2.0;
	aminoAcidData[GLN].composition[S] = 0.0;
	aminoAcidData[GLN].vol = 146.91;	
	
	aminoAcidData[ARG].formula[0] = 'R';  /* assumed to be in basic form pK=12 */
	aminoAcidData[ARG].exH = 6.0;
	aminoAcidData[ARG].composition[H] = 13.0;
	aminoAcidData[ARG].composition[C] = 6.0;
	aminoAcidData[ARG].composition[N] = 4.0;
	aminoAcidData[ARG].composition[O] = 1.0;
	aminoAcidData[ARG].composition[S] = 0.0;
	aminoAcidData[ARG].vol = 190.33;
	
	aminoAcidData[SER].formula[0] = 'S';
	aminoAcidData[SER].exH = 2.0;
	aminoAcidData[SER].composition[H] = 5.0;
	aminoAcidData[SER].composition[C] = 3.0;
	aminoAcidData[SER].composition[N] = 1.0;
	aminoAcidData[SER].composition[O] = 2.0;
	aminoAcidData[SER].composition[S] = 0.0;
	aminoAcidData[SER].vol = 93.50;
	
	aminoAcidData[THR].formula[0] = 'T';
	aminoAcidData[THR].exH = 2.0;
	aminoAcidData[THR].composition[H] = 7.0;
	aminoAcidData[THR].composition[C] = 4.0;
	aminoAcidData[THR].composition[N] = 1.0;
	aminoAcidData[THR].composition[O] = 2.0;
	aminoAcidData[THR].composition[S] = 0.0;
	aminoAcidData[THR].vol = 119.61;
	
	aminoAcidData[VAL].formula[0] = 'V';
	aminoAcidData[VAL].exH = 1;
	aminoAcidData[VAL].composition[H] = 9.0;
	aminoAcidData[VAL].composition[C] = 5.0;
	aminoAcidData[VAL].composition[N] = 1.0;
	aminoAcidData[VAL].composition[O] = 1.0;
	aminoAcidData[VAL].composition[S] = 0.0;
	aminoAcidData[VAL].vol = 138.16;
	
	aminoAcidData[TRP].formula[0] = 'W';
	aminoAcidData[TRP].exH = 2.0;
	aminoAcidData[TRP].composition[H] = 10.0;
	aminoAcidData[TRP].composition[C] = 11.0;
	aminoAcidData[TRP].composition[N] = 2.0;
	aminoAcidData[TRP].composition[O] = 1.0;
	aminoAcidData[TRP].composition[S] = 0.0;
	aminoAcidData[TRP].vol = 226.38;
	
	aminoAcidData[TYR].formula[0] = 'Y';
	aminoAcidData[TYR].exH = 2.0;
	aminoAcidData[TYR].composition[H] = 9.0;
	aminoAcidData[TYR].composition[C] = 9.0;
	aminoAcidData[TYR].composition[N] = 1.0;
	aminoAcidData[TYR].composition[O] = 2.0;
	aminoAcidData[TYR].composition[S] = 0.0;
	aminoAcidData[TYR].vol = 194.63;
	
	return(1);
		
}

int initDNAData(void){   /* initialise the DNA data */

	int i, j;

	/* i and j are counters */

	MOLECULE(DNAData, 0, NUCLEOTIDES, "DNAData", "initDNAData", i, j);
		
	/* Volumes taken from Nadassy, Nuc. Acid Res., 29, 3362-3376 (2001) */
	
	DNAData[DADENINE].formula[0] = 'A';   /* quantities have been broken into nucleotide + backbone contributions */
	DNAData[DADENINE].exH = 2.0;
	DNAData[DADENINE].composition[H] = 4.0+7.0;
	DNAData[DADENINE].composition[C] = 5.0+5.0;
	DNAData[DADENINE].composition[N] = 5.0;
	DNAData[DADENINE].composition[O] = 0.0+5.0;
	DNAData[DADENINE].composition[P] = 0.0+1.0;
	DNAData[DADENINE].vol = 133.9+168.1;

	DNAData[DCYTOSINE].formula[0] = 'C';
	DNAData[DCYTOSINE].exH = 2.0;
	DNAData[DCYTOSINE].composition[H] = 4.0+7.0;
	DNAData[DCYTOSINE].composition[C] = 4.0+5.0;
	DNAData[DCYTOSINE].composition[N] = 3.0;
	DNAData[DCYTOSINE].composition[O] = 1.0+5.0;
	DNAData[DCYTOSINE].composition[P] = 0.0+1.0;
	DNAData[DCYTOSINE].vol = 115.6+168.1;
	
	DNAData[DGUANINE].formula[0] = 'G';
	DNAData[DGUANINE].exH = 3.0;
	DNAData[DGUANINE].composition[H] = 4.0+7.0;
	DNAData[DGUANINE].composition[C] = 5.0+5.0;
	DNAData[DGUANINE].composition[N] = 5.0;
	DNAData[DGUANINE].composition[O] = 1.0+5.0;
	DNAData[DGUANINE].composition[P] = 0.0+1.0;
	DNAData[DGUANINE].vol = 146.6+168.1;
	
	DNAData[DTHYMINE].formula[0] = 'T';
	DNAData[DTHYMINE].exH = 1.0;
	DNAData[DTHYMINE].composition[H] = 5.0+7.0;
	DNAData[DTHYMINE].composition[C] = 5.0+5.0;
	DNAData[DTHYMINE].composition[N] = 2.0;
	DNAData[DTHYMINE].composition[O] = 2.0+5.0;
	DNAData[DTHYMINE].composition[P] = 0.0+1.0;
	DNAData[DTHYMINE].vol = 133.6+168.1;

	return(1);
		
}

int initRNAData(void){   /* initialise the RNA data */

	int i, j;

	/* i and j are counters */

	MOLECULE(RNAData, 0, NUCLEOTIDES, "RNAData", "initRNAData", i, j);
		
	/* Volumes taken from Voss, JMB., 346, 477-492 (2005) */

	RNAData[RADENINE].formula[0] = 'A';   /* quantities have been broken into nucleotide + backbone contributions */
	RNAData[RADENINE].exH = 2.0+1.0;
	RNAData[RADENINE].composition[H] = 4.0+7.0;
	RNAData[RADENINE].composition[C] = 5.0+5.0;
	RNAData[RADENINE].composition[N] = 5.0;
	RNAData[RADENINE].composition[O] = 0.0+6.0;
	RNAData[RADENINE].composition[P] = 0.0+1.0;
	RNAData[RADENINE].vol = 139.2+176.1;

	RNAData[RCYTOSINE].formula[0] = 'C';
	RNAData[RCYTOSINE].exH = 2.0+1.0;
	RNAData[RCYTOSINE].composition[H] = 4.0+7.0;
	RNAData[RCYTOSINE].composition[C] = 4.0+5.0;
	RNAData[RCYTOSINE].composition[N] = 3.0;
	RNAData[RCYTOSINE].composition[O] = 1.0+6.0;
	RNAData[RCYTOSINE].composition[P] = 0.0+1.0;
	RNAData[RCYTOSINE].vol = 115.0+176.1;
	
	RNAData[RGUANINE].formula[0] = 'G';
	RNAData[RGUANINE].exH = 3.0+1.0;
	RNAData[RGUANINE].composition[H] = 4.0+7.0;
	RNAData[RGUANINE].composition[C] = 5.0+5.0;
	RNAData[RGUANINE].composition[N] = 5.0;
	RNAData[RGUANINE].composition[O] = 1.0+6.0;
	RNAData[RGUANINE].composition[P] = 0.0+1.0;
	RNAData[RGUANINE].vol = 145.9+176.1;
	
	RNAData[RURACIL].formula[0] = 'U';
	RNAData[RURACIL].exH = 1.0+1.0;
	RNAData[RURACIL].composition[H] = 3.0+7.0;
	RNAData[RURACIL].composition[C] = 4.0+5.0;
	RNAData[RURACIL].composition[N] = 2.0;
	RNAData[RURACIL].composition[O] = 2.0+6.0;
	RNAData[RURACIL].composition[P] = 0.0+1.0;
	RNAData[RURACIL].vol = 110.8+176.1;

	return(1);
		
}

int getSequenceData(molecule *mol){   /* Function that coordinates the extraction of formulas from sequences or molecular formulas */

	/* mol is where the molecule data is to be stored */

	int i, j, nProcessed = 0, nIncorrect = 0;

	/* i and j are counters; nProcessed is a counter to keep track of the number of residues/nucleotides
	 * nIncorrect counts the number of residues/nucleotides that do not have a valid code */

	double volume;
	
	/* volume is for temporary storage of volume information */
	
	char buffer[9999], tempChar;

	/* buffer is used for reading information from the input file; tempChar is a temporary character variable */

	fgets(buffer, 9999, stdin);
	sscanf(buffer,"%lf%s%s%lf", &mol->n, &tempChar, mol->formula, &mol->vol);

	if(mol->vol > 0.0) mol->status = 1;   /* if the molecular volume has been specified do not calculate the volume */
	else mol->status = 0;

	if(tempChar == 'P') mol->type = PROTEIN;
	else if(tempChar == 'R') mol->type = RNA;
	else if(tempChar == 'D') mol->type = DNA;
	else mol->type = MOLECULE;

	//fprintf(stdout,"\n<BR>tempChar : %c<BR>\n", tempChar);

	switch(mol->type){
		
		case PROTEIN:
			if(mol->status == 1) volume = mol->vol;
			for(j = 0; mol->formula[j] != '\0'; j ++){

				for(i = 0; i < AMINOACIDS; i ++)  /* check for lowercase or uppercase */
					if(mol->formula[j] == aminoAcidData[i].formula[0] || mol->formula[j] == aminoAcidData[i].formula[0] + 32){

						nProcessed ++;
						addMoleculeComponents(mol, &aminoAcidData[i]);   /* update the molecule information */
					
					}

				if(j != nProcessed + nIncorrect - 1){

					nIncorrect ++;
					fprintf(stdout,"\n%c (#%d) does not correspond to a valid protein residue.<BR><BR>\n", mol->formula[j], j + 1);

				}

			}

			mol->composition[H] += 2.0;   /* because all residues are assumed to be peptide bonded, the C-terminus */
			mol->composition[O] += 1.0;   /* has a hydroxyl group attached and the N-terminus has another hydrogen */
			mol->exH += 2.0;
			if(mol->status == 1) mol->vol = volume;
			else mol->vol += 2.0*getAtomicVol(H) + getAtomicVol(O);	
			
			fprintf(stdout,"Protein(%d residues)", nProcessed);
			break;
			
		case RNA:
			if(mol->status == 1) volume = mol->vol;
			for(j = 0; mol->formula[j] != '\0'; j ++){
				
				for(i = 0; i < NUCLEOTIDES; i ++)
					if(mol->formula[j] == RNAData[i].formula[0] || mol->formula[j] == RNAData[i].formula[0] + 32){

						addMoleculeComponents(mol, &RNAData[i]);
						nProcessed ++;

					}

				if(j != nProcessed + nIncorrect - 1){

					nIncorrect ++;
					fprintf(stdout,"\n%c (#%d) does not correspond to a valid RNA nucleotide.<BR><BR>\n", mol->formula[j], j + 1);

				}

			}
			
			mol->composition[P] -= 1.0;   /* this assumes that the 5' end of the RNA or DNA is not phosphorylated */
			mol->composition[O] -= 3.0;
			mol->composition[H] += 2.0;   /* these cap the ends of the protein (or DNA or RNA), OH on one end H on the other */
			mol->composition[O] += 1.0;
			mol->exH += 2.0;
			if(mol->status == 1) mol->vol = volume;
			else mol->vol += 2.0*getAtomicVol(H) - getAtomicVol(O)- getAtomicVol(P);
			
			fprintf(stdout,"RNA(%d nucleotides)", nProcessed);
			break;

		case DNA:
			if(mol->status == 1) volume = mol->vol;
			for(j = 0; mol->formula[j] != '\0'; j ++){

				for(i = 0; i < NUCLEOTIDES; i ++)
					if(mol->formula[j] == DNAData[i].formula[0] || mol->formula[j] == DNAData[i].formula[0] + 32){

						addMoleculeComponents(mol, &DNAData[i]);
						nProcessed ++;

					}

				if(j != nProcessed + nIncorrect - 1){

					nIncorrect ++;
					fprintf(stdout,"\n%c (#%d) does not correspond to a valid DNA nucleotide.<BR><BR>\n", mol->formula[j], j + 1);

				}

			}

			mol->composition[P] -= 1.0;   /* this assumes that the 5' end of the RNA or DNA is not phosphorylated */
			mol->composition[O] -= 3.0;
			mol->composition[H] += 2.0;   /* these cap the ends of the protein (or DNA or RNA), OH on one end H on the other */
			mol->composition[O] += 1.0;
			mol->exH += 2.0;
			if(mol->status == 1) mol->vol = volume;
			else mol->vol += 2.0*getAtomicVol(H) - getAtomicVol(O)- getAtomicVol(P);
			
			fprintf(stdout,"DNA(%d nucleotides)", j);
			break;
		
		case MOLECULE:
			getComposition(mol);
			//fprintf(stdout,"getSequenceData, MOLECULE, after getComposition, mol->status = %d", mol->status );
			if(mol->status == 0) getMolecularVol(mol);
			j = 1;
			
			fprintf(stdout,"<P>Molecule(");
		
			for(i = 0; i < TABLESIZE; i ++){

				if(mol->composition[i] == 1.0) fprintf(stdout,"%s",getElementLabel(i));
				else if(mol->composition[i] > 0.0) fprintf(stdout,"%s<sub>%.0lf</sub>",getElementLabel(i), mol->composition[i]);

			}

			fprintf(stdout,")");
			break;
		
		default:
			fprintf(stderr, "<P>The type of molecule specified was invalid!<BR>");
			return(0);
	
	}
	
	return(j);   /* returns the number of residues, or molecules processed */
	
}

int getI0Data(double *X, double *Y, double *sigY, int n){   /* reads I0 data from input and approximately determines the matchpoint */

	/* X is fD2O; Y is I(0) on entry, and is +/- sqrt(I(0)/c) on exit; sigY is the ESD of Y; n is the number of contrast points */
	
	char buffer[256];

	/* buffer is used for reading data from the input file */

	int i, j;
	
	/* i and j  are counters */

	double *abc, **corr_abc, conc[n], approxMatch;

	/* abc are the polynomial coefficients; corr_abc is the correlation matrix (unused);
	 * conc is the concentration at each contrast point; approxMatch is the approximate match-point */

	for(i = 0; i < n; i ++){

		fgets(buffer, 256, stdin);
		sscanf(buffer,"%lf%lf%lf%lf", &X[i], &Y[i], &sigY[i], &conc[i]);
		Y[i] = Y[i]/conc[i];
		sigY[i] = sigY[i]/conc[i];

	}
	
	DOUBLE(abc, 3, "abc", "getI0Data", i);
	DOUBLE2D(corr_abc, 3, 3, "corr_abc", "getI0Data", i, j);
	
	polynomialFit(2, X, Y, sigY, n, abc, corr_abc, 1);

	approxMatch = -abc[1]/2.0/abc[0];

	fprintf(stdout,"\t\t<P>Approximate <i>X</i>-intercept = %.3lf (Determined from 2<sup>nd</sup> order polynomial fit to <i>I</i>(0)/<i>c</i> <i>vs.</i> <i>f</i><sub>D<sub>2</sub>O</sub>)<BR>\n", approxMatch);

	for(i = 0; i < n; i ++){
		
		sigY[i] = sigY[i]/2.0/sqrt(Y[i]);

		if(X[i] < approxMatch) Y[i] = sqrt(Y[i]);
		else Y[i] = -sqrt(Y[i]);   /* for those points after the matchpoint swap the sign */

	}
	
	free(abc);
	free(corr_abc);

	return(1);
	
}

int convertMolecularData(molecule *mol, contrast_t *p){   /* converts molecular data into contrast data */

	/* mol is the molecular data; p is the contrast data */

	int i;

	/* i is a counter */

	double H, D;

	/* H and D refer to the amount of hydrogen and deuterium in the molecule */
	
	for(i = 0; i < mol->size; i ++){
		
		H = mol[i].composition[1];
		D = mol[i].composition[0];

		mol[i].composition[1] = H + D;
		mol[i].composition[0] = 0.0;
		
		p->nH += mol[i].n*mol[i].composition[1];


		p->bN[0] += mol[i].n*getNFScat(&mol[i]);   /* these are the normal total scattering lengths and masses, free of isotope effects */
		p->bX += mol[i].n*getXFScat(&mol[i]);
		p->mass[0] += mol[i].n*getMolecularMass(&mol[i])/1000.0;
		
		mol[i].composition[1] = H;
		mol[i].composition[0] = D;

		p->XchH += mol[i].n*mol[i].exH;
		p->V += mol[i].n*mol[i].vol;

	}

	p->bN[1] = (p->nH - p->XchH)*10.412;   /* 10.412 = bD - bH */
	p->bN[2] = p->XchH*10.412;

	p->mass[1] = (p->nH - p->XchH)*1.006/1000.0;   /* 1.006 = mD - mH */
	p->mass[2] = p->XchH*p->fAccessH*1.006/1000.0;

	return(1);
}

int printMoleculeData(molecule mol){   /* prints the contrast data for each molecule */

	double vol, bX, bN, rhoX, rhoN, mass;

	/* vol is the volume; bX is the total X-ray scattering length; bN is the total neutron scattering length;
	 * rhoX is the X-ray scattering length density; rhoN is the neutron scattering length density; mass is the mass */

	vol = mol.vol;
	bX = getXFScat(&mol);
	bN = getNFScat(&mol);
	rhoX = bX/vol*10.0;
	rhoN = bN/vol*10.0;
	mass = getMolecularMass(&mol)/1000.0;
	
	fprintf(stdout,"\t\t\t<TR align=center><TH>%.2lf<TD>%.1lf<TD>%.2lf<TD>%.2lf<TD>%.3lf<TD>%.3lf<TD>%.3lf\n", mol.n, vol, bX, bN, rhoX, rhoN, mass);

	return(1);

}
int main(void){
	
	int nMol1, nMol2, nSolv = 0, i, j, nI0;

	/* nMol1 is the number of fragments in molecule 1; nMol2 is the number of fragments in molecule 2
	 * nSolv is the number of components in the solvent (in addition to water); i and j are counters
	 * nI0 is the number of contrast points used to determine the match-point  */

	char buffer[999], calcTitle[256], tempChar;

	/* buffer is used to read information from the input file, calcTitle is the title of the calculation
	 * tempchar is a temporary character */

	double soluteVol = 0.0, fD1, fD2, rho1, rho2, rhoS, rho, mass1, mass2, slope[3], intercept[3], densityS = 0.0, psv = 0.0;

	/* soluteVol is the volume of the solutes in solution; fD1 is the deuteration level of molecule 1; 
	 * fD2 is the deuteration level of molecule 2; rho1 is the scattering length density (SLD) of molecule 1 
	 * rho2 is the SLD of molecule 2; rhoS is the SLD of the solvent; rho is the total SLD of the particle 
	 * mass1 and mass2 are the masses of molecule 1 and molecule 2 respectively; intercept and slope are
	 * used as temporary storage when printing out data; densityS is the density of the solvent; psv is the 
	 * partial specific volume */

	molecule *m1, *m2, *mS;

	/* m1 is molecule 1; molecule 2 is molecule 2; mS is the solvent */

	contrast_t *p1, *p2, *s;

	/*  p1 is the SL information for molecule 1; p2 and s as above */

	double *X, *Y, *sigY, *ab, **corr_ab, chi2;

	/* X is fD2O; Y is +/- sqrt(I(0))/conc; sigY is the ESD of Y; ab is the slope and y-intercept of a linear line of best fit
	 * corr_ab is the correlation matrix; chi2 is chi^2*/

	double errXint = 0.0, dFdP[2];

	/* errXint is the ESD of the X-intercept; dFdP are derivates used to calculate the error in the X-intercept */

	initAminoAcidData();   /* initialise the amino acid data */
	initDNAData();   /* initialise the DNA data */
	initRNAData();   /* initialise the RNA data */
	
	CONTRAST(p1, 0, 1, "p1", "main", i);
	CONTRAST(p2, 0, 1, "p2", "main", i);
	CONTRAST(s, 0, 1, "s", "main", i);
	
	fgets(calcTitle, 256, stdin);   /* read in the title of the calculation */

//	fprintf(stdout,"<!DOCTYPE html PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">\n<HTML>\n");
//	fprintf(stdout,"\t<HEAD>\n\t\t<TITLE>\n\t\t\tProject: %s\t\t</TITLE>\n", calcTitle);
//	fprintf(stdout,"\t\t<META HTTP-EQUIV=\"Content-Type\" CONTENT=\"text/html; charset=iso-8859-1\">\n");
//	fprintf(stdout,"\t</HEAD>\n");
//	fprintf(stdout,"\t<BODY>\n");
	fprintf(stdout,"Project: %s\n", calcTitle);

	/* do the I(0) analysis */

	fprintf(stdout,"\t\t<P><B>Analysis of the dependence of <i>I</i>(0) upon <i>f</i><sub>D<sub>2</sub>O</sub></B><BR>\n");
	fgets(buffer, 256, stdin);
	sscanf(buffer,"%d", &nI0);
	
	if(nI0 < 2) fprintf(stdout,"\t\t<P><i>I</i>(0) analysis not performed...<BR>\n");

	else{

		DOUBLE(X, nI0, "X", "main", i);
		DOUBLE(Y, nI0, "Y", "main", i);
		DOUBLE(sigY, nI0, "sigY", "main", i);
	
		getI0Data(X, Y, sigY, nI0);
	
		DOUBLE(ab, 2, "ab", "main", i)
		DOUBLE2D(corr_ab, 2, 2, "corr_ab", "main", i, j);
	
		chi2 = polynomialFit(1, X, Y, sigY, nI0, ab, corr_ab, 1);   /* determine the matchpoint */

		dFdP[0] = ab[1]/ab[0]/ab[0];
		dFdP[1] = -1.0/ab[0];

		for(i = 0; i < 2; i ++)
			for(j = 0; j < 2; j ++)
				errXint += chi2*dFdP[i]*corr_ab[i][j]*dFdP[j];


		fprintf(stdout,"\t\t<P>Linear fit to &radic;<i>I</i>(0)/<i>c</i> <i>vs.</i> <i>f</i><sub>D<sub>2</sub>O</sub> (sign of &radic;<i>I</i>(0)/<i>c</i> for <i>f</i><sub>D<sub>2</sub>O</sub> &gt; approximate <i>X</i>-intercept adjusted)<BR><BR>\n");
		
		fprintf(stdout,"\t\t<TABLE align=center width=\"30%%\" border=1>\n");
		fprintf(stdout,"\t\t\t<CAPTION><EM>I(0) analysis</EM></CAPTION>\n");
		fprintf(stdout,"\t\t\t<TR>");
		fprintf(stdout,"<TH width=\"33%%\">Parameter");
		fprintf(stdout,"<TH width=\"33%%\">Value");
		fprintf(stdout,"<TH width=\"33%%\">ESD\n");

		fprintf(stdout,"\t\t\t<TR align=center><TH><i>&chi;</i><sup>2</sup><TD>%.2lf<TD>-\n", chi2);
		fprintf(stdout,"\t\t\t<TR align=center><TH>Slope<TD>%.3lf<TD>%.3lf\n", ab[0], sqrt(chi2*corr_ab[0][0]));
		fprintf(stdout,"\t\t\t<TR align=center><TH><i>Y</i>-intercept<TD>%.3lf<TD>%.3lf\n", ab[1], sqrt(chi2*corr_ab[1][1]));
		fprintf(stdout,"\t\t\t<TR align=center><TH><i>X</i>-intercept<TD>%.3lf<TD>%.3lf\n", -ab[1]/ab[0], sqrt(errXint));
		fprintf(stdout,"\t\t</TABLE><BR>\n");

	}

	fgets(buffer, 256, stdin);
	sscanf(buffer, "%d", &nSolv);

	nSolv ++;   /* increment nSolvent so that the number of molecules now includes water */

	MOLECULE(mS, 0, nSolv, "mS", "main", i, j);
	
	fprintf(stdout,"\t\t<P><B>Contrast calculation parameters</B><BR>\n");
	
	fprintf(stdout,"\t\t<P>In addition to water the solvent contains %d other substance(s)<BR>\n", nSolv - 1);

	fprintf(stdout,"\t\t<OL>\n");
	
	mS[0].n = MOLESOFWATERPERLITRE;
	mS[0].exH = 2.0;
	mS[0].composition[1] = 2.0;
	mS[0].composition[8] = 1.0;
	mS[0].vol = 1.0E-3/NAVAGADRO/MOLESOFWATERPERLITRE*1.0E+30;  /* 1.0E-3 is m^3 in a litre, and 1E+30 converts m^3 to angstroms */
		
	s->fAccessH = 1.0;
	
	for(i = 1; i < mS->size; i ++){   /* Read in and process the solvent information */
		
		fprintf(stdout,"\t\t\t<LI>");
		getSequenceData(&mS[i]);   /* Any exchange between water and small molecules in the solvent is not considered */

		fprintf(stdout,", Conc. = %.3lf mol/L", mS[i].n);
		if(mS[i].status == 1) fprintf(stdout,", Volume of each formula unit = %.2lf &Aring;<sup>3</sup> (Entered)\n", mS[i].vol);
		else fprintf(stdout,", Volume of each formula unit to be calculated<BR>\n");

	}
	
	//fprintf(stdout,"<br>\n");
	for(i = 1; i < mS->size; i ++) {
		soluteVol += mS[i].n*mS[i].vol*NAVAGADRO;
		//fprintf(stdout,"calculate solute volume : i : %i, n : %.13f, volume : %.13f, soluteVol += n * vol * Avogadro = %.13fE23<br>\n", i, mS[i].n, mS[i].vol, soluteVol/1.0E23);
	}

	//fprintf(stdout,"calculate water for solute volume : n : %.13f, n = n * (1E27 - soluteVol %.13f) / 1E27", mS[0].n, soluteVol);
	mS[0].n = mS[0].n*(1.0E+27 - soluteVol)/1.0E+27;  /* 1.0E+27 is the number of cubic angstroms per litre */ 
	//fprintf(stdout,"= %.13f<br>\n", mS[0].n);
	
	/* If the concentration of solute is large, then the adjustments must be made for the volume of water
	 * they displace. i.e. if the concentration of solute is large, the concentration of water is no longer 55.5 mol/L */

	for(i = 0; i < mS->size; i ++) {
		densityS += mS[i].n*getMolecularMass(&mS[i])/1000.0;
		//fprintf(stdout,"calculate densityS : i : %i, n : %.13f, MolecularMass : %.13f, densityS += n * MolecularMass / 1000 = %.13f<br>\n", i, mS[i].n, getMolecularMass(&mS[i]), densityS);
	}

	convertMolecularData(mS, s);

	fprintf(stdout,"\t\t</OL>\n");

	fgets(buffer, 256, stdin);
	sscanf(buffer,"%lf", &p1->fAccessH);
		
	fgets(buffer, 256, stdin);
	sscanf(buffer,"%lf", &fD1);
	
	fgets(buffer, 256,stdin);
	sscanf(buffer, "%d", &nMol1);	
	
	MOLECULE(m1, 0, nMol1, "m1", "main", i, j);
	
	fprintf(stdout,"\t\t<P>Subunit (1) - possesses %d component(s), with a deuteration level of %.0lf%% and %.0lf%% of the total exchangeable protons accessible by the solvent<BR>\n", nMol1, fD1*100.0, p1->fAccessH*100.0);
	fprintf(stdout,"\t\t<OL>\n");

	for(i = 0; i < m1->size; i ++){
		
		fprintf(stdout,"\t\t\t<LI>");
		getSequenceData(&m1[i]);   /* transform a sequence of formula into a molecular compostion */
		
		m1[i].composition[0] = fD1*(m1[i].composition[1] - m1[i].exH);
		m1[i].composition[1] -= m1[i].composition[0];
		
		fprintf(stdout,", %.0lf formula  unit(s)", m1[i].n);
		
		if(m1[i].status == 1) fprintf(stdout,", Volume of each formula unit = %.2lf &Aring;<sup>3</sup> (Entered)\n", m1[i].vol);
		else fprintf(stdout,", Volume of each formula unit to be calculated<BR>\n");
		

	}

	convertMolecularData(m1, p1);

	fprintf(stdout,"\t\t</OL>\n");

	fgets(buffer, 256, stdin);
	sscanf(buffer,"%lf", &p2->fAccessH);
		
	fgets(buffer, 256, stdin);
	sscanf(buffer,"%lf", &fD2);
	
	fgets(buffer, 256,stdin);
	sscanf(buffer, "%d", &nMol2);	

	MOLECULE(m2, 0, nMol2, "m2", "main", i, j);
	
	fprintf(stdout,"\t\t<P>Subunit (2) - possesses %d components(s), with a deuteration level of %.0lf%% and %.0lf%% of the total exchangeable protons accessible by the solvent<BR>\n", nMol2, fD2*100.0, p2->fAccessH*100.0);

	fprintf(stdout,"\t\t<OL>\n");

	for(i = 0; i < m2->size; i ++){
		
		fprintf(stdout,"\t\t\t<LI>");
		getSequenceData(&m2[i]);   /* transform a sequence of formula into a molecular compostion */
		
		m2[i].composition[0] = fD2*(m2[i].composition[1] - m2[i].exH);
		m2[i].composition[1] -= m2[i].composition[0];

		fprintf(stdout,", %.0lf formula unit(s)", m2[i].n);
		
		if(m2[i].status == 1) fprintf(stdout,", Volume of each formula unit = %.2lf &Aring;<sup>3</sup> (Entered)\n", m2[i].vol);
		else fprintf(stdout,", Volume of each formula unit to be calculated<BR>\n");
		
	}

	convertMolecularData(m2, p2);

	fprintf(stdout,"\t\t</OL>\n");
	
	fprintf(stdout,"\t\t<P><BR><B>Contrast calculation results</B><BR><BR>\n");
	
	/* Print out this data mainly for the puposes of checking the calculation */

	fprintf(stdout,"\t\t<TABLE align=center width=\"80%%\" border=1>\n");
	fprintf(stdout,"\t\t\t<CAPTION><EM>Solvent and particle information in H<sub>2</sub>O</EM></CAPTION>\n");
	fprintf(stdout,"\t\t\t<TR>");
	fprintf(stdout,"<TH>");
	fprintf(stdout,"<TH>Volume<BR>(&Aring;<sup>3</sup>)");
	fprintf(stdout,"<TH>&Sigma;<i>b<sub>X</sub><BR></i>(10<sup>-13</sup>cm)");
	fprintf(stdout,"<TH>&Sigma;<i>b<sub>N</sub><BR></i>(10<sup>-13</sup>cm)");
	fprintf(stdout,"<TH><i>&rho;<sub>X</sub></i><BR>(10<sup>10</sup>cm<sup>-2</sup>)");
	fprintf(stdout,"<TH><i>&rho;<sub>N</sub></i><BR>(10<sup>10</sup>cm<sup>-2</sup>)");
	fprintf(stdout,"<TH>Mass<BR>(kDa)\n");

	fprintf(stdout,"\t\t\t<TR align=center><TH colspan=\"7\">Solvent\n");

	for(i = 0; i < mS->size; i ++)
		printMoleculeData(mS[i]);
	
	fprintf(stdout,"\t\t\t<TR align=center><TH>Total<TD><B>-</B><TD><B>%.2lf</B><TD><B>%.2lf</B><TD><B>%.3lf</B><TD><B>%.3lf</B><TD><B>-</B>\n", s->bX, s->bN[0], s->bX/s->V*10.0, s->bN[0]/s->V*10.0);

	fprintf(stdout,"\t\t\t<TR align=center><TH colspan=\"7\">Subunit (1) - Number of hydrogen atoms = %.1lf, Number of exchangeable hydrogen atoms = %.1lf\n", p1->nH + p1->nD, p1->XchH);

	for(i = 0; i < m1->size; i ++)
		printMoleculeData(m1[i]);
	
	fprintf(stdout,"\t\t\t<TR align=center><TH>Total<TD><B>%.1lf</B><TD><B>%.2lf</B><TD><B>%.2lf</B><TD><B>%.3lf</B><TD><B>%.3lf</B><TD><B>%.3lf</B>\n", p1->V, p1->bX, p1->bN[0]+p1->bN[1]*fD1, p1->bX/p1->V*10.0, (p1->bN[0]+p1->bN[1]*fD1)/p1->V*10.0, p1->mass[0] + p1->mass[1]*fD1);

	fprintf(stdout,"\t\t\t<TR align=center><TH colspan=\"7\">Subunit (2) - Number of hydrogen atoms = %.1lf, Number of exchangeable hydrogen atoms = %.1lf\n", p2->nH + p2->nD, p2->XchH);

	for(i = 0; i < m2->size; i ++)
		printMoleculeData(m2[i]);

	fprintf(stdout,"\t\t\t<TR align=center><TH>Total<TD><B>%.1lf</B><TD><B>%.2lf</B><TD><B>%.2lf</B><TD><B>%.3lf</B><TD><B>%.3lf</B><TD><B>%.3lf</B>\n", p2->V, p2->bX, p2->bN[0]+p2->bN[1]*fD2, p2->bX/p2->V*10.0, (p2->bN[0]+p2->bN[1]*fD2)/p2->V*10.0, p2->mass[0] + p2->mass[1]*fD2);
	
	fprintf(stdout,"\t\t</TABLE>\n");   /* end of the table */
	
	fprintf(stdout,"\t\t<P>The partial specific volume of the complex (fully proteated) is: %.3lf cm<sup>3</sup>g<sup>-1</sup>\n",(p1->V+p2->V)/(p1->mass[0]+p2->mass[0])*NAVAGADRO*1.0E-27);

	/* scattering length density information */

	/* The scattering lengths (and masses) are broken into 3 components:
	 * 1: That for an object with no deuterium - (X->bN[0])
	 * 2: A correction for the deuteration level (fD) - (X->bN[1])
	 * 3: A correction for the deuterium content of the solvent (fD2O) and the accessibility of exchangeable protons (fAccessH) - (X->bN[2])
	 * Conversion to scattering length densities is then achieved via dividing the scattering lengths by the volume (X->V)
	 * Multiplication by a factor of 10 brings the SLDs onto the correct scale (i.e. 10^10 cm^-2) */

	fprintf(stdout,"\t\t<P>The neutron scattering length density (10<sup>10</sup>cm<sup>-2</sup>) at any D<sub>2</sub>O fraction is:<BR>\n");
	fprintf(stdout,"\t\t<UL>\n");
	fprintf(stdout,"\t\t\t<LI><i>&rho;</i><sub>1</sub> =  %.3lf +  %.3lf<i>f</i><sub>D,1</sub> + %.3lf<i>f</i><sub>AccessH,1</sub><i>f</i><sub>D<sub>2</sub>O</sub>\n", p1->bN[0]/p1->V*10.0, p1->bN[1]/p1->V*10.0, p1->bN[2]/p1->V*10.0);
	fprintf(stdout," =  %.3lf + %.3lf<i>f</i><sub>D<sub>2</sub>O</sub>\n", p1->bN[0]/p1->V*10.0 + p1->bN[1]/p1->V*10.0*fD1, p1->bN[2]*p1->fAccessH/p1->V*10.0);
	fprintf(stdout,"\t\t\t<LI><i>&rho;</i><sub>2</sub> =  %.3lf +  %.3lf<i>f</i><sub>D,2</sub> + %.3lf<i>f</i><sub>AccessH,2</sub><i>f</i><sub>D<sub>2</sub>O</sub>\n", p2->bN[0]/p2->V*10.0, p2->bN[1]/p2->V*10.0, p2->bN[2]/p2->V*10.0);
	fprintf(stdout," =  %.3lf + %.3lf<i>f</i><sub>D<sub>2</sub>O</sub>\n", p2->bN[0]/p2->V*10.0 + p2->bN[1]/p2->V*10.0*fD2, p2->bN[2]*p2->fAccessH/p2->V*10.0);
//	fprintf(stdout,"\t\t\t<LI><i>&rho;</i>  =  %.3lf +  %.3lf<i>f</i><sub>D,1</sub> + %.3lf<i>f</i><sub>D,2</sub> + %.3lf<i>f</i><sub>AccessH,1</sub><i>f</i><sub>D<sub>2</sub>O</sub> + %.3lf<i>f</i><sub>AccessH,2</sub><i>f</i><sub>D<sub>2</sub>O</sub>\n", (p1->bN[0] + p2->bN[0])/(p1->V + p2->V)*10.0, p1->bN[1]/(p1->V + p2->V)*10.0, p2->bN[1]/(p1->V + p2->V)*10.0, p1->bN[2]/(p1->V + p2->V)*10.0, p2->bN[2]/(p1->V + p2->V)*10.0);
//	fprintf(stdout," =  %.3lf + %.3lf<i>f</i><sub>D<sub>2</sub>O</sub>\n", (p1->bN[0] + p2->bN[0] + p1->bN[1]*fD1 + p2->bN[1]*fD2)/(p1->V + p2->V)*10.0, (p1->bN[2]*p1->fAccessH + p2->bN[2]*p2->fAccessH)/(p1->V + p2->V)*10.0);
	fprintf(stdout,"\t\t\t<LI><i>&rho;</i><sub>s</sub> =  %.3lf +  %.3lf<i>f</i><sub>D<sub>2</sub>O</sub>\n", s->bN[0]/s->V*10.0, s->bN[2]/s->V*10.0);
	fprintf(stdout,"\t\t</UL>\n");

	/* contrast information */

	fprintf(stdout,"\t\t<P>The neutron contrast (10<sup>10</sup>cm<sup>-2</sup>) at any D<sub>2</sub>O fraction is:<BR>\n");
	fprintf(stdout,"\t\t<UL>\n");
	fprintf(stdout,"\t\t\t<LI>&Delta;<i>&rho;</i><sub>1</sub> = %.3lf +  %.3lf<i>f</i><sub>D,1</sub> - %.3lf<i>f</i><sub>AccessH,1</sub><i>f</i><sub>D<sub>2</sub>O</sub>\n", (p1->bN[0]/p1->V - s->bN[0]/s->V)*10.0, p1->bN[1]/p1->V*10.0,  -(p1->bN[2]/p1->V - s->bN[2]/s->V)*10.0);
	fprintf(stdout," = %.3lf - %.3lf<i>f</i><sub>D<sub>2</sub>O</sub>\n", (p1->bN[0]/p1->V - s->bN[0]/s->V +  p1->bN[1]/p1->V*fD1)*10.0, -(p1->bN[2]/p1->V*p1->fAccessH - s->bN[2]/s->V)*10.0);
	fprintf(stdout,"\t\t\t<LI>&Delta;<i>&rho;</i><sub>2</sub> = %.3lf +  %.3lf<i>f</i><sub>D,2</sub> - %.3lf<i>f</i><sub>AccessH,2</sub><i>f</i><sub>D<sub>2</sub>O</sub>\n", (p2->bN[0]/p2->V - s->bN[0]/s->V)*10.0, p2->bN[1]/p2->V*10.0,  -(p2->bN[2]/p2->V - s->bN[2]/s->V)*10.0);
	fprintf(stdout," = %.3lf - %.3lf<i>f</i><sub>D<sub>2</sub>O</sub>\n", (p2->bN[0]/p2->V - s->bN[0]/s->V +  p2->bN[1]/p2->V*fD2)*10.0, -(p2->bN[2]/p2->V*p2->fAccessH - s->bN[2]/s->V)*10.0);
	fprintf(stdout,"\t\t\t<LI>&Delta;<i>&rho;</i> =  %.3lf +  %.3lf<i>f</i><sub>D,1</sub> +  %.3lf<i>f</i><sub>D,2</sub> + %.3lf<i>f</i><sub>AccessH,1</sub></i><i>f</i><sub>D<sub>2</sub>O</sub> + %.3lf<i>f</i><sub>AccessH,2</sub></i><i>f</i><sub>D<sub>2</sub>O</sub> - %.3lf<i>f</i><sub>D<sub>2</sub>O</sub>\n", ((p1->bN[0] +  p2->bN[0])/(p1->V + p2->V) - s->bN[0]/s->V)*10.0, p1->bN[1]/(p1->V + p2->V)*10.0, p2->bN[1]/(p1->V + p2->V)*10.0, p1->bN[2]/(p1->V + p2->V)*10.0, p2->bN[2]/(p1->V + p2->V)*10.0, s->bN[2]/s->V*10.0);
	fprintf(stdout," =  %.3lf - %.3lf<i>f</i><sub>D<sub>2</sub>O</sub>\n", ((p1->bN[0] +  p2->bN[0] + p1->bN[1]*fD1 +  p2->bN[1]*fD2)/(p1->V + p2->V) - s->bN[0]/s->V)*10.0, -((p1->bN[2]*p1->fAccessH +  p2->bN[2]*p2->fAccessH)/(p1->V + p2->V) - s->bN[2]/s->V)*10.0);
	fprintf(stdout,"\t\t</UL>\n");

	/* table of representative values (fD2O = 0, 0.1, ..., 1.0) of the SLD, contrast, including X-ray values */

	fprintf(stdout,"\t\t<BR><TABLE align=center width=\"80%%\" border=1>\n");
	fprintf(stdout,"\t\t\t<CAPTION><EM>Tabulated scattering length densities and contrasts</EM></CAPTION>\n");
	fprintf(stdout,"\t\t\t<TR>");
	fprintf(stdout,"<TH rowspan=2>");
	fprintf(stdout,"<TH colspan=3><i>&rho;</i> (10<sup>10</sup>cm<sup>-2</sup>)");
	fprintf(stdout,"<TH colspan=3>&Delta;<i>&rho;</i> (10<sup>10</sup>cm<sup>-2</sup>)");
//	fprintf(stdout,"<TH colspan=2>Mass (kDa)\n");
	fprintf(stdout,"\t\t\t<TR>");
	fprintf(stdout,"<TH>1");
	fprintf(stdout,"<TH>2");
	fprintf(stdout,"<TH>Solvent");
	fprintf(stdout,"<TH>1");
	fprintf(stdout,"<TH>2");
	fprintf(stdout,"<TH>Total\n");
//	fprintf(stdout,"<TH>1");
//	fprintf(stdout,"<TH>2\n");
//	fprintf(stdout,"\t\t\t<TR align=center><TH>X-RAY<TD>%.3lf<TD>%.3lf<TD>%.3lf<TD>%.3lf<TD>%.3lf<TD>%.3lf<TD>%.3lf<TD>%.3lf\n",p1->bX/p1->V*10.0,p2->bX/p2->V*10.0,s->bX/s->V*10.0,(p1->bX/p1->V - s->bX/s->V)*10.0, (p2->bX/p2->V - s->bX/s->V)*10.0,((p1->bX+p2->bX)/(p1->V+p2->V) - s->bX/s->V)*10.0,p1->mass[0],p2->mass[0]);
	fprintf(stdout,"\t\t\t<TR align=center><TH>X-RAY<TD>%.3lf<TD>%.3lf<TD>%.3lf<TD>%.3lf<TD>%.3lf<TD>%.3lf\n",p1->bX/p1->V*10.0,p2->bX/p2->V*10.0,s->bX/s->V*10.0,(p1->bX/p1->V - s->bX/s->V)*10.0, (p2->bX/p2->V - s->bX/s->V)*10.0,((p1->bX+p2->bX)/(p1->V+p2->V) - s->bX/s->V)*10.0);
	fprintf(stdout,"\t\t\t<TR><TH>NEUTRON\n");
	for(i = 0; i <= 100; i += 10){

		rho1 = (p1->bN[0] + p1->bN[1]*fD1 + p1->bN[2]*p1->fAccessH*(double)i/100.0)/p1->V*10.0;
		rho2 = (p2->bN[0] + p2->bN[1]*fD2 + p2->bN[2]*p2->fAccessH*(double)i/100.0)/p2->V*10.0;
		rho = (p1->bN[0] + p1->bN[1]*fD1 + p1->bN[2]*p1->fAccessH*(double)i/100.0 + p2->bN[0] + p2->bN[1]*fD2 + p2->bN[2]*p2->fAccessH*(double)i/100.0)/(p1->V + p2->V)*10.0;
		rhoS = (s->bN[0] + s->bN[2]*(double)i/100.0)/s->V*10.0;
//		mass1 = p1->mass[0] + p1->mass[1]*fD1 + p1->mass[2]*p1->fAccessH*(double)i/100.0; 
//		mass2 = p2->mass[0] + p2->mass[1]*fD2 + p2->mass[2]*p2->fAccessH*(double)i/100.0; 

//		fprintf(stdout,"\t\t\t<TR align=center><TH>%.1lf<TD>%.3lf<TD>%.3lf<TD>%.3lf<TD>%.3lf<TD>%.3lf<TD>%.3lf<TD>%.3lf<TD>%.3lf\n",(double)i/100.0, rho1, rho2, rhoS, rho1 - rhoS, rho2 - rhoS, rho - rhoS, mass1, mass2);
		fprintf(stdout,"\t\t\t<TR align=center><TH>%.1lf<TD>%.3lf<TD>%.3lf<TD>%.3lf<TD>%.3lf<TD>%.3lf<TD>%.3lf\n",(double)i/100.0, rho1, rho2, rhoS, rho1 - rhoS, rho2 - rhoS, rho - rhoS);
	       
	}

	intercept[0] = (p1->bN[0] + p1->bN[1]*fD1)/p1->V - s->bN[0]/s->V;
	slope[0] = p1->bN[2]*p1->fAccessH/p1->V - s->bN[2]/s->V;
	intercept[1] = (p2->bN[0] + p2->bN[1]*fD2)/p2->V - s->bN[0]/s->V;
	slope[1] = p2->bN[2]*p2->fAccessH/p2->V - s->bN[2]/s->V;
	intercept[2] = (p1->bN[0] + p1->bN[1]*fD1 + p2->bN[0] + p2->bN[1]*fD2)/(p1->V + p2->V) - s->bN[0]/s->V;
	slope[2] = (p1->bN[2]*p1->fAccessH + p2->bN[2]*p2->fAccessH)/(p1->V + p2->V) - s->bN[2]/s->V;
	
	fprintf(stdout,"\t\t\t<TR align=center><TH align=right colspan=4>Calculated match-point (<i>f</i><sub>D<sub>2</sub>O</sub>)<TD>%.3lf<TD>%.3lf<TD>%.3lf\n", -intercept[0]/slope[0], -intercept[1]/slope[1], -intercept[2]/slope[2]);	
	if(nI0 > 1)
		fprintf(stdout,"\t\t\t<TR align=center><TH align=right colspan=4>Experimental match-point (<i>f</i><sub>D<sub>2</sub>O</sub>)<TD>-<TD>-<TD>%.3lf(%.0lf)\n", -ab[1]/ab[0], 1000.0*sqrt(errXint));	
	fprintf(stdout,"\t\t</TABLE>\n");
	
	fprintf(stdout,"\t\t<P><BR><B>Other related quantities</B><BR>\n");
	/* print data relating to the density of the solvent and the complex */

	fprintf(stdout,"\t\t<P>The density (in g.cm<sup>-3</sup> @ 20&deg;C) of the solvent at any D<sub>2</sub>O fraction is:<BR>\n");
	fprintf(stdout,"\t\t<UL>\n");
	fprintf(stdout,"\t\t\t<LI>Density = %.3lf + %.3lf<i>f</i><sub>D<sub>2</sub>O</sub>\n", densityS, mS[0].n/1000.0*1.928);
	fprintf(stdout,"\t\t</UL>\n");
		/* 2.012 = mD2O - mH2O, however I am using an adjusted figure that incorporates a slight volume difference between H2O and D2O of 
		 * 1.928 to ensure the density of pure D2O comes out as 1.1050 @ 20 degrees C */

	/* mass information */

	fprintf(stdout,"\t\t<P>The molecular mass (kDa) at any D<sub>2</sub>O fraction is:<BR>\n");
	fprintf(stdout,"\t\t<UL>\n");
	fprintf(stdout,"\t\t\t<LI>Mass<sub>1</sub> = %.3lf + %.3lf<i>f</i><sub>D,1</sub> + %.3lf<i>f</i><sub>AccessH,1</sub><i>f</i><sub>D<sub>2</sub>O</sub>\n",p1->mass[0], p1->mass[1], p1->mass[2]);
	fprintf(stdout," = %.3lf + %.3lf<i>f</i><sub>D<sub>2</sub>O</sub>\n",p1->mass[0] + p1->mass[1]*fD1, p1->mass[2]*p1->fAccessH);
	fprintf(stdout,"\t\t\t<LI>Mass<sub>2</sub> = %.3lf + %.3lf<i>f</i><sub>D,2</sub> + %.3lf<i>f</i><sub>AccessH,2</sub><i>f</i><sub>D<sub>2</sub>O</sub>\n",p2->mass[0], p2->mass[1], p2->mass[2]);
	fprintf(stdout," = %.3lf + %.3lf<i>f</i><sub>D<sub>2</sub>O</sub>\n",p2->mass[0] + p2->mass[1]*fD2, p2->mass[2]*p2->fAccessH);
	fprintf(stdout,"\t\t</UL>\n");

//	fprintf(stdout,"\t</BODY>\n");
//	fprintf(stdout,"</HTML>\n");
	
	/* print a comment at the end of the .html file, that can be read in by other modules */

	fprintf(stdout,"\t\t\n<!-- %lf %lf %lf %lf %lf %lf -->\n", slope[0]*10.0, intercept[0]*10.0, p1->V, slope[1]*10.0, intercept[1]*10.0, p2->V);

	free(m1); free(m2); free(mS);
	free(p1); free(p2); free(s);

	if(nI0 > 2){
		
		free(X); free(Y); free(sigY); free(ab); free(corr_ab);

	}

	return(1);

}

