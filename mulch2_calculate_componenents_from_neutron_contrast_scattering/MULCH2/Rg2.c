/****************************************************************** 
 ***                             Rg2.c                          *** 
 ******************************************************************
 *
 * Additional Info:
 *
 * Date: 9/4/06 last modified 24/9/08
 * Author: A.E. Whitten
 * Description:
 * 
 * Program that analyses the dependence of Rg's from a NCV series 
 * to determine information about the complex
 *
 * Version 1.0: File created
 *
 * Version 1.1: Changes made to array allocation, commenting and
 *		typos corrected
 *
 * Oct 2016 : Create Rg2.c from Rg.c
 *
 * Oct 2016 : Changes include file name and parameter changes
 *		Compost was changed to allow more than 2 components
 *		and to calculate form factor scattering for more than 2 components.
 *		Rg was not changed to allow more than 2 components
 *		due to lack of algorithm for parallel axis theorum involving more than 2 components.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mathFunctions2.h"

#define RGT(p, columns, name, func_name){									\
														\
	if (( p = (Rg_t *)malloc((columns)*sizeof(Rg_t))) == NULL){						\
														\
		fprintf(stderr, "\nError: Can't allocate memory for *%s in function %s.\n", name, func_name);	\
		return(0);											\
														\
	}													\
														\
	p->n = columns;												\
														\
}													
/* Rg array memory allocation macro */

typedef struct{

	int n;

	/* n is the number of contrast points */

	double Rg, RgErr, Drho1, fV1, Drho2, fV2, Drho;

	/* Rg is the radius of gyration; RgErr is the ESD of Rg; Drho1 is the contrast of component 1;
	 * fV1 is the volume fraction of component1; Drho2 and fV2 are analogous quantities for component 2;
	 * Drho is the total contrast */

}Rg_t;   /* data structure for radius of gyration information */

int getData(Rg_t *data){   /* retrives data from the input file */

	/* data stores the information from this file */

	char buffer[256];

	/* buffer stores information as it is read line by line */

	int i;

	/* i is a counter */

	for(i = 0; i < data->n; i ++){

		fgets(buffer, 256, stdin);
		sscanf(buffer,"%lf%lf%lf%lf%lf", &data[i].Rg, &data[i].RgErr, &data[i].Drho1, &data[i].Drho2, &data[i].fV1);

		if(data[i].RgErr == 0.0) data[i].RgErr = 0.01;   /* prevent any divide by zero problems */
		data[i].fV2 = 1.0 - data[i].fV1;
		data[i].Drho = data[i].fV1*data[i].Drho1 + data[i].fV2*data[i].Drho2;

	}
	
	return(1);
	
}

int doParallelAxisAnalysis(Rg_t *data){   /* Perform the parallel axis analysis */
	
	/* data is the Rg data */
	
	int i, j, k, l;

	/* i, j, k and l are all counters */

	double **X, **Xi, *Y, *I, chi2 = 0.0, paaVector[3], w; 

	/* X is the Hessian; Xi is the inverse Hessian; Y is the LSQ vector; I is the result
	 * of the LSQ; chi2 is chi^2; paaVector contains products of contrasts and volumes, and is
	 * the coefficients in the parallel axis theorem; w is the weight */
	
	DOUBLE(Y, 3, "Y", "doParallelAxisAnalysis", i);
	DOUBLE(I, 3, "I", "doParallelAxisAnalysis", i);
	DOUBLE2D(X, 3, 3, "X", "doParallelAxisAnalysis", i, j);
	DOUBLE2D(Xi, 3, 3, "Xi", "doParallelAxisAnalysis", i, j);
			
	for(i = 0; i < data->n; i ++){
		
		paaVector[0] = data[i].Drho1*data[i].fV1/(data[i].Drho1*data[i].fV1+data[i].Drho2*data[i].fV2);
		paaVector[1] = 1.0 - paaVector[0];   /* calculate the coefficients */
		paaVector[2] = paaVector[0]*paaVector[1];
		
		w = 1.0/(2.0*data[i].Rg*data[i].RgErr);   /* convert 1/sig(Rg) to 1/sig(R^2) */ 
		
		for(k = 0; k < 3; k ++){
			
			Y[k] += w*w*paaVector[k]*data[i].Rg*data[i].Rg;
			
			for(l = 0; l < 3; l ++)
				X[k][l] += w*w*paaVector[k]*paaVector[l];
				
		}
	
	}

	invert3x3Matrix(X, Xi);   /* solve LSQ problem */

	multiplyMatrixAndVector(Xi, 3, 3, Y, I);
	
	for(i = 0; i < data->n; i ++){
		
		paaVector[0] = data[i].Drho1*data[i].fV1/(data[i].Drho1*data[i].fV1+data[i].Drho2*data[i].fV2);
		paaVector[1] = 1.0 - paaVector[0];
		paaVector[2] = paaVector[0]*paaVector[1];

		w = 1.0/(2.0*data[i].Rg*data[i].RgErr);
		
		chi2 += w*w*(data[i].Rg*data[i].Rg - paaVector[0]*I[0] - paaVector[1]*I[1] - paaVector[2]*I[2])*(data[i].Rg*data[i].Rg - paaVector[0]*I[0] - paaVector[1]*I[1] - paaVector[2]*I[2])/(double)(data->n-3);
		
	}

        fprintf(stdout,"\t\t<TABLE align=center width=\"30%%\" border=1>\n");   /* construct table */
        fprintf(stdout,"\t\t\t<CAPTION><EM>Parallel Axis Analysis</EM></CAPTION>\n");
        fprintf(stdout,"\t\t\t<TR>");
        fprintf(stdout,"<TH width=\"33%%\">Parameter");
        fprintf(stdout,"<TH width=\"33%%\">Value");
        fprintf(stdout,"<TH width=\"33%%\">ESD\n");

	fprintf(stdout,"\t\t\t<TR align=center><TH><i>&chi;</i><sup>2</sup><TD>%.2lf<TD>-\n", chi2);
	fprintf(stdout,"\t\t\t<TR align=center><TH><i>R</i><sub>1</sub><sup>2</sup><TD>%.2lf<TD>%.2lf\n", I[0], sqrt(chi2*Xi[0][0]));
	fprintf(stdout,"\t\t\t<TR align=center><TH><i>R</i><sub>2</sub><sup>2</sup><TD>%.2lf<TD>%.2lf\n", I[1], sqrt(chi2*Xi[1][1]));
	fprintf(stdout,"\t\t\t<TR align=center><TH><i>D</i><sup>2</sup><TD>%.2lf<TD>%.2lf\n", I[2], sqrt(chi2*Xi[2][2]));
	fprintf(stdout,"\t\t\t<TR align=center><TH><i>R</i><sub>1</sub><TD>%.2lf<TD>%.2lf\n", sqrt(I[0]), 0.5/sqrt(I[0])*sqrt(chi2*Xi[0][0]));
	fprintf(stdout,"\t\t\t<TR align=center><TH><i>R</i><sub>2</sub><TD>%.2lf<TD>%.2lf\n", sqrt(I[1]), 0.5/sqrt(I[1])*sqrt(chi2*Xi[1][1]));
	fprintf(stdout,"\t\t\t<TR align=center><TH><i>D</i><TD>%.2lf<TD>%.2lf\n", sqrt(I[2]), 0.5/sqrt(I[2])*sqrt(chi2*Xi[2][2]));
	
	fprintf(stdout,"\t\t</TABLE>\n");	

	free(X);
	free(Xi);
	free(Y);
	free(I);
	
	return(1);	

}

int doStuhrmannAnalysis(Rg_t *data){   /* Performs the Stuhrmann analysis */

	/* data stores the Rg data */
	
	double X[20], Y[20], sigY[20], *fit, **corr, chi2 = 0.0;

	/* X contains the x-coordinates (1/Drho); Y contains the y-coordinates (Rg^2)
	 * sigY contains the error in the y-coordinates (sig(Rg^2)); fit are the cooeficients
	 * from the polynomial fit; corr is the correlation matrix; chi2 is chi^2*/

	double *DR1R2D, **B, **Bi, *C;

	/* DR1R2D are the error estimates (accounting for parameter correlations of the derived parameters
	 * R1, R2 and D; B is the RHS of the equations 5a-c given by Olah 1994
	 * Bi is the inverse of B; C contains the solution (i.e. R1, R2 and D) */

	int i, j, k;

	/* i, j and k are counters */

	DOUBLE(fit, 3, "fit", "doStuhrmannAnalysis", i);
	DOUBLE2D(corr, 3, 3, "corr", "doStuhrmannAnalysis", i, j);

	for(i = 0; i < data->n; i ++){

		X[i] = 1.0/data[i].Drho;
		Y[i] = data[i].Rg*data[i].Rg;
		sigY[i] = 2.0*data[i].Rg*data[i].RgErr;	

	}

	chi2 = polynomialFit(2, X, Y, sigY, data->n, fit, corr, 1);

	DOUBLE(DR1R2D, 3, "DR1R2D", "doStuhrmannAnalysis", i);
	
	/* Now solve linear equations to determin RH^2, RD^2 and D^2 */
	
	DOUBLE2D(Bi, 3, 3, "Bi", "doStuhrmannAnalysis", i, j);
	DOUBLE2D(B, 3, 3, "B", "doStuhrmannAnalysis", i, j);
	DOUBLE(C, 3, "C", "doStuhrmannAnalysis", i);

	B[2][0] = data[0].fV1;
	B[2][1] = data[0].fV2;
	B[2][2] = data[0].fV1*data[0].fV2;
	B[1][0] = (data[0].Drho1 - data[0].Drho2)*data[0].fV1*data[0].fV2;
	B[1][1] = -(data[0].Drho1 - data[0].Drho2)*data[0].fV1*data[0].fV2;
	B[1][2] = (data[0].Drho1 - data[0].Drho2)*data[0].fV1*data[0].fV2*(data[0].fV2*data[0].fV2 - data[0].fV1*data[0].fV1);
	B[0][0] = 0.0;
	B[0][1] = 0.0;
	B[0][2] = -(data[0].Drho1 - data[0].Drho2)*(data[0].Drho1 - data[0].Drho2)*data[0].fV1*data[0].fV1*data[0].fV2*data[0].fV2;

	invert3x3Matrix(B, Bi);

	multiplyMatrixAndVector(Bi, 3, 3, fit, C);
	
	fprintf(stdout,"\n\t<TABLE align=center width=\"40%%\" border=1>\n");   /* construct table */
	fprintf(stdout,"\t\t\t<CAPTION><EM>Stuhrmann plot data</EM></CAPTION>\n");
	fprintf(stdout,"\t\t\t<TR>");
	fprintf(stdout,"<TH width=\"25%%\">&Delta;<i>&rho;</i><sup>-1</sup>");
	fprintf(stdout,"<TH width=\"25%%\"><i>R<sub>g</sub></i><sup>2</sup> (exp)");
	fprintf(stdout,"<TH width=\"25%%\"><i>&sigma;</i>(<i>R<sub>g</sub></i><sup>2</sup>)\n");
	fprintf(stdout,"<TH width=\"25%%\"><i>R<sub>g</sub></i><sup>2</sup> (calc)\n");
	
	for(i = 0; i < data->n; i ++)
		fprintf(stdout,"\t\t\t<TR align=center><TD>%.3lf<TD>%.1lf<TD>%.1lf<TD>%.1lf\n", X[i], Y[i], sigY[i], fit[2] + fit[1]*X[i] + fit[0]*X[i]*X[i]);
		
	fprintf(stdout,"\t\t</TABLE><BR><BR>\n");
	
	for(i = 0; i < 3; i ++)
		for(j = 0; j < 3; j ++)
			for(k = 0; k < 3; k ++)
				DR1R2D[i] += chi2*Bi[i][j]*corr[k][j]*Bi[i][k];
			
		
	fprintf(stdout,"\n\t<TABLE align=center width=\"30%%\" border=1>\n");
	fprintf(stdout,"\t\t\t<CAPTION><EM>Stuhrmann Analysis</EM></CAPTION>\n");
	fprintf(stdout,"\t\t\t<TR>");
	fprintf(stdout,"<TH width=\"33%%\">Parameter");
	fprintf(stdout,"<TH width=\"33%%\">Value");
	fprintf(stdout,"<TH width=\"33%%\">ESD\n");

	fprintf(stdout,"\t\t\t<TR align=center><TH><i>&chi;</i><sup>2</sup><TD>%.2lf<TD>-\n", chi2);
	fprintf(stdout,"\t\t\t<TR><TH colspan=3 align=left>Stuhrmann plot coefficients\n");
	fprintf(stdout,"\t\t\t<TR align=center><TH><i>R<sub>m</sub></i><sup>2</sup><TD>%.2lf<TD>%.2lf\n", fit[2], sqrt(chi2*corr[2][2]));
	fprintf(stdout,"\t\t\t<TR align=center><TH><i>&alpha;</i><sup>&dagger;</sup></i><TD>%.2lf<TD>%.2lf\n", fit[1], sqrt(chi2*corr[1][1]));
	fprintf(stdout,"\t\t\t<TR align=center><TH><i>&beta;</i><TD>%.2lf<TD>%.2lf\n", -fit[0], sqrt(chi2*corr[0][0]));
	fprintf(stdout,"\t\t\t<TR><TH colspan=3 align=left>Extracted parameters\n");
	fprintf(stdout,"\t\t\t<TR align=center><TH><i>R</i><sub>1</sub><sup>2</sup><TD>%.2lf<TD>%.2lf\n", C[0], sqrt(DR1R2D[0]));
	fprintf(stdout,"\t\t\t<TR align=center><TH><i>R</i><sub>2</sub><sup>2</sup><TD>%.2lf<TD>%.2lf\n", C[1], sqrt(DR1R2D[1]));
	fprintf(stdout,"\t\t\t<TR align=center><TH><i>D</i><sup>2</sup><TD>%.2lf<TD>%.2lf\n",  C[2], sqrt(DR1R2D[2]));
	fprintf(stdout,"\t\t\t<TR align=center><TH><i>R</i><sub>1</sub><TD>%.2lf<TD>%.2lf\n", sqrt(C[0]), 0.5/sqrt(C[0])*sqrt(DR1R2D[0]));
	fprintf(stdout,"\t\t\t<TR align=center><TH><i>R</i><sub>2</sub><TD>%.2lf<TD>%.2lf\n", sqrt(C[1]), 0.5/sqrt(C[1])*sqrt(DR1R2D[1]));
	fprintf(stdout,"\t\t\t<TR align=center><TH><i>D</i><TD>%.2lf<TD>%.2lf\n", sqrt(C[2]), 0.5/sqrt(C[2])*sqrt(DR1R2D[2]));
	fprintf(stdout,"\t\t\t<TR align=center><TH><i>R</i><sub>m</sub><TD>%.2lf<TD>%.2lf\n", sqrt(fit[2]), 0.5/sqrt(fit[2])*sqrt(chi2*corr[2][2]));
	
	fprintf(stdout,"\t\t</TABLE>\n");	
	
	if(fit[1] > 0) fprintf(stdout,"\t\t<P align=center><sup>&dagger;</sup>Because the sign of alpha is positive, the component with a higher <BR>scattering density lies toward the periphery of the complex.<BR><BR><BR>\n");
	if(fit[1] < 0) fprintf(stdout,"\t\t<P align=center><sup>&dagger;</sup>Because the sign of alpha is negative, the component with a lower <BR>scattering density lies toward the periphery of the complex.<BR><BR><BR>\n");

	free(DR1R2D);
	free(B);
	free(Bi);
	free(C);
	free(fit);
	free(corr);

	return(1);
		
}

int main(void){
	
	char buffer[256];

	/* buffer stores data read from stdin */

	int n;

	/* n is the number of contrast points */

	Rg_t *data;

	/* data is the radius of gyration data */

	fgets(buffer, 256, stdin);

//	fprintf(stdout,"<!DOCTYPE html PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">\n<HTML>\n");
//	fprintf(stdout,"\t<HEAD>\n\t\t<TITLE>\n\t\t\tProject: %s\t\t</TITLE>\n", buffer);
//	fprintf(stdout,"\t\t<META HTTP-EQUIV=\"Content-Type\" CONTENT=\"text/html; charset=iso-8859-1\">\n");
//	fprintf(stdout,"\t</HEAD>\n");
//	fprintf(stdout,"\t<BODY>\n");
	fprintf(stdout,"\n\n<P>Project: %s<BR><BR>", buffer);
	
	fgets(buffer, 256, stdin);
	sscanf(buffer,"%d", &n);
	
	RGT(data, n, "data", "main");
	data->n = n;
	getData(data);

	doStuhrmannAnalysis(data);
	doParallelAxisAnalysis(data);
	
//	fprintf(stdout,"\t</BODY>\n");
//	fprintf(stdout,"</HTML>\n");

	free(data);

	return(1);
		
}
