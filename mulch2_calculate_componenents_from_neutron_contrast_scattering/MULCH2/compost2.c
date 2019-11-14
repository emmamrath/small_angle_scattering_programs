/****************************************************************** 
 ***                           compost2.c                       *** 
 ******************************************************************
 *
 * Additional Info:
 *
 * Date: 9/4/06 last modified 23/9/08
 * Author: A.E. Whitten
 * Description:
 * 
 * Program that decomposes a NCV series into three composite profiles
 *
 * Version 1.0: File created
 *
 * Version 1.1: Changes made to array allocation, commenting and
 * 		typos corrected
 *
 * Oct 2016 : Create compost2.c from compost.c
 *            Coding changes by E.M.Rath.
 *	      The Cholesky and invertNxNMatrix in mathFunctions2.c written by G. Doherty.
 *
 * Oct 2016 : Changes to compost.c to create compost2.c include :
 *            Input file has an extra p1 parameter
 *            so that p1 for regression can be different to the p1 for Rg analysis.
 *            Call multiple regression and invert matrix with changed parameters
 *            to allow more than 2 components in the scattering system.
 *            Add input parameters for : calc_structure_factors, divide_intensity_by_error, SVC_threshold.
 *            Compost was changed to allow more than 2 components
 *            and to calculate form factor scattering for more than 2 components.
 *            Rg was not changed to allow more than 2 components
 *            due to lack of algorithm for parallel axis theorum involving more than 2 components.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mathFunctions2.h"

#define MAX_FORM_FACTORS 10
#define MAX_FACTORS MAX_FORM_FACTORS * (MAX_FORM_FACTORS + 1) / 2
	/* eg. for 2 form factors 11 & 22, the factors are 11, 12 & 22, and the structure factor is 12 */
	/* eg. for 3 form factors 11, 22, 33, the factors are 11, 12, 13, 22, 23, 33, and the structure factors are 12, 13, 23 */
	/* eg. for 4 form factors 11, 22, 33, 44, the factors are 11, 12, 13, 14, 22, 23, 24, 33, 34, 44, and the structure factors are 12, 13, 14, 23, 24, 34 */
#define MAX_CONTRAST_POINTS 100
#define MAX_Q_POINTS 1000

#define INTENSITY(p, columns, name, func_name){									\
														\
	if (( p = (intensity *)malloc((columns)*sizeof(intensity))) == NULL){					\
														\
		fprintf(stderr, "\nError: Can't allocate memory for *%s in function %s.\n", name, func_name);	\
		return(0);											\
														\
	}													\
														\
	p->num_contrast_points = columns;											\
														\
}													
/* intensity array memory allocation macro */

#define OPEN_FILE(fpr, name, mode){								\
												\
	if((fpr = fopen(name, mode)) == NULL){							\
												\
		fprintf(stderr,"\nError: Can't open file %s in %s mode!!!\n", name, mode);	\
		return(0);									\
												\
	}											\
												\
}												\
/* file opening macro */

int num_form_factors, num_factors, calc_structure_factors, divide_intensity_by_error, order_of_factors_1[MAX_FORM_FACTORS], order_of_factors_2[MAX_FORM_FACTORS];
double SVC_threshold;
int max_num_read_Q_points_in_data_files = 0;
int min_start_Q_point = MAX_Q_POINTS;

	/* num_form_factors are the number of different scattering components, 
	   and can include solvent when buffer subtraction has not been done. 
	   num_factors = num_form_factors + number of structure factors.
	   number of structure factors : there is one structure factor for each pair of form factors,
	and represents the scattering caused by electron pair where one electron is in one component and the other electron in the other component */



typedef struct{

        int num_contrast_points, n, p1, Rgp1, num_form_factors;
	
	/* num_contrast_points is the number of contrast points; n is the number of q-values in each file; p1 is the first point to be included in the fitting; */
	/* Rgp1 is the first point to be included in the Guinier analysis */

        double q[MAX_Q_POINTS], I[MAX_Q_POINTS], IErr[MAX_Q_POINTS], SLD[MAX_FORM_FACTORS], Volume[MAX_FORM_FACTORS], Drho;

	/* q is the momentum transfer (scattering angle); I is the intensity; Ierr is the ESD in I; 
	 * SLD is the scattering length density of each of the component; 
	 * Volume is the volume of each component; 
	 * Drho is the contrast of the complex. */

	double I0, I0Err, Rg, RgErr, chi2, m;

	/* I0 is the extropolated zero angle scattering; I0Err is the ESD of I0; Rg is the radius of gyration; RgErr is the ESD of Rg;
	 * chi2 is the quality of the Guinier fit; m is the scale factor that used to bring all datasets onto the correct scale */

	int from_data_point, to_data_point;

	/* from_data_point and to_data_point are the q points (from 1 to num_q_points) used to derive Rg and I(0) using Guinier region analysis */
	
	char correct;

	/* correct is a flag indicating whether a scale factor correction should be applied to each data set */

	int qRg;

	/* qRg limits the number of points included in the Guinier fit */

}intensity;   /* intensity data structure */



int getData(char *filename, intensity *data){   /* gets the data from the input file */
	
	/* filename is the filename of the file from which data is read; data stores the intensity data */

	FILE *in;

	/* in is a pointer to the scattering data */

	char buffer[256], buffer2[256];
	
	/* buffer store data line by line as it is read from the input file */

	int i, status = 0;

	/* i is a counter; status is a flag that switches to 1 when header information has been skipped over */
	
	OPEN_FILE(in, filename, "r");
	//fprintf(stdout,"filename = %s<BR>\n", filename );

	int num_Q_points_in_data_file = 0;

	while(status == 0){   /* scan over the title */

		fgets(buffer, 256, in);
		sscanf(buffer,"%lf%lf%lf", &data->q[0], &data->I[0], &data->IErr[0]);
		if ((data->q[0] != 0.0) && (data->I[0] != 0.0) && (data->IErr[0] != 0.0)) {
			status = 1;
			num_Q_points_in_data_file = 1;
		}
		//fprintf(stdout,"data->q[0] = %lf, data->I[0] = %lf, data->IErr[0] = %lf<BR>\n", data->q[0], data->I[0], data->IErr[0] );

	}

	//if(data->p1 - 1 > 0) data->IErr[0] = 0.0; /* set error = 0.0 so that it is not used in further fitting */ 
	// Don't do this preceding statement. Error of 0.0 will be treated as an excellent accurate point.
	// Points not to be used will be flagged a different way.

	status = 0;

	for(i = 1; fgets(buffer,256,in) != (char *)NULL && status == 0; i ++){

		sscanf(buffer,"%lf%lf%lf", &data->q[i], &data->I[i], &data->IErr[i]);
		//fprintf(stdout,"data->q[i] = %lf, data->I[i] = %lf, data->IErr[i] = %lf<BR>\n", data->q[i], data->I[i], data->IErr[i] );

		if ((data->q[i] != 0.0) && (data->I[i] != 0.0) && (data->IErr[i] != 0.0)) {
			num_Q_points_in_data_file = num_Q_points_in_data_file + 1;
		} else {
			status = 1; /* must be comments at the end of the file - data is finished */
		}
		
		//if(data->p1 - i - 1 > 0) data->IErr[i] = 0.0; /* set error = 0.0 so that it is not used in further fitting */
		// Don't do this preceding statement. Error of 0.0 will be treated as an excellent accurate point.
		// Points not to be used will be flagged a different way.
	}

	data->n = num_Q_points_in_data_file;   /* set the number of q-values in this file */

	if (num_Q_points_in_data_file > max_num_read_Q_points_in_data_files) {
		max_num_read_Q_points_in_data_files = num_Q_points_in_data_file;
	}

	fclose(in);

	return(1);
	
}



double getChi2(intensity *data, intensity *model, int parameter){   /* calculates chi^2 */

	/* data stores the intensity data; model stores the model data; parameter refers to the q-value being fitted */

	int i, ix1, ix2, ix3, num_parameters_in_fit = 0;

	/* i is a counter; num_parameters_in_fit is the number of parameters included in the fit */

	double weight, D, subtract_from_D, chi2 = 0.0, DrhoXY[MAX_FACTORS], IExp, IExpErr;

	/* weight is the weight; D is the difference between the data and the model; chi2 is  a measure of the quality of the fit;
	 * DrhoXY is DrhoX x DrhoY; Iexp is the experimental intensity; IExpErr is the ESD of IExp */

	fprintf(stdout,"\n\t\t<TR align=center><TH><nobr>&nbsp;&nbsp;&nbsp;%d (Q = %lf)&nbsp;&nbsp;&nbsp;</nobr>", parameter + 1, data->q[parameter]);;
	
	for( i = 0; i < data->num_contrast_points; i++ ) {

		for( ix3 = 0; ix3 < num_factors; ix3++ ) {

			ix1 = order_of_factors_1[ix3];
			ix2 = order_of_factors_2[ix3];

			DrhoXY[ix3] = data[i].SLD[ix1] * data[i].SLD[ix2];
		}
		/* eg. for 2 form_factors :
			DrhoXY[0] = data[i].SLD[0] * data[i].SLD[0];
			DrhoXY[1] = data[i].SLD[1] * data[i].SLD[1];
			DrhoXY[2] = data[i].SLD[0] * data[i].SLD[1];
		*/
		/* eg. for 4 form_factors :
			DrhoXY[0] = data[i].SLD[0] * data[i].SLD[0];
			DrhoXY[1] = data[i].SLD[1] * data[i].SLD[1];
			DrhoXY[2] = data[i].SLD[2] * data[i].SLD[2];
			DrhoXY[3] = data[i].SLD[3] * data[i].SLD[3];
			DrhoXY[4] = data[i].SLD[0] * data[i].SLD[1];
			DrhoXY[5] = data[i].SLD[0] * data[i].SLD[2];
			DrhoXY[6] = data[i].SLD[0] * data[i].SLD[3];
			DrhoXY[7] = data[i].SLD[1] * data[i].SLD[2];
			DrhoXY[8] = data[i].SLD[1] * data[i].SLD[3];
			DrhoXY[9] = data[i].SLD[2] * data[i].SLD[3];
		*/

		IExp = data[i].m*(data[i].I[parameter]);   /* rescale as appropriate */
		IExpErr = data[i].m*data[i].IErr[parameter];

		weight = 1.0;
		if(IExpErr > 0.0){
			weight = 1.0/IExpErr;
			if (divide_intensity_by_error == 0) {
				weight = 1.0;
			}
		}

		if (IExpErr > 0.0) {

			num_parameters_in_fit ++;   /* increase the number of points included in the fit */

			/* D = IExp - DrhoXY[0]*model[0].I[parameter] - DrhoXY[1]*model[1].I[parameter] - DrhoXY[2]*model[2].I[parameter]; */
			subtract_from_D = 0;
			for( ix1 = 0; ix1 < num_factors; ix1++ ) {
				subtract_from_D += DrhoXY[ix1]*model[ix1].I[parameter];
			}
			D = IExp - subtract_from_D;
		
			if(weight*D*weight*D < 9.0) fprintf(stdout,"<TD>%.2lf", weight*D*weight*D);   /* print the weighted mean square difference (MSD) */
			else fprintf(stdout,"<TD>%.2lf*", weight*D*weight*D);   /* if the MSD exceeds 9, then put an asterisk on it */
			chi2 += weight*weight*D*D;

		} else {
			fprintf(stdout,"<TD>-\n");
		}
	}

	if(num_parameters_in_fit <= 3) chi2 = 0.0;   /* I don't know how exactly to handle this situation at present, this is the best I can think of at the moment */

	else chi2 = chi2 / ((double)num_parameters_in_fit - 3.0);   /* calculate chi^2 */
	
	fprintf(stdout, "<TD>%.2lf",chi2);
	
	return(chi2);

}



int doMultipleLinearRegression(intensity *data, intensity *model, int parameter) {   /* does the LSQ fit */
	
	/* data stores the intensity data; model stores the model data; parameter refers to the q-value being fitted */
	/* This doMultipleLinearRegression is called for each q data point */

	int i, j, k, l, n_num_contrast_points, ix1, ix2, ix3, ix_n;

	/* i, j, k and l are counters; n is equivalent to data->num_contrast_points */

	double **X, **Xi, *Y, *I, weight, chi2, DrhoXY[MAX_FACTORS], IExp, IExpErr; 

	/* X is the hessian; Xi is the inverse hessian; Y is the corresponding vector; I are the calculated composite intensities;
	 * weight is the weight; chi2 is chi^2; DrhoXY is DrhoX x DrhoY; Iexp is the experimental intensity; IExpErr is the ESD of IExp */

	DOUBLE2D(X, MAX_FACTORS, MAX_FACTORS, "X", "doMultipleLinearRegression", i, j);   /* allocate memory for the arrays */
	DOUBLE2D(Xi, MAX_FACTORS, MAX_FACTORS, "Xi", "doMultipleLinearRegression", i, j);
	DOUBLE(Y, MAX_FACTORS, "Y", "doMultipleLinearRegression", i);
	DOUBLE(I, MAX_FACTORS, "I", "doMultipleLinearRegression", i);

	n_num_contrast_points = data->num_contrast_points;

	for(k = 0; k < num_factors; k ++) {
		Y[k] = 0;
		for(l = 0; l < num_factors; l ++) {
			X[k][l] = 0;
		}
	}

	ix_n = 0;
	for (i = 0; i < n_num_contrast_points; i ++) { /* for each contrast point (within this Q data point)... */

		if ((parameter + 1) >= data[i].p1) { /* include this Q data point for this contrast file only if we are up to or past the first Q point to use in this particular contrast data file */

			for( ix3 = 0; ix3 < num_factors; ix3++ ) {

				ix1 = order_of_factors_1[ix3];
				ix2 = order_of_factors_2[ix3];

				DrhoXY[ix3] = data[i].SLD[ix1] * data[i].SLD[ix2]; /* ...work out the delta-contrast value */
			}
			/* eg. for 2 form factors :
				DrhoXY[0] = data[i].SLD[0] * data[i].SLD[0];
				DrhoXY[1] = data[i].SLD[1] * data[i].SLD[1];
				DrhoXY[2] = data[i].SLD[0] * data[i].SLD[1];
			*/

			/* rescale as appropriate the experimental intensity and experimental error for this Q data point for this contrast point */
			IExp = data[i].m * data[i].I[parameter];
			IExpErr = data[i].m * data[i].IErr[parameter];

			weight = 1;
			if (IExpErr > 0.0) {
				weight = 1.0/IExpErr;
				if (divide_intensity_by_error == 0) {
					weight = 1.0;
				}
			}

			if (IExpErr > 0.0) {

				for (k = 0; k < num_factors; k ++) {
		
					Y[k] += weight*weight*DrhoXY[k]*IExp;
		
					for (l = 0; l < num_factors; l ++) {
						X[k][l] += weight*weight*DrhoXY[k]*DrhoXY[l];
					}
				}
			}

			ix_n = ix_n + 1;
		}
	}

	int matrix_n = num_factors;
	int num_columns_A = num_factors;

	// ix_n is the number of contrast points minus any contrast points for which we are not yet up to the first Q point to use

	if (ix_n > 0) {

		invertNxNMatrix(matrix_n, X, Xi, num_columns_A);

		multiplyMatrixAndVector(Xi, num_factors, num_factors, Y, I);

		for(i = 0; i < num_factors; i ++) {
			model[i].I[parameter] = I[i];  /* set the model values */ 
		}
	
		chi2 = getChi2(data, model, parameter);
		
		for(i = 0; i < num_factors; i ++) {
			model[i].IErr[parameter] = sqrt(chi2*Xi[i][i]);
		}
	}

	free(Y);
	free(I);
	free(X);
	free(Xi);	

	return(1);	

}



int doSingularValueDecomposition( intensity *data, intensity *model, int parameter, double *SVD_values, double *SVD_y_div_sigma ) {

	/*

	Solve simultaneous equations of intensity to get the individual form factors and structure factors.

	For example, for 4 contrast points and 2 form factors (11 and 22) and 1 structure factor (12) :

	IExp1 = prot1.prot1.I11 + poil1.poil1.I22 + prot1.poil1.I12
	IExp2 = prot2.prot2.I11 + poil2.poil2.I22 + prot2.poil2.I12
	IExp3 = prot3.prot3.I11 + poil3.poil3.I22 + prot3.poil3.I12
	IExp4 = prot4.prot4.I11 + poil4.poil4.I22 + prot4.poil4.I12

	where :	IExp1,IExp2,IExp3,IExp4 are the experimental intensities for contrast points 1,2,3,4, of the entire molecule (having 2 subcomponents),
		prot1,prot2,prot3,prot4 are the contrasts (scattering length density (SLD) x volume) of the first subcomponent (having form factor 11),
		poil1,poil2,poil3,poil4 are the contrasts (scattering length density (SLD) x volume) of the second subcomponent (having form factor 22),
		I11 is the form factor for the first subcomponent (scattering from 2 electrons within this one subcomponent),
		I22 is the form factor for the second subcomponent (scattering from 2 electrons within this one subcomponent),
		I12 is the structure factor for the scattering from the first and second subcomponents (an electron from each).

	This gives the matrix equation to solve :

	[  prot1.prot1  poil1.poil1  prot1.poil1  ]  .  [  I11  ]  =  [  IExp1  ]
	[  prot2.prot2  poil2.poil2  prot2.poil2  ]     [  I22  ]     [  IExp2  ]
	[  prot3.prot3  poil3.poil3  prot3.poil3  ]     [  112  ]     [  IExp3  ]
	[  prot4.prot4  poil4.poil4  prot4.poil4  ]                   [  IExp4  ]

	                                       A     .      x      =      b

	To solve [A] . [x] = [b] 
	do Singular Value Decomposition
	from "Numerical Recipes in C++" Press, Teukolsky, Vetterling, Flannery

	[x] = [V] . [diag(1/wj)] . [U^T] . [b]

	1x3   3x3  4x4(3x3filled)   4x4    4x1		???

	*/

	/* data stores the intensity data; model stores the model data; parameter refers to the q-value being fitted */
	/* This doMultipleLinearRegression is called for each q data point */

	int i, j, k, l, n_num_contrast_points, ix1, ix2, ix3, ix_n;

	/* i, j, k and l are counters; n is equivalent to data->num_contrast_points */

	double weight, chi2, IExp, IExpErr; 

	/* weight is the weight; chi2 is chi^2; DrhoXY is DrhoX x DrhoY; Iexp is the experimental intensity; IExpErr is the ESD of IExp */

	double **A_array, *w_array, **v_array, *b_array, *x_array, **v_transpose, *vector_for_IErr, *vector_for_IErr_2;

	DOUBLE2D( A_array, MAX_CONTRAST_POINTS, MAX_FACTORS, "A_array", "doSingularValueDecomposition", i, j );   /* allocate memory for the arrays */
	DOUBLE2D( v_array, MAX_FACTORS, MAX_FACTORS, "v_array", "doSingularValueDecomposition", i, j );
	DOUBLE( w_array, MAX_CONTRAST_POINTS, "w_array", "doSingularValueDecomposition", i );
	DOUBLE( b_array, MAX_CONTRAST_POINTS, "b_array", "doSingularValueDecomposition", i );
	DOUBLE( x_array, MAX_FACTORS, "x_array", "doSingularValueDecomposition", i );
	DOUBLE2D( v_transpose, MAX_FACTORS, MAX_FACTORS, "v_transpose", "doSingularValueDecomposition", i, j );
	DOUBLE( vector_for_IErr, MAX_FACTORS, "vector_for_IErr", "doSingularValueDecomposition", i );
	DOUBLE( vector_for_IErr_2, MAX_FACTORS, "vector_for_IErr_2", "doSingularValueDecomposition", i );
	//fprintf(stdout,"MAX_CONTRAST_POINTS %i, n_num_contrast_points %i <br>\n", MAX_CONTRAST_POINTS, n_num_contrast_points );
	for ( i=0; i<MAX_CONTRAST_POINTS; i++ ) {
		w_array[i] = 0.0;
		b_array[i] = 0.0;
		for ( j=0; j<MAX_FACTORS; j++ ) {
			A_array[i][j] = 0.0;
		}
	}
	for ( i=0; i<MAX_FACTORS; i++ ) {
		for ( j=0; j<MAX_FACTORS; j++ ) {
			v_array[i][j] = 0.0;
		}
	}
	for ( j=0; j<MAX_FACTORS; j++ ) {
		x_array[j] = 0.0;
	}

	n_num_contrast_points = data->num_contrast_points;

	ix_n = 0;
	for (i = 0; i < n_num_contrast_points; i ++) { /* for each contrast point (within this Q data point)... */

		if ((parameter + 1) >= data[i].p1) { /* include this Q data point for this contrast file only if we are up to or past the first Q point to use in this particular contrast data file */

			/* rescale as appropriate the experimental intensity and experimental error for this Q data point for this contrast point */
			IExp = data[i].m * data[i].I[parameter];
			IExpErr = data[i].m * data[i].IErr[parameter];

			weight = 1.0;
			if (IExpErr > 0.0) {
				weight = 1.0/IExpErr;
				if (divide_intensity_by_error == 0) {
					weight = 1.0;
					fprintf(stdout,"HERE");
				}
			}

			for( ix3 = 0; ix3 < num_factors; ix3++ ) {

				ix1 = order_of_factors_1[ix3];
				ix2 = order_of_factors_2[ix3];

				A_array[ix_n][ix3] = (data[i].SLD[ix1] * data[i].SLD[ix2]) * weight; // ...work out the delta-contrast value

			}

			b_array[ix_n] = IExp * weight;

			ix_n = ix_n + 1;
		}
	}

	// ix_n is the number of contrast points minus any contrast points for which we are not yet up to the first Q point to use

	/* Given a matrix a[1..m][1..n], this routine computes its singular value decomposition, A = U ·W ·V^T. 
	The matrix U replaces a on output. The diagonal matrix of singular values W is output as a vector w[1..n]. 
	The matrix V (not the transpose V^T ) is output as v[1..n][1..n]. */

	if (ix_n > 0) {

		//One can test the SVD routine separately with these three matrices:
		//1 1 1          1 1 1        1 1 1
		//2 2 2          1 2 3        1 2 4
		//3 3 3          2 3 4        1 3 9
		//4 4 4          5 6 7        8 9 10
		//You should get 1 non-zero singular value from the first matrix, 2 from the second, and 3 from the third.

		svdcmp( A_array, ix_n, num_factors, w_array, v_array ); //SVD the non-square matrix A_array.

		for ( j=0; j<num_factors; j++ ) {
			SVD_values[j] = w_array[j];
		}

		svbksb( A_array, w_array, v_array, ix_n, num_factors, b_array, x_array, SVC_threshold, SVD_y_div_sigma ); //Now we can backsubstitute.

		for(i = 0; i < num_factors; i ++) {
			model[i].I[parameter] = x_array[i];  /* set the model values */ 
		}
	
		chi2 = getChi2(data, model, parameter);	

		// multiply arrays : V x 1/SVD^2 x V^-T
		// model[i].IErr[parameter] = chi2 * 1/SVD^2[i][i]

		for ( i=0; i<num_factors; i++ ) {
			for ( j=0; j<num_factors; j++ ) {
				v_transpose[i][j] = v_array[j][i];
			}
		}
		multiplyMatrixAndVector(v_array, num_factors, num_factors, w_array, vector_for_IErr);
		multiplyMatrixAndVector(v_transpose, num_factors, num_factors, vector_for_IErr, vector_for_IErr_2);

		for(i = 0; i < num_factors; i ++) {
			//model[i].IErr[parameter] = sqrt(chi2*Xi[i][i]);
			//model[i].IErr[parameter] = 0.0001;
			model[i].IErr[parameter] = chi2 * vector_for_IErr_2[i];
		}

	}

	free(A_array);
	free(v_array);
	free(w_array);
	free(b_array);
	free(x_array);	

	return(1);	
}



int main(int argc, char *argv[]){
	
	/* argc is the command line argument counter; argv contains the command line arguments,
	 * argv[0] is the executable
	 * argv[1] is the directory in which the file is to be created
	 * argv[2] is a timestamp or something to distinguish one file from another
	 * argv[3] is a flag - if is is set to 'n' or 'N', then the program will not calculate structure factors and will only calculate the form factors.
	 *			and the structure factor scattering will be in the form factor scattering results.
				Do this when the scattering results are so noisy that teasing out the structure factors makes the results worse.
	 * argv[4] is a flag - if is is set to 'n' or 'N', then the program will not divide each intensity value by its error value when doing multiple regression.
	 * argv[5] is the threshold for the SVD Single Value Decomposition.
				A higher value gives smoother data that does not have as good a fit.
				A smaller value gives noiser data that has a better fit.
				The default is 1.0e-6 */

	/* For 2 components, the form factors are i11, i22, and the structure factors are i12.
	 * For 3 components, the form factors are i11, i22, i33, and the structure factors are i12, i13, i23.
	 * For 4 components, the form factors are i11, i22, i33, i44, and the structure factors are i12, i13, i14, i23, i24, i34. */

	FILE *output_file_ptr, *output_SVD_file_ptr, *output_YdivSigma_file_ptr;

	/* i11, i12 and i22 are file pointers to the 3 column output files */

	char filename[256], buffer[256], buffer2[256], output_filename[256], output_SVD_filename[256], output_YdivSigma_filename[256], command[256], SVD_threshold_string[256];

	/* filename is the filename of the file containing the experimental data; buffer is used to store 
	 * information as it is read; i11_, i12_ and i22_fname are the filenames where the 3 column output is written */

	intensity *data, *model;

	/* data is where the experimental intensity data is stored; model is where the model intensity data is stored */

	double q2[MAX_Q_POINTS], lnI[MAX_Q_POINTS], errlnI[MAX_Q_POINTS], *fit, **corr, chi2, qRgLimit;

	/* q2 is q^2; lnI is the natural log of I; errlnI is the ESD of the natural log of I; fit contains slope and intercept
	 * for the Guinier fit, corr is the correlation matrix for the fit; chi2 is chi^2 and qRgLimit is maximum number of points
	 * to use in the Guinier fit */
	
	int num_contrast_points_input_param, i, j, ix1, ix2, ix3, status;

	/* num_contrast_points_input_param is the number of contrast points; i, j, ix1 and ix2 are counters;
	 * status is a flag */

	calc_structure_factors = 1; /* default is to calculate form factors and structure factors. */
	if (argc > 3) {
		char struc_flag;
		struc_flag = *argv[3];
		if ((struc_flag == 'n') || (struc_flag == 'N')) {
			calc_structure_factors = 0; /* don't calculate the structure factors, only calculate the form factors */
		}
		//fprintf(stdout,"\n<BR>struc_flag=%c<BR><BR>\n", struc_flag );
	}

	divide_intensity_by_error = 1; /* default is to divide each intensity value by its error value when doing multiple regression */
	if (argc > 4) {
		char div_flag;
		div_flag = *argv[4];
		if ((div_flag == 'n') || (div_flag == 'N')) {
			divide_intensity_by_error = 0; /* don't divide each intensity value by its error value when doing multiple regression */
		}
		//fprintf(stdout,"\n<BR>div_flag=%c<BR><BR>\n", div_flag );
	}

	SVC_threshold = 1.0e-6;
	if (argc > 5) {
		strcpy(SVD_threshold_string, argv[5]);
		SVC_threshold = atof( SVD_threshold_string );
	}

	fgets(buffer, 256, stdin);

	fprintf(stdout,"<html>\n<head>\n<style type=\"text/css\">\n<!--\n");
	fprintf(stdout,".smallblueitalic {\n\tfont-size: 10px;\n\tfont-style: italic;\n\tfont-family: Arial, Helvetica, Verdana, sans-serif;\n\tcolor: #4343e3;\n}\n");
	fprintf(stdout,"-->\n</style>\n</head>\n<body>\n");

	fprintf(stdout,"\n\n<P>Project: %s<BR><BR>", buffer);

	fgets(buffer, 256, stdin);
	sscanf(buffer,"%d", &num_contrast_points_input_param);
	
	fgets(buffer, 256, stdin);
	sscanf(buffer,"%lf", &qRgLimit);

	fgets(buffer, 256, stdin);
	sscanf(buffer,"%d", &num_form_factors);

	num_factors = 0;
	for ( ix1 = 0; ix1 < num_form_factors; ix1++ ) {
		for ( ix2 = ix1; ix2 < num_form_factors; ix2++ ) {
			int process_this_factor = 1;
			if (calc_structure_factors == 0) {
				if (ix1 != ix2) {
					process_this_factor = 0;
				}
			}
			if (process_this_factor == 1) {
				num_factors++;
			}
		}
	}	

	// The order of filling and outputing the factors is to do the form factors first, then the structure/interference factors.
	// Ie, 11,22,33,12,13,23 and not 11,12,13,22,23,33
	// Thus, if singular value decomposition does not find enough signals, it will fill the molecue form factor signals first.

	ix3 = 0;
	for( ix1 = 0; ix1 < num_form_factors; ix1++ ) {
		for( ix2 = ix1; ix2 < num_form_factors; ix2++ ) {

			int process_this_factor = 0;
			if (ix1 == ix2) {
				process_this_factor = 1;
			}
			if (process_this_factor == 1) {
				order_of_factors_1[ix3] = ix1;
				order_of_factors_2[ix3] = ix2;
				ix3++;
			}
		}
	}
	if (calc_structure_factors == 1) {
		for( ix1 = 0; ix1 < num_form_factors; ix1++ ) {
			for( ix2 = ix1; ix2 < num_form_factors; ix2++ ) {

				int process_this_factor = 0;
				if (ix1 != ix2) {
					process_this_factor = 1;
				}
				if (process_this_factor == 1) {
					order_of_factors_1[ix3] = ix1;
					order_of_factors_2[ix3] = ix2;
					ix3++;
				}
			}
		}
	}

	INTENSITY(data, num_contrast_points_input_param, "data", "main");
	INTENSITY(model, num_factors, "model", "main");
	
	for (i = 0; i < data->num_contrast_points; i ++) { /* for each contrast point... */
		
		fgets(buffer, 256, stdin);

		/* read in the following : correct, scale, 1st-point, 1st-point-for-Rg-calculations, filename */

		sscanf( buffer, "%c%lf%d%d%s", &data[i].correct, &data[i].m, &data[i].p1, &data[i].Rgp1, filename );

		if (data[i].p1 < min_start_Q_point) {
			min_start_Q_point = data[i].p1;
		}

		/* for each form factor, read in the following : scattering-length-density-SLD, Volume */

		data[i].Drho = 0;

		for (ix1 = 0; ix1 < num_form_factors; ix1++) {
			fgets(buffer, 256, stdin);
			sscanf( buffer, "%lf%lf", &data[i].SLD[ix1], &data[i].Volume[ix1] );
			data[i].Drho += data[i].SLD[ix1] * data[i].Volume[ix1];
		}

		getData(filename, &data[i]); /* ...read in the Q data points */
	
	}

	for(i = 0; i < data->n; i ++) {   /* Determine whether the binning for all data sets is the same */
		for(j = 1; j < data->num_contrast_points; j ++) {

			//fprintf(stdout,"i = %i, j = %i, data[0].q[i] = %0.3f, data[j].q[i] = %0.3f<br>", i, j, data[0].q[i], data[j].q[i]);

			if (data[0].q[i] != data[j].q[i]) {

				//fprintf(stdout,"<P>The data binning for point %d in dataset %d differs from the reference<br>\n", i + 1, j + 1);
				fprintf(stdout,"<P>The data binning for point %d (%lf) in dataset %d differs from the reference (%lf)<br>\n", i + 1, data[0].q[i], j + 1, data[j].q[i]);
				return(0);

			}
		}
	}

	DOUBLE(fit, num_form_factors, "fit", "main", i);
	DOUBLE2D(corr, num_form_factors, num_form_factors, "corr", "main", i, j);

	/*Guinier analysis to determine scale factors*/	

	for (j = 0; j < data->num_contrast_points; j ++) {

		//fprintf(stdout,"j = %i<br>\n", j );
			
		fit[0] = 0.0;
		fit[1] = 0.0;

		for(i = 0; i < data->n; i ++){

			//fprintf(stdout,"calc q2, lnI, errlnI, i = %i<br>\n", i );

			q2[i] = data[j].q[i]*data[j].q[i];
			lnI[i] = log(data[j].I[i]);
			errlnI[i] = data[j].IErr[i] / (data[j].I[i]);

			//fprintf(stdout,"calc q2, lnI, errlnI, i = %i, q2[i] = %lf, lnI[i] = %lf, errlnI[i] = %lf<br>\n", i, q2[i], lnI[i], errlnI[i] );
		}

		status = 0;

		//fprintf(stdout,"j = %i, data[j].Rgp1 = %i, 5 + data[j].Rgp1 = %i<br>\n", j, data[j].Rgp1, 5 + data[j].Rgp1 );

		data[j].from_data_point = 0;
		data[j].to_data_point = 0;

		for (i = 5 + data[j].Rgp1; ((status == 0) && (i < data->n)); i ++) {  /* assumes that we need at least 5 points in the guinier region */

			if (i < data->n) {

				if (data[j].from_data_point == 0)
					data[j].from_data_point = i + 1;

				chi2 = polynomialFit2(1, q2, lnI, errlnI, i, fit, corr, data[j].Rgp1);

				//fprintf(stdout,"calc chi2, i = %i, fit[0] slope = %0.3f, fit[1] intercept = %0.3f<br>\n", i, fit[0], fit[1] );

				//fprintf(stdout,"i = %i, fit[0] = %lf, data[j].q[i - 1] = %lf, sqrt(-3.0*fit[0])*data[j].q[i - 1] = %lf, qRgLimit = %lf<br>\n", i, fit[0], data[j].q[i - 1], sqrt(-3.0*fit[0])*data[j].q[i - 1], qRgLimit );
			
				if(sqrt(-3.0*fit[0])*data[j].q[i - 1] < qRgLimit || fit[0] > 0.0){   /* keep repeating the calculation until qRg is as large as it can be */

					data[j].qRg = i;
					data[j].chi2 = chi2;
					data[j].Rg = sqrt(-3.0*fit[0]);
					data[j].RgErr = 3.0*sqrt(data[j].chi2*corr[0][0])/2.0/sqrt(-3.0*fit[0]);
					data[j].I0 = exp(fit[1]);
					data[j].I0Err = exp(fit[1])*sqrt(data[j].chi2*corr[1][1]);
					data[j].to_data_point = i + 1;
					//fprintf(stdout,"j=%i, i=%i, data[j].I0=%lf, data[j].I0Err=%lf, to_data_point=%i, data->n=%i<br>\n", j, i, data[j].I0, data[j].I0Err, data->n );

				} else {
					status = 1;
				}
			} else {
				status = 0;
			}
		}
		
		if (data[j].correct == 'R') { /* if required, calculate the scale factor */
			data[j].m = data[0].I0 / data[0].Drho / data[0].Drho / (data[j].I0 / data[j].Drho / data[j].Drho);
			//fprintf(stdout,"data[%i].m=%lf, data[0].I0=%lf, data[0].Drho=%lf, data[%i].I0=%lf, data[%i].Drho=%lf\n", j, data[j].m, data[0].I0, data[0].Drho, j, data[j].I0, j, data[j].Drho );
		}

	}
	
	fprintf(stdout,"\t\t<TABLE align=center width=\"80%%\" border=1>\n");   /* construct a table for the Guinier analysis results */
	fprintf(stdout,"\t\t<CAPTION><EM>Guinier analysis and scale factor calculation</EM></CAPTION>\n");
	fprintf(stdout,"\t\t\t<TR>");
	fprintf(stdout,"<TH rowspan=2 width=\"10%%\">");
	fprintf(stdout,"<TH colspan=2><i>R<sub>g</sub></i>");
	fprintf(stdout,"<TH colspan=2><i>I</i>(0)");
	fprintf(stdout,"<TH colspan=2><i>qR<sub>g</sub></i>");
	fprintf(stdout,"<TH colspan=2><i>q points used</i>");
	fprintf(stdout,"<TH rowspan=2 width=\"10%%\"><i>&chi;</i><sup>2</sup>");
	fprintf(stdout,"<TH rowspan=2>&Delta;<i>&rho;</i></i><sup>2</sup><BR>(10<sup>20</sup> cm<sup>-4</sup>)");
	fprintf(stdout,"<TH rowspan=2>Scale Factor\n");
	
	fprintf(stdout,"\t\t\t<TR>");
	fprintf(stdout,"<TH width=\"10%%\">Value");
	fprintf(stdout,"<TH width=\"10%%\">ESD");
	fprintf(stdout,"<TH width=\"10%%\">Value");
	fprintf(stdout,"<TH width=\"10%%\">ESD");
	fprintf(stdout,"<TH width=\"10%%\">Min");
	fprintf(stdout,"<TH width=\"10%%\">Max");
	fprintf(stdout,"<TH width=\"10%%\">From");
	fprintf(stdout,"<TH width=\"10%%\">To");

	for (j = 0; j < data->num_contrast_points; j ++) { /* print the results of the Guinier analysis */
		fprintf(stdout,"\t\t\t<TR align=center><TH>%d<TD>%.2lf<TD>%.2lf<TD>%.3lf<TD>%.3lf<TD>%.2lf<TD>%.2lf<TD>%d<TD>%d<TD>%.2lf<TD>%.3lf<TD>%.2lf\n",j + 1, data[j].Rg, data[j].RgErr, data[j].I0, data[j].I0Err, data[j].Rg*data[j].q[data[j].Rgp1 - 1], data[j].Rg*data[j].q[data[j].qRg - 1], data[j].from_data_point, data[j].to_data_point, data[j].chi2, data[j].Drho*data[j].Drho, data[j].m);
	}

	fprintf(stdout,"\t\t</TABLE>\n");
	
	fprintf(stdout,"<BR><BR>\n");

	fprintf(stdout,"\t\t<TABLE align=center width=\"80%%\" border=1>\n");   /* construct a table for the mean-square deviations */
	fprintf(stdout,"\t\t\t<CAPTION><EM>Deviations between the composite scattering functions and the contrast variation series</EM></CAPTION>\n");
	fprintf(stdout,"\t\t\t<TR>");
	fprintf(stdout,"<TH rowspan=\"2\" width=\"%d%%\">", (int)(100/(data->num_contrast_points+2)));
	fprintf(stdout,"<TH colspan=\"%d\">[&Delta;<i>I</i>(<i>q</i>)/<i>&sigma;</i>(<i>I</i>(<i>q</i>))]<sup>2</sup>", data->num_contrast_points);
	fprintf(stdout,"<TH width=\"%d%%\" rowspan=\"2\"><i>&chi;</i><sup>2</sup>", (int)(100/(data->num_contrast_points+2)));

	fprintf(stdout,"\t\t\t<TR>");
	for(i = 0; i < data->num_contrast_points; i ++) {
		fprintf(stdout,"<TH width=\"%d%%\">%d",(int)(100/(data->num_contrast_points+2)), i+1);
	}

	/* do the CSF analysis - MSDs are tabulated from within this function */

	double **SVD_data, **SVD_y_div_sigma;
	DOUBLE2D(SVD_data, MAX_Q_POINTS, num_factors, "SVD_data", "main", i, j);
	DOUBLE2D(SVD_y_div_sigma, MAX_Q_POINTS, num_factors, "SVD_y_div_sigma", "main", i, j);

	for (i = 0; i < data->n; i ++) { /* for each q data point... */

		doSingularValueDecomposition(data, model, i, SVD_data[i], SVD_y_div_sigma[i]); /* ...will work out what the intensity is for each form factor at this Q data point */

	}

	fprintf(stdout,"\t\t</TABLE>\n");
	fprintf(stdout,"<BR><BR>\n");

	fprintf(stdout,"\t\t<TABLE align=center border=1>\n");   /* construct a table for the mean-square deviations */
	fprintf(stdout,"\t\t\t<CAPTION><EM>Singular Value Decomposition (SVD) values</EM><br><br>\n");
	fprintf(stdout,"For each row, if some values are large and others are small, then at this Q point the number of signals seen in the data is only as many signals having large values,<br>\n");
	fprintf(stdout,"and the SVD values of the signals having small SVD values should be set to zero (for no signal) by setting the SVD threshold higher than the small SVD values at this Q point.\n<br>\n");
	fprintf(stdout,"The SVD threshold value used is %17.30lf. SVD values less than this are marked with * and are set to zero.</CAPTION>\n<br>\n", SVC_threshold );

	fprintf(stdout,"\t\t\t<TR>\n");
	fprintf( stdout,"<TD align=center width=10><nobr><B>&nbsp;&nbsp;&nbsp;Q point&nbsp;&nbsp;&nbsp;</B></nobr></TD>" );
	for( ix3 = 0; ix3 < num_factors; ix3++ ) {

		ix1 = order_of_factors_1[ix3] + 1;
		ix2 = order_of_factors_2[ix3] + 1;

		fprintf( stdout,"<TD align=center width=20><B>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<nobr>SVD %d,%d</nobr>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</B></TD>\n", ix1, ix2 );

	}
	fprintf(stdout,"\t\t\t</TR><TR>\n");
	fprintf( stdout,"<TD align=center width=10><nobr><B>&nbsp;</B></nobr></TD>" );
	for( ix3 = 0; ix3 < num_factors; ix3++ ) {

		fprintf( stdout,"<TD align=right width=20><span class='smallblueitalic'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<nobr>y%d/sigma%d</nobr>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span></TD>\n", ix3+1, ix3+1 );

	}
	for(i = (min_start_Q_point - 1); i < max_num_read_Q_points_in_data_files; i ++) {
		fprintf( stdout,"<TR><TD align=center>%i</TD>\n", i+1 );

		for( ix3 = 0; ix3 < num_factors; ix3++ ) {

			if (SVD_data[i][ix3] < SVC_threshold) {
				fprintf( stdout,"<TD align=center>&nbsp;&nbsp;&nbsp;*%17.30lf&nbsp;&nbsp;&nbsp;</TD>", SVD_data[i][ix3] );
			} else {
				fprintf( stdout,"<TD align=center>&nbsp;&nbsp;&nbsp;&nbsp;%17.30lf&nbsp;&nbsp;&nbsp;</TD>", SVD_data[i][ix3] );
			}

		}
		fprintf( stdout,"\n</TR>\n" );
		fprintf( stdout,"<TR><TD align=center>&nbsp;</TD>\n" );

		for( ix3 = 0; ix3 < num_factors; ix3++ ) {

			fprintf( stdout,"<TD align=right><span class='smallblueitalic'>&nbsp;&nbsp;&nbsp;&nbsp;%17.30lf&nbsp;&nbsp;&nbsp;</span></TD>", SVD_y_div_sigma[i][ix3] );
		}
		fprintf( stdout,"\n</TR>\n" );
	}
	fprintf(stdout,"\t\t\t</TR>\n");
	fprintf(stdout,"\t\t</TABLE>\n");
	fprintf(stdout,"<BR><BR>\n");


	/* output the scattering file for each factor, including the form factors and structure factors. */
	/* these scattering files are decomposed from the original scattering file. */
	/* adding up all these component scattering files should give the original scattering file. */

	char c1[2], c2[2];

	for( ix3 = 0; ix3 < num_factors; ix3++ ) {

		ix1 = order_of_factors_1[ix3] + 1;
		ix2 = order_of_factors_2[ix3] + 1;

		/* assign the output filename for this particular scattering file */

		if (num_form_factors > 9) {
			sprintf( c1, "%02i", ix1 );
			sprintf( c2, "%02i", ix2 );
		} else {
			sprintf( c1, "%d", ix1 );
			sprintf( c2, "%d", ix2 );
		}

		if (argc > 2) { /* construct the output filename */
			strcpy(output_filename, argv[1]);
			strcat(output_filename, "i");
			strcat(output_filename, c1);
			strcat(output_filename, c2);
			strcat(output_filename, "_");
			strcat(output_filename, argv[2]);
			strcat(output_filename, ".dat");
		} else { /* if there are no cmd line arguments just go for simple filenames */
			strcpy(output_filename, "i");
			strcat(output_filename, c1);
			strcat(output_filename, c2);
			strcat(output_filename, "_");
			strcat(output_filename, ".dat");
		}

		/* open the output file for this particular scattering file */

		OPEN_FILE(output_file_ptr, output_filename, "w");

		/* output file header for this particular scattering file */

		if (ix1 == ix2) {
			fprintf( output_file_ptr, "# Subunit-%d composite scattering function", ix1 );
		} else {
			fprintf( output_file_ptr, "# Inter-subunit %d-%d composite scattering function", ix1, ix2 );
		}
		fprintf( output_file_ptr, "\n#     q           I(q)         sigma(I(q))");

		/* output the scattering angles, intensities and errors for this particular scattering file */

		for(i = (min_start_Q_point - 1); i < max_num_read_Q_points_in_data_files; i ++) {

			if (model[ix3].IErr[i] > 0.0) {
				fprintf( output_file_ptr, "\n%10.5lf%16.5E%16.5E", data->q[i], model[ix3].I[i], model[ix3].IErr[i] );
			}
		}

		fclose(output_file_ptr);
	}


	/* output the singular value decomposition (SVD) values for each Q point */

	if (argc > 2) { /* construct the output filename */
		strcpy(output_SVD_filename, argv[1]);
		strcat(output_SVD_filename, "SVD_");
		strcat(output_SVD_filename, argv[2]);
		strcat(output_SVD_filename, ".dat");

		strcpy(output_YdivSigma_filename, argv[1]);
		strcat(output_YdivSigma_filename, "YdivSigma_");
		strcat(output_YdivSigma_filename, argv[2]);
		strcat(output_YdivSigma_filename, ".dat");
	} else { /* if there are no cmd line arguments just go for simple filenames */
		strcpy(output_SVD_filename, "SVD");
		strcat(output_SVD_filename, ".dat");

		strcpy(output_YdivSigma_filename, "YdivSigma");
		strcat(output_YdivSigma_filename, ".dat");
	}

	OPEN_FILE(output_SVD_file_ptr, output_SVD_filename, "w");
	fprintf( output_SVD_file_ptr, "# Singular Value Decomposition (SVD) values for the %d factors", num_factors );
	fprintf( output_SVD_file_ptr, "\n#\t\tQ" );

	OPEN_FILE(output_YdivSigma_file_ptr, output_YdivSigma_filename, "w");
	fprintf( output_YdivSigma_file_ptr, "# The %d intermediate Y/Sigma values (Y = experimental intensity, Sigma = the Singular Value Decomposition (SVD) values)", num_factors );
	fprintf( output_YdivSigma_file_ptr, "\n#\t\tQ" );

	for( ix3 = 0; ix3 < num_factors; ix3++ ) {

		ix1 = order_of_factors_1[ix3] + 1;
		ix2 = order_of_factors_2[ix3] + 1;

		fprintf( output_SVD_file_ptr, "\t\t\tSVD-%d,%d", ix1, ix2 );
		fprintf( output_YdivSigma_file_ptr, "\t\t\tY%d/Sigma%d", ix3+1, ix3+1 );
	}

	for (i = (min_start_Q_point - 1); i < max_num_read_Q_points_in_data_files; i ++) {

		fprintf( output_SVD_file_ptr, "\n\t%10.5lf", data->q[i] );
		fprintf( output_YdivSigma_file_ptr, "\n\t%10.5lf", data->q[i] );

		for( ix3 = 0; ix3 < num_factors; ix3++ ) {

			ix1 = order_of_factors_1[ix3] + 1;
			ix2 = order_of_factors_2[ix3] + 1;

			fprintf( output_SVD_file_ptr, "\t%17.30lf", SVD_data[i][ix3] );
			fprintf( output_YdivSigma_file_ptr, "\t%17.30lf", SVD_y_div_sigma[i][ix3] );
		}
	}
	fclose(output_SVD_file_ptr);
	fclose(output_YdivSigma_file_ptr);

	// Output the information for each contrast point, eg. p1.p1.i11 + p2.p2.i22 + p1.p2.i12 = Iexp

	fprintf(stdout,"\t\t<TABLE align=center width=\"80%%\" border=1>\n");   /* construct a table for the mean-square deviations */
	fprintf(stdout,"\t\t\t<CAPTION><EM>Information for each contrast point to derive formula (eg. SLD1 x SLD1 x i11 + SLD2 x SLD2 x i22 + SLD1 x SLD2 x i12 = Iexp x scale-correction)<br>\n");
	fprintf(stdout,"\t\t\tand to derive scale correction for contrast points marked 'R'<br>(eg. scale-correction for contrast-point-X = I0[contrast-point-1] / Rho[contrast-point-1] / Rho[contrast-point-1] / (I0[contrast-point-X] / Rho[contrast-point-X] / Rho[contrast-point-X])<br>where Rho = volume1 x SLD1 + volume2 x SLD2 and I0 is derived by the Guinier analysis. )</EM></CAPTION>\n");
	fprintf(stdout,"\t\t\t<TR>");
	fprintf(stdout,"<td align=center>Contrast Point</td>" );
	for (ix1 = 0; ix1 < num_form_factors; ix1++) {
		fprintf(stdout, "<td align=center>SLD-%i</td><td align=center>Volume-%i</td>", ix1+1, ix1+1 );
	}
	fprintf(stdout,"<td align=center>Rho</td><td align=center>Scale Correction</td>\n" );
	fprintf(stdout,"\t\t\t</TR>\n");
	for (i = 0; i < data->num_contrast_points; i ++) { /* for each contrast point... */
		fprintf(stdout,"\t\t\t<TR>");
		fprintf(stdout,"<td align=center>%d</td>", i+1 );
		for (ix1 = 0; ix1 < num_form_factors; ix1++) {
			fprintf(stdout,"<td align=center>%lf</td><td align=center>%lf</td>", data[i].SLD[ix1], data[i].Volume[ix1] );
		}
		fprintf(stdout,"<td align=center>%lf</td><td align=center>%lf</td>\n", data[i].Drho, data[i].m );
		fprintf(stdout,"\t\t\t</TR>\n\n");
	}
	fprintf(stdout,"\t\t</TABLE>\n");
	fprintf(stdout,"<BR><BR>\n");

	free(data);
	free(model);
	free(fit);
	free(corr);

//	fprintf(stdout,"\t</BODY>\n");
//	fprintf(stdout,"</HTML>\n");
	fprintf(stdout,"</body>\n</html>\n");

	return(1);
		
}

