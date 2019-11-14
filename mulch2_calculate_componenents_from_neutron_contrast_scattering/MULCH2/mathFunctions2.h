/****************************************************************** 
 ***                        mathFunctions2.h                    *** 
 ******************************************************************
 *
 * Additional Info:
 *
 * Date: 9/4/06 last modfied on 25/9/08
 * Author: A.E. Whitten
 * Description:
 * 
 * Function prototypes and definitions from mathFunctions.c
 *
 * Version 1.0: File created
 *
 * Version 1.1: memory allocation functions changed to macros
 *
 * Oct 2016 : Create mathFunctions2.h from mathFunctions.h
 *
 * Oct 2016 : Changes to mathFunctions.c to create mathFunctions2.c include :
 *            Create invertNxNMatrix, using Cholesky instead of matrix determinant.
 *            This Cholesky and invertNxNMatrix written by G. Doherty.
 */
 

#define DOUBLE(p, columns, name, func_name, i){									\
														\
	if (( p = (double *)malloc((columns)*sizeof(double))) == NULL){						\
														\
		fprintf(stderr, "\nError: Can't allocate memory for *%s in function %s.\n", name, func_name);	\
		return(0);											\
														\
	}													\
														\
	for(i = 0; i < columns; i ++)										\
		p[i] = 0.0;											\
														\
}														\
/* double memory allocation macro */

#define DOUBLE2D(p, rows, columns, name, func_name, i, j){							\
														\
	if (( p = (double **)malloc((rows)*sizeof(double))) == NULL){						\
														\
		fprintf(stderr, "\nError: Can't allocate memory for **%s in function %s.\n", name, func_name);	\
		return(0);											\
														\
	}													\
														\
	if (( p[0] = (double *)malloc((rows)*(columns)*sizeof(double))) == NULL){				\
														\
		fprintf(stderr, "\nError: Can't allocate memory for *%s in function %s.\n", name, func_name);	\
		return(0);											\
														\
	}													\
														\
	for(i = 1; i < rows; i ++)										\
		p[i] = p[i - 1] + columns;  									\
														\
	for(i = 0; i < rows; i ++)										\
		for(j = 0; j < columns; j ++)									\
			p[i][j] = 0.0;										\
														\
}														\
/* 2D double memory allocation macro */

#define LONG(p, columns, name, func_name){									\
														\
	if (( p = (long *)malloc((columns)*sizeof(long))) == NULL){						\
														\
		fprintf(stderr, "\nError: Can't allocate memory for *%s in function %s.\n", name, func_name);	\
		return(0);											\
														\
	}													\
														\
	for(i = 0; i < columns; i ++)										\
		p[i] = 0.0;											\
														\
}														\
/* long memory allocation macro */

int invert2x2Matrix(double **A, double **Ai);

int invert3x3Matrix(double **A, double **Ai);

void Cholesky (int n, double a[], int num_columns_A, double *b);

int invertNxNMatrix(int n, double **A, double **Ai, int num_columns_A);

int multiplyMatrixAndVector(double **A, int rowA, int colA, double *B, double *AxB);

int multiplyMatrixAndMatrix(double **A, int rowA, int colA, double **B, double **AxB);

double polynomialFit(int order, double *X, double *Y, double *sigY, int n, double *fit, double **corr, int start_point);

double polynomialFit2(int order, double *X, double *Y, double *sigY, int n, double *fit, double **corr, int start_point);

void svdcmp(double **a, int m, int n, double w[], double **v);

void svbksb(double **u, double w[], double **v, int m, int n, double b[], double x[], double SVC_threshold, double *SVD_y_div_sigma);

