/****************************************************************** 
 ***                        mathFunctions2.c                    *** 
 ******************************************************************
 *
 * Additional Info:
 *
 * Date: 9/4/06
 * Author: A.E. Whitten
 * Description:
 * 
 * Function prototypes and definitions for matrix operations
 *
 * Version 1.0: File created
 *
 * Version 1.1: memory allocation functions changed to macros
 *
 * Oct 2016 : Create mathFunctions2.c from mathFunctions.c
 *
 * Oct 2016 : Changes to mathFunctions.c to create mathFunctions2.c include :
 *            Create invertNxNMatrix, using Cholesky instead of matrix determinant.
 *            This Cholesky and invertNxNMatrix written by G. Doherty.
 *            Add svdcmp and svbksb (Singular Value Decomposition) from
 *            "Numerical Recipes in C++ - The Art of Scientific Computing" Second Edition
 *            Press, Teukolsky, Vetterling and Flannery,
 *            Cambridge University Press 2002.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "mathFunctions2.h"

int invert2x2Matrix(double **A, double **Ai){

	Ai[0][0] =  1.0/(A[0][0]*A[1][1] - A[0][1]*A[1][0])*A[1][1];
	Ai[0][1] = -1.0/(A[0][0]*A[1][1] - A[0][1]*A[1][0])*A[0][1];
	Ai[1][0] = -1.0/(A[0][0]*A[1][1] - A[0][1]*A[1][0])*A[1][0];
	Ai[1][1] =  1.0/(A[0][0]*A[1][1] - A[0][1]*A[1][0])*A[0][0];

	return(1);

}

int invert3x3Matrix(double **A, double **Ai){   /* Inverts a given 3x3 matrix */

	/* A is the matrix; Ai is the inverse of matrix A */

        double det;

	/* det is the determinant of matrix A */

        det = A[0][0]*A[1][1]*A[2][2]+ 
              A[1][0]*A[2][1]*A[0][2]+
              A[2][0]*A[1][2]*A[0][1]-
              A[0][2]*A[1][1]*A[2][0]-
              A[0][0]*A[2][1]*A[1][2]-
              A[0][1]*A[1][0]*A[2][2];

        /* http://www.dr-lex.be/random/matrix_inv.html */

        Ai[0][0] = (A[1][1]*A[2][2] - A[1][2]*A[2][1])/det;
        Ai[0][1] = (A[0][2]*A[2][1] - A[0][1]*A[2][2])/det;
        Ai[0][2] = (A[0][1]*A[1][2] - A[0][2]*A[1][1])/det;
        Ai[1][0] = (A[1][2]*A[2][0] - A[1][0]*A[2][2])/det;
        Ai[1][1] = (A[0][0]*A[2][2] - A[0][2]*A[2][0])/det;
        Ai[1][2] = (A[0][2]*A[1][0] - A[0][0]*A[1][2])/det;
        Ai[2][0] = (A[1][0]*A[2][1] - A[1][1]*A[2][0])/det;
        Ai[2][1] = (A[0][1]*A[2][0] - A[0][0]*A[2][1])/det;
        Ai[2][2] = (A[0][0]*A[1][1] - A[0][1]*A[1][0])/det;

        return(1);

}

void Cholesky (int n, double a[], int num_columns_A, double *b){

	/* Cholesky decomposition A = L D LT					*/
	/* Computes the inverse of a real, symmetric, positive definite A	*/
	/* stored in an array dimensioned a[][num_columns_A]			*/
	/* on entry, B contains the identity I, on exit inverse of A		*/

	int	i, j, k, nn, kcol;
	int	ii, jj, kk, kkk;
	double	fac, fac1, fac2, *bb, *bbi, *bbj;

	nn = n - 1;
	ii = 0;
	for ( i=0; i< nn; i++) {
		fac2 = a[ii];
		jj = ii;
		fac1 = 1.0 / fac2;
		for ( j= i+1; j< n; j++) {
			jj += num_columns_A;
			a[jj] *= fac1;
		}
		jj = ii;
		for ( j= i+1; j< n; j++) {
			jj += num_columns_A;
			fac = a[jj] * fac2;
			kk = jj;
			kkk = j*num_columns_A + j;
			for ( k= j; k< n; k++) {
				a[kkk] -= a[kk] * fac;
				kk += num_columns_A;
				kkk += num_columns_A;
			}
		}
		a[ii] = fac1;
		ii += (num_columns_A+1);
	}
	a[ii] = 1.0 / a[ii];

	/* invert steps */

	bb = b;
	for (kcol=0; kcol< n; kcol++) {
		bbi = bb;
		for ( i= 0; i< nn; i++) {
			jj = i * (num_columns_A + 1);
			bbj = bbi;
			for ( j=i+1; j< n; j++) {
				jj += num_columns_A;
				bbj += num_columns_A;
				(*bbj) -= a[jj] * (*bbi);
			}
			bbi += num_columns_A;
		}
		ii = 0;
		bbi = bb;
		for ( i= 0; i< n; i++) {
			(*bbi) *= a[ii];
			bbi += num_columns_A;
			ii += (num_columns_A+1);
		}
		bbi = bb + nn*num_columns_A;
		for ( i= nn; i> 0; i--) {
			jj = i*num_columns_A;
			bbj = bb;
			for ( j=0; j< i; j++) {
				(*bbj) -= a[jj] * (*bbi);
				jj++;
				bbj += num_columns_A;
			}
			bbi -= num_columns_A;
		}
		bb++;
	}

	return;

} /* Cholesky */


int invertNxNMatrix(int n, double **A, double **Ai, int num_columns_A){   /* Inverts a given NxN matrix using Cholesky */

	double *b, *bb;
	int	i, j, ii, ij, k, kk;

	//fprintf(stdout,"</table>Start invertNxNMatrix<BR>\n" );

	/* the matrix */

	b = (double *) malloc(n*n*sizeof(double));

	/* the identity matrix */	

	bb = (double *) malloc(n*n*sizeof(double));

	kk = 0;
	ii = 0;
	k = 0;
	for (i=0; i<n; i++) {
		ij = ii;
		//fprintf(stdout,"i=%i, n=%i, ii=%i, ij=%i<BR>\n", i, n, ii, ij );
		for (j=0; j<n; j++) {
			//fprintf(stdout,"   j=%i, n=%i, k=%i, A[i][j]=%0.3f<BR>\n", j, n, k, A[i][j] );
			b[k] = A[i][j];
			ij++;
			bb[k] = 0.0;
			k++;
		}
		ii += num_columns_A;
		bb[kk] = 1.0;
		kk += (n+1); 
	}

	Cholesky (n, b, n, bb);
	kk = 0;
	ii = 0;
	for (i=0; i<n; i++) {
		ij = ii;
		for(j=0; j<n; j++) {
			Ai[i][j]= bb[kk];
			ij++;
			kk++;
		}
		ii += num_columns_A;
	}
	free(b);
	free(bb);

	return(1);
}


int multiplyMatrixAndVector(double **A, int rowA, int colA, double *B, double *AxB){   /* Multiplies a matrix by a vector */

	/* A is a matrix; rowA is the number of rows in matrix A; colA is the number of columns in matrix A;
	 * B is a vector; AxB is A x B */
        int i, k;
        
	for(i = 0; i < rowA; i ++){

                AxB[i] = 0.0;

                for(k = 0; k < colA; k ++) {
                        AxB[i] += A[i][k]*B[k];
		}

        }

        return(1);

}


int multiplyMatrixAndMatrix(double **A, int rowA, int colA, double **B, double **AxB){   /* Multiplies a matrix by a vector */

	/* A is a matrix; rowA is the number of rows in matrix A; colA is the number of columns in matrix A;
	 * B is a matrix; AxB is A x B */

        int i, j, k;

        for(i = 0; i < rowA; i ++)
		for(j = 0; j < colA; j ++){
               
			AxB[i][j] = 0.0;

               		for(k = 0; k < colA; k ++) {
                        	AxB[i][j] += A[i][k]*B[k][j];
			}
        	}

        return(1);

}

double polynomialFit(int order, double *X, double *Y, double *sigY, int n, double *fit, double **corr, int start_point){

	int i, j, k, nDisregarded = 0;
	double *A, **B, **Bi, *M, chi2, w, delta, this_sigY;

	if(order > 2){

		fprintf(stderr, "\npolynomial fit (polynomialFit) can only fit a 2nd order polynomial or less/\n");
		return(0);

	}

        DOUBLE(A, order + 1, "A", "polynomialFit", i);
        DOUBLE(M, order + 1, "M", "polynomialFit", i);
        DOUBLE2D(B, order + 1, order + 1, "B", "polynomialFit", i, j);
        DOUBLE2D(Bi, order + 1, order + 1, "Bi", "polynomialFit", i, j);

	for(i = 0; i < n; i ++){

		this_sigY = sigY[i];
		if (this_sigY == 0) 
			this_sigY = 1;

		if ((i + 1) >= start_point) {
			
			w = 1.0/(this_sigY*this_sigY);

			for(j = 0; j <= order; j ++){

				A[j] += w*Y[i]*pow(X[i], order - j);

				for(k = 0; k <= order; k ++)
					B[j][k] += w*pow(X[i], order - j)*pow(X[i], order - k);


			}

		}

	}

	//if(order == 1) invert2x2Matrix(B, Bi);
	//else invert3x3Matrix(B, Bi);
	int matrix_n = order + 1;
	int num_columns_A = order + 1;
	invertNxNMatrix(matrix_n, B, Bi, num_columns_A);

	multiplyMatrixAndVector(Bi, order + 1, order + 1, A, M);

	chi2 = 0.0;

	if(n > order + 1)
		for(i = 0; i < n; i ++){

			this_sigY = sigY[i];
			if (this_sigY == 0) 
				this_sigY = 1;

			if ((i + 1) >= start_point) {
				
				delta = 0.0;

				for(j = 0; j <= order; j ++)
					delta += M[j]*pow(X[i], order - j);

				chi2 += (delta - Y[i])*(delta - Y[i])/this_sigY/this_sigY;
			
			}

			else nDisregarded ++;

		}

	chi2 = chi2/(double)((n - nDisregarded) - (order + 1));

	for(j = 0; j <= order; j ++){

		fit[j] = M[j];

		for(k = 0; k <= order; k ++)
			 corr[j][k] = Bi[j][k];

	}
	
	free(A);
	free(M);
	free(B);
	free(Bi);

	return(chi2);

}

double polynomialFit2(int order, double *X, double *Y, double *sigY, int n, double *fit, double **corr, int start_point){

	int i, j, k, nDisregarded = 0;
	double *A, **B, **Bi, *M, chi2, w, delta, this_sigY;

	if(order > 2){

		fprintf(stderr, "\npolynomial fit (polynomialFit2) can only fit a 2nd order polynomial or less.\n");
		return(0);

	}

        DOUBLE(A, order + 1, "A", "polynomialFit", i);
        DOUBLE(M, order + 1, "M", "polynomialFit", i);
        DOUBLE2D(B, order + 1, order + 1, "B", "polynomialFit", i, j);
        DOUBLE2D(Bi, order + 1, order + 1, "Bi", "polynomialFit", i, j);

	for (i = 0; i < n; i ++) {

		this_sigY = sigY[i];
		if (this_sigY == 0) 
			this_sigY = 1;

		if ((i + 1) >= start_point) {
			
			w = 1.0/(this_sigY*this_sigY);

			for (j = 0; j <= order; j ++) {

				A[j] += w*Y[i]*pow(X[i], order - j);

				for (k = 0; k <= order; k ++) {
					B[j][k] += w*pow(X[i], order - j)*pow(X[i], order - k);
				}
			}

		}

	}

	//if(order == 1) invert2x2Matrix(B, Bi);
	//else invert3x3Matrix(B, Bi);
	int matrix_n = order + 1;
	int num_columns_A = order + 1;
	invertNxNMatrix(matrix_n, B, Bi, num_columns_A);

	multiplyMatrixAndVector(Bi, order + 1, order + 1, A, M);

	chi2 = 0.0;

	if (n > order + 1) {
		for (i = 0; i < n; i ++) {

			this_sigY = sigY[i];
			if (this_sigY == 0) 
				this_sigY = 1;

			if ((i + 1) >= start_point) {
				
				delta = 0.0;

				for (j = 0; j <= order; j ++)
					delta += M[j]*pow(X[i], order - j);

				chi2 += (delta - Y[i])*(delta - Y[i])/this_sigY/this_sigY;
			
			}

			else nDisregarded ++;

		}
	}

	chi2 = chi2/(double)((n - nDisregarded) - (order + 1));

	for (j = 0; j <= order; j ++) {

		fit[j] = M[j];

		for(k = 0; k <= order; k ++)
			 corr[j][k] = Bi[j][k];

	}
	
	free(A);
	free(M);
	free(B);
	free(Bi);

	return(chi2);

}


static double sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -1 * fabs(a))
static double maxarg1,maxarg2;
#define DMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))
static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ? (iminarg1) : (iminarg2))


double pythag(const double a, const double b) {
// Computes .a2 C b 2 /1=2 without destructive underflow or overflow.
//	double absa=abs(a);
//	double absb=abs(b);
//	return (absa > absb ? absa*sqrt(1.0+SQR(absb/absa)) :
//		(absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb))));
	double return_value = 0.0;
	double absa=fabs(a);
	double absb=fabs(b);
	if (absa > absb) {
		if ((absa != 0.0) && (absa != -0.0)) {
			return_value = absa*sqrt(1.0+SQR(absb/absa));
		} else {
			return_value = absa*sqrt(1.0+SQR(absb));
			fprintf(stdout,"PYTHAG-0001 : DIVIDE BY ZERO<br>\n" );// This is not defined if denominator is zero
		}
	} else {
		if (absb == 0.0) {
			return_value = 0.0;
		} else {
			if ((absb != 0.0) && (absb != -0.0)) {
				return_value = absb*sqrt(1.0+SQR(absa/absb));
			} else {
				return_value = absb*sqrt(1.0+SQR(absa));
				fprintf(stdout,"PYTHAG-0002 : DIVIDE BY ZERO<br>\n" );// This is not defined if denominator is zero
			}
		}
	}
	return return_value;
}



void svdcmp(double **a, int m, int n, double *w, double **v)
/* Given a matrix a[1..m][1..n], this routine computes its singular value decomposition, A = U ·W ·V^T. 
The matrix U replaces a on output. The diagonal matrix of singular values W is output as a vector w[1..n]. 
The matrix V (not the transpose V^T ) is output as v[1..n][1..n]. */
{
	double pythag(double a, double b);
	int flag,i,its,j,jj,k,l,nm;
	double anorm,c,f,g,h,s,scale,x,y,z,*rv1;
	int max_iterations = 1000;

	DOUBLE(rv1, 1000, "rv1", "svdcmp", i);
	//DOUBLE(rv1, n, "rv1", "svdcmp", i);
	for ( i=0; i<1000; i++ ) {
		rv1[i] = 0.0;
	}

	g = 0.0; //Householder reduction to bidiagonal form.
	scale = 0.0;
	anorm = 0.0;
	for (i=0;i<n;i++) {
		l = i + 1;
		rv1[i] = scale * g;
		g=0.0;
		s=0.0;
		scale=0.0;
		if (i < m) {
			for (k=i;k<m;k++) {
				scale = scale + fabs(a[k][i]);
			}
			if ((scale != 0.0) && (scale != -0.0)) {
				for (k=i;k<m;k++) {
					a[k][i] = a[k][i] / scale;
					s = s + (a[k][i] * a[k][i]);
				}
				f=a[i][i];
				g = -1 * SIGN(sqrt(s),f);
				h = (f * g) - s;
				a[i][i] = f - g;
				for (j=l;j<n;j++) {
					for (s=0.0,k=i;k<m;k++) {
						s = s + a[k][i]*a[k][j];
					}
					if ((h != 0.0) && (h != -0.0)) {
						f = s / h;
					} else {
						f = s;
						fprintf(stdout,"SVDCMP-0001 : DIVIDE BY ZERO<br>\n" );// f=s/h; This is not defined if denominator is zero
					}
					for (k=i;k<m;k++) {
						a[k][j] = a[k][j] + (f * a[k][i]);
					}
				}
				for (k=i;k<m;k++) {
					a[k][i] = a[k][i] * scale;
				}
			}
		}
		w[i]=scale *g;
		g=0.0;
		s=0.0;
		scale=0.0;
		if (i < m && i != (n-1)) {
			for (k=l;k<n;k++) {
				scale = scale + fabs(a[i][k]);
			}
			if ((scale != 0.0) && (scale != -0.0)) {
				for (k=l;k<n;k++) {
					a[i][k] = a[i][k] / scale;
					s = s + (a[i][k] * a[i][k]);
				}
				f=a[i][l];
				g = -1 * SIGN(sqrt(s),f);
				h = (f * g) - s;
				a[i][l] = f - g;
				for (k=l;k<n;k++) {
					if ((h != 0.0) && (h != -0.0)) {
						rv1[k] = a[i][k] / h;
					} else {
						rv1[k] = a[i][k];
						fprintf(stdout,"SVDCMP-0002 : DIVIDE BY ZERO<br>\n" );// rv1[k]=a[i][k]/h; This is not defined if denominator is zero
					}
				}
				for (j=l;j<m;j++) {
					for (s=0.0,k=l;k<n;k++) {
						s = s + (a[j][k] * a[i][k]);
					}
					for (k=l;k<n;k++) {
						a[j][k] = a[j][k] + (s * rv1[k]);
					}
				}
				for (k=l;k<n;k++) {
					a[i][k] = a[i][k] * scale;
				}
			}
		}
		anorm=DMAX(anorm,(fabs(w[i])+fabs(rv1[i])));
	}
	for (i=(n-1);i>=0;i--) { //Accumulation of right-hand transformations.
		if (i < (n - 1)) {
			if ((g != 0.0) && (g != -0.0)) {

				for (j=l;j<n;j++) {//Double division to avoid possible underflow.
					v[j][i]= a[i][j];
					if ((a[i][l] != 0.0) && (a[i][l] != -0.0)) {
						v[j][i]= a[i][j] / a[i][l]; 
					} else {
						v[j][i]= a[i][j];
						fprintf(stdout,"SVDCMP-0003 : DIVIDE BY ZERO<br>\n" );//This is not defined if denominator is zero
					}
					if ((g != 0.0) && (g != -0.0)) {
						v[j][i] = v[j][i] / g;
					} else {
						fprintf(stdout,"SVDCMP-0004 : DIVIDE BY ZERO<br>\n" );//This is not defined if denominator is zero
					}
				}
				for (j=l;j<n;j++) {
					for (s=0.0,k=l;k<n;k++) {
						s = s + a[i][k] * v[k][j];
					}
					for (k=l;k<n;k++) {
						v[k][j] = v[k][j] + (s * v[k][i]);
					}
				}
			}
			for (j=l;j<n;j++) {
				v[i][j] = 0.0;
				v[j][i] = 0.0;
			}
		}
		v[i][i] = 1.0;
		g = rv1[i];
		l = i;
	}
	for ( i=(IMIN(m,n)-1); i>=0; i-- ) { //Accumulation of left-hand transformations.
		l = i + 1;
		g = w[i];
		for (j=l;j<n;j++) {
			a[i][j] = 0.0;
		}
		if ((g != 0.0) && (g != -0.0)) {
			if ((g != 0.0) && (g != -0.0)) {
				g = 1.0 / g;
			} else {
				g = 1.0;
				fprintf(stdout,"SVDCMP-0005 : DIVIDE BY ZERO<br>\n" );//This is not defined if denominator is zero
			}
			for (j=l;j<n;j++) {
				for (s=0.0,k=l;k<m;k++) {
					s = s + (a[k][i] * a[k][j]);
				}
				//f = ( s / a[i][i] ) * g; //??? WHAT SHOULD THIS BE IF DENOMINATOR IS ZERO?
				if ((a[i][i] != 0.0) && (a[i][i] != -0.0)) {
					f = (s / a[i][i]) * g;
				} else {
					f = s * g;
					fprintf(stdout,"SVDCMP-0006 : DIVIDE BY ZERO<br>\n" );//This is not defined if denominator is zero
				}
				for (k=i;k<m;k++) {
					a[k][j] = a[k][j] + (f * a[k][i]);
				}
			}
			for (j=i;j<m;j++) {
				a[j][i] = a[j][i] * g;
			}
		} else {
			for (j=i;j<m;j++) {
				a[j][i] = 0.0;
			}
		}
		a[i][i] = a[i][i] + 1;
	}
	for ( k=(n-1); k>=0; k-- ) { //Diagonalization of the bidiagonal form: Loop over
		for (its=0;its<max_iterations;its++) { //singular values, and over allowed iterations.
			flag = 1;
			nm = k; // "nm" needs to be initialised, otherwise "y = a[j][nm];" FURTHER BELOW IN "for (j=1;j<m;j++) {" WILL CAUSE "Segmentation fault".
			for (l=k;l>=0;l--) { //Test for splitting.
				nm = l - 1; //Note that rv1[1] is always zero.
				if ((fabs(rv1[l])+anorm) == anorm) {
					flag = 0;
					break;
				}
				if ((fabs(w[nm])+anorm) == anorm) {
					break;
				}
			}
			if (flag == 1) {
				c = 0.0; //Cancellation of rv1[l], if l > 1.
				s = 1.0;
				for (i=l;i<=k;i++) {
					f = s * rv1[i];
					rv1[i] = c * rv1[i];
					if ((fabs(f)+anorm) == anorm) {
						break;
					}
					g = w[i];
					h = pythag(f, g);
					w[i] = h;
					if ((h != 0.0) && (h != -0.0)) {
						h = 1.0 / h;
					} else {
						h = 1.0;
						fprintf(stdout,"SVDCMP-0007 : DIVIDE BY ZERO<br>\n" );//h = 1.0 / h; This is not defined if denominator is zero
					}
					c = g * h;
					s = -1 * f * h;
					for (j=1;j<m;j++) {
						y = a[j][nm];
						z = a[j][i];
						a[j][nm] = (y * c) + (z * s);
						a[j][i] = (z * c) - (y * s);
					}
				}
			}
			z = w[k];
			if (l == k) { //Convergence.
				if (z < 0.0) { //Singular value is made nonnegative.
					w[k] = -1 * z;
					for (j=0;j<n;j++) {
						v[j][k] = -1 * v[j][k];
					}
				}
				break;
			}
			if (its == max_iterations) {
				nrerror("no convergence in max_iterations svdcmp iterations");
			}
			x = w[l]; //Shift from bottom 2-by-2 minor.
			nm = k - 1;
			y = w[nm];
			g = rv1[nm];
			h = rv1[k];
			f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
			g = pythag(f,1.0);
			f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;
			c = 1.0;
			s = 1.0;
			//Next QR transformation:
			for (j=l;j<=nm;j++) {
				i = j+1;
				g = rv1[i];
				y = w[i];
				h = s * g;
				g = c * g;
				z = pythag(f, h);
				rv1[j] = z;
				if ((z != 0.0) & (z != -0.0)) { //This is not defined if denominator is zero
					c = f / z;
					s = h / z;
				} else {
					c = f;
					s = h;
					fprintf(stdout,"SVDCMP-0008 : DIVIDE BY ZERO. f=%lf, h=%lf, z=%lf<br>\n", f, h, z );//This is not defined if denominator is zero
				}
				f = (x * c) + (g * s);
				g = (g * c) - (x * s);
				h = y * s;
				y = y * c;
				for (jj=0;jj<n;jj++) {
					x = v[jj][j];
					z = v[jj][i];
					v[jj][j] = (x * c) + (z * s);
					v[jj][i] = (z * c) - (x * s);
				}
				z = pythag(f, h);
				w[j] = z;
				//Rotation can be arbitrary if z = 0.
				if ((z != 0.0) && (z != -0.0)) {
					z = 1.0 / z;
					c = f * z;
					s = h * z;
				}
				f = (c * g) + (s * y);
				x = (c * y) - (s * g);
				for (jj=0;jj<m;jj++) {
					y = a[jj][j];
					z = a[jj][i];
					a[jj][j] = (y * c) + (z * s);
					a[jj][i] = (z * c) - (y * s);
				}
			}
			rv1[l] = 0.0;
			rv1[k] = f;
			w[k] = x;
		}
	}
	free(rv1);
}


void svbksb(double **u, double *w, double **v, int m, int n, double *b, double *x, double SVC_threshold, double *SVD_y_div_sigma)

/* Solves A·X = B for a vector X, where A is specified by the arrays u[1..m][1..n] , w[1..n],
v[1..n][1..n] as returned by svdcmp. m and n are the dimensions of a, and will be equal for
square matrices. b[1..m] is the input right-hand side. x[1..n] is the output solution vector.
No input quantities are destroyed, so the routine may be called sequentially with different b’s. */

/* For evaluating equation (2.6.7) and obtaining a solution vector x from a right-hand side b, 
given that the SVD of a matrix A has already been calculated by a call to svdcmp. 
Note that this original routine presumed that you have already zeroed the small wj’s. 
It does not do this for you. If you hadn’t zeroed the small wj’s, 
then this routine would have been just as ill-conditioned as any direct method, and you are misusing SVD. 
However, the lines of code containing wmax and wmin do zero out the small wj's. */

/* To help understanding of the SVD values, pass back the tmp array (called SVD_y_div_sigma) to be printed.
As an example, the last step of the SVD algorithm for 6 factors of which 3 SVDs are zero,
where x are the unknowns, y are the experimental intensities, and sigma are the SVD values :

        x  =  V[6x6]  x  [  y1/sigma1  ]
                         [  y2/sigma2  ]
                         [  y3/sigma2  ]
                         [      0      ]
                         [      0      ]
                         [      0      ]

        x  =  y1/sigma1 . (1st col of V) + y2/sigma2 . (2nd col of V) + y3/sigma3 . (3rd col of V)
*/

{
	int jj,j,i;
	double s,*tmp;

	double wmax,wmin;
	wmax=0.0; 	//Will be the maximum singular value obtained.
	for(j=0;j<n;j++) {
		if (w[j] > wmax) {
			wmax=w[j];
		}
	}
	for(j=0;j<n;j++) {
		SVD_y_div_sigma[j] = 999999999; //if see this initial value, then something has gone wrong
	}
	//This is where we set the threshold for singular values allowed to be nonzero. 
	//The constant is typical, but not universal. You have to experiment with your own application.
	//w are the singular values of the 
	//These sigma values are the square roots of the eigenvalues of the Whitten/Trewhella matrix.
	//wmin = wmax * 1.0e-6;
	wmin = wmax * SVC_threshold;
	for(j=0;j<n;j++) {
		if (w[j] < wmin) {
			w[j]=0.0;
		}
	}

	//tmp=vector(1,n);
	DOUBLE(tmp, 1000, "tmp", "svbksb", i);
	for (j=0;j<n;j++) { //Calculate U^T . B
		s=0.0;
		if (w[j] != 0.0) { //Nonzero result only if wj is nonzero.
			for (i=0;i<m;i++) {
				s = s + (u[i][j]*b[i]);
			}
			s = s / w[j]; //This is the divide by wj .
		}
		tmp[j]=s;
		SVD_y_div_sigma[j]=s;
	}
	for (j=0;j<n;j++) { //Matrix multiply by V to get answer.
		s=0.0;
		for (jj=0;jj<n;jj++) {
			s = s + (v[j][jj]*tmp[jj]);
		}
		x[j]=s;
	}
	//free_vector(tmp,1,n);
	free(tmp);
}


void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}


