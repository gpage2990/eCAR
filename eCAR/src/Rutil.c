/*************************************************************
 * Copyright (c) 2009 Steven McKay Curtis and Garritt L. Page
 *
 * I give permission for anyone to use and/or modify these
 * programs provided the following conditions are met:
 *
 * 1) This copyright message is retained.
 * 2) Any changes to the code are noted.
 *
 * These programs are provided WITHOUT ANY WARRNTY.
 *
 *************************************************************/

#include <R.h>
#include <Rmath.h>
#include "matrix.h"
#include <math.h>

#include <R_ext/Lapack.h>
#include <time.h>


//============================================================
// Allocating and Freeing Vectors and Matrices
//
// Note that the Matrix allocation functions allocate the
// memory for the entire matrix in a contiguous block.

//==================================================
// Vectors
double* R_Vector(int n){
    double* v;
    v = (double *) R_alloc(n,sizeof(double));
    return(v);
}

double* R_VectorInit(int n,double init){
    double* v;
    v = (double *) R_alloc(n,sizeof(double));
    for(int i=0; i<n; i++)
	v[i]=init;
    return(v);
}

/* For possible Calloc version
void R_FreeVector(double *v, int n){
    Free(v);
}
*/

//==================================================
// Matrices
double** R_Matrix(int nr, int nc){
    double** m;
    m = (double **) R_alloc(nr,sizeof(double*));
    *m = (double *) R_alloc(nr*nc,sizeof(double));
    for (int i=1; i<nr; i++)
	*(m+i) = *m + i*nc;
    return(m);
}

double** R_MatrixInit(int nr, int nc, double init){
    double** m;
    int i, j;
    m = (double **) R_alloc(nr,sizeof(double*));
    *m = (double *) R_alloc(nr*nc,sizeof(double));
    for (i=1; i<nr; i++)
	*(m+i) = *m + i*nc;
    for (i=0; i<nr; i++)
	for (j=0; j<nc; j++)
	    m[i][j] = init;
    return(m);
}

/* for possible Calloc versions
void R_FreeMatrix(double** m, int nr, int nc){
   Free(*m); // because it was allocated in a contiguous block
   Free(m);
//for (int i=nr-1; i>=nr; i--)
//Free(m[i]);
//Free(m);
}
*/

double** R_Data2Matrix(double* d, int nr, int nc){
    double** m;
    m = (double **) R_alloc(nr,sizeof(double*));
    for (int i=0; i<nr; i++)
	*(m+i) = d + i*nc;
    return(m);
}

int** R_Data2iMatrix(int* d, int nr, int nc){
    int** m;
    m = (int **) R_alloc(nr,sizeof(int*));
    for (int i=0; i<nr; i++)
	*(m+i) = d + i*nc;
    return(m);
}

//============================================================
// Matrix Operations
//


double quform(double *x, double* A, int dim){
    int i,j;
    double sm=0.0;
    for (i=1; i<dim; i++)
	for (j=0; j<i; j++)
	    sm += x[i]*x[j]*A[i*dim+j];
    sm*=2;
    for(i=0; i<dim; i++)
	sm += x[i]*x[i]*A[i*dim+i];
    return(sm);
}

double quform2(double *x, double *y, double* A, int dim){
    int i,j;
    double sm=0.0;
    for (i=0; i<dim; i++)
	for (j=0; j<dim; j++)
	    sm = sm + x[i]*y[j]*A[i*dim+j];
    return(sm);
}



//************************************************************
// Printing Matrices and Vectors
//
// The following functions
//
// Rprintvec, Rprintmat, RprintIvec, RprintImat,
//
// are modified versions of functions
// provided by Howard Seltman at the following web page:
// http://www.stat.cmu.edu/~hseltman/files/vmr.c
//
// I have modified the functions to work with R and to
// provide slightly modified output to suit my tastes.
//

// for doubles
void Rprintvec(char* title, double *v, int l){
    if (title!=NULL)
	Rprintf("%s\n",title);
    for (int i=0; i<l; i++)
	Rprintf("%f\n", v[i]);
    Rprintf("\n");
    return;
}

void Rprintmat(char* title, double **m, int nr, int nc){
    if (title!=NULL)
	Rprintf("%s\n", title);
    for (int i=0; i<nr; i++) {
	for (int j=0; j<nc; j++)
	    Rprintf("%f ", m[i][j]);
	Rprintf("\n");
    }
    Rprintf("\n");
    return;
}



// for integers
void RprintIvec(char* title, int* v, int n){
    if (title!=NULL)
	Rprintf("%s\n", title);
    for (int i=0; i<n; i++)
	Rprintf("%i\n", v[i]);
    Rprintf("\n");
    return;
}

void RprintImat(char* title, int** m, int nr, int nc){
    if (title!=NULL)
	Rprintf("%s\n", title);
    for (int i=0; i<nr; i++) {
	for (int j=0; j<nc; j++)
	    Rprintf("%i ", m[i][j]);
	Rprintf("\n");
    }
    Rprintf("\n");
    return;
}


// Print a vector as a matrix
void RprintVecAsMat(char* title, double *v, int nr, int nc){
    if (title!=NULL)
	Rprintf("%s\n", title);
    for (int i=0; i<nr; i++) {
	for (int j=0; j<nc; j++)
	    Rprintf("%f ", v[i*nc + j]);
	Rprintf("\n");
    }
    Rprintf("\n");
    return;
}


// Print a vector as a matrix
void RprintIVecAsMat(char* title, int *v, int nr, int nc){
    if (title!=NULL)
	Rprintf("%s\n", title);
    for (int i=0; i<nr; i++) {
	for (int j=0; j<nc; j++)
	    Rprintf("%d ", v[i*nc + j]);
	Rprintf("\n");
    }
    Rprintf("\n");
    return;
}


// Matrix tranpose
void mat_transpose(double *mat, double *tmat, int nr, int nc){

	int i, j;
	for(i = 0; i < nr; i++){
		for(j = 0; j < nc; j++){
			tmat[j*nr + i] = mat[i*nc + j];
		}
	}
}


// Matrix tranpose that destroys original matrix.
void mat_transpose2(double *m, int w, int h){
	int start, next, i;
	double tmp;

	for (start = 0; start <= w * h - 1; start++) {
		next = start;
		i = 0;
		do {	i++;
			next = (next % h) * w + next / h;
		} while (next > start);
		if (next < start || i == 1) continue;

		tmp = m[next = start];
		do {
			i = (next % h) * w + next / h;
			m[next] = (i == start) ? tmp : m[i];
			next = i;
		} while (next > start);
	}
}




//============================================================
// Random Number Generation
//


/* Unequal probability sampling; with-replacement case
 * n are the lengths of p and perm. p contains probabilities, perm
 * contains the actual outcomes, and ans contains an array of values
 * that were sampled.
 */

int ran_discrete_unif(int L, int U){
/*************************************************************
 * PURPOSE:
 * Generate a draw from discrete distribution.
 *
 * INPUTS:
 * L = Lower bound (inclusive)
 * U = Upper bound (inclusive) *
 *
 *************************************************************/

    int i, out=L;
    int n = U - L+1;
	double uu, cp;
	double p = 1.0/(double) n;

	uu = runif(0,1);

	cp = 0.0;
    for(i = 0; i < n; i++){
        cp = cp + p;
		if(uu < cp) {
			break;
		}
		out = out + 1;
	}
	return(out);
}


void ran_mvnorm(double* m, double* cholV, int dim, double* z, double* out){
/*************************************************************
 * PURPOSE:
 * Generate a draw from a multivariate normal distribution.
 *
 * INPUTS:
 * m     = an array for the mean
 * cholV = an array for the cholesky decomposition
 *         of the variance (note this is a 1D array that
 *         that will be treated as 2D).
 * dim   = dimension of the multivariate normal distribution
 * z     = a scratch array of length dim
 *         to be filled with N(0,1) draws
 *
 * OUTPUTS:
 * out   = final output array to be filled with a random
 *         multivariate normal draw
 *
 *************************************************************/
    int i,j;
    for (i=0; i<dim; i++){
	z[i] = rnorm(0,1);
	out[i] = m[i];
	for (j=0; j<=i; j++)
	    out[i] += cholV[i*dim + j]*z[j];
    }
}


void ran_wish(int nu, double* cholS, int dim, double* z, double* x, double* zeros, double* out){
/*************************************************************
 * PURPOSE:
 * Generate a random draw from Wishart distribution with
 * degrees of freedom nu and parameter matrix S
 *
 * INPUTS:
 * nu    = degrees of freedom
 * cholS = cholesky decomposition of matrix parameter S
 *         This is a 1D array but is accessed
 *         as a 2D array as in cholS[i*dim + j]
 * dim   = the dimension of the Wishart distribution
 * z     = scratch vector of length dim to be passed
 *         to the function ran_mvnorm
 * x     = scratch vector of length dim to be passed
 *         to the function ran_mvnorm
 * zeros = vector of zeros of length dim to be passed
 *         to the function ran_mvnorm as the mean
 * out   = matrix to contain the random draw from the
 *         Wishart distribution
 *
 *************************************************************/

    int i, j, k;

    /* Zero the "out" matrix */
    for (j=0; j<dim; j++)
	for (k=0; k<dim; k++)
	    out[j*dim + k] = 0.0;

    for (i=0; i<nu; i++){
	ran_mvnorm(zeros,cholS,dim,z,x);
	for (j=0; j<dim; j++)
	    for (k=0; k<=j; k++)
		out[j*dim + k] += x[j]*x[k];
    }

    /* fill the upper triangular part with lower triangular part */
    for (j=0; j<dim; j++)
	for (k=0; k<j; k++)
	    out[k*dim + j] = out[j*dim + k];
}


void ran_dirich(double *alpha, int dim, double* scratch, double* out){
/*************************************************************
 * PURPOSE:
 * Generate a random draw from a Dirichlet distribution with
 * parameters alpha
 *
 * INPUTS:
 * alpha   = parameter vector associated with dirichlet
 * dim     = dimension of alpha vector
 * scratch = vector of length dim to be passed
 * out     = array that holds the random values
 *
 *************************************************************/

    int h;
	double sg = 0;
    /* Zero the "out" matrix */
//	RprintVecAsMat("alpha", alpha, 1, dim);
	for(h = 0; h < dim; h++)
		{

			scratch[h] = rgamma(alpha[h], 1);
//			Rprintf("rgamma = %f\n", scratch[h]);
			sg = sg + scratch[h];
//			Rprintf("sg = %f\n", sg);
		}

//	Rprintf("sg = %f\n", sg);
	for(h = 0; h < dim; h++) out[h] = scratch[h]/sg;

//	RprintVecAsMat("out", out, 1, dim);
}





/* The following provides a function that allows me to sample from a
   truncated normal distribution with the following arguments.  This
   function relies heavly on the Rmath library.  As an example here
   are the arguments for the pnorm function
   double pnorm(double x, double mu, double sigma, int lower_tail,
				int give_log );


  	m - mean of the normal prior to truncation;
  	s - standard deviation of normal distribution prior to truncation;
  	a - Lower bound of the truncation
  	b - Upper bound of the truncation */

double rtnorm(double m, double s, double a, double b){

  int jj;
  double a_term, b_term, Un, p, rtdraw, tmp, z;

//	Rprintf("a = %f\n", a);
//	Rprintf("b = %f\n", b);

//	Rprintf("(a-m)/s = %f\n", (a-m)/s);
	a_term = pnorm((a-m)/s, 0.0, 1.0, 1, 0);

//	Rprintf("(b-m)/s = %f\n", (b-m)/s);
	b_term = pnorm((b-m)/s, 0.0, 1.0, 1, 0);

//	Rprintf("a_term = %f\n", a_term);
//	Rprintf("b_term = %f\n", b_term);

	Un = runif(0,1);

	p = a_term + Un*(b_term - a_term);

	rtdraw = m + s*qnorm(p, 0.0, 1.0, 1, 0);

//	Rprintf("rtdraw = %f\n", rtdraw);
//	Rprintf("p = %f\n", p);

	if(p == 1.0){
		jj = 1;

		while(jj != 0){

			a = (a - m)/s;
    		b = (b - m)/s;

//			Rprintf("a = %f\n", a);
			tmp = (a + sqrt(a*a + 4))/2;
//			Rprintf("tmp = %f\n", tmp);
        	z = rexp(1/tmp) + a;
//			Rprintf("z = %f\n", z);
        	Un = runif(0,1);
//			Rprintf("Un = %f\n", Un);
//			Rprintf("exp(-(z - tmp)*(z - tmp)/2) = %f\n", exp(-(z - tmp)*(z - tmp)/2));

        	if((Un <= exp(-(z - tmp)*(z - tmp)/2)) & (z <= b)) jj = 0;
//			Rprintf("jj = %d\n", jj);

		}

//		Rprintf("z = %f\n", z);

		rtdraw = z*s + m;


	}


//	Rprintf("rtdraw = %f\n", rtdraw);

	return(rtdraw);

}


//============================================================
// Multivariate density functions
//
/* The following provides a density function for the multivariate
   normal distribution.  This function relies heavily on matrix.c
   To use it I must create an matrix.o(object file)

	*y -  is the observation vector for which the density will be computed
  	*mu - mean vector
  	*iSig - the inverse covariance matrix as a contiguos vector in row-major form
  	dim - is the dimension of the multivariate distribution
	ld - the log of the determinant of Sigma
	scr - scratch vector to hold
	logout - logical determines if log density is returned
*/

//Multivariate normal
double dmvnorm(double* y, double* mu, double* iSig, int dim, double ld, double* scr, int logout){
    int i;
    double qf, out;
    for(i=0; i<dim; i++) scr[i] = y[i] - mu[i];
    qf = quform(scr,iSig,dim);
    out = -(double) dim*M_LN_SQRT_2PI - 0.5*(ld + qf);
    if (logout)	return out;
    return exp(out);
}



//============================================================
// univariate scaled and shifted t-distribution
//

//
double dsst(double y, double mu, double s, double nu, int logout){

    double qf, out, lc;

	lc = lgamma(0.5*(nu + 1)) - (0.5*log(nu*M_PI) + lgamma(0.5*nu));
	qf = ((y-mu)/s)*((y-mu)/s);

    out = lc - 0.5*(nu+1)*log(1 + (1/nu)*qf);

    if (logout)	return out;
    return exp(out);
}



/* The following provides a density function for the Inverse Wishart
   distribution.  This function relies heavily on matrix.c
   To use it I must create an matrix.o(object file)

	*ASigInv - is ASigma^{-1} that is found in the trace of the density.  It is a
	           contiguos vector of memory that is in row-major form
	detSig - is the determinant of the argument matrix
	*A - is the dimxdim scale matrix that is in contiguos vector row-major form
	detA - is the determinant of A
  	df - degrees of freedom of the inverse-wishart function
  	dim - is the dimension of the scale matrix distribution */

double dinvwish(double *ASigInv, double detSig, double detA, int df, int dim)
	{

		int i;
		double C, density;
		double gamprod;
		double trace;
		double p1, p2, p3, p4;


		gamprod = 1;
		for(i = 0; i < dim; i++)
			{

				gamprod = gamprod*gammafn(0.5*(df + 1 - (i+1)));

			}
//		Rprintf("gamprod = %f\n", gamprod);
		trace = 0.0;
		for(i = 0; i < dim*dim; i++)
			{

				if(i % (dim+1) == 0) trace = trace + ASigInv[i];

			}
//		Rprintf("trace = %f\n", trace);
		p1 = 0.5*df*dim;
		p2 = 0.25*dim*(dim - 1);
		p3 = 0.5*df;
		p4 = -0.5*(df + dim + 1);
//		Rprintf("p1 = %f\n", p1);Rprintf("p2 = %f\n", p2);Rprintf("p3 = %f\n", p3);Rprintf("p4 = %f\n", p4);

		C = 1/(pow(2, p1)*pow(M_PI, p2)*gamprod);
//		Rprintf("detA = %f\n", detA);	Rprintf("detSig = %f\n", detSig);
//		Rprintf("C = %f\n", C);
		density = C*pow(detA, p3)*pow(detSig, p4)*exp(-0.5*trace);
//		Rprintf("density = %f\n", density);
		return(density);


	}



// Inverse Gamma density
// the parameterization here provides a mean of beta/(alpha - 1)
double dinvgamma(double y, double alpha, double beta, int logout){

//	Rprintf("alpha = %f\n", alpha);
//	Rprintf("beta = %f\n", beta);

	double ldens;

	ldens = alpha*log(beta) - lgamma(alpha) - (alpha + 1)*log(y) - (beta/y);

	if(logout) return ldens;
	return exp(ldens);

}


// Normal-Inverse Gamma density
// the parameterization here is such that mu0 is prior mean value with k0 prior "observations"
// and s20 is the prior guess for sig2 with nu0 prior "observations"
double dN_IG(double mu, double sig2, double mu0, double k0, double a0, double b0, int logout){

//	Rprintf("alpha = %f\n", alpha);
//	Rprintf("beta = %f\n", beta);

	double ldens;

//	ldens =  0.5*(log(k0) - log(2*M_PI*sig2)) + a0*log(b0) -
//	         lgammafn(a0) - (a0 + 1)*log(sig2) -
//	         0.5*(1/sig2)*(k0*(mu-mu0)*(mu-mu0) + 2*b0);

	ldens = dnorm(mu, mu0, sqrt(sig2/k0),logout) +
	        dinvgamma(sig2, a0, b0, logout);
//	Rprintf("ldens = %f\n", ldens);
	if(logout){ return ldens;
	}else{return exp(ldens);}

}



double ddirich(double *pi, double *alpha, int C, int logout){
/*************************************************************
 * PURPOSE:
 * Evaluate Dirichlet density with parameters alpha
 *
 * INPUTS:
 * pi   = C - dimensional simplex values
 * alpha     = C - dimensional parameter vector
 * C = dimension of probability vector
 * logout = logical indicating whether log density should be returned
 *
 *************************************************************/

    int ii;
	double suma = 0.0, sumla = 0.0, ldensity=0.0;
    /* Zero the "out" matrix */
//	RprintVecAsMat("alpha", alpha, 1, dim);
	for(ii = 0; ii < C; ii++){

		suma = suma + alpha[ii];
		sumla = sumla + lgammafn(alpha[ii]);

	}

	for(ii = 0; ii < C; ii++){

		ldensity = ldensity + (alpha[ii]-1)*log(pi[ii]);

	}

	ldensity = ldensity + lgammafn(suma) - sumla;

	if(logout) return ldensity;
	return exp(ldensity);

}


// density for univariate truncated normal
double dtnorm(double x, double mu, double sigma, double l, double u, int logout){

	double den, num, ldensity;

	den = pnorm5(u,mu,sigma, 1, 0) - pnorm5(l,mu,sigma, 1, 0);
	num = dnorm(x, mu, sigma,1);

//	Rprintf("den = %f\n", den);

	ldensity = num - log(den);
	if(logout) return ldensity;
	return exp(ldensity);
}




/*************************************************************
 * PURPOSE:
 * Use the singular value decomposition to obtain
 * the square root of square a matrix (or inverse)
 *
 * INPUTS:
 * matrix = a nrxnr matrix in row-major form.  This will have
 *          to be transposed to use the LAPACK routine and
 *          then the U matrix will have to be transposed to be
 *          in row-major form.
 *
 * inverse = a scalar indicating of inverse square matrix
 *           should be returned
 *
 *
 * output = a nrxnr matrix that is the inverse square root
 *
 *
 * svs     = scratch memory to create diagonal matrix whose
 *           diagonal values are the singular values.
 *
 * nr     = number of rows and columns
 *
 * VT    = a scratch array of length n*n.  This is were the
 *         elements of the orthogonal matrix V' will be stored
 *
 * S     = a scratch array of length n.  This is were the positive
 *         values of the diagonal will be stored
 *
 *************************************************************/


void sqrtmat(double *matrix, int inverse, int nr,  double *svs, double *Up, double *VT,
                  double *S, double *scr1){

  int i, j;
  int LDA=nr, LDU=nr,  LDVT = nr;
//  int LDwork=1;
  int INFO=0, LWORK=10*nr;

  double *work = R_Vector(LWORK);
  double one = 1.0, zero = 0.0;

  F77_CALL(dgesvd)( "All", "All", &nr, &nr, matrix, &LDA, S, Up, &LDU, VT, &LDVT, work, &LWORK, &INFO);

//	RprintVecAsMat("S", S, 1, nr);
//	RprintVecAsMat("Up", Up, nr, nr);
//	RprintVecAsMat("VT", VT, nr, nr);



//	RprintVecAsMat("U", scratch2, nr, nr);

	for(i = 0; i < nr; i++){

		for(j = 0; j < nr; j++){
      svs[i*nr + j] = 0.0;
			if(i == j){
				if(inverse) svs[i*nr + j] = 1/sqrt(S[i]);
				if(!inverse) svs[i*nr + j] = sqrt(S[i]);
			}

		}
	}

//  RprintVecAsMat("svs", svs, nr, nr);
  F77_CALL(dgemm)("N", "N", &nr, &nr, &nr, &one, Up, &nr, svs, &nr, &zero, scr1, &nr);

//  RprintVecAsMat("UpS", scr1, nr, nr);

  F77_CALL(dgemm)("N", "N", &nr, &nr, &nr, &one, scr1, &nr, VT, &nr, &zero, matrix, &nr);

//  RprintVecAsMat("UpS^(-1/2)Vt", matrix, nr, nr);


}
