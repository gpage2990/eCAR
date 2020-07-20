/*************************************************************
 * Copyright (c) 2012 Garritt Leland Page
 *
 * C-code that fits the joint eCAR Leroux model
 *
 *************************************************************/

#include "matrix.h"
#include "Rutil.h"

#include <R_ext/Lapack.h>
#include <R.h>
#include <Rmath.h>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/*****************************************************************************************
* The following are the inputs of the function that are read from R
*
* draws = total number of MCMC draws
* burn = number of MCMC draws discarded as burn-in
* thin = indicates how much MCMC chain should be thinned
*
* nobs = vector whose entries indicate number of observations per subject
* y = nobs x 1 vector containing nobs transformed response
* x = nobs x 1 vector containing nobs transformed explanatory variable
* W = nobs x nobs neighbor matrix
* evals = nobs x 1 eigenvalues of R=M-W
*
* modelPriors = vector containing prior values of model
*
* verbose = integer determining whether to print to screen information regarding run.
*
* Output:
* beta = nout x 1 scratch vector that holds beta MCMC iterates
* alpha = nout x 1 scratch vector that holds alpha MCMC iterates
* tau = nout x 1 scratch vector that holds tau MCMC iterates
* sig2x = nout x 1 scratch vector that holds sig2x MCMC iterates
* lamx = nout x 1 scratch vector that holds lamx MCMC iterates
* lamz = nout x 1 scratch vector that holds lamz MCMC iterates
* sig2 = nout x 1 scratch vectro that holds sig2 MCMC iterates
*****************************************************************************************/


void mcmcloop_leroux(int *draws, int *burn, int *thin, int *nobs, double *y, double *x,
			                double *evals,  double *modelPriors, int *verbose,
			                double *beta, double *alpha, double *tau, double *sig2x, double *lamx,
			                double *lamz, double *sig2){



	// i - MCMC iterate
	// ii - MCMC iterate that is saved
	// j - observation iterate
	// jj - second observation iterate
	// k - cluster iterate
	// kk - knot iterate

	int i, j, jj;
	int ii = 0;


	int nout = (*draws - *burn)/(*thin);

	if(*verbose){
		Rprintf("nobs = %d\n", *nobs);
		Rprintf("nout = %d\n", nout);
	}




	// ===================================================================================
	//
	// Variables to hold MCMC iterates for non cluster specific parameters
	//
	// ===================================================================================

	double sig2x_iter=1, beta_iter=0, alpha_iter=1;
  double lamx_iter=0.5, lamz_iter=0.9, tau_iter=1.0;
  double sig2_iter = 0.5;
	// ===================================================================================
	//
	// scratch vectors of memory needed to update parameters
	//
	// ===================================================================================

	// These are made particularly big to make sure there is enough memory
	double *scr1 = R_Vector(*nobs);

	// stuff needed to update lam.z and lam.x and tau;
  double llo, lln, llr, uu, Ho, Hn, vo, vn, v, s2xo, s2xn, Hx, ld;
	double lamo, lamn, tauo, taun, s2o, s2n, ssq, astar, bstar, Det;

	// stuff I need to update beta and alpha simultaneously
	double *Mstar = R_Vector(2);
	double *Sstar = R_Vector(2*2);
	double *outrmvnorm = R_Vector(2);


	// ===================================================================================
	//
	// Prior parameter values
	//
	// ===================================================================================

	// prior values
  double m=modelPriors[0], s2=modelPriors[1]; // prior values for (alpha, beta) ~ N(m, s2I)
	double alamx = modelPriors[2], blamx=modelPriors[3]; // prior values for for lamx ~ Beta(a,b)
	double alamz = modelPriors[4], blamz=modelPriors[5]; // prior values for for lamz ~ Beta(a,b)
	double asig = modelPriors[6], bsig=modelPriors[7]; // shape and scale for sig2 ~ Gamma(a,b)
	double atau = modelPriors[8], btau=modelPriors[9]; // shape and scale for tau ~ Gamma(a,b)
	double asigx = modelPriors[10], bsigx=modelPriors[11]; // shape and scale for sig2x ~ InvGamma(a,b)

	if(*verbose){
    Rprintf("Prior values are: m=%f, s2=%f, alamx=%f, blamx=%f,\n alamz=%f, blamz=%f, asig=%f, bsig=%f, \n atau=%f, btau=%f, asigx=%f, bsigx=%f\n",
              m, s2, alamx,blamx,alamz,blamz,asig, bsig, atau, btau,asigx, bsigx);
	}
	GetRNGstate();


	// ===================================================================================
	//
	// start of the mcmc algorithm;
	//
	// ===================================================================================

	for(i = 0; i < *draws; i++){

		if(*verbose){
			if((i+1) % 10000 == 0){
				time_t now;
				time(&now);

				Rprintf("mcmc iter = %d ===================== \n", i+1);
				Rprintf("%s", ctime(&now));

			}
		}


		//////////////////////////////////////////////////////////////////////////////////
		//
		// begin by updating lam.z.
		//
		//////////////////////////////////////////////////////////////////////////////////
		lamo = lamz_iter;
    lamn = rnorm(lamo, 0.05);

//    Rprintf("lamo = %f\n", lamo);
//    Rprintf("lamn = %f\n", lamn);
    if(lamn > 0 & lamn < 1){

      // Calcuate inverse sqrt matrix
      llo = 0.0; lln = 0.0;
      for(j=0; j<*nobs; j++){

        vo = tau_iter/(1.0 - lamo + lamo*evals[j]) + sig2_iter;
        vn = tau_iter/(1.0 - lamn + lamn*evals[j]) + sig2_iter;

        Ho = sqrt((1.0 - lamx_iter + lamx_iter*evals[j])/(1.0 - lamo + lamo*evals[j]));
        Hn = sqrt((1.0 - lamx_iter + lamx_iter*evals[j])/(1.0 - lamn + lamn*evals[j]));

        llo = llo + dnorm(y[j], beta_iter*x[j] + alpha_iter*Ho*x[j], sqrt(vo), 1);
        lln = lln + dnorm(y[j], beta_iter*x[j] + alpha_iter*Hn*x[j], sqrt(vn), 1);
      }

      llo = llo + dbeta(lamo, alamz, blamz, 1);
      lln = lln + dbeta(lamn, alamz, blamz, 1);

//      Rprintf("llo = %f\n", llo);
//      Rprintf("lln = %f\n", lln);

      llr = lln - llo;
      uu = runif(0,1);

      if(llr > log(uu)) lamz_iter = lamn;

    }
//    lamz_iter=0.9;
//    Rprintf("lamz_iter = %f\n", lamz_iter);


		//////////////////////////////////////////////////////////////////////////////////
		//
		// updating lam.x
		//
		//////////////////////////////////////////////////////////////////////////////////
    lamo = lamx_iter;
    lamn = rnorm(lamo, 0.1);
    if(lamn > 0 & lamn < 1){

      llo = 0.0; lln = 0.0;
      for(j=0; j<*nobs; j++){

        v = tau_iter/(1.0 - lamz_iter + lamz_iter*evals[j]) + sig2_iter;

        Ho = sqrt((1.0 - lamo + lamo*evals[j])/(1.0 - lamz_iter + lamz_iter*evals[j]));
        Hn = sqrt((1.0 - lamn + lamn*evals[j])/(1.0 - lamz_iter + lamz_iter*evals[j]));

        s2xo = sig2x_iter/(1.0 - lamo + lamo*evals[j]);
        s2xn = sig2x_iter/(1.0 - lamn + lamn*evals[j]);

        llo = llo + dnorm(y[j], beta_iter*x[j] + alpha_iter*Ho*x[j], sqrt(v), 1) +
                    dnorm(x[j], 0, sqrt(s2xo), 1);

        lln = lln + dnorm(y[j], beta_iter*x[j] + alpha_iter*Hn*x[j], sqrt(v), 1) +
                    dnorm(x[j], 0, sqrt(s2xn), 1);

      }

      llo = llo + dbeta(lamo, alamx, blamx, 1);
      lln = lln + dbeta(lamn, alamx, blamx, 1);

      llr = lln - llo;
      uu = runif(0,1);

      if(llr > log(uu)) lamx_iter = lamn;

    }
//    lamx_iter=0.1;
//    Rprintf("lamx_iter = %f\n", lamx_iter);


		//////////////////////////////////////////////////////////////////////////////////
		//
		// updating tau
		//
		//////////////////////////////////////////////////////////////////////////////////
    tauo = tau_iter;
    taun = rnorm(tauo, 0.1);
    if(taun > 0){

      lln=0.0, llo=0.0;
      for(j=0; j < *nobs; j++){

        Hx = sqrt((1.0 - lamx_iter + lamx_iter*evals[j])/(1.0 - lamz_iter + lamz_iter*evals[j]))*x[j];

        vo = tauo/(1.0 - lamz_iter + lamz_iter*evals[j]) + sig2_iter;
        vn = taun/(1.0 - lamz_iter + lamz_iter*evals[j]) + sig2_iter;

        llo = llo + dnorm(y[j], beta_iter*x[j] + alpha_iter*Hx, sqrt(vo), 1);
        lln = lln + dnorm(y[j], beta_iter*x[j] + alpha_iter*Hx, sqrt(vn), 1);

      }

      llo = llo + dgamma(tauo, atau, btau, 1);
      lln = lln + dgamma(taun, atau, btau, 1);

      llr = lln - llo;
      uu = runif(0,1);

      if(llr > log(uu)) tau_iter = taun;

    }
//    tau_iter = 1.1475;
//    Rprintf("tau = %f\n", tau_iter);

		//////////////////////////////////////////////////////////////////////////////////
		//
		// updating sigma2
		//
		//////////////////////////////////////////////////////////////////////////////////
    s2o = sig2_iter;
    s2n = rnorm(s2o, 0.1);
    if(s2n > 0){

      lln=0.0, llo=0.0;
      for(j=0; j < *nobs; j++){

        Hx = sqrt((1.0 - lamx_iter + lamx_iter*evals[j])/(1.0 - lamz_iter + lamz_iter*evals[j]))*x[j];

        vo = tau_iter/(1.0 - lamz_iter + lamz_iter*evals[j]) + s2o;
        vn = tau_iter/(1.0 - lamz_iter + lamz_iter*evals[j]) + s2n;

        llo = llo + dnorm(y[j], beta_iter*x[j] + alpha_iter*Hx, sqrt(vo), 1);
        lln = lln + dnorm(y[j], beta_iter*x[j] + alpha_iter*Hx, sqrt(vn), 1);

      }

      llo = llo + dgamma(s2o, asig, bsig, 1);
      lln = lln + dgamma(s2n, asig, bsig, 1);

      llr = lln - llo;
      uu = runif(0,1);

      if(llr > log(uu)) sig2_iter = s2n;

    }
//    sig2_iter = 0.1;
//    Rprintf("sig2 = %f\n", sig2_iter);
		//////////////////////////////////////////////////////////////////////////////////
		//
		// updating sig2x
		//
		//////////////////////////////////////////////////////////////////////////////////
    ssq = 0.0;
		for(j=0; j<*nobs; j++){

		  ssq = ssq + (1 - lamx_iter + lamx_iter*evals[j])*x[j]*x[j];
		}

    astar = 0.5*(*nobs) + asigx;
		bstar = 0.5*ssq + bsigx;
		sig2x_iter = 1/rgamma(astar, 1/bstar);

//    Rprintf("sig2x = %f\n", sig2x_iter);


		//////////////////////////////////////////
		//									                    //
		// udate beta and alpha simultaneously  //
		//									                    //
		/////////////////////////////////////////
		// In an attempt to avoid matrix inversion and
		// multiplication, I am creating the vectors
		// X.star %*% Sigma.inv and
		// HX.star %*% Sigma.inv
    // Note that since Sigma.inv is a diagonal, then this
    // matrix multiplication results in a mutliplying each
    // entry of X.star or HX.star with the corresponding
    // diagonal value of Sigma.inv.  Then I can use
    // inner and cross products to get Sstar
    for(j = 0; j < 2; j++){
      scr1[j] = 0.0;
      for(jj = 0; jj < 2; jj++){
        Sstar[j*2 + jj] = 0.0;
      }
    }

//    RprintVecAsMat("Sstar", Sstar, 2, 2);


		for(j = 0; j < *nobs; j++){

      Sstar[0] = Sstar[0] + sqrt((1.0 - lamx_iter + lamx_iter*evals[j])/(1.0 - lamz_iter + lamz_iter*evals[j]))*x[j]*
		                          1/(tau_iter/(1.0 - lamz_iter + lamz_iter*evals[j]) + sig2_iter)*
		                        sqrt((1.0 - lamx_iter + lamx_iter*evals[j])/(1.0 - lamz_iter + lamz_iter*evals[j]))*x[j];


      Sstar[1] = Sstar[1] - x[j]*1/(tau_iter/(1.0 - lamz_iter + lamz_iter*evals[j]) + sig2_iter)*
                            sqrt((1.0 - lamx_iter + lamx_iter*evals[j])/(1.0 - lamz_iter + lamz_iter*evals[j]))*x[j];

      Sstar[2] = Sstar[1];


      Sstar[3] = Sstar[3] + x[j]*1/(tau_iter/(1.0 - lamz_iter + lamz_iter*evals[j]) + sig2_iter)*x[j];

      scr1[0] = scr1[0] + x[j]*1/(tau_iter/(1.0 - lamz_iter + lamz_iter*evals[j]) + sig2_iter)*y[j] ;

      scr1[1] = scr1[1] + sqrt((1.0 - lamx_iter + lamx_iter*evals[j])/(1.0 - lamz_iter + lamz_iter*evals[j]))*x[j]*
                              1/(tau_iter/(1.0 - lamz_iter + lamz_iter*evals[j]) + sig2_iter)*y[j] ;

		}

		scr1[0] = scr1[0] + (1/s2)*m; scr1[1] = scr1[1] + (1/s2)*m; // add prior in

    Sstar[0] = Sstar[0] + 1/s2; Sstar[3] = Sstar[3] + 1/s2; // add prior

//    RprintVecAsMat("Sstar", Sstar, 2, 2);

    // Get inverse of the 2 x 2 Sstar matrix
		Det = (Sstar[0]*Sstar[3] - Sstar[1]*Sstar[2]);

		Sstar[0] = Sstar[0]/Det;
		Sstar[1] = Sstar[1]/Det;
		Sstar[2] = Sstar[2]/Det;
		Sstar[3] = Sstar[3]/Det;


//    RprintVecAsMat("Sstar", Sstar, 2, 2);

    Mstar[0] = Sstar[0]*scr1[0] + Sstar[1]*scr1[1];
    Mstar[1] = Sstar[2]*scr1[0] + Sstar[3]*scr1[1];

//    RprintVecAsMat("Mstar", Mstar, 1, 2);


		cholesky(Sstar, (2) , &ld);

		ran_mvnorm(Mstar, Sstar, (2), scr1, outrmvnorm);


//    RprintVecAsMat("outrmvnorm", outrmvnorm, 1, 2);

    beta_iter = outrmvnorm[0];
    alpha_iter = outrmvnorm[1];


		//////////////////////////////////////////////////////////////////////////////////
		//																				                                      //
		// Save MCMC iterates															                              //
		//																				                                      //
		//////////////////////////////////////////////////////////////////////////////////
		if((i > (*burn-1)) & ((i+1) % *thin == 0)){

			lamz[ii] = lamz_iter;
			lamx[ii] = lamx_iter;
      tau[ii] = tau_iter;
      sig2[ii] = sig2_iter;
      sig2x[ii] = sig2x_iter;
      beta[ii] = beta_iter;
      alpha[ii] = alpha_iter;

			ii = ii+1;
		}

	}

	PutRNGstate();




}


