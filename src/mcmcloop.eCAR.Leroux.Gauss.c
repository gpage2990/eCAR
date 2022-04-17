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
* evals = nobs x 1 eigenvalues of R=M-W
*
* modelPriors = vector containing prior values of model
* MHsd = vector containing candidate density sd values in MH algorithm for tau and sig2.
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
* eta = nout x ncov scratch vector that holds eta MCMC iterates
*****************************************************************************************/


void mcmcloop_leroux_gauss(int *draws, int *burn, int *thin, int *nobs, double *y, double *x,
			                double *evals,  double *Cstar, int *ncov, double *modelPriors,
			                double *MHsd, int *verbose, int* joint_prior_lamx_lamz, 
			                int *updateXparms, double *lamx_fix, double *sig2x_fix,
			                double *beta, double *alpha, double *tau, double *sig2x, double *lamx,
			                double *lamz, double *sig2, double *eta){



	// i - MCMC iterate
	// ii - MCMC iterate that is saved
	// j - observation iterate
	// jj - second observation iterate
	// k - cluster iterate
	// kk - knot iterate
	// b - control covariate counter

	int i, j, jj, b, bb;
	int ii = 0;


	int nout = (*draws - *burn)/(*thin);

	if(*verbose){
		Rprintf("nobs = %d, nout = %d\n", *nobs, nout);
		double ystar_mean = 0.0;
		double xstar_mean = 0.0;
		for(j=0; j<*nobs; j++){
		  ystar_mean = ystar_mean + y[j]/(*nobs);
		  xstar_mean = xstar_mean + x[j]/(*nobs);
		}
		Rprintf("ystar_mean = %f, xstar_mean = %f\n", ystar_mean, xstar_mean);
	}

  	int maxcov_2 = *ncov;
	if(*ncov < 2) maxcov_2 = 2;


	// ===================================================================================
	//
	// Variables to hold MCMC iterates for non cluster specific parameters
	//
	// ===================================================================================

	double sig2x_iter=*sig2x_fix, lamx_iter=*lamx_fix;
	double beta_iter=0, alpha_iter=1;
  	double lamz_iter=0.35, tau_iter=1.0;
  	double sig2_iter = 0.5;
	double *eta_iter = R_VectorInit(*ncov, 0.0);


	// ===================================================================================
	//
	// scratch vectors of memory needed to update parameters
	//
	// ===================================================================================

	// These are made particularly big to make sure there is enough memory
	double *scr1 = R_Vector(*nobs);
	double *scr2 = R_Vector(*nobs*(*ncov));

	// stuff needed to update lam.z and lam.x and tau;
  	double llo, lln, llr, uu, Ho, Hn, vo, vn, v, s2xo, s2xn, Hx, ld;
	double lamo, lamn, tauo, taun, s2o, s2n, ssq, astar, bstar, Det;

	// stuff I need to update beta and alpha simultaneously
	double *Mstar = R_Vector(maxcov_2);
	double *Mstar0 = R_Vector(maxcov_2);
	double *Sstar = R_Vector(maxcov_2*maxcov_2);
	double *SstarInv = R_Vector(maxcov_2*maxcov_2);
	double *outrmvnorm = R_Vector(maxcov_2);

	// stuff I need to update eta
  	double sum;
  	double *Ce = R_VectorInit(*nobs,0.0);
  	if(*ncov > 0){
    	for(j=0; j<*nobs; j++){
      		for(b=0; b< *ncov; b++){
        		Ce[j] = Ce[j] + Cstar[j*(*ncov) + b]*eta_iter[b];
      		}
    	}
  	}



	// Proposal distribution standard deviation for MH steps
	double tau_cand_sd = MHsd[0];
	double sig2_cand_sd = MHsd[1];
	// ===================================================================================
	//
	// Prior parameter values
	//
	// ===================================================================================

	// prior values
  	double mb=modelPriors[0], s2b=modelPriors[1]; // prior values for beta ~ N(mb, s2b)
	double alamx = modelPriors[2], blamx=modelPriors[3]; // prior values for for lamx ~ Beta(a,b)
	double alamz = modelPriors[4], blamz=modelPriors[5]; // prior values for for lamz ~ Beta(a,b)
	double asig = modelPriors[6], bsig=modelPriors[7]; // shape and scale for sig2 ~ Gamma(a,b)
	double atau = modelPriors[8], btau=modelPriors[9]; // shape and scale for tau ~ Gamma(a,b)
	double asigx = modelPriors[10], bsigx=modelPriors[11]; // shape and scale for sig2x ~ InvGamma(a,b)
	double me = modelPriors[14], s2e=modelPriors[15]; // mean and variance for eta ~ N(me*j,se*I)
    double ma = modelPriors[18], s2a = modelPriors[19]; // prior values for alpha ~ N(ma, s2a)

	if(*verbose){
    	Rprintf("Prior values are: mb=%0.1f, s2b=%0.1f, ma=%0.1f, s2a=%0.1f, alamx=%0.1f, blamx=%0.1f,\n alamz=%0.1f, blamz=%0.1f, asig=%0.1f, bsig=%0.1f, \n atau=%0.1f, btau=%0.1f, asigx=%0.1f, bsigx=%0.1f\n",
              mb, s2b, ma, s2a,  alamx,blamx,alamz,blamz,asig, bsig, atau, btau,asigx, bsigx);
//	  if(*joint_prior_lamx_lamz) Rprinf("Prior Pr(lam_z > lam_x) = 1 \n");
	}
	GetRNGstate();


	// ===================================================================================
	//
	// start of the mcmc algorithm;
	//
	// ===================================================================================
	double calc_time = 0.0;
	clock_t  begin = clock();

	for(i = 0; i < *draws; i++){

		if(*verbose){
			if((i+1) % 1000 == 0){
//				time_t now;
//				time(&now);

//				Rprintf("mcmc iter = %d ===================== \n", i+1);
//				Rprintf("%s", ctime(&now));

			}
		}
		if(*verbose){
      		clock_t ith_iterate = clock();
		  	calc_time = (ith_iterate - begin)/CLOCKS_PER_SEC;

      		Rprintf("  Progress:%.1f%%, Time:%.1f seconds\r", ((double) (i+1) / (double) (*draws))*100.0, calc_time);
//      	fflush(stdout);
		}


		//////////////////////////////////////////////////////////////////////////////////
		//
		// begin by updating lam.z.
		//
		//////////////////////////////////////////////////////////////////////////////////
		lamo = lamz_iter;
    	lamn = rnorm(lamo, 0.1);

    	if((lamn > 0) & (lamn < 1)){

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

      		if(*joint_prior_lamx_lamz){
        		llo = llo + log(0.5); // Employs prior that forces lamz > lamx
        		if(lamx_iter > lamn) lln = lln + log(0);
        		if(lamx_iter < lamn) lln = lln + log(0.5);
     	 	} else {
        		llo = llo + dbeta(lamo, alamz, blamz, 1);
        		lln = lln + dbeta(lamn, alamz, blamz, 1);
      		}


      		llr = lln - llo;
      		uu = runif(0,1);

            if(llr > log(uu)) lamz_iter = lamn;

    	}


		//////////////////////////////////////////////////////////////////////////////////
		//
		// updating lam.x
		//
		//////////////////////////////////////////////////////////////////////////////////
		if(*updateXparms){
    		lamo = lamx_iter;
    		lamn = rnorm(lamo, 0.1);
    		if((lamn > 0) & (lamn < 1)){

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

      			if(*joint_prior_lamx_lamz){
        			llo = llo + log(0.5); // Employs prior that forces lamz > lamx
        			if(lamn > lamz_iter) lln = lln + log(0);
        			if(lamn < lamz_iter) lln = lln + log(0.5);
      			} else {
        			llo = llo + dbeta(lamo, alamx, blamx, 1);
        			lln = lln + dbeta(lamn, alamx, blamx, 1);
      			}

      			llr = lln - llo;
      			uu = runif(0,1);

                if(llr > log(uu)) lamx_iter = lamn;
    		}
		}

		//////////////////////////////////////////////////////////////////////////////////
		//
		// updating tau
		//
		//////////////////////////////////////////////////////////////////////////////////
    	tauo = tau_iter;
    	taun = rnorm(tauo, tau_cand_sd);
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

		//////////////////////////////////////////////////////////////////////////////////
		//
		// updating sigma2
		//
		//////////////////////////////////////////////////////////////////////////////////
    	s2o = sig2_iter;
    	s2n = rnorm(s2o, sig2_cand_sd);
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

		//////////////////////////////////////////////////////////////////////////////////
		//
		// updating sig2x
		//
		//////////////////////////////////////////////////////////////////////////////////
		if(*updateXparms){
    		ssq = 0.0;
			for(j=0; j<*nobs; j++){
				ssq = ssq + (1 - lamx_iter + lamx_iter*evals[j])*x[j]*x[j];
			}

    		astar = 0.5*(*nobs) + asigx;
			bstar = 0.5*ssq + bsigx;
			sig2x_iter = 1/rgamma(astar, 1/bstar);
		}

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


		for(j = 0; j < *nobs; j++){

      		Sstar[0] = Sstar[0] + x[j]*1/(tau_iter/(1.0 - lamz_iter + lamz_iter*evals[j]) + sig2_iter)*x[j];


      		Sstar[1] = Sstar[1] + x[j]*1/(tau_iter/(1.0 - lamz_iter + lamz_iter*evals[j]) + sig2_iter)*
                            sqrt((1.0 - lamx_iter + lamx_iter*evals[j])/(1.0 - lamz_iter + lamz_iter*evals[j]))*x[j];

      		Sstar[2] = Sstar[1];

      		Sstar[3] = Sstar[3] + sqrt((1.0 - lamx_iter + lamx_iter*evals[j])/(1.0 - lamz_iter + lamz_iter*evals[j]))*x[j]*
		                          1/(tau_iter/(1.0 - lamz_iter + lamz_iter*evals[j]) + sig2_iter)*
		                          sqrt((1.0 - lamx_iter + lamx_iter*evals[j])/(1.0 - lamz_iter + lamz_iter*evals[j]))*x[j];



      		scr1[0] = scr1[0] + x[j]*1/(tau_iter/(1.0 - lamz_iter + lamz_iter*evals[j]) + sig2_iter)*y[j] ;

      		scr1[1] = scr1[1] + sqrt((1.0 - lamx_iter + lamx_iter*evals[j])/(1.0 - lamz_iter + lamz_iter*evals[j]))*x[j]*
                              1/(tau_iter/(1.0 - lamz_iter + lamz_iter*evals[j]) + sig2_iter)*y[j] ;

		}

		scr1[0] = scr1[0] + (1/s2b)*mb; scr1[1] = scr1[1] + (1/s2a)*ma; // add prior in

    	Sstar[0] = Sstar[0] + 1/s2b; Sstar[3] = Sstar[3] + 1/s2a; // add prior

    	// Get inverse of the 2 x 2 Sstar matrix
		Det = (Sstar[0]*Sstar[3] - Sstar[1]*Sstar[2]);

		SstarInv[0] = Sstar[3]/Det;
		SstarInv[1] = -Sstar[1]/Det;
		SstarInv[2] = -Sstar[2]/Det;
		SstarInv[3] = Sstar[0]/Det;


    	Mstar[0] = SstarInv[0]*scr1[0] + SstarInv[1]*scr1[1];
    	Mstar[1] = SstarInv[2]*scr1[0] + SstarInv[3]*scr1[1];

		cholesky(SstarInv, (2) , &ld);

		ran_mvnorm(Mstar, SstarInv, (2), scr1, outrmvnorm);

    	beta_iter = outrmvnorm[0];
        alpha_iter = outrmvnorm[1];


		////////////////////////
		//					  //
		// Update eta         //
		//					  //
		////////////////////////
    	if(*ncov>0){

  			for(b = 0; b < *ncov; b++){

		  		for(bb = 0; bb < *ncov; bb++){

          			sum=0.0;
		  	  		for(j=0;j<*nobs;j++){
            			vo = tau_iter/(1.0 - lamz_iter + lamz_iter*evals[j]) + sig2_iter;
		  	    		sum = sum + Cstar[j*(*ncov)+b]*Cstar[j*(*ncov)+bb]*(1.0/vo);
		  	  		}

			  		Sstar[b*(*ncov) + bb] = sum;

  					if(b==bb) Sstar[b*(*ncov) + bb] = Sstar[b*(*ncov) + bb] + 1/s2e;

	  			}

        		sum=0.0;
		    	for(j=0; j<*nobs; j++){
           			vo = tau_iter/(1.0 - lamz_iter + lamz_iter*evals[j]) + sig2_iter;
		  	   		sum = sum + Cstar[j*(*ncov)+b]*(1.0/vo)*(y[j] -
		                                            beta_iter*x[j] -
		                                            alpha_iter*sqrt((1.0 - lamx_iter + lamx_iter*evals[j])/(1.0 - lamz_iter + lamz_iter*evals[j]))*x[j]);
		    	}

		    	Mstar0[b] = sum + me*(1/s2e);

		  	}

		  	cholesky(Sstar, *ncov, &ld);
		  	inverse_from_cholesky(Sstar, scr1, scr2, *ncov); //Sstar is now an inverse;

		  	matrix_product(Sstar, Mstar0, Mstar, *ncov, 1, *ncov);

		  	cholesky(Sstar, *ncov , &ld);

		  	ran_mvnorm(Mstar, Sstar, *ncov, scr1, scr2);

		  	for(b = 0; b < *ncov; b++){
		    	eta_iter[b]=scr2[b];
		  	}

      		for(j=0; j<*nobs; j++){
        		Ce[j] = 0.0;
        		for(b=0; b< *ncov; b++){
          			Ce[j] = Ce[j] + Cstar[j*(*ncov) + b]*eta_iter[b];
        		}
      		}
    	}



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
      		for(b=0; b<*ncov; b++){
        		eta[ii*(*ncov) + b] = eta_iter[b];
      		}

			ii = ii+1;
		}

	}
	PutRNGstate();




}


