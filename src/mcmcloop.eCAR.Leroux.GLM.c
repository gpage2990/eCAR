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
* E = nobs x 1 vector containing number of trials in each county
* evals = nobs x 1 eigenvalues of R=M-W
* evecs = nobs x nobs matrix of eigenvectors of R=M-W
* W = nobs x nobs neighborhood matrix
* C = nobs x q design matrix associated with covariates that act as controls (i.e., causal inference is not of interest)
*
* modelPriors = vector containing prior values of model
* MHsd = vector containing candidate density sd values in MH algorithm for tau and sig2.
* verbose = integer determining whether to print to screen information regarding run.
*
* Output:
* beta0 = nout x 1 scratch vector that holds beta0 MCMC iterates
* beta = nout x 1 scratch vector that holds beta MCMC iterates
* alpha = nout x 1 scratch vector that holds alpha MCMC iterates
* tau = nout x 1 scratch vector that holds tau MCMC iterates
* sig2x = nout x 1 scratch vector that holds sig2x MCMC iterates
* lamx = nout x 1 scratch vector that holds lamx MCMC iterates
* lamz = nout x 1 scratch vector that holds lamz MCMC iterates
* theta = nout x nobs scratch vector that holds sig2 MCMC iterates
* eta = nout x ncov scratch vector that holds eta MCMC iterates
*****************************************************************************************/


void mcmcloop_leroux_GLM(int *draws, int *burn, int *thin, int *nobs, double *y, double *x, double *E,
			                double *evals,  double *evecs, double *W, double *C, int *ncov, int* modelnum,
			                double *modelPriors, double *MHsd, int *verbose, int *joint_prior_lamx_lamz,
			                int *updateXparms,
			                double *beta, double *alpha, double *tau, double *sig2x,
			                double *lamx, double *lamz, double *theta, double *beta0,
			                double *eta, double *nb_r){



	// i - MCMC iterate
	// ii - MCMC iterate that is saved
	// j - observation iterate
	// jj - second observation iterate
	// b - control covariate counter

	int i, j, jj, b, bb;
	int ii = 0;


	int nout = (*draws - *burn)/(*thin);

	if(*verbose){
		Rprintf("nobs = %d\n", *nobs);
		Rprintf("ncov = %d\n", *ncov);
		Rprintf("nout = %d\n", nout);
	}
  	int maxcov_2 = *ncov;
	if(*ncov < 2) maxcov_2 = 2;


  // Create a vector containing the number of neighbors for each county
	double *neighbor_vec = R_VectorInit(*nobs, 0.0);

  	for(j=0; j<*nobs; j++){
    	for(jj=0; jj<*nobs; jj++){
      		neighbor_vec[j] = neighbor_vec[j] + W[j*(*nobs) + jj];
    	}
  	}

  	// Create the x-star vector which requires making Gamma transpose
  	double *evecsT = R_Vector((*nobs)*(*nobs));
	double *xstar = R_VectorInit(*nobs, 0.0);
	double *CC = R_VectorInit((*ncov)*(*ncov), 0.0);
	double *CMC = R_VectorInit((*ncov)*(*ncov), 0.0);
	double *CWC = R_VectorInit((*ncov)*(*ncov), 0.0);
	double *Wmn = R_VectorInit((*ncov)*(*nobs), 0.0);

	mat_transpose(evecs, evecsT, *nobs, *nobs);
  	matrix_product(evecsT, x, xstar, *nobs, 1, *nobs);


	// ===================================================================================
	//
	// Variables to hold MCMC iterates for non cluster specific parameters
	//
	// ===================================================================================

	double sig2x_iter=sig2x[0], lamx_iter=lamx[0];
	double beta_iter=0.5, alpha_iter=0.5752237, beta0_iter=0.0;
  	double lamz_iter=0.96, tau_iter=0.4375;
	double *theta_iter = R_VectorInit(*nobs, 0.0);
	double *eta_iter = R_VectorInit(*ncov, 0.0);
	double *xi_iter = R_VectorInit(*nobs, 0.0);
	double *likelihood_iter = R_VectorInit(*nobs, 0.0);

	// ===================================================================================
	//
	// scratch vectors of memory needed to update parameters
	//
	// ===================================================================================

	// These are made particularly big to make sure there is enough memory
	double *scr1 = R_Vector(*nobs);
	double *scr2 = R_Vector(*nobs*(*ncov));

	// stuff needed to store likelihood values

	// stuff needed to update lam.z and lam.x and tau;
  	double llo, lln, llr, uu, ld, ldeto, ldetn, to, tn, qf1, qf2, qf3, qf4, qf5, xio, xin;
	double lamo, lamn, ssq, astar, bstar, Det, ssqn, ssqo, summn, mstar, s2star;

	// stuff I need to update beta and alpha simultaneously
	double *Mstar = R_Vector(maxcov_2);
	double *Mstar0 = R_Vector(maxcov_2);
	double *Sstar = R_Vector(maxcov_2*maxcov_2);
	double *SstarInv = R_Vector(maxcov_2*maxcov_2);
	double *outrmvnorm = R_Vector(maxcov_2);

  	//Stuff I need to update theta and lamz, lamx, sigma2x, sigma2z,
  	double *mnvec = R_Vector((*nobs));
  	double *mnveco = R_Vector((*nobs));
  	double *mnvecn = R_Vector((*nobs));
  	double *GAxstar = R_Vector((*nobs));
  	double *GAxstaro = R_Vector((*nobs));
  	double *GAxstarn = R_Vector((*nobs));
  	double *Sigma_inv = R_VectorInit((*nobs)*(*nobs),0.0);
  	double *Sigma_invx = R_VectorInit((*nobs)*(*nobs),0.0);
  	double *Sigma_invz = R_VectorInit((*nobs)*(*nobs),0.0);
  	double *Sigma_invo = R_VectorInit((*nobs)*(*nobs),0.0);
  	double *Sigma_invn = R_VectorInit((*nobs)*(*nobs),0.0);

  	//Stuff I need to update eta
  	double *Ce = R_VectorInit(*nobs,0.0);

  	if(*ncov > 0){
    	for(b = 0; b < *ncov; b++){
      		for(bb = 0; bb < *ncov; bb++){
        		for(j = 0; j < *nobs; j++){
          			CC[b*(*ncov) + bb] = CC[b*(*ncov) + bb] + C[j*(*ncov) + b]*C[j*(*ncov) + bb];
          			CMC[b*(*ncov) + bb] = CMC[b*(*ncov) + bb] + C[j*(*ncov) + b]*C[j*(*ncov) + bb]*neighbor_vec[j];
          			for(jj = 0; jj < *nobs; jj++){
           	 			CWC[b*(*ncov) + bb] = CWC[b*(*ncov) + bb] + C[j*(*ncov) + b]*C[jj*(*ncov) + bb]*W[j*(*nobs) + jj];
          			}
        		}
      		}
    	}
    	for(j=0; j<*nobs; j++){
      		for(b=0; b< *ncov; b++){
        		Ce[j] = Ce[j] + C[j*(*ncov) + b]*eta_iter[b];
      		}
    	}
  	}

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
	double mb0 = modelPriors[12], s2b0=modelPriors[13]; // mean and variance for beta0 ~ N(mb0,s2b0)
	double me = modelPriors[14], s2e=modelPriors[15]; // mean and variance for eta ~ N(me*j,s2e*I)
	double mx = modelPriors[16], s2x=modelPriors[17]; // mean and variance for xi ~ N(mx,s2x)

	if(*verbose){
    Rprintf("Prior values being used are:\n m=%.1f, s2=%.1f, alamx=%.1f, blamx=%.1f, \n alamz=%.1f, blamz=%.1f, asig=%.1f, bsig=%.1f,  atau=%.1f, btau=%.1f, \n asigx=%.1f, bsigx=%.1f, mb0=%.1f, s2b0=%.1f, me=%.1f, s2e=%.1f\n\n",
              m, s2, alamx,blamx,alamz,blamz,asig, bsig, atau, btau,asigx, bsigx, mb0, s2b0, me, s2e);
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

		if(*verbose & ((i+1) % 1000 == 0)){
//			time_t now;
//			time(&now);

//		  Rprintf("mcmc iter = %d ===================================================== \n", i+1);
//      Rprintf("%s", ctime(&now));

		}

		if(*verbose){
      		clock_t ith_iterate = clock();
			calc_time = (ith_iterate - begin)/CLOCKS_PER_SEC;

      		Rprintf("  Progress:%.1f%%, Time:%.1f seconds\r", ((double) (i+1) / (double) (*draws))*100.0, calc_time);
//      fflush(stdout);
		}




		//////////////////////////////////////////////////////////////////////////////////
		//
		// begin by updating sig2x and tau
		//
		//////////////////////////////////////////////////////////////////////////////////
		ssq = 0.0;
		for(j=0; j<*nobs; j++){
			GAxstar[j] = 0.0;
		  	for(jj=0; jj<*nobs; jj++){

				GAxstar[j] = GAxstar[j] +
                      		evecs[j*(*nobs)+jj]*
                      		sqrt((1.0 - lamx_iter + lamx_iter*evals[jj])/(1.0 - lamz_iter + lamz_iter*evals[jj]))*xstar[jj];

		    	if(W[j*(*nobs) + jj] == 1){
		      		Sigma_invx[j*(*nobs) + jj] = -lamx_iter;
		      		Sigma_invz[j*(*nobs) + jj] = -lamz_iter;
		    	}
        		if(j == jj){
          			Sigma_invx[j*(*nobs) + jj] = ((1-lamx_iter) + lamx_iter*neighbor_vec[j]);
          			Sigma_invz[j*(*nobs) + jj] = ((1-lamz_iter) + lamz_iter*neighbor_vec[j]);
        		}

        		ssq = ssq + x[j] *Sigma_invx[j*(*nobs) + jj]*x[jj];

		  	}
		  	mnvec[j] = theta_iter[j] - beta0_iter - beta_iter*x[j] - alpha_iter*GAxstar[j] - Ce[j];
		}

        if(*updateXparms){
    		astar = 0.5*(*nobs) + asigx;
			bstar = 0.5*ssq + bsigx;
			sig2x_iter = 1/rgamma(astar, 1/bstar);
		}


		ssq = quform(mnvec, Sigma_invz, *nobs);

    	astar = 0.5*(*nobs) + atau;
		bstar = 0.5*ssq + btau;
		tau_iter = 1/rgamma(astar, 1/bstar);




		//////////////////////////////////////////////////////////////////////////////////
		//
		// Next update lam.z.
		//
		//////////////////////////////////////////////////////////////////////////////////
		lamo = lamz_iter;
    	lamn = rnorm(lamo, 0.05);

		if((lamn > 0) & (lamn < 1)){

      		ldetn=0.0; ldeto=0.0;
			for(j=0; j<*nobs; j++){
		    	GAxstaro[j] = 0.0;
		    	GAxstarn[j] = 0.0;

		    	for(jj=0; jj<*nobs; jj++){
          			GAxstaro[j] = GAxstaro[j] +
                					evecs[j*(*nobs)+jj]*
                  					sqrt((1.0 - lamx_iter + lamx_iter*evals[jj])/(1.0 - lamo + lamo*evals[jj]))*xstar[jj];

          			GAxstarn[j] = GAxstarn[j] +
                					evecs[j*(*nobs)+jj]*
                  					sqrt((1.0 - lamx_iter + lamx_iter*evals[jj])/(1.0 - lamn + lamn*evals[jj]))*xstar[jj];


		      		if(W[j*(*nobs) + jj] == 1){
		        		Sigma_invo[j*(*nobs) + jj] = -lamo/tau_iter;
		        		Sigma_invn[j*(*nobs) + jj] = -lamn/tau_iter;
		      		}
          			if(j == jj){
            			Sigma_invo[j*(*nobs) + jj] = ((1-lamo) + lamo*neighbor_vec[j])/tau_iter;
            			Sigma_invn[j*(*nobs) + jj] = ((1-lamn) + lamn*neighbor_vec[j])/tau_iter;
          			}
		    	}

		    	mnveco[j] = theta_iter[j] - beta0_iter - beta_iter*x[j] - alpha_iter*GAxstaro[j] - Ce[j];
		    	mnvecn[j] = theta_iter[j] - beta0_iter - beta_iter*x[j] - alpha_iter*GAxstarn[j] - Ce[j];
        		ldeto = ldeto + log(1-lamo+lamo*evals[j]);
        		ldetn = ldetn + log(1-lamn+lamn*evals[j]);
		  	}


      		ssqo = quform(mnveco, Sigma_invo, *nobs);
      		ssqn = quform(mnvecn, Sigma_invn, *nobs);

      		if(*joint_prior_lamx_lamz){
        		llo = 0.5*ldeto + -0.5*ssqo + log(0.5); // Employs prior that forces lamz > lamx
        		if(lamx_iter > lamn) lln = 0.5*ldetn + -0.5*ssqn + log(0);
        		if(lamx_iter < lamn) lln = 0.5*ldetn + -0.5*ssqn + log(0.5);
      		} else {
        		llo = 0.5*ldeto + -0.5*ssqo + dbeta(lamo, alamz, blamz, 1);
        		lln = 0.5*ldetn + -0.5*ssqn + dbeta(lamn, alamz, blamz, 1);
      		}


      		llr = lln - llo;
      		uu = runif(0,1);

      		if(llr > log(uu)) lamz_iter = lamn;

    	}

//    	Rprintf("lamz_iter = %f\n", lamz_iter);


		//////////////////////////////////////////////////////////////////////////////////
		//
		// updating lam.x
		//
		//////////////////////////////////////////////////////////////////////////////////
		if(*updateXparms){
    		lamo = lamx_iter;
    		lamn = rnorm(lamo, 0.05);

    		if((lamn > 0) & (lamn < 1)){

      			ldetn=0.0; ldeto=0.0;
		  		for(j=0; j<*nobs; j++){

		    		GAxstaro[j] = 0.0;
		    		GAxstarn[j] = 0.0;

		    		for(jj=0; jj<*nobs; jj++){

          				GAxstaro[j] = GAxstaro[j] +
                						evecs[j*(*nobs)+jj]*
                  						sqrt((1.0 - lamo + lamo*evals[jj])/(1.0 - lamz_iter + lamz_iter*evals[jj]))*xstar[jj];

          				GAxstarn[j] = GAxstarn[j] +
                						evecs[j*(*nobs)+jj]*
                  						sqrt((1.0 - lamn + lamn*evals[jj])/(1.0 - lamz_iter + lamz_iter*evals[jj]))*xstar[jj];

		      			if(W[j*(*nobs) + jj] == 1){
		        			Sigma_inv[j*(*nobs) + jj] = -lamz_iter/tau_iter;
		        			Sigma_invo[j*(*nobs) + jj] = -lamo/sig2x_iter;
		        			Sigma_invn[j*(*nobs) + jj] = -lamn/sig2x_iter;
		      			}
          				if(j == jj){
            				Sigma_inv[j*(*nobs) + jj] = ((1-lamz_iter) + lamz_iter*neighbor_vec[j])/tau_iter;
            				Sigma_invo[j*(*nobs) + jj] = ((1-lamo) + lamo*neighbor_vec[j])/sig2x_iter;
            				Sigma_invn[j*(*nobs) + jj] = ((1-lamn) + lamn*neighbor_vec[j])/sig2x_iter;
          				}
		    		}

		    		mnveco[j] = theta_iter[j] - beta0_iter - beta_iter*x[j] - alpha_iter*GAxstaro[j] - Ce[j];
		    		mnvecn[j] = theta_iter[j] - beta0_iter - beta_iter*x[j] - alpha_iter*GAxstarn[j] - Ce[j];

        			ldeto = ldeto + log(1-lamo+lamo*evals[j]);
        			ldetn = ldetn + log(1-lamn+lamn*evals[j]);
 		  		}


      			qf1=0.0, qf2=0.0, qf3=0.0, qf4=0.0;
      			for(j=0; j<*nobs; j++){
        			for(jj=0; jj<*nobs; jj++){
          				qf1 = qf1 + mnveco[j]*mnveco[jj]*Sigma_inv[j*(*nobs)+jj];
          				qf2 = qf2 + mnvecn[j]*mnvecn[jj]*Sigma_inv[j*(*nobs)+jj];
          				qf3 = qf3 + x[j]*x[jj]*Sigma_invo[j*(*nobs)+jj];
          				qf4 = qf4 + x[j]*x[jj]*Sigma_invn[j*(*nobs)+jj];
        			}
      			}


      			if(*joint_prior_lamx_lamz){
        			llo = -0.5*qf1 + 0.5*ldeto + -0.5*qf3 + log(0.5); // Uniform prior over unit square where lamz > lamx
        			if(lamn > lamz_iter) lln = -0.5*qf2 + 0.5*ldetn + -0.5*qf4 + log(0);
        			if(lamn < lamz_iter) lln = -0.5*qf2 + 0.5*ldetn + -0.5*qf4 + log(0.5);
      			} else {
        			llo = -0.5*qf1 + 0.5*ldeto + -0.5*qf3 + dbeta(lamo, alamx, blamx, 1);
        			lln = -0.5*qf2 + 0.5*ldetn + -0.5*qf4 + dbeta(lamn, alamx, blamx, 1);
      			}

      			llr = lln - llo;
      			uu = runif(0,1);

      			if(llr > log(uu)) lamx_iter = lamn;

    		}
		}

//    	Rprintf("lamx_iter = %f\n", lamx_iter);




		//////////////////////////////////////////
		//									    //
		// udate beta and alpha simultaneously  //
		//									    //
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

		for(j=0; j<*nobs; j++){
			GAxstar[j] = 0.0;
			for(jj=0; jj<*nobs; jj++){
        		GAxstar[j] = GAxstar[j] +
                				evecs[j*(*nobs)+jj]*
                  				sqrt((1.0 - lamx_iter + lamx_iter*evals[jj])/(1.0 - lamz_iter + lamz_iter*evals[jj]))*xstar[jj];
		    	if(W[j*(*nobs) + jj] == 1){
		      		Sigma_inv[j*(*nobs) + jj] = -lamz_iter/tau_iter;
		    	}
        		if(j == jj){
          			Sigma_inv[j*(*nobs) + jj] = ((1-lamz_iter) + lamz_iter*neighbor_vec[j])/tau_iter;
        		}
		  	}
		}


    	qf1=0.0, qf2=0.0, qf3=0.0, qf4=0.0, qf5=0.0;
    	for(j=0; j<*nobs; j++){
      		for(jj=0; jj<*nobs; jj++){
        		qf1 = qf1 + x[j]*x[jj]*Sigma_inv[j*(*nobs)+jj];
        		qf2 = qf2 + x[j]*GAxstar[jj]*Sigma_inv[j*(*nobs)+jj];
        		qf3 = qf3 + GAxstar[j]*GAxstar[jj]*Sigma_inv[j*(*nobs)+jj];
        		qf4 = qf4 + x[j]*(theta_iter[jj] - beta0_iter - Ce[jj])*Sigma_inv[j*(*nobs)+jj];
        		qf5 = qf5 + GAxstar[j]*(theta_iter[jj] - beta0_iter - Ce[jj])*Sigma_inv[j*(*nobs)+jj];
      		}
    	}
    	Sstar[0] = qf1 + (1/s2); // Note I am adding prior
    	Sstar[1] = qf2; // since prior is diagonal I don't add anything here.
    	Sstar[2] = qf2;
    	Sstar[3] = qf3 + (1/s2); // Note that I am adding prior


    	scr1[0] = qf4 + (1/s2)*m; // I am adding the prior
    	scr1[1] = qf5 + (1/s2)*m; // I am adding the prior

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



		//////////////////////////////////////////////////////////////////////////////////
		//																				                                      //
		// Update the theta value for each county (or areal unit)	                      //
		//																				                                      //
		//////////////////////////////////////////////////////////////////////////////////

		for(j=0; j<*nobs; j++){
			mnvec[j] = theta_iter[j] - beta0_iter - beta_iter*x[j] - alpha_iter*GAxstar[j] - Ce[j];
		}

    	for(j=0; j < *nobs; j++){

      		to = theta_iter[j];
      		tn = rnorm(to, 1);


      		ssq=0.0;
      		for(jj=0; jj < *nobs; jj++){

        		if(jj != j){
          			ssq = ssq + (mnvec[jj])*Sigma_inv[j*(*nobs) + jj];
        		}

      		}

      		ssqo =  (to - beta0_iter - beta_iter*x[j] - alpha_iter*GAxstar[j] - Ce[j])*
             		(to - beta0_iter - beta_iter*x[j] - alpha_iter*GAxstar[j] - Ce[j])*Sigma_inv[j*(*nobs) + j] +
             		(to - beta0_iter - beta_iter*x[j] - alpha_iter*GAxstar[j] - Ce[j])*ssq;

      		ssqn =  (tn - beta0_iter - beta_iter*x[j] - alpha_iter*GAxstar[j] - Ce[j])*
             		(tn - beta0_iter - beta_iter*x[j] - alpha_iter*GAxstar[j] - Ce[j])*Sigma_inv[j*(*nobs) + j] +
             		(tn - beta0_iter - beta_iter*x[j] - alpha_iter*GAxstar[j] - Ce[j])*ssq;



      		if(*modelnum == 1){ // Poisson likelihood
        		llo = -E[j]*exp(to) + y[j]*to -0.5*ssqo;
        		lln = -E[j]*exp(tn) + y[j]*tn -0.5*ssqn;
      		}

      		if(*modelnum == 2){ // Binomial likelihood
        		llo = y[j]*to - E[j]*log(exp(to) + 1) - 0.5*ssqo;
        		lln = y[j]*tn - E[j]*log(exp(tn) + 1) - 0.5*ssqn;
      		}

      		if(*modelnum == 3){ // Negative Binomial likelihood
        		llo = y[j]*(log(E[j]) + to) - (y[j] + exp(xi_iter[j]))*log(E[j]*exp(to) + exp(xi_iter[j])) - 0.5*ssqo;
        		lln = y[j]*(log(E[j]) + tn) - (y[j] + exp(xi_iter[j]))*log(E[j]*exp(tn) + exp(xi_iter[j])) - 0.5*ssqn;
      		}



      		llr = lln - llo;
      		uu = runif(0,1);


      		if(llr > log(uu)) theta_iter[j] = tn;

      		mnvec[j] = theta_iter[j] - beta0_iter - beta_iter*x[j] - alpha_iter*GAxstar[j] - Ce[j];



 		  	////////////////////////////////////////////////////////////////////////
		  	//										                              //
		  	// For the NB model within in the same for-loop, update xi_iter where //
		  	// the number of success until the Y_i failures is r_i = exp(xi_i)    //
		  	//										                              //
		 	 ////////////////////////////////////////////////////////////////////////
      		if(*modelnum == 3){
  		  		xio = xi_iter[j];
	  	  		xin = rnorm(xio, 0.5);

  		  		llo = lgammafn(y[j] + exp(xio)) - lgammafn(exp(xio)) +
	  	        		xio*exp(xio) - (y[j] + exp(xio))*log(exp(xio) + E[j]*exp(theta_iter[j])) -
		          		0.5*(1/s2x)*(xio*xio - 2*xio*mx);

  		  		lln = lgammafn(y[j] + exp(xin)) - lgammafn(exp(xin)) +
	  	        		xin*exp(xin) - (y[j] + exp(xin))*log(exp(xin) + E[j]*exp(theta_iter[j])) -
		          		0.5*(1/s2x)*(xin*xin - 2*xin*mx);

  		  		llr = lln - llo;
	  	  		uu = runif(0,1);
        		if(llr > log(uu)) xi_iter[j] = xin;

      		}

    	}


		////////////////////////
		//					  //
		// Update beta0	      //
		//					  //
		////////////////////////
    	summn=0.0;
		for(j=0; j<*nobs;j++){
			summn = summn + theta_iter[j] - beta_iter*x[j] - alpha_iter*GAxstar[j] - Ce[j];
		}

    	s2star = 1/((*nobs/tau_iter)*(1-lamz_iter) + 1/s2b0);
    	mstar = s2star*(((1-lamz_iter)/tau_iter)*summn + mb0*(1/s2b0));

    	beta0_iter = rnorm(mstar, sqrt(s2star));


		////////////////////////
		//					  //
		// Update eta         //
		//					  //
		////////////////////////
    	if(*ncov>0){

  			for(b = 0; b < *ncov; b++){
	  	  		Mstar0[b] = 0.0;
		  		for(bb = 0; bb < *ncov; bb++){
			  		Sstar[b*(*ncov) + bb] = (1/tau_iter)*((1-lamz_iter)*CC[b*(*ncov)+bb] + lamz_iter*(CMC[b*(*ncov)+bb] - CWC[b*(*ncov)+bb]));

  					if(b==bb) Sstar[b*(*ncov) + bb] = Sstar[b*(*ncov) + bb] + 1/s2e;

	  			}

		    	for(j=0; j<*nobs; j++){

  	  	  			Wmn[j] = 0.0;
	  	    		for(jj=0;jj<*nobs;jj++){
		        		Wmn[j] = Wmn[j] + W[j*(*nobs) + jj]*(theta_iter[jj] - beta0_iter - beta_iter*x[jj] - alpha_iter*GAxstar[jj]);
		      		}
          			mnvec[j] = theta_iter[j] - beta0_iter - beta_iter*x[j] - alpha_iter*GAxstar[j];

          			Mstar0[b] =  Mstar0[b] + (1/tau_iter)*((1-lamz_iter)*C[j*(*ncov) + b]*mnvec[j] +
                                            lamz_iter*(C[j*(*ncov) + b]*mnvec[j]*neighbor_vec[j] -
                                                     C[j*(*ncov) + b]*Wmn[j]));

		    	}

		    	Mstar0[b] = Mstar0[b] + me*(1/s2e);

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
          			Ce[j] = Ce[j] + C[j*(*ncov) + b]*eta_iter[b];
        		}
      		}

    	}




    	/////////////////////////////////////////////////
    	//                                            //
		// Evaluate likelihood to compute DIC         //
		//										      //
		////////////////////////////////////////////////
    	for(j=0; j < *nobs; j++){
      		if(*modelnum == 1){ // Poisson likelihood
        		likelihood_iter[j] = dpois(y[j], E[j]*exp(theta_iter[j]), 1);
      		}

      		if(*modelnum == 2){ // Binomial likelihood
        		likelihood_iter[j] = dbinom(y[j], E[j],   exp(theta_iter[j]), 1);
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
      		sig2x[ii] = sig2x_iter;
      		beta[ii] = beta_iter;
      		beta0[ii] = beta0_iter;
      		alpha[ii] = alpha_iter;
      		for(j=0; j<*nobs; j++){
        		theta[ii*(*nobs) + j] = theta_iter[j];
        		if(*modelnum==3) nb_r[ii*(*nobs) + j] = exp(xi_iter[j]);

      		}
      		for(b=0; b<*ncov; b++){
        		eta[ii*(*ncov) + b] = eta_iter[b];
      		}

			ii = ii+1;
		}

	}

	PutRNGstate();




}


