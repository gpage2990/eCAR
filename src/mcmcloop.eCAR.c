/*************************************************************
 * Copyright (c) 2012 Garritt Leland Page
 *
 * C-code that fits the joint eCAR model
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
* D = nobs x 1 vector containing the number of neighbors for each location
* evals = nobs x 1 eigenvalues of M^{-1/2} W M^{-1/2}
* P = nobs x nobs projection matrix Gamma'M^{1/2}
* Pinv = nobs x nobs inverse projection matrix Gamma'M^{1/2}
*
* modelPriors = vector containing prior values from the hierarchical model
*
* verbose = integer determining whether to print to screen information regarding run.
*
* Output:
* beta = nout x 1 scratch vector that holds beta MCMC iterates
* alpha = nout x 1 scratch vector that holds beta MCMC iterates
* tau = nout x 1 scratch vector that holds beta MCMC iterates
* sig2x = nout x 1 scratch vector that holds beta MCMC iterates
* rhox = nout x 1 scratch vector that holds beta MCMC iterates
* rhoz = nout x 1 scratch vector that holds beta MCMC iterates
*
*****************************************************************************************/


void mcmcloop(int *draws, int *burn, int *thin, int *nobs, double *y, double *x,
			          double *W, double *D, double *evals, double *P, double *Pinv,
			          double *modelPriors, int *verbose,
			          double *beta, double *alpha, double *tau, double *sig2x, double *rhox, double *rhoz){




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

//	double sig2x_iter=1, beta_iter=0, alpha_iter=1;
  double rhox_iter=0.9, rhoz_iter=0.5, tau_iter=1.0;

	// ===================================================================================
	//
	// scratch vectors of memory needed to update parameters
	//
	// ===================================================================================

	// These are made particularly big to make sure there is enough memory
	double *scr1 = R_Vector(((*nobs)*(*nobs)));
	double *scr2 = R_Vector(((*nobs)*(*nobs)));
	double *scr3 = R_Vector(((*nobs)*(*nobs)));
	double *scr4 = R_Vector(((*nobs)*(*nobs)));
	double *scr5 = R_Vector(((*nobs)*(*nobs)));

	// stuff needed to update rho.z;
	double rhoo, rhon;
	double *Rx_sqrt_inv = R_VectorInit((*nobs)*(*nobs),0.0);
	double *Rz_sqrto = R_VectorInit((*nobs)*(*nobs),0.0);
	double *Rz_sqrtn = R_VectorInit((*nobs)*(*nobs),0.0);
	double *Sigovec = R_VectorInit((*nobs),0.0);
	double *Signvec = R_VectorInit((*nobs),0.0);


	// ===================================================================================
	//
	// Prior parameter values
	//
	// ===================================================================================

	// prior values for sig2
//	double asig = modelPriors[0];
//	double bsig = modelPriors[1];


	// IG parameters for tau2
	double at = modelPriors[2]; double bt = modelPriors[3];



	if(*verbose){

		Rprintf("at = %f\n", at);
		Rprintf("bt = %f\n", bt);


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
//				time_t now;
//				time(&now);

//				Rprintf("mcmc iter = %d =========================================== \n", i+1);
//				Rprintf("%s", ctime(&now));

			}
		}


		//////////////////////////////////////////////////////////////////////////////////
		//
		// begin by updating rho.z.
		//
		//////////////////////////////////////////////////////////////////////////////////
		rhoo = rhoz_iter;
    rhon = rnorm(rhoo, 0.05);

//    Rprintf("rhoo = %f\n", rhoo);
//    Rprintf("rhon = %f\n", rhon);
    if(rhon > 0 & rhon < 1){

      // Calcuate inverse sqrt matrix
      for(j=0; j<*nobs; j++){

        Sigovec[j] = (1/tau_iter)*(1/(1-rhoo*evals[j]));
        Signvec[j] = (1/tau_iter)*(1/(1-rhon*evals[j]));

        for(jj=0; jj<*nobs; jj++){

          Rx_sqrt_inv[j*(*nobs)+jj] = -rhox_iter*W[j*(*nobs)+jj];
          Rz_sqrto[j*(*nobs)+jj] = -rhoo*W[j*(*nobs)+jj];
          Rz_sqrtn[j*(*nobs)+jj] = -rhon*W[j*(*nobs)+jj];

          if(j == jj){
            Rx_sqrt_inv[j*(*nobs)+jj] = D[j] - rhox_iter*W[j*(*nobs)+jj];
            Rz_sqrto[j*(*nobs)+jj] = D[j] - rhoo*W[j*(*nobs)+jj];
            Rz_sqrtn[j*(*nobs)+jj] = D[j] - rhon*W[j*(*nobs)+jj];
          }
        }
      }

//      RprintVecAsMat("M-rhoxW", scr1, *nobs, *nobs);
      sqrtmat(Rx_sqrt_inv, 1, *nobs,  scr1, scr2, scr3, scr4, scr5);

      sqrtmat(Rz_sqrto, 0, *nobs,  scr1, scr2, scr3, scr4, scr5);
      sqrtmat(Rz_sqrtn, 0, *nobs,  scr1, scr2, scr3, scr4, scr5);




//      RprintVecAsMat("(M-rhoxW)^{-1/2}", scr1, *nobs, *nobs);


    }




/*
			/////////////////////////////////////////
			//									   //
			// udate beta within the same loop.    //
			//									   //
			/////////////////////////////////////////
			for(t = 0; t < nobs[j]; t++){
				z_b0[t] = z_tmp[t];
			}




			mat_transpose(H, tH, nobs[j], (*nb));

			matrix_product(tH, H, HtH, (*nb), (*nb), nobs[j]);
			matrix_product(tH, z_tmp, Htz, (*nb), 1, nobs[j]);



			for(b = 0; b < (*nb); b++){
				for(bb = 0; bb < (*nb); bb++){

					if(*nsubject > 1){

						Sstar[b*(*nb)+bb] = (1/sig2_iter[j])*HtH[b*(*nb)+bb];

						if(b == bb){
							Sstar[b*(*nb)+bb] = (1/sig2_iter[j])*HtH[b*(*nb)+bb] +
											 (1/(lamh[gvec[j]-1]*lamh[gvec[j]-1]));
						}
					}

					if(*nsubject==1){

						Sstar[b*(*nb)+bb] = (1/sig2_iter[j])*HtH[b*(*nb)+bb] +
						                    (1/tau2h[0])*K[b*(*nb)+bb];
					}



				}

			}

			cholesky(Sstar, (*nb), &ld);
			inverse_from_cholesky(Sstar, scr1, scr2, (*nb));


			for(b = 0; b < (*nb); b++){
				if(*nsubject > 1){
					scr3[b] = (1/sig2_iter[j])*Htz[b] +
							 (1/(lamh[gvec[j]-1]*lamh[gvec[j]-1]))*
							 thetah[(gvec[j]-1)*((*nb)) + b];
				}

				if(*nsubject == 1){
					scr3[b] = (1/sig2_iter[j])*Htz[b];

				}

			}


			matrix_product(Sstar, scr3, Mstar, (*nb), 1, (*nb));


			cholesky(Sstar, (*nb) , &ld);

			ran_mvnorm(Mstar, Sstar, (*nb), scr1, outrmvnorm);





			for(b = 0; b < (*nb); b++){

				beta_iter[b*(*nsubject) + j] = outrmvnorm[b];
				btmp[b] = beta_iter[b*(*nsubject) + j];

			}


			matrix_product(H, btmp, Hb, nobs[j], 1, (*nb));



			/////////////////////////////////////////
			//									   //
			// udate sigma2 within the same loop.  //
			//									   //
			/////////////////////////////////////////


			sumz_Hb = 0.0;

			for(jj = 0; jj < nobs[j]; jj++){

				scr3[jj] = z_b0[jj] - Hb[jj];
				sumz_Hb = sumz_Hb + (z_tmp[jj] - Hb[jj]);
			}

			sumsq = inner_product(scr3, 1, scr3, 1, nobs[j]);

			astar = 0.5*nobs[j] + asig;
			bstar = 0.5*sumsq + 1/bsig;

			//bstar is rate and rgamma requires scale hence inverse

			sig2_iter[j] = 1/rgamma(astar, 1/bstar);

			for(t = 0; t < nobs[j]; t++){
				fprime_iter[j*max_nobs + t] = Hb[t];
			}


		}







		//////////////////////////////////////////////////////////////////////////////////
		//
		// Update lam (using MH-step) this is the standard deviation not the variance
		//
		//////////////////////////////////////////////////////////////////////////////////
		if(*nsubject > 1){

			for(k = 0; k < *ng; k++){


				olam = lamh[k];
				nlam = rnorm(olam,csigLAM);

				if((nlam > 0) & (nlam < A)){

					for(b = 0; b < (*nb); b++){

						thtmp[b] = thetah[k*(*nb) + b];

						for(bb = 0; bb < (*nb); bb++){

							oV[b*(*nb)+bb] = 0.0;
							nV[b*(*nb)+bb] = 0.0;

							if(b == bb){

								oV[b*(*nb)+bb] = 1/(olam*olam);
								nV[b*(*nb)+bb] = 1/(nlam*nlam);
							}
						}
					}

					ldo = 2.0*(*nb)*log(olam);
					ldn = 2.0*(*nb)*log(nlam);

					lln = 0.0;
					llo = 0.0;
					for(j = 0; j < *nsubject; j++){

						if(gvec[j] == k+1){

							for(b = 0; b < (*nb); b++){

								btmp[b] = beta_iter[b*(*nsubject) + j];

							}

							llo = llo + dmvnorm(btmp, thtmp, oV, (*nb), ldo, scr1, 1);
							lln = lln + dmvnorm(btmp, thtmp, nV, (*nb), ldn, scr1, 1);

						}

					}


					llo = llo + dunif(olam, 0.0, A, 1);
					lln = lln + dunif(nlam, 0.0, A, 1);


					llr = lln - llo;
					uu = runif(0.0,1.0);

					if(log(uu) < llr) lamh[k] = nlam;
				}





				//////////////////////////////////////////////////////////////////////////////
				//																			//
				// udpate thetah each of the cluster specific coefficients;					//
				//																			//
				//////////////////////////////////////////////////////////////////////////////

				for(b = 0; b < (*nb); b++){
					for(bb = 0; bb < (*nb); bb++){

				 		Sstar[b*(*nb)+bb] = (1/tau2h[k])*K[b*(*nb)+bb];

				 		if(b == bb){ Sstar[b*(*nb)+bb] = ((double) n_group[k]/(lamh[k]*lamh[k])) +
				 							  (1/tau2h[k])*K[b*(*nb)+bb];}

					}

					sumbeta[b] = 0.0;

				}

				cholesky(Sstar, (*nb), &ld);
				inverse_from_cholesky(Sstar, scr1, scr2, (*nb));

				for(j = 0; j < *nsubject; j++){

					if(gvec[j] == k+1){

						for(b = 0; b < (*nb); b++){

							sumbeta[b] = sumbeta[b] + (1/(lamh[k]*lamh[k]))*
														beta_iter[b*(*nsubject) + j];
						}
					}
				}


				matrix_product(Sstar, sumbeta, Mstar, (*nb), 1, (*nb));


				cholesky(Sstar, (*nb) , &ld);



				ran_mvnorm(Mstar, Sstar, (*nb), scr1, outrmvnorm);


				for(b = 0; b < (*nb); b++){

					thetah[k*(*nb) + b] = outrmvnorm[b];

				}


 				matrix_product(Htheta, outrmvnorm, Ht, *nrHt, 1, (*nb));


				for(t = 0; t < *nrHt; t++){
					fgprime_iter[k*(*nrHt) + t] = Ht[t];
				}

			}
		}

		//////////////////////////////////////////////////////////////////////////////
		//
		// Update tau2 for each of the clusters (P-spline smoothing parameter)
		//
		//////////////////////////////////////////////////////////////////////////////
		for(k = 0; k < *ng; k++){


			for(b = 0; b < (*nb); b++){

				if(*nsubject > 1){
					thtmp[b] = thetah[k*(*nb) + b];
				}
				if(*nsubject == 1){
					thtmp[b] = beta_iter[b*(*nsubject) + 0];
				}
			}


			sumsq = quform(thtmp,K,(*nb));


			astar = 0.5*((*nb)-2) + at;
			bstar = 1/bt + 0.5*sumsq;


			tau2h[k] = 1/rgamma(astar, 1/bstar);// E(tau2) = astarbstar for gamma.  bstar is scale




		}

*/

		//////////////////////////////////////////////////////////////////////////////////
		//																				//
		// Save MCMC iterates															//
		//																				//
		//////////////////////////////////////////////////////////////////////////////////
		if((i > (*burn-1)) & ((i+1) % *thin == 0)){

      rhox[ii] = rhox_iter;

			ii = ii+1;
		}

	}

	PutRNGstate();




}


