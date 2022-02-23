/*
	Filtered-X ECAP-Like algorithm for Active Noise Control
	Single channel version
	Copyright (C) 2016-2018 Alejandro Rodriguez Silva. All rights reserved.
	Visit https://www.hindawi.com/journals/sv/2017/3864951/ for more details.
*/

#include "DSK6713_AIC23.h"
#include <math.h>
#include "shift.h"
#include "fir.h"
Uint32 fs = DSK6713_AIC23_FREQ_16KHZ;		// let's start with 8kHz sample frequency
#define LEFT 0
#define RIGHT 1
#define MAX	64							// maximum order of arrays
#define LS 	64							// order for S^(z)
#define eps  1e-20						// small value to avoid divisions between zeros
#define sf   1e-6;						// scaling factor (worked with 1e-3 in Comsol)
#define N	   20							// number of filter's coefficients
#define L	   5							// Projection order (we should try to use less than 10)
#define b    8							// number of bits to codify the error (worked with 2 in Comsol Simulations)

int aux_threshold[L];		  // auxiliary variable for setting the threshold

float w[MAX] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
float sh[MAX]= {0.000298444, 0.002760473, 0.005002137, 0.007336578, 0.00843871, 0.005184635, 0.001288979, 0.000524694, 0.003400356, 0.006647305, 0.007597509, 0.007114698, 0.00473923,
		0.002065019, 0.002522527, 0.003002069, 0.0033328, 0.004270365, 0.006414848, 0.007061463, 0.005742801, 0.002548069, 0.000160405, -0.000827542, 0.000204436, 0.00323687,
		0.007519222, 0.009524123, 0.006984652, 0.002856652, 0.001770516, 0.001485647, 0.003230608, 0.005948947, 0.008637768, 0.007863206, 0.006591664, 0.00107428, 0.000268603,
		0.002275554, 0.00572296, 0.002577792, 0.002433581, -0.002645024, 0.02212714, 0.098646, 0.1086554, 0.01307926, -0.09262124, -0.07705191, 1.85E-05, 0.03597173, 0.02012938,
		-0.01215009, -0.02718293, -0.02116016, -0.02309607, -0.0326416, -0.03286942, -0.01541614, 0.004441432, 0.006768666, -0.001856694,-0.01383966
};
float x[MAX] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};                 // signal vector of x(n)
float yp[MAX];                // signal vector of y'(n)
float d[MAX];                 // signal vector of d(n)
float xp[MAX];                // signal vector of x'(n)
float error = 0;			  			// error signal
float Xp[N][L] = {{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},
				  {0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}};			  	  // reuse matrix of x'(n)					-
float Xp_trans[L][N];		  		// transpose of the reuse matrix
float XX[L][L];               // (Xp_trans)(Xp)
float Xp_ec[N];			     	 		// (Xp)(ec)
float error_vector[L] = {0,0,0,0,0};		  // error vector
float coded_error[L];		  		// coded error
float up[N];				  				// (Xp)(error_vector)
float down[L];				  			// (Xp_trans)(Xp)(error_vector)

volatile union{unsigned int unit;short channel[2];} AIC23_data;

interrupt void c_int11(){
	int Ls = 64;                  // order of S^(z)
	int add = 0;									// auxiliary variable for setting the threshold
	int aux = 0;									// auxiliary value for numerator/denominator calculations
	int counter1 = 0;							// counter increments if the algorithm needs to be calculated
	int counter2 = 0;							// counter increments if the algorithm doesn't need to be calculated
	float mu;											// convergence factor for the ECAP-Like
	float em = 0.9;								// maximum probable error
	float res = em/(pow(2,b)-1);	// coder resolution
	float up_norm = 0;						// norm of numerator
	float down_norm = 0; 					// norm of denominator

	/*  Following variables are used in on-line ANC mode  */

	float xx = 0;                 // new data of x(n)
	float xpn = 0;                // new data of x'(n)

	DSK6713_LED_init();
	DSK6713_LED_on(0);											// turn on led #0
	AIC23_data.unit = input_sample();				// input 32-bit from both channel
	xx = (AIC23_data.channel[LEFT]);				// input left channel

	shift(x, MAX, xx);        							// update signal vector of x(n)
	yn = fir(x, MAX, w, N);  								// FIR filtering of x(n) by W(z) to get y(n)

	output_sample(yn);						// send output to the speaker

	error = (AIC23_data.channel[RIGHT]);		// input right channel (error)
	error *= 5;										// gain for the error

	/*
		Update the error vector
	*/
	for (i = L-1; i > 0; i--){
		error_vector[i] = error_vector[i-1];
	}
	error_vector[0] = error;

	/*
		Coded error
	*/
	for (i = 0; i < L; i++){
		coded_error[i] = roundf(error_vector[i]/res);
	}

	/*
		Threshold rule
	*/
	for (i = 0; i < L; i++){
		if (coded_error[i] == 0 || coded_error[i] == 0){
			aux_threshold[i] = 0;
			}
		else{
			aux_threshold[i] = 1;
			}
			add = add + aux_threshold[i];
	}

	xpn = fir(x, MAX, sh, Ls); 	// FIR filtering of x(n) by S^(z) to get x'(n)
	shift(xp, MAX, xpn);      	// update signal vector of x'(n)

	/*
		Fill up the reuse matrix
	*/
	for (i = N-1; i >= 0; i--){
		for (j = L-1; j > 0; j--){
				Xp[i][j] = Xp[i][j-1];
				if (i==0){
					Xp[j][i]=Xp[j-1][i];
				}
			}
	}
	Xp[0][0] = xpn; 				// set the first row and first column of the matrix to x'(n)


	/*
		Threshold
	*/
	if (add == L){

	/*
		Transpose the matrix
	*/
	for (i = 0;i < L; i++){
		for (j =0; j < N; j++){
			Xp_trans[i][j] = Xp[j][i];
		}
	}
	//

	/*
		(X')(X)
	*/
	for (i = 0; i < L; i++){
		for (j = 0; j < L; j++){
			XX[i][j] = 0;
			for (k = 0; k < N; k++){
				XX[i][j]=(XX[i][j]+(Xp_trans[i][k] * Xp[k][j]));
			}
		}
	}

	/*
		Compute the numerator
	*/
	for (i = 0; i < N; i++){
		up[i] = 0;
		aux = 0;
		for (j = 0; j < L; j++){
		aux = (Xp[i][j] * error_vector[j]) + aux;
		up[i] = aux;
		}
		up_norm = up_norm + (up[i] * up[i]);
	}
	up_norm = sqrt(up_norm);

	/*
		Compute the denominator
	*/
	for (i = 0; i < L; i++){
		down[i] = 0;
		aux = 0;
	for (j = 0; j < L; j++){
		aux = (XX[i][j]*error_vector[j]) + aux;
		down[i] = aux;
		}
		down_norm = down_norm + (down[i] * down[i]);
	}
	down_norm = sqrt(down_norm);

	/*
		Compute the adaptive convergence factor
	*/
	mu = ((up_norm * up_norm) / (eps + (down_norm *down_norm))) * sf * res;

	/*
		(Xp)(ec)
	*/
	for (i = 0; i < N; i++){
		Xp_ec[i] = 0;
		aux = 0;
		for (j = 0; j < L; j++){
		aux = (Xp[i][j]*coded_error[j]) + aux;
		Xp_ec[i] = aux;
		}
	}

	/*
		update coefficients of W(z)
		using the FXLMS algorithm
	*/
	for (i = 0; i < N; i++) {
		w[i] = w[i] + (mu * Xp_ec[i]);
	}
	counter1++;
	} // threshold
	counter2++;

	DSK6713_LED_off(0);								 // turn off led #0
	return;
}
