/*
 * fir.h
 *
 *  Created on: Aug 26, 2014
 *      Author: Sina
 */
/**************************************************************************
*  FIR.H - This function performs FIR filtering                           *
*                                                                         *
*          ntap-1                                                         *
*    y(n) = sum  wi * x(n-i)                                              *
*           i=0                                                           *
*                                                                         *
**************************************************************************/

#ifndef FIR_H_
#define FIR_H_

float  fir(float *x, int N, float *w, int ntap)

{
  float yn;                   // output of FIR filter
  int i;                      // index


  yn = (float) 0.0;           // y(n) = 0.
  for ( i=0; i<=ntap-1; ++i)
  {
  yn = yn + w[i] * x[N-1-i];  // FIR filtering of x(n) to get y(n)
  }
  return(yn);                 // return y(n) to main function
}


float  fir_realtime(float *pLeft, float *xLeft, float *w, int N)

{

  int i;
  float output, *p;

  output = 0;								// set up for LEFT channel
  p = pLeft;								// save current sample pointer

  if(++pLeft > &xLeft[N]) {					// update pointer, wrap if necessary
  		pLeft = xLeft;						// and store
  }

	for (i = 0; i <= N; i++) {				// do LEFT channel FIR
	        output += *p-- * w[i];  		// multiply and accumulate
	        if(p < &xLeft[0])       		// check for pointer wrap around
      	    p = &xLeft[N];
	}

  return(output);                 // return output to main function
}




#endif /* FIR_H_ */
