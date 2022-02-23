/*
 * shift.h
 *
 *  Created on: Aug 24, 2014
 *      Author: Sina
 */

#ifndef SHIFT_H_
#define SHIFT_H_


/**************************************************************************
*  SHIFT.H - This function updates signal vector                          *
*                                                                         *
**************************************************************************/

void  shift(float *x, int N, float new)
{
  int i;                      // index
  for (i=0; i<=N-2; ++i)
  {
    x[i] = x[i+1];            // shift old data x(n-i), i = 1, 2, ... N-1
  }
  x[N-1] = (float) new;       // insert new data x(n)
}


#endif /* SHIFT_H_ */
