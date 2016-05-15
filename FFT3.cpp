#include <iostream>
#include <stdio.h>
#include <math.h>
#define PI 3.141592653589793238462643383279502884197169
using namespace std;

struct COMPLEX
{
	double re;
	double im;
};


void scramble(int numpoints, COMPLEX *f)
{
    int i, j, m;          /* loop variables */
    double temp;          /* temporary storage during swapping */

    j = 0;
    for (i=0;i<numpoints;i++)
	{
	if (i > j)
	    {                    /* swap f[j] and f[i] */
	    temp = f[j].re;      /* swap real */
	    f[j].re = f[i].re;
	    f[i].re = temp;
	    temp = f[j].im;      /* swap imaginary */
	    f[j].im = f[i].im;
	    f[i].im = temp;
	    }
	m = numpoints>>1;
	while ((j >= m) & (m >= 2))
	    {
	    j -= m;
	    m = m >> 1;
	    }
	j += m;
	}
}

void butterfly(int numpoints, int logN, int dir, COMPLEX *f)
{
    double angle;                  /* polar angle */
    COMPLEX w, wp, temp;           /* intermediate complex numbers */
    int i, j, k, offset;           /* loop variables */
    int N, half_N;                 /* storage for powers of 2 */
    double wtemp;                  /* temporary storage for w.re */

    N = 1;
    for(k=0; k<logN; k++)
	{
	half_N = N;
	N <<= 1;                       /* multiply N by 2 */
	angle = -2.0 * PI / N * dir;
	wtemp = sin(0.5 * angle);
	wp.re = -2.0 * wtemp * wtemp;
	wp.im = sin(angle);
	w.re = 1.0;
	w.im = 0.0;
	for(offset=0; offset<half_N; offset++)
	    {
	    for(i=offset; i<numpoints; i+=N)
		{
		j = i + half_N;
		temp.re = (w.re * f[j].re) - (w.im * f[j].im);
		temp.im = (w.im * f[j].re) + (w.re * f[j].im);
		f[j].re = f[i].re - temp.re;
		f[i].re += temp.re;
		f[j].im = f[i].im - temp.im;
		f[i].im += temp.im;
		}
	    wtemp = w.re;
	    w.re = wtemp * wp.re - w.im * wp.im + w.re;
	    w.im = w.im * wp.re + wtemp * wp.im + w.im;
	    }
	}
    if (dir == -1)                 /* if inverse fft, divide by numpoints */
	for (i=0;i<numpoints;i++)
	    {
	    f[i].re /= numpoints;
	    f[i].im /= numpoints;
	   }
}


void m_fft(COMPLEX *f, int numpoints, int dir)
{
   int logN;
   logN=0;
   while(numpoints != (1<<logN)) logN++;
   scramble(numpoints, f);
   butterfly(numpoints, logN, dir, f);
}



/*
int main()
{
	m_fft(ff,8,1);
	for(int i=0;i<8;i++)
	{
		cout<<ff[i].re<<"+"<<ff[i].im<<endl;
	}
	return 0;
}

*/

