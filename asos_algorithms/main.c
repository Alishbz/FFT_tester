#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h> 
#include <time.h>


#define q	10		/* for 2^3 points */
#define N	(1<<q)		/* N-point FFT, iFFT */

typedef float real;
typedef struct { real Re; real Im; } complex;

#ifndef PI
# define PI	3.14159265358979323846264338327950288
#endif


/* Print a vector of complexes as ordered pairs. */
static void
print_vector(
    const char* title,
    complex* x,
    int n)
{
    int i;
    printf("%s (dim=%d):", title, n);
    for (i = 0; i < n; i++) printf(" %5.2f,%5.2f ", x[i].Re, x[i].Im);
    putchar('\n');
    return;
}
 
void
fft(complex* v, int n, complex* tmp)
{
    if (n > 1) {			/* otherwise, do nothing and return */
        int k, m;    complex z, w, * vo, * ve;
        ve = tmp; vo = tmp + n / 2;
        for (k = 0; k < n / 2; k++) {
            ve[k] = v[2 * k];
            vo[k] = v[2 * k + 1];
        }
        fft(ve, n / 2, v);		/* FFT on even-indexed elements of v[] */
        fft(vo, n / 2, v);		/* FFT on odd-indexed elements of v[] */
        for (m = 0; m < n / 2; m++) {
            w.Re = cos(2 * PI * m / (double)n);
            w.Im = -sin(2 * PI * m / (double)n);
            z.Re = w.Re * vo[m].Re - w.Im * vo[m].Im;	/* Re(w*vo[m]) */
            z.Im = w.Re * vo[m].Im + w.Im * vo[m].Re;	/* Im(w*vo[m]) */
            v[m].Re = ve[m].Re + z.Re;
            v[m].Im = ve[m].Im + z.Im;
            v[m + n / 2].Re = ve[m].Re - z.Re;
            v[m + n / 2].Im = ve[m].Im - z.Im;
        }
    }
    return;
}
 
void
ifft(complex* v, int n, complex* tmp)
{
    if (n > 1) {			/* otherwise, do nothing and return */
        int k, m;    complex z, w, * vo, * ve;
        ve = tmp; vo = tmp + n / 2;
        for (k = 0; k < n / 2; k++) {
            ve[k] = v[2 * k];
            vo[k] = v[2 * k + 1];
        }
        ifft(ve, n / 2, v);		/* FFT on even-indexed elements of v[] */
        ifft(vo, n / 2, v);		/* FFT on odd-indexed elements of v[] */
        for (m = 0; m < n / 2; m++) {
            w.Re = cos(2 * PI * m / (double)n);
            w.Im = sin(2 * PI * m / (double)n);
            z.Re = w.Re * vo[m].Re - w.Im * vo[m].Im;	/* Re(w*vo[m]) */
            z.Im = w.Re * vo[m].Im + w.Im * vo[m].Re;	/* Im(w*vo[m]) */
            v[m].Re = ve[m].Re + z.Re;
            v[m].Im = ve[m].Im + z.Im;
            v[m + n / 2].Re = ve[m].Re - z.Re;
            v[m + n / 2].Im = ve[m].Im - z.Im;
        }
    }
    return;
}


int
main(void)
{
    complex v[N], v1[N], scratch[N];
    int k;

    printf("const int test_size = %d;\n", N);
    printf(" double[] test_x = new double[] { \n");

    srand(time(0));
      
  
    #define SAMPLE_RATE 44100.0 //  (  44.1 kHz)
     
    double frequency_la = 440.0; //  "la" nota
    double frequency_fa = 349.23; // "fa" nota

    for (k = 0; k < N; k++) {

        double noise = (double)rand() / RAND_MAX / 2;

        v[k].Re = 0.5 * sin(2.0 * PI * frequency_la * k / SAMPLE_RATE)  ;
        v[k].Re += 0.5 * sin(2.0 * PI * frequency_fa * k / SAMPLE_RATE) + noise;

        //v[k].Re = 2 * cos(20 * PI * k / (double)N) + noise;
        //v[k].Re += 2 * sin(30) + noise;

        printf("%lf, " , v[k].Re); 
    }

    printf(" };\n\n");
 
    fft(v, N, scratch);
     
    printf("const int fft_arr_size = %d;\n", N);
    printf(" double[] fft_arr = new double[] { \n");
     
    for (k = 0; k < N; k++) { 
        printf("%lf, ", scratch[k].Re);
    }

    printf(" };\n\n");

    ifft(v, N, scratch);
     
    printf("const int ifft_arr_size = %d;\n", N);
    printf(" double[] ifft_arr = new double[] { \n");
     
    for (k = 0; k < N; k++) {
        printf("%lf, ", scratch[k].Re);
    }

    printf(" };\n\n");
     
    exit(EXIT_SUCCESS);
}