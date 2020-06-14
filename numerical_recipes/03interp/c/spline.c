/**
 * @file    spline.c
 * @brief   Cubic spline interpolation.
 * @author  Stanley A. Baronett <stanley.a.baronett@gmail.com>
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>

void spline(float x[], float y[], int n, float yp1, float ypn, float y2[])
/**
 * Given arrays x[0..(n-1)] and y[0..(n-1)] containing a tabulated function,
 * i.e., yᵢ = f(xᵢ), with x₀ < x₁ < ... < xₙ₋₁, and given values of yp1 and
 * ypn for the first derivative of the interpolating function at points
 * 1 and n, respectively, this routine returns an array y2[0..(n-1)] that
 * contains the second derivatives of the interpolating function at the
 * tabulated points xᵢ. If yp1 and/or ypn are equal to 1 × 10³⁰ or larger,
 * the routine is signaled to set the corresponding boundary condition for
 * a natural spline, with zero second derivative on that boundary.
 * 
 * —Numerican Recipes for C, 2nd Ed., §3.3, p. 115
 */
{
    int i, k;
    float p, qn, sig, un, u[n];

    if (yp1 > 0.99e30) // lower boundary condition is set to "natural"
        y2[0] = u[0] = 0.0;
    else {             // or else, have a specified first derivative
        y2[0] = -0.5;
        u[0] = (3. / (x[1] - x[0])) * ((y[1] - y[0]) / (x[1] - x[0]) - yp1);
    }
    for (i = 1; i < n-1; i++) {
        // the decomposition loop of the tridiagonal algorithm.
        // y2 and u are used for temporary storage of the decompsed factors.
        sig = (x[i] - x[i-1]) / (x[i+1] - x[i-1]);
        p = sig*y2[i-1] + 2.;
        y2[i] = (sig - 1.)/p;
        u[i] = (y[i+1]-y[i]) / (x[i+1]-x[i]) - (y[i]-y[i-1]) / (x[i]-x[i-1]);
        u[i] = (6.*u[i] / (x[i+1] - x[i-1]) - sig*u[i-1]) / p;
    }
    if (ypn > 0.99e30) // upper boundary condition is set to "natural"
        qn = un = 0.0;
    else {             // or else, have a specified first derivative
        qn = 0.5;
        un = (3. / (x[n-1]-x[n-2]))*(ypn - (y[n-1]-y[n-2]) / (x[n-1]-x[n-2]));
    }
    y2[n-1] = (un - qn*u[n-2]) / (qn*y2[n-2] + 1.);
    for (k = n-2; k >= 0; k--) // backsubstitution loop of tridiagonal alg.
        y2[k] = y2[k]*y2[k+1] + u[k];
}

float splint(float xa[], float ya[], float y2a[], int n, float x, int klo)
/**
 * Given the arrays xa[0..(n-1)] and ya[0..(n-1)], which tabulate a function
 * (with the xaᵢ's in order), and given the array y2a[0..(n-1)], which is the
 * output from spline above, and given a value of x, this routine returns a
 * cubic-spline interpolated value y.
 * 
 * —Numerican Recipes for C, 2nd Ed., §3.3, p. 116
 */
{
    float h, b, a;

    // since sequential calls are in increasing order,
    // find and update place for current and future calls
    while (xa[klo+1] < x) klo++;
    h = xa[klo+1] - xa[klo];
    if (h == 0.0) { // xa's must be distinct
        fprintf(stderr,"Cubic spline run-time error...\n");
        fprintf(stderr,"Bad xa input to routine splint\n");
        fprintf(stderr,"...now exiting to system...\n");
        exit(1);
    }
    a = (xa[klo+1] - x) / h;
    b = (x - xa[klo]) / h;
    // evaluate cubic spline
    return a*ya[klo]+b*ya[klo+1]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[klo+1])*(h*h)/6.;
}

void output(float x[], float y[], int n, int interp)
/* output xy values to readable text file */
{
    FILE* file;
    if (interp) {
        system("rm -f interp.txt");
        file = fopen("interp.txt","a");
    }
    else {
        system("rm -f dataset.txt");
        file = fopen("dataset.txt","a");
    }
    fprintf(file, "x\t\ty\n");
    for (int i = 0; i < n; i++)
    fprintf(file,"%e\t\t%e\n",x[i], y[i]);
    fclose(file);
}

void linspace(float lo, float hi, int n, float *arr)
{
    float step = (hi - lo) / n;

    for (int i = 0; i < n; i++) {
        arr[i] = lo;
        lo += step;
    }
}

int main()
{
    int nints = 100;
    float x[nints];
    float y_int[nints];
    int klo = 0;

    // DATASET #1
    // int n = 4;
    // float xa[] = {0, 1, 5, 10};
    // float ya[] = {10, 10, 5, 0};
    // float y2a[n];
    // spline(xa, ya, n, 0, 0, y2a);

    // DATASET #2
    int n = 6;
    float xa[] = {0, 2, 5, 6, 8, 10};
    float ya[] = {0, 1, 3, 6, 7, 10};
    float y2a[n];
    spline(xa, ya, n, 1e30, 1e30, y2a);

    linspace(0, 10, nints, x);
    for (int i = 0; i < nints; i++) {
        y_int[i] = splint(xa, ya, y2a, n, x[i], klo);
    }
    output(xa, ya, n, 0);
    output(x, y_int, nints, 1);

    return 0;
}
