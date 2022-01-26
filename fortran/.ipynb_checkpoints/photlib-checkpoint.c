#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define PI 3.14159265358979323846

void flux(double *c1, double *c2, double *rp, double *rm, double bp2[], double bm2[], double cth[], double sth[], double **lc, int *j);

int main(){
    return 0;
}

void LC(double c1, double c2, double rp, double rm, double *bp, double *bpm, double *cth, double*sth, double **lc, int j){
    flux(&c1, &c2, &rp, &rm, bp, bpm, cth, sth, lc, &j);
}