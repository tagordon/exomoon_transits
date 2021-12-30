#include <math.h>
#include <stdlib.h>
#include <stdio.h>

void impacts(double t[], double *ms, double *t0p, double *ep, double *Pp, double *Op, double *wp, double *ip, double *mp, \
    double *t0m, double *em, double *Pm, double *Om, double *wm, double *im, double *mm, int *j, \
    double x[], double y[], double xbc[], double ybc[], double bp2[], double bm2[], double bpm2[]);

void system_impacts(double *t, double ms, double t0p, double ep, double Pp, double Op, double wp, double ip, double mp, \
    double t0m, double em, double Pm, double Om, double wm, double im, double mm, int j, \
    double *x, double *y, double *xbc, double*ybc, double *bp2, double *bm2, double *bpm2){
    impacts(t, &ms, &t0p, &ep, &Pp, &Op, &wp, &ip, &mp, &t0m, &em, &Pm, &Om, &wm, &im, &mm, &j, x, y, xbc, ybc, bp2, bm2, bpm2);
}