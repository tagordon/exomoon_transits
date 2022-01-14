#include <math.h>
#include <stdlib.h>
#include <stdio.h>

void kepler_solve(double M[], double *e, double cosf[], double sinf[], double f_e[], double f_M[], int *j);

void grad_impacts(double t[], double *ms, double *t0p, double *ep, double *Pp, double *Op, double *wp, double *ip, double *mp, \
    double *t0m, double *em, double *Pm, double *Om, double *wm, double *im, double *mm, int *j, \
    double bp[], double bpm[], double theta[], double dbp[], double dbpm[], double dtheta[]);

void impacts(double t[], double *ms, double *t0p, double *ep, double *Pp, double *Op, double *wp, double *ip, double *mp, \
    double *t0m, double *em, double *Pm, double *Om, double *wm, double *im, double *mm, int *j, \
    double bp[], double bpm[], double theta[]);

void coords(double t[], double *ms, double *t0p, double *ep, double *Pp, double *Op, double *wp, double *ip, double *mp, \
    double *t0m, double *em, double *Pm, double *Om, double *wm, double *im, double *mm, int *j, \
    double xp[], double yp[], double zp[], double xm[], double ym[], double zm[]);

void impacts_(double *t, double ms, double t0p, double ep, double Pp, double Op, double wp, double ip, double mp, \
    double t0m, double em, double Pm, double Om, double wm, double im, double mm, int j, \
    double *bp, double *bpm, double *theta){
    impacts(t, &ms, &t0p, &ep, &Pp, &Op, &wp, &ip, &mp, &t0m, &em, &Pm, &Om, &wm, &im, &mm, &j, bp, bpm, theta);
}

void grad_impacts_(double *t, double ms, double t0p, double ep, double Pp, double Op, double wp, double ip, double mp, \
    double t0m, double em, double Pm, double Om, double wm, double im, double mm, int j, double *bp, double *bpm, double *theta, \
    double *dbp, double *dbpm, double *dtheta){
    grad_impacts(t, &ms, &t0p, &ep, &Pp, &Op, &wp, &ip, &mp, &t0m, &em, &Pm, &Om, &wm, &im, &mm, &j, bp, bpm, theta, dbp, dbpm, dtheta);
}

void coords_(double *t, double ms, double t0p, double ep, double Pp, double Op, double wp, double ip, double mp, \
    double t0m, double em, double Pm, double Om, double wm, double im, double mm, int j, \
    double *xp, double *yp, double *zp, double *xm, double *ym, double *zm){
    coords(t, &ms, &t0p, &ep, &Pp, &Op, &wp, &ip, &mp, &t0m, &em, &Pm, &Om, &wm, &im, &mm, &j, xp, yp, zp, xm, ym, zm);
}

void kepler_solve_(double *M, double e, double *cosf, double *sinf, double *f_e, double *f_M, int j){
    kepler_solve(M, &e, cosf, sinf, f_e, f_M, &j);
}