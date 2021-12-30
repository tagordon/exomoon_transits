#include <math.h>
#include <stdlib.h>
#include <stdio.h>

void input_coords(double t[], double *ms, double *t0p, double *ep, double *Pp, double *Op, double *wp, double *ip, double *mp, \
    double *t0m, double *em, double *Pm, double *Om, double *wm, double *im, double *mm, int *j, \
    double xp[], double yp[], double zp[], double xm[], double ym[], double zm[], double bp[], double bpm[], double theta[]);

void coords(double t[], double *ms, double *t0p, double *ep, double *Pp, double *Op, double *wp, double *ip, double *mp, \
    double *t0m, double *em, double *Pm, double *Om, double *wm, double *im, double *mm, int *j, \
    double xs[], double ys[], double zs[], double xp[], double yp[], double zp[], double xm[], double ym[], double zm[]);

void input_coords_(double *t, double ms, double t0p, double ep, double Pp, double Op, double wp, double ip, double mp, \
    double t0m, double em, double Pm, double Om, double wm, double im, double mm, int j, \
    double *xp, double *yp, double *zp, double *xm, double *ym, double *zm, double *bp, double *bpm, double *theta){
    input_coords(t, &ms, &t0p, &ep, &Pp, &Op, &wp, &ip, &mp, &t0m, &em, &Pm, &Om, &wm, &im, &mm, &j, xp, yp, zp, xm, ym, zm, bp, bpm, theta);
}

void coords_(double *t, double ms, double t0p, double ep, double Pp, double Op, double wp, double ip, double mp, \
    double t0m, double em, double Pm, double Om, double wm, double im, double mm, int j, \
    double *xs, double *ys, double *zs, double *xp, double *yp, double *zp, double *xm, double *ym, double *zm){
    coords(t, &ms, &t0p, &ep, &Pp, &Op, &wp, &ip, &mp, &t0m, &em, &Pm, &Om, &wm, &im, &mm, &j, xs, ys, zs, xp, yp, zp, xm, ym, zm);
}