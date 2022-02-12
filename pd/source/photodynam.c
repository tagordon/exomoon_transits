#include <stdio.h> 
#include <string.h> 
#include <memory.h> 
#include <stdlib.h> 
#include <math.h>

#include "icirc.h"
#include "scpolyint.h"
#include "elliptic.h"

double mttr_flux_general(circle *circles,int ncircle,double c0,double c1,double c2);

int main(){
    return 0;
}

void flux(double *f, double *xp, double *yp, double *xm, double *ym, double rp, double rm, double u1, double u2, int j){
    
    int i;
    circle c[3];
    double c1, c2, c3;
    
    c[0].x0 = 0.0;
    c[0].y0 = 0.0;
    c[0].r = 1.0;
    
    c[1].r = rp;
    c[2].r = rm;
    
    c1 = 1.0 - u1 - 2 * u2;
    c2 = u1 + 2 * u2;
    c3 = u2;
    
    for (i = 0; i < j; i++){
        c[1].x0 = xp[i];
        c[1].y0 = yp[i];
        c[2].x0 = xm[i];
        c[2].y0 = ym[i];
        f[i] = mttr_flux_general(c, 3, c1, c2, c3);
    }
}