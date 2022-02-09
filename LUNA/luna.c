#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

#define PI 3.14159265358979323846
#define G 8.887641315935487e-10


// compute the mean intensity in the annulus overlapped by the occulting body 
// according to the small-body approximation from Mandel & Agol 
double intensity(double z, double p, double c1, double c2, double c3, double c4){
    
    if (p == 0){
        return 0.0;
    }
    else if (z > 1 + p){
        return 0.0;
    }
    else{
        // partial occultation 
        if ((z < 1 + p) & (z > 1 - p)){
            double alpha = - pow(p, 2.) - 2 * p * z - pow(z, 2.) + 1;
            double beta = alpha + 4 * p * z;

            double I = c1 * ((4./5.) * pow(beta, 5./4.) - beta)
                + c2 * ((2./3.) * pow(beta, 3./2.) - beta)
                + c3 * ((4./7.) * pow(beta, 7./4.) - beta)
                + c4 * (pow(z - p, 4.) - 1) / 2. + beta;
            
            double denom = 1 - pow(z - p, 2.);
            return I / denom;
        }
        // full occultation 
        else {
            
            double alpha = - pow(p, 2.) - 2 * p * z - pow(z, 2.) + 1;
            double beta = alpha + 4 * p * z;
            double I = (4./5.) * c1 * (pow(beta, 5./4.) - pow(alpha, 5./4.))
                + (2./3.) * c2 * (pow(beta, 3./2.) - pow(alpha, 3./2.))
                + (4./7.) * c3 * (pow(beta, 7./4.) - pow(alpha, (7./4.))) 
                + 4 * p * z * (1 - c1 - c2 - c3 - c4 * (pow(p, 2.) + pow(z, 2.)));
            
            double denom = 4. * z * p;
            return I / denom;
        }
    }
}

// overlap of two circles 
double overlap(double bigr, double r, double d, bool computebigr){
        
    double alpha = (pow(d, 2.) + pow(r, 2.) - pow(bigr, 2.)) / (2. * d * r);
    double beta = (pow(d, 2.) + pow(bigr, 2.) - pow(r, 2.)) / (2. * d * bigr);
            
    double k0 = acos(alpha);
    double k1 = acos(beta);
    double k2 = sqrt(pow(2. * d * bigr, 2.) - pow(pow(d, 2.) + pow(bigr, 2.) - pow(r, 2.), 2.)) / 2.;
    return pow(r, 2.) * k0 + pow(bigr, 2.) * k1 - k2;
}

// compute the area occulted by the planet + moon according to Kipping, 2012 
// //casenums are for debugging and do not correspond to the case labels from Kipping 
void area(double * a, double zp, double zm, double zpm, double pp, double pm){
    
    double ap, am;
    double sumfact;
    double c14_1a_overlap;
    double c14_1b_overlap;
    double x12, y12, x13, y13, x23, y23;
    double x13_p, y13_p, x23_pp, y23_pp;
    double cos_theta_pp, sin_theta_pp;
    double d1, d2, d3;
    double sqrtarg;
    double cos_theta, sin_theta;
    double alpha, alpha1, alpha2;
    
    double casenum = 0;
    
    if (zp >= 1 + pp){
        if (zm >= 1 + pm){
            ap = 0;
            am = 0;
        }
        else if ((1 - pm < zm) & (zm < 1 + pm)){
            alpha = overlap(1, pm, zm, false);
            ap = 0;
            am = alpha;
        }
        else{
            ap = 0;
            am = pow(pm, 2) * PI;
        }
    }
    
    else if (((1 - pp) < zp) & (zp < (1 + pp))){
        casenum = 1;
        if (zm >= 1 + pm){
            casenum = 2;
            alpha = overlap(1, pp, zp, false);
            ap = alpha;
            am = 0;
        }

        else if (((1 - pm) < zm) & (zm < (1 + pm))){
            casenum = 3;
            if (zpm >= (pm + pp)){
                casenum = 50;
                alpha1 = overlap(1, pp, zp, false);
                alpha2 = overlap(1, pm, zm, false);
                
                ap = alpha1;
                am = alpha2;
            }
            else if (((pp - pm) < zpm) & (zpm < (pp + pm))){
                casenum = 4;
                alpha = overlap(1., pp, zp, false);
                x12 = (1. - pow(pp, 2.) + pow(zp, 2.)) / (2 * zp);
                y12 = sqrt(2. * pow(zp, 2.) * (1. + pow(pp, 2.)) - 
                                  pow((1. - pow(pp, 2.)), 2.) - pow(zp, 4.)) / (2. * zp);
                cos_theta = (pow(zp, 2.) + pow(zm, 2.) - pow(zpm, 2.)) / (2. * zp * zm);
                if(cos_theta > 1){
                    sin_theta = 0.0;
                }
                else{
                    sin_theta = sqrt(1. - pow(cos_theta, 2.));
                }
                if ((pow((x12 - zm * cos_theta), 2.) + pow((y12 - zm * sin_theta), 2.)) < pow(pm, 2.)){
                    if ((pow((x12 - zm*cos_theta), 2.) + pow((y12 + zm * sin_theta), 2.)) < pow(pm, 2.)){
                        if (zp > 1){
                            casenum = 5;
                            ap = alpha;
                            
                            alpha1 = overlap(1, pm, zm, false);
                            alpha2 = overlap(1, pp, zp, false);
                            am = (alpha1 - alpha2);
                        }
                        else{
                            casenum = 6;
                            alpha1 = overlap(1, pm, zm, false);
                            alpha2 = overlap(pp, pm, zpm, true);
                            
                            ap = alpha;
                            am = pow(pp, 2.) * PI + alpha1 - ap - alpha2;
                        }
                    }
                    else{
                        casenum = 7;
                        x13_p = (1. - pow(pm, 2.) + pow(zm, 2.)) / (2. * zm);
                        y13_p = - sqrt(2. * pow(zm, 2.) * (1. + pow(pm, 2.)) - 
                                              pow(1. - pow(pm, 2.), 2.) - pow(zm, 4.)) / (2. * zm);
                        x13 = x13_p * cos_theta - y13_p * sin_theta;
                        y13 = x13_p * sin_theta + y13_p * cos_theta;
                        x23_pp = (pow(pp, 2.) - pow(pm, 2.) + pow(zpm, 2.)) / (2. * zpm);
                        y23_pp = sqrt(2. * pow(zpm, 2.) * (pow(pp, 2.) + pow(pm, 2.)) - 
                                             pow(pow(pp, 2.) - pow(pm, 2.), 2.) - pow(zpm, 4.)) / (2. * zpm);
                        cos_theta_pp = - (pow(zp, 2.) + pow(zpm, 2.) - pow(zm, 2.)) / (2. * zp * zpm);
                        if (cos_theta_pp > 1){
                            sin_theta_pp = 0.0;
                        }
                        else{
                            sin_theta_pp = sqrt(1. - pow(cos_theta_pp, 2.));
                        }
                        x23 = x23_pp * cos_theta_pp - y23_pp * sin_theta_pp + zp;
                        y23 = x23_pp * sin_theta_pp + y23_pp * cos_theta_pp;
                        d1 = sqrt(pow((x12 - x13), 2.) + pow((y12 - y13), 2.));
                        d2 = sqrt(pow((x12 - x23), 2.) + pow((y12 - y23), 2.));
                        d3 = sqrt(pow((x13 - x23), 2.) + pow((y13 - y23), 2.));
                        if ((zm * sin_theta) > (y13 + (zm * cos_theta - x13) * (y23 - y13) / (x23 - x13))){
                            casenum = 8;
                            sumfact = 0;
                            sumfact += asin(d1 / 2.) - sqrt(4 - pow(d1, 2.)) * d1 / 4.;
                            sumfact += pow(pp, 2) * asin(d2 / (2. * pp)) - 
                                sqrt(4 * pow(pp, 2) - pow(d2, 2.)) * d2 / 4.;
                            sumfact += pow(pm, 2) * asin(d3 / (2. * pm)) - 
                                sqrt(4 * pow(pm, 2) - pow(d3, 2.)) * d3 / 4.;
                            sqrtarg = (d1 + d2 + d3) * 
                                                      (d2 + d3 - d1) * 
                                                      (d1 + d3 - d2) * 
                                                      (d1 + d2 - d3);
                            c14_1a_overlap = (sqrt(sqrtarg) / 4. + sumfact);
                            
                            alpha1 = overlap(1, pm, zm, false);
                            
                            ap = alpha;
                            am = (alpha1 - c14_1a_overlap);
                        }
                        else{
                            casenum = 9;
                            sumfact = asin(d1 / 2.);
                            sumfact += pow(pp, 2.) * asin(d2 / (2. * pp));
                            sumfact -= pow(pm, 2.) * asin(d3 / (2. * pm));
                            c14_1b_overlap = sqrt((d1 + d2 + d3) * 
                                                         (d2 + d3 - d1) * 
                                                         (d1 + d3 - d2) * 
                                                         (d1 + d2 - d3)) / 4. + sumfact - 
                                d1 * sqrt(4. - pow(d1, 2.)) / 4. - 
                                d2 * sqrt(4. * pow(pp, 2.) - pow(d2, 2.)) / 4. + 
                                d3 * sqrt(4. * pow(pm, 2.) - pow(d3, 2.)) / 4. + PI * pow(pm, 2.);
                            alpha1 = overlap(1, pm, zm, false);
                            
                            ap = alpha;
                            am = (alpha1 - c14_1b_overlap);
                        }
                    }
                }
                else{
                    casenum = 10;
                    x13_p = (1 - pow(pm, 2) + pow(zm, 2)) / (2 * zm);
                    y13_p = - sqrt(2 * pow(zm, 2) * (1 + pow(pm, 2)) - 
                                          pow((1 - pow(pm, 2)), 2) - pow(zm, 4)) / (2 * zm);
                    x13 = x13_p * cos_theta - y13_p * sin_theta;
                    y13 = x13_p * sin_theta + y13_p * cos_theta;
                    if ((pow((x13 - zp), 2) + pow(y13, 2)) < (pow(pp, 2))){
                        // this condition needs fixed at some point. 
                        // I think it should be "is either of the planet-moon intersections
                        // inside of the star? If not, only count planet-star overlap."
                        if ((zpm > (pp - pm)) & (zm < zp)){
                            alpha1 = overlap(pp, pm, zpm, true);
                            
                            ap = alpha;
                            am = pow(pm, 2) * PI - alpha1;
                        }
                        else{
                            casenum = 11;
                            ap = alpha;
                            am = 0;
                        }
                    }
                    else{
                        casenum = 12;
                        x23_pp = (pow(pp, 2) - pow(pm, 2) + pow(zpm, 2)) / (2 * zpm);
                        y23_pp = sqrt(2 * pow(zpm, 2) * (pow(pp, 2) + pow(pm, 2)) - 
                                         pow((pow(pp, 2) - pow(pm, 2)), 2) - pow(zpm, 4)) / (2 * zpm);
                        cos_theta_pp = - (pow(zp, 2) + pow(zpm, 2) - pow(zm, 2)) / (2 * zp * zpm);
                        sin_theta_pp = sqrt(1 - pow(cos_theta_pp, 2));
                        x23 = x23_pp * cos_theta_pp - y23_pp * sin_theta_pp + zp;
                        y23 = x23_pp * sin_theta_pp + y23_pp * cos_theta_pp;
                        if ((pow(x23, 2) + pow(y23, 2)) < 1){
                            casenum = 13;
                            alpha1 = overlap(1, pm, zm, false);
                            alpha2 = overlap(pp, pm, zpm, true);
                            
                            ap = alpha;
                            am = (alpha1 - alpha2);
                        }
                        else{
                            casenum = 14;
                            alpha1 = overlap(1, pm, zm, false);
                            
                            ap = alpha;
                            am = alpha1;
                        }
                    }
                }
            }
            else{
                casenum = 15;
                alpha = overlap(1, pp, zp, false);
                
                ap = alpha;
                am = 0;
            }
        }
        else{
            casenum = 16;
            if (zpm >= pm + pp){
                casenum = 17;
                alpha = overlap(1, pp, zp, false);
                
                ap = alpha;
                am = pow(pm, 2) * PI;
            }
            else if ((pp - pm < zpm) & (zpm < pp + pm)){
                casenum = 18;
                alpha = overlap(1, pp, zp, false);
                alpha1 = overlap(pp, pm, zpm, true);
                ap = alpha;
                am = pow(pm, 2) * PI - alpha1;
            }
            else{
                casenum = 19;
                alpha = overlap(1, pp, zp, false);
                ap = alpha;
                am = 0;
            }
        }
    }
            
    else{
        casenum = 20;
        if (zm >= (1 + pm)){
            casenum = 21;
            ap = pow(pp, 2) * PI;
            am = 0;
        }
        else if ((zm > (1 - pm)) & (zm < (1 + pm))){
            casenum = 22;
            if (zpm >= pm + pp){
                casenum = 23;
                alpha = overlap(1, pm, zm, false);
                ap = pow(pp, 2) * PI;
                am = alpha;
            }
            else{
                casenum = 24;
                alpha1 = overlap(1, pm, zm, false);
                alpha1 = overlap(pp, pm, zpm, true);
                
                ap = pow(pp, 2) * PI;
                am = (alpha1 - alpha2);
            }
        }

        else{
            casenum = 25;
            if (zpm >= pm + pp){
                ap = pow(pp, 2) * PI;
                am = pow(pm, 2) * PI;
            }
            else if ((pp - pm < zpm) & (zpm < pp + pm)){
                casenum = 26;
                alpha = overlap(pp, pm, zpm, true);
                ap = pow(pp, 2) * PI;
                am = pow(pm, 2) * PI - alpha;
            }
            else{
                casenum = 27;
                ap = pow(pp, 2) * PI;
                am = 0;
            }
        }
    }
    
    a[0] = ap;
    a[1] = am;
}

// compute the light curve
void flux(double *flux, double *zp, double *zm, double *zpm, double pp, double pm, double c1, double c2, double c3, double c4, int m){
    
    int j;
    double a[2];
    double Ip = 0;
    double Im = 0;
    
    double c0 = 1. - c1 - c2 - c3 - c4;
    double Sigma = c0 / 4.0;
    Sigma += c1 / 5.0;
    Sigma += c2 / 6.0;
    Sigma += c3 / 7.0;
    Sigma += c4 / 8.0;
    
    for (j = 0; j < m; j++){
                            
        area(a, zp[j], zm[j], zpm[j], pp, pm);
        //Ip = intensity(zp[j], pp, c1, c2, c3, c4);
        Im = intensity(zm[j], pm, c1, c2, c3, c4);
        flux[j] = -(a[1] * Im) / (4. * Sigma * PI);
    }
}