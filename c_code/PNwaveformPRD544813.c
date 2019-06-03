//$Id: PNwaveformPRD544813.C,v 1.1.1.1 2016/12/30 06:03:09 zjcao Exp $
#ifdef newc
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <complex>
using namespace std;
#else
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#endif

// PN waveform, PRD 54, 4813, Eq.(6.10)
// we have assumed M = R = 1
// so dm := m1-m2 = sqrt(m(m-4mu))=sqrt(1-4 eta)
void PNwaveformPRD544813rdotc_22mode(double *hr,double *hi,
			 const double x1,const double x2,const double x3, // test particle position
			 const double v1,const double v2,const double v3, // test particle velocity
			 const double eta) // symmetric mass ratio
{
    const double dm = sqrt(1-4*eta);
    double r = sqrt(x1*x1+x2*x2+x3*x3);
    double n1 = x1/r,n2 = x2/r,n3 = x3/r;
    double rdot = v1*n1+v2*n2+v3*n3;
    double vsqr = v1*v1+v2*v2+v3*v3;

    double vn1,vn2,vn3;
    vn1 = rdot*n1;
    vn2 = rdot*n2;
    vn3 = rdot*n3;
    double lambda1,lambda2,lambda3;
    lambda1 = v1-vn1;
    lambda2 = v2-vn2;
    lambda3 = v3-vn3;
    double ln = sqrt(lambda1*lambda1+lambda2*lambda2+lambda3*lambda3);
    lambda1 = lambda1/ln;
    lambda2 = lambda2/ln;
    lambda3 = lambda3/ln;

    double complex hm22;

    double complex h11,h12,h13,h22,h23,h33;
    double complex Q11,Q12,Q13,Q22,Q23,Q33;

    ///////////////////////////////////////00 part///////////////////////////////////////
    Q11 = sqrt(M_PI/5)/3;
    Q12 = -sqrt(M_PI/5)/3*I;
    Q13 = 0;
    Q22 = -sqrt(M_PI/5)/3;
    Q23 = 8/sqrt(5*M_PI)/3;
    Q33 = 0;

    h11 = 0;
    h12 = 0;
    h13 = 0;
    h22 = 0;
    h23 = 0;
    h33 = 0;
   
    h11 += 1.0/3*(2/r*rdot*(5+3*eta)*(n1*v1+v1*n1)+(3*(1-3*eta)*rdot*rdot)/r*n1*n1);
    h12 += 1.0/3*(2/r*rdot*(5+3*eta)*(n1*v2+v1*n2)+(3*(1-3*eta)*rdot*rdot)/r*n1*n2);
    h13 += 1.0/3*(2/r*rdot*(5+3*eta)*(n1*v3+v1*n3)+(3*(1-3*eta)*rdot*rdot)/r*n1*n3);
    h22 += 1.0/3*(2/r*rdot*(5+3*eta)*(n2*v2+v2*n2)+(3*(1-3*eta)*rdot*rdot)/r*n2*n2);
    h23 += 1.0/3*(2/r*rdot*(5+3*eta)*(n2*v3+v2*n3)+(3*(1-3*eta)*rdot*rdot)/r*n2*n3);
    h33 += 1.0/3*(2/r*rdot*(5+3*eta)*(n3*v3+v3*n3)+(3*(1-3*eta)*rdot*rdot)/r*n3*n3);

    hm22 = Q11*h11+Q22*h22+Q33*h33+2.*(Q12*h12+Q13*h13+Q23*h23);

    ///////////////////////////////////////10 part///////////////////////////////////////
    Q11 = 0;
    Q12 = 0;
    Q13 = -sqrt(M_PI/5)*(x1 - I*x2)/(3.*r);
    Q22 = 0;
    Q23 = sqrt(M_PI/5)*(I*x1 + x2)/(3.*r);
    Q33 = 0;

    h11 = dm*3/r*(-rdot*n1*n1);
    h12 = dm*3/r*(-rdot*n1*n2);
    h13 = dm*3/r*(-rdot*n1*n2);
    h22 = dm*3/r*(-rdot*n2*n2);
    h23 = dm*3/r*(-rdot*n2*n3);
    h33 = dm*3/r*(-rdot*n3*n3);
    
    h11 += dm/12/r*((n1*v1+v1*n1)*(rdot*rdot*(63+54*eta))
			 +n1*n1*rdot*(rdot*rdot*(15-90*eta)-vsqr*(63-54*eta)+(242-24*eta)/r)-rdot*v1*v1*(186+24*eta));
    h12 += dm/12/r*((n1*v2+v1*n2)*(rdot*rdot*(63+54*eta))
			 +n1*n2*rdot*(rdot*rdot*(15-90*eta)-vsqr*(63-54*eta)+(242-24*eta)/r)-rdot*v1*v2*(186+24*eta));
    h13 += dm/12/r*((n1*v3+v1*n3)*(rdot*rdot*(63+54*eta))
			 +n1*n3*rdot*(rdot*rdot*(15-90*eta)-vsqr*(63-54*eta)+(242-24*eta)/r)-rdot*v1*v3*(186+24*eta));
    h22 += dm/12/r*((n2*v2+v2*n2)*(rdot*rdot*(63+54*eta))
			 +n2*n2*rdot*(rdot*rdot*(15-90*eta)-vsqr*(63-54*eta)+(242-24*eta)/r)-rdot*v2*v2*(186+24*eta));
    h23 += dm/12/r*((n2*v3+v2*n3)*(rdot*rdot*(63+54*eta))
			 +n2*n3*rdot*(rdot*rdot*(15-90*eta)-vsqr*(63-54*eta)+(242-24*eta)/r)-rdot*v2*v3*(186+24*eta));
    h33 += dm/12/r*((n3*v3+v3*n3)*(rdot*rdot*(63+54*eta))
			 +n3*n3*rdot*(rdot*rdot*(15-90*eta)-vsqr*(63-54*eta)+(242-24*eta)/r)-rdot*v3*v3*(186+24*eta));
    
    hm22 += Q11*h11+Q22*h22+Q33*h33+2.*(Q12*h12+Q13*h13+Q23*h23);
    
    ///////////////////////////////////////01 part///////////////////////////////////////
    Q11 = 0;
    Q12 = 0;
    Q13 = -sqrt(M_PI/5)*(v1 - I*v2)/3.;
    Q22 = 0;
    Q23 = sqrt(M_PI/5)*(I*v1 + v2)/3.;
    Q33 = 0;

    h11 = 0;
    h12 = 0;
    h13 = 0;
    h22 = 0;
    h23 = 0;
    h33 = 0;

    h11 += dm*(-(n1*v1+v1*n1)/2/r*rdot*(7+4*eta)
	            -n1*n1/r*(0.75*(1-2*eta)*rdot*rdot));
    h12 += dm*(-(n1*v2+v1*n2)/2/r*rdot*(7+4*eta)
	            -n1*n2/r*(0.75*(1-2*eta)*rdot*rdot));
    h13 += dm*(-(n1*v3+v1*n3)/2/r*rdot*(7+4*eta)
	            -n1*n3/r*(0.75*(1-2*eta)*rdot*rdot));
    h22 += dm*(-(n2*v2+v2*n2)/2/r*rdot*(7+4*eta)
	            -n2*n2/r*(0.75*(1-2*eta)*rdot*rdot));
    h23 += dm*(-(n2*v3+v2*n3)/2/r*rdot*(7+4*eta)
	            -n2*n3/r*(0.75*(1-2*eta)*rdot*rdot));
    h33 += dm*(-(n3*v3+v3*n3)/2/r*rdot*(7+4*eta)
	            -n3*n3/r*(0.75*(1-2*eta)*rdot*rdot));

    hm22 += Q11*h11+Q22*h22+Q33*h33+2.*(Q12*h12+Q13*h13+Q23*h23);

    ///////////////////////////////////////20 part///////////////////////////////////////
    Q11 = -sqrt(M_PI/5)*(x1*x1-8.*I*x1*x2-7*x2*x2-x3*x3)/21./r/r;
    Q12 = -I*sqrt(M_PI/5)*(3*(x1*x1+x2*x2)+x3*x3)/21./r/r;
    Q13 = 2*sqrt(M_PI/5)*(x1-I*x2)*x3/21./r/r;
    Q22 = sqrt(M_PI/5)*(-7*x1*x1+8.*I*x1*x2+x2*x2-x3*x3)/21./r/r;
    Q23 = -2.*I*sqrt(M_PI/5)*(x1-I*x2)*x3/21./r/r;
    Q33 = 8*sqrt(M_PI/5)*pow(x1-I*x2,2)/21./r/r;

    h11 = (1-3*eta)/3/r*((-15*rdot*rdot)*n1*n1+15*rdot*(n1*v1+v1*n1));
    h12 = (1-3*eta)/3/r*((-15*rdot*rdot)*n1*n2+15*rdot*(n1*v2+v1*n2));
    h13 = (1-3*eta)/3/r*((-15*rdot*rdot)*n1*n3+15*rdot*(n1*v3+v1*n3));
    h22 = (1-3*eta)/3/r*((-15*rdot*rdot)*n2*n2+15*rdot*(n2*v2+v2*n2));
    h23 = (1-3*eta)/3/r*((-15*rdot*rdot)*n2*n3+15*rdot*(n2*v3+v2*n3));
    h33 = (1-3*eta)/3/r*((-15*rdot*rdot)*n3*n3+15*rdot*(n3*v3+v3*n3));

    hm22 += Q11*h11+Q22*h22+Q33*h33+2.*(Q12*h12+Q13*h13+Q23*h23);

    ///////////////////////////////////////11 part///////////////////////////////////////
    Q11 = sqrt(M_PI/5)*(-v1*x1+4.*I*v2*x1+4.*I*v1*x2+7*v2*x2+v3*x3)/21./r;
    Q12 = -I*sqrt(M_PI/5)*(3*v1*x1+3*v2*x2+v3*x3)/21./r;
    Q13 = sqrt(M_PI/5)*(v3*(x1-I*x2)+x3*(v1-I*v2))/21./r;
    Q22 = sqrt(M_PI/5)*(-7*v1*x1+4.*I*v2*x1+4.*I*v1*x2+v2*x2-v3*x3)/21./r;
    Q23 = -sqrt(M_PI/5)*(v3*(I*x1+x2)+x3*(I*v1+v2))/21./r;
    Q33 = 8*sqrt(M_PI/5)*(v1-I*v2)*(x1-I*x2)/21./r;

    h11 = (1-3*eta)/3/r*(12*rdot*n1*n1);
    h12 = (1-3*eta)/3/r*(12*rdot*n1*n2);
    h13 = (1-3*eta)/3/r*(12*rdot*n1*n3);
    h22 = (1-3*eta)/3/r*(12*rdot*n2*n2);
    h23 = (1-3*eta)/3/r*(12*rdot*n2*n3);
    h33 = (1-3*eta)/3/r*(12*rdot*n3*n3);

    hm22 += Q11*h11+Q22*h22+Q33*h33+2.*(Q12*h12+Q13*h13+Q23*h23);

    ///////////////////////////////////////02 part///////////////////////////////////////
    Q11 = sqrt(M_PI/5)*(-v1*v1+8.*I*v1*v2+7*v2*v2+v3*v3)/21.;
    Q12 = -I*sqrt(M_PI/5)*(3*(v1*v1+v2*v2)+v3*v3)/21.;
    Q13 = 2*sqrt(M_PI/5)*(v1-I*v2)*v3/21.;
    Q22 = sqrt(M_PI/5)*(-7*v1*v1+8.*I*v1*v2+v2*v2-v3*v3)/21.;
    Q23 = -2.*I*sqrt(M_PI/5)*(v1-I*v2)*v3/21.;
    Q33 = 8*sqrt(M_PI/5)*pow(v1-I*v2,2)/21.;

    h11 = 0;
    h12 = 0;
    h13 = 0;
    h22 = 0;
    h23 = 0;
    h33 = 0;

    hm22 += Q11*h11+Q22*h22+Q33*h33+2.*(Q12*h12+Q13*h13+Q23*h23);

    ///////////////////////////////////////30 part///////////////////////////////////////
    Q11 = sqrt(M_PI/5)*pow(x1-I*x2,2)*x3/7./r/r/r;
    Q12 = 0;
    Q13 = -sqrt(M_PI/5)*(x1-I*x2)*(2.*x1*x1+5.*I*x1*x2+7*x2*x2+3*x3*x3)/21./r/r/r;
    Q22 = sqrt(M_PI/5)*pow(x1-I*x2,2)*x3/7./r/r/r;
    Q23 = sqrt(M_PI/5)*(I*x1+x2)*(7*x1*x1-5.*I*x1*x2+2*x2*x2+3*x3*x3)/21./r/r/r;
    Q33 = -2*sqrt(M_PI/5)*pow(x1-I*x2,2)*x3/7./r/r/r;

    h11 = dm*(1-2*eta)/r*(1.25*(3*vsqr-7*rdot*rdot+6/r)*rdot*n1*n1-8.5*rdot*v1*v1-(-105*rdot*rdot)/12*(n1*v1+v1*n1));
    h12 = dm*(1-2*eta)/r*(1.25*(3*vsqr-7*rdot*rdot+6/r)*rdot*n1*n2-8.5*rdot*v1*v2-(-105*rdot*rdot)/12*(n1*v2+v1*n2));
    h13 = dm*(1-2*eta)/r*(1.25*(3*vsqr-7*rdot*rdot+6/r)*rdot*n1*n3-8.5*rdot*v1*v3-(-105*rdot*rdot)/12*(n1*v3+v1*n3));
    h22 = dm*(1-2*eta)/r*(1.25*(3*vsqr-7*rdot*rdot+6/r)*rdot*n2*n2-8.5*rdot*v2*v2-(-105*rdot*rdot)/12*(n2*v2+v2*n2));
    h23 = dm*(1-2*eta)/r*(1.25*(3*vsqr-7*rdot*rdot+6/r)*rdot*n2*n3-8.5*rdot*v2*v3-(-105*rdot*rdot)/12*(n2*v3+v2*n3));
    h33 = dm*(1-2*eta)/r*(1.25*(3*vsqr-7*rdot*rdot+6/r)*rdot*n3*n3-8.5*rdot*v3*v3-(-105*rdot*rdot)/12*(n3*v3+v3*n3));

    hm22 += Q11*h11+Q22*h22+Q33*h33+2.*(Q12*h12+Q13*h13+Q23*h23);

    ///////////////////////////////////////21 part///////////////////////////////////////
    Q11 = sqrt(M_PI/5)*(x1-I*x2)*(v3*(x1-I*x2)+2.*(v1-I*v2)*x3)/21./r/r;
    Q12 = 0;
    Q13 = -sqrt(M_PI/5)*((x1-I*x2)*(2*v1*x1+I*x1*v2+4.*I*v1*x2+7*x2*v2)+2*v3*(x1-I*x2)*x3+(v1-I*v2)*x3*x3)/21./r/r;
    Q22 = sqrt(M_PI/5)*(x1-I*x2)*(v3*(x1-I*x2)+2.*(v1-I*v2)*x3)/21./r/r;
    Q23 = sqrt(M_PI/5)*((x1-I*x2)*(v2*(4*x1+2.*I*x2)+v1*(7.*I*x1+x2))+2*v3*(I*x1+x2)*x3+(I*v1+v2)*x3*x3)/21./r/r;
    Q33 = -2*sqrt(M_PI/5)*(x1-I*x2)*(v3*(x1-I*x2)+2.*(v1-I*v2)*x3)/21./r/r;

    h11 = dm*(1-2*eta)/4/r*((45*rdot*rdot)*n1*n1-54*rdot*(n1*v1+v1*n1));
    h12 = dm*(1-2*eta)/4/r*((45*rdot*rdot)*n1*n2-54*rdot*(n1*v2+v1*n2));
    h13 = dm*(1-2*eta)/4/r*((45*rdot*rdot)*n1*n3-54*rdot*(n1*v3+v1*n3));
    h22 = dm*(1-2*eta)/4/r*((45*rdot*rdot)*n2*n2-54*rdot*(n2*v2+v2*n2));
    h23 = dm*(1-2*eta)/4/r*((45*rdot*rdot)*n2*n3-54*rdot*(n2*v3+v2*n3));
    h33 = dm*(1-2*eta)/4/r*((45*rdot*rdot)*n3*n3-54*rdot*(n3*v3+v3*n3));

    hm22 += Q11*h11+Q22*h22+Q33*h33+2.*(Q12*h12+Q13*h13+Q23*h23);

    ///////////////////////////////////////12 part///////////////////////////////////////
    Q11 = sqrt(M_PI/5)*(v1-I*v2)*(2*v3*(x1-I*x2)+(v1-I*v2)*x3)/21./r;
    Q12 = 0;
    Q13 = -sqrt(M_PI/5)*(v3*v3*(x1-I*x2)+v1*v1*(2*x1+I*x2)+v2*v2*(4*x1-7.*I*x2)-2.*I*v2*v3*x3+2*v1*(I*v2*x1+4*v2*x2+v3*x3))/21./r;
    Q22 = sqrt(M_PI/5)*(v1-I*v2)*(2*v3*(x1-I*x2)+(v1-I*v2)*x3)/21./r;
    Q23 = sqrt(M_PI/5)*(v3*v3*(I*x1+x2)+v2*v2*(-I*x1+2*x2)+v1*v1*(7.*I*x1+4*x2)+2*v2*v3*x3+v1*(8*v2*x1-2.*I*v2*x2+2.*I*v3*x3))/21./r;
    Q33 = -2*sqrt(M_PI/5)*(v1-I*v2)*(2*v3*(x1-I*x2)+(v1-I*v2)*x3)/21./r;

    h11 = dm*(1-2*eta)*1.5/r*(-3*rdot*n1*n1);
    h12 = dm*(1-2*eta)*1.5/r*(-3*rdot*n1*n2);
    h13 = dm*(1-2*eta)*1.5/r*(-3*rdot*n1*n3);
    h22 = dm*(1-2*eta)*1.5/r*(-3*rdot*n2*n2);
    h23 = dm*(1-2*eta)*1.5/r*(-3*rdot*n2*n3);
    h33 = dm*(1-2*eta)*1.5/r*(-3*rdot*n3*n3);

    hm22 += Q11*h11+Q22*h22+Q33*h33+2.*(Q12*h12+Q13*h13+Q23*h23);

    ///////////////////////////////////////03 part///////////////////////////////////////
    Q11 = sqrt(M_PI/5)*pow(v1-I*v2,2)*v3/7.;
    Q12 = 0;
    Q13 = -sqrt(M_PI/5)*(v1-I*v2)*(2*v1*v1+5.*I*v1*v2+7*v2*v2+3*v3*v3)/21.;
    Q22 = sqrt(M_PI/5)*pow(v1-I*v2,2)*v3/7.;
    Q23 = sqrt(M_PI/5)*(I*v1+v2)*(7*v1*v1-5.*I*v1*v2+2*v2*v2+3*v3*v3)/21.;
    Q33 = -2*sqrt(M_PI/5)*pow(v1-I*v2,2)*v3/7.;

    h11 = 0;
    h12 = 0;
    h13 = 0;
    h22 = 0;
    h23 = 0;
    h33 = 0;

    hm22 += Q11*h11+Q22*h22+Q33*h33+2.*(Q12*h12+Q13*h13+Q23*h23);
    //-----------------------------------complete-----------------------------------------

    hm22 = hm22*2.*eta;
    *hr = creal(hm22);
    *hi = cimag(hm22);
    //exit(0);
}
