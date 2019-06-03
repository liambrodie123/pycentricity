//$Id: PNwaveformPRD544813.h,v 1.1.1.1 2016/12/30 06:03:09 zjcao Exp $
#ifdef newc
#include <cstdlib>
#include <cstring>
#include <cmath>
using namespace std;
#else
#include <stdlib.h>
#include <string.h>
#include <math.h>
#endif

// PN waveform, PRD 54, 4813, Eq.(6.10)
// we have assumed M = R = 1
// so dm := m1-m2 = sqrt(m(m-4mu))=sqrt(1-4 eta)
void PNwaveformPRD544813rdotc_22mode(double *hr,double *hi,
			 const double x1,const double x2,const double x3, // test particle position
			 const double v1,const double v2,const double v3, // test particle velocity
			 const double eta); // symmetric mass ratio
