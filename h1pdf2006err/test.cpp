extern "C" {
    void qcd_2006_(double *z,double *q2, int *ifit, double *xPq, double *f2, double *fl, double *c2, double *cl);
    void h12006flux_(double *xpom, double *t, int *Int, int *ifit, int *ipom, double *flux);
    void h12006pdf_(double *z, double *q2, int *ifit, int *ipdf, double *xpq, double *f2, double *fl, double *c2, double *cl);
}

#include <iostream>

using namespace std;

int main()
{

    double z = 0.3;
    for(double q2 = 1.75; q2 <= 1000; q2 *= 1.3) {
        int ifit = 1, ierr=0 ;

        double xPq[13];
        double f2[2]={}, fl[2]={}, c2[2]={}, cl[2]={};

        double xPqE[13];
        double f2E[2]={}, flE[2]={}, c2E[2]={}, clE[2]={};

        qcd_2006_(&z,&q2, &ifit, xPq, f2, fl, c2, cl);
        h12006pdf_(&z,&q2,&ifit,&ierr,xPqE,f2E,f2E,c2E,clE);

        cout <<"R  : " <<  z <<" "<< q2 <<" "<< xPq[6] <<" "<< xPq[7] <<" : "<< f2[0] <<" "<< fl[0] << " "<<  f2[1] <<" "<< fl[1]  << endl;
        cout <<"Re : " <<  z <<" "<< q2 <<" "<< xPqE[6] <<" "<< xPqE[7] <<" :  "<< f2E[0] <<" "<< flE[0] <<" "<<  f2E[1] <<" "<< flE[1] << endl;

    }

    return 0;
}
