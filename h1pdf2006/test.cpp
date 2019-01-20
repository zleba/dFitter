extern "C" {
    void qcd_2006_(double *z,double *q2, int *ifit, double *xPq, double *f2, double *fl, double *c2, double *cl);
    void h12006flux_(double *xpom, double *t, int *Int, int *ifit, int *ipom, double *flux);
}

#include <iostream>

using namespace std;

int main()
{

    double z = 0.3;
    for(double q2 = 1.75; q2 <= 1000; q2 *= 1.3) {
        int ifit = 1;
        double xPq[13];
        double f2[4]={}, fl[4]={}, c2[4]={}, cl[4]={};
        qcd_2006_(&z,&q2, &ifit, xPq, f2, fl, c2, cl);

        cout << z <<" "<< q2 <<" "<< xPq[6] <<" "<< xPq[7] <<" "<<
            f2[0] <<" "<< fl[0] << " " <<
            f2[1] <<" "<< fl[1] << " " <<
            c2[0] <<" "<< cl[0] << " " <<
            c2[1] <<" "<< cl[1] << " " <<
            
            endl;

    }

    return 0;
}
