#include "pdf.h"

#include <tuple>
#include <iostream>
#include "QCDNUM/QCDNUM.h"

using namespace std;

PDF *fitB;

double func(int* ipdf, double* x) {
  int i = *ipdf;
  double xb = *x;

  double xg, xq;
  tie(xg, xq) = fitB->eval0(xb);

  if(i == 0)
      return xg;
  else if(4 <= i  && i <= 6)
      return xq;
  else
      return 0;
}




void PDF::evolve()
{
    int    ityp = 1, iord = 2, nfin = 3;                //unpol, NLO, VFNS
    double as0 = 0.364, r20 = 2.0;                          //input alphas
    double xmin[] = {1.e-3};                                      //x-grid
    int    iwt[] = {1}, ng = 1, nxin = 280, iosp = 3;             //x-grid
    int    nqin = 120;                                           //mu2-grid
    double qq[] = { 2.5e0, 3000 }, wt[] = { 1e0, 1e0};              //mu2-grid
    double q2c = 3, q2b = 25, q0 = 2.5;                   //thresholds, mu20
    double x = 1e-3, q = 1e3, qmz2 = 8315.25, pdf[13];            //output

    double def[] =                             //input flavour composition
        // tb  bb  cb  sb  ub  db   g   d   u   s   c   b   t
    { 0., 0., 0., 0., 0.,-1., 0., 1., 0., 0., 0., 0., 0.,      // 1=dval
        0., 0., 0., 0.,-1., 0., 0., 0., 1., 0., 0., 0., 0.,      // 2=uval
        0., 0., 0.,-1., 0., 0., 0., 0., 0., 1., 0., 0., 0.,      // 3=sval
        0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0.,      // 4=dbar
        0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0.,      // 5=ubar
        0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.,      // 6=sbar
        0., 0.,-1., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0.,      // 7=cval
        0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,      // 8=cbar
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,      // 9=zero
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,      //10=zero
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,      //11=zero
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};     //12=zero

    int nx, nq, id1, id2, nw, nfout ; double eps;
    int lun = 6 ; string outfile = " ";

    QCDNUM::qcinit(lun,outfile);                              //initialize
    QCDNUM::gxmake(xmin,iwt,ng,nxin,nx,iosp);                     //x-grid
    QCDNUM::gqmake(qq,wt,2,nqin,nq);                            //mu2-grid
    QCDNUM::fillwt(1,id1,id2,nw);                      //calculate weights
    QCDNUM::setord(iord);                                  //LO, NLO, NNLO
    QCDNUM::setalf(as0,r20);                                //input alphas
    int iqc  = QCDNUM::iqfrmq(q2c);                      //charm threshold
    int iqb  = QCDNUM::iqfrmq(q2b);                     //bottom threshold
    QCDNUM::setcbt(nfin,iqc,iqb,999);             //thresholds in the VFNS

    int iq0  = QCDNUM::iqfrmq(q0);                           //start scale
    QCDNUM::evolfg(1,func,def,iq0,eps);                 //evolve all pdf's
}

pair<double,double> PDF::eval(double z, double q2) const
{
    double pdf[13];
    QCDNUM::allfxq(1,z,q2,pdf,0,1);                 //interpolate all pdf's
    return {pdf[6], pdf[7]};
}


void PDF::initConv()
{
    //Structure function init
    // Try to read the weight file and create one if that fails
    string wname = "unpol.wgt";
    int idmin = 0 ;
    int idmax = 0;
    int nwpdf = 0;
    int ierr  = 0;
    QCDNUM::readwt(22,wname,idmin,idmax,nwpdf,ierr);      //read weights
    if(ierr != 0) {
        QCDNUM::fillwt(1,idmin,idmax,nwpdf);                //calculate weights
        QCDNUM::dmpwgt(1,22,wname);                         //dump weights
    }

    int nztot, nzuse;
    QCDNUM::zmwords(nztot,nzuse);
    cout << " nztot, nzuse = " <<  nztot << " " << nzuse << endl;

    // Try to read the weight file and create one if that fails
    string zmname = "zmstf.wgt";
    int nwords = 0;
    QCDNUM::zmreadw(22,zmname,nwords,ierr);
    if(ierr != 0) {
        QCDNUM::zmfillw(nwords);
        QCDNUM::zmdumpw(22,zmname);
    }
    cout << " QCDNUM: words used = " << nwords << endl;
    cout << endl;

    QCDNUM::zmwords(nztot,nzuse);
    cout << " nztot, nzuse = " <<  nztot << "  " << nzuse << endl;




}


extern "C" {
    void qcd_2006_(double *z,double *q2, int *ifit, double *xPq, double *f2, double *fl, double *c2, double *cl);
    void h12006flux_(double *xpom, double *t, int *Int, int *ifit, int *ipom, double *flux);
}

std::pair<double,double> PDF::evalFitB(double z, double q2) const
{
    int ifit = 2;
    double xPq[13];
    double f2, fl, c2, cl;
    qcd_2006_(&z,&q2, &ifit, xPq, &f2, &fl, &c2, &cl);
    return {xPq[6], xPq[7]};
}

double PDF::evalFitBred(double xpom, double z, double q2) const
{
    //Get F2 and FL
    int ifit = 2;
    double xPq[13];
    double f2, fl, c2, cl;
    qcd_2006_(&z,&q2, &ifit, xPq, &f2, &fl, &c2, &cl);


    //Multiply by flux
    double t = -1, flux;
    int Int = 1, ipom = 1;
    h12006flux_(&xpom, &t, &Int, &ifit, &ipom, &flux);
    f2 *= flux;
    fl *= flux;


    //Get the reduced xSec

    const double mp2 = pow(0.92, 2);
    double Ep = (q2 < 120) ? 820 : 920;
    double Ee = 27.5;
    const double s = 4*Ep * Ee;

    double x = z*xpom;
    double y = q2/(s-mp2)/x;

    double sRed = f2 - y*y/(1 + pow(1-y,2)) * fl;
    return xpom*sRed;
    return sRed;
}
