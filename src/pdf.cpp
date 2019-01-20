#include "pdf.h"

#include <tuple>
#include <iostream>
#include "QCDNUM/QCDNUM.h"

using namespace std;

PDF *pdfNow;

double func(int* ipdf, double* x) {
  int i = *ipdf;
  double xb = *x;

  double xg, xq;
  tie(xg, xq) = pdfNow->eval0(xb);

  if(i == 0)
      return xg;
  else if(4 <= i  && i <= 6)
      return xq;
  //else if(i == 1)
      //return xq*6;
  else
      return 0;
}




void PDF::evolve()
{
    //lamda QCD = 0.399(37) GeV  | http://www-h1.desy.de/psfiles/papers/desy06-049.pdf (page 18)
    int    ityp = 1, iord = 2, nfin = 3;                //unpol, NLO, VFNS
    double as0 = 0.368915, r20 = 1.75;                       //input alphas (corresponds to 0.399
    //double as0 = 0.331497, r20 = 2.5;                       //input alphas

    double xmin[] = {1.e-4, 0.02, 0.1, 0.4, 0.7};                                      //x-grid
    //double xmin[] = {1.e-4};
    //double xmin[] = {1.e-5, 0.2, 0.4, 0.6, 0.75};                                      //x-grid
    int    iwt[] = {1,2,4,8,16}, ng = 1, nxin = 1600, iosp = 2;             //x-grid
    int    nqin = 480;                                           //mu2-grid
    double qq[] = { 1.75e0, 1600 }, wt[] = { 1e0, 1e0};              //mu2-grid
    double q2c = pow(1.4,2), q2b = pow(4.5,2), q0 = 1.75;                   //thresholds, mu20
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
    QCDNUM::setcbt(nfin,iqc,iqb,9999);             //thresholds in the VFNS

    int iq0  = QCDNUM::iqfrmq(q0);                           //start scale
    pdfNow = this;
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

std::pair<double,double> PDF::evalFitA(double z, double q2)
{
    int ifit = 1;
    double xPq[13];
    double f2[2], fl[2], c2[2], cl[2];
    qcd_2006_(&z,&q2, &ifit, xPq, f2, fl, c2, cl);
    return {xPq[6], xPq[7]};
}

double PDF::evalFitRed(double xpom, double z, double q2) const
{
    //Get F2 and FL
    int ifit = 1;//FitA
    double xPq[13];
    double f2[2], fl[2], c2[2], cl[2]; //for pomeron & regeon
    qcd_2006_(&z,&q2, &ifit, xPq, f2, fl, c2, cl);

    //f2 = 2*(2*1./9 + 1.5*4./9) * xPq[7];
    //fl = 0;
    //cout << "Radek " <<z<<" "<< q2 <<" | "<<  2*(2*1./9 + 4./9) * xPq[7] <<" "<< f2 << endl;

    //Multiply by flux
    double t = -1, flux;

    //Pomeron
    int Int = 1, ipom = 1;
    h12006flux_(&xpom, &t, &Int, &ifit, &ipom, &flux);
    f2[0] *= flux;
    fl[0] *= flux;
    c2[0] *= flux;
    cl[0] *= flux;

    //Regeon
    ipom = 2;
    h12006flux_(&xpom, &t, &Int, &ifit, &ipom, &flux);
    f2[1] *= flux;
    fl[1] *= flux;
    c2[1] *= flux;
    cl[1] *= flux;



    //Get the reduced xSec

    const double mp2 = pow(0.92, 2);
    double Ep = (q2 < 120) ? 820 : 920;
    double Ee = 27.5;
    const double s = 4*Ep * Ee;

    double x = z*xpom;
    double y = q2/(s-mp2)/x;

    double sRedP = (f2[0]) - y*y/(1 + pow(1-y,2)) * (fl[0]);
    double sRedR = (f2[1]) - y*y/(1 + pow(1-y,2)) * (fl[1]);
    double sRedC = (c2[0]) - y*y/(1 + pow(1-y,2)) * (cl[0]); //should not be included

    return xpom*(sRedP+sRedR);
   // return sRed;
}

double PDF::evalFitRedMy(double xpom, double z, double q2) const
{
    //Get F2 and FL
    double pro[] = { 4., 1., 4., 1., 4., 1., 0., 1., 4., 1., 4., 1., 4. };
    for(int i=0; i<13; i++) pro[i] /= 9;
    int ichk = 1;
    double f2, fl;
    double xAr[] = {z};
    double qAr[] = {q2};
    QCDNUM::zmstfun(2,pro,xAr,qAr,&f2, 1,ichk);
    QCDNUM::zmstfun(1,pro,xAr,qAr,&fl, 1,ichk);


    //Multiply by flux
    double t = -1, flux;
    int ifit = 1, Int = 1, ipom = 1;
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

    double sRed = (f2) - y*y/(1 + pow(1-y,2)) * (fl);
    return xpom*sRed;// TODO why is it wrong?
   // return sRed;
}
