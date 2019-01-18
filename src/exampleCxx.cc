/*
 * ---------------------------------------------------------------------
 * 23-01-17  Basic QCDNUM example job in C++
 * ---------------------------------------------------------------------
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include "QCDNUM/QCDNUM.h"
using namespace std;

extern "C" {
    void qcd_2006_(double *z,double *q2, int *ifit, double *xPq, double *f2, double *fl, double *c2, double *cl);
}

double ag = 0.36781;
double bg = 0;
double cg = 0;

double aq = 0.69876;
double bq = 1.5024;
double cq = 0.44690;

double xglu(double x) {
  double pd = ag * pow(x,bg) * pow((1-x),cg);
  return pd;
}

double xq(double x) {
  double pd = aq/6. * pow(x,bq) * pow((1-x),cq);
  return pd;
}


double sup(double z)
{
    return exp(-0.01/(1-z));
}

double func(int* ipdf, double* x) {
  int i = *ipdf;
  double xb = *x;
  double f = 0;
  double r = exp(-0.01/(1-xb));
  if(i ==  0) f = r * xglu(xb);
  if(i ==  1) f = 0;//xdnv(xb); //xdnv
  if(i ==  2) f = 0;//xupv(xb); //xupv
  if(i ==  3) f = 0;
  if(i ==  4) f = r * xq(xb);
  if(i ==  5) f = r * xq(xb);
  if(i ==  6) f = r * xq(xb);
  if(i ==  7) f = 0;
  if(i ==  8) f = 0;
  if(i ==  9) f = 0;
  if(i == 10) f = 0;
  if(i == 11) f = 0;
  if(i == 12) f = 0;
  return f;
}



//----------------------------------------------------------------------
int main() {
    double z  = 0.3;
    double q2 = 2.5;
    int ifit = 2;
    double xPq[13];
    double f2, fl, c2, cl;
    qcd_2006_(&z,&q2, &ifit, xPq, &f2, &fl, &c2, &cl);

    //return 0;



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














  int iq0  = QCDNUM::iqfrmq(q0);                           //start scale
  QCDNUM::evolfg(1,func,def,iq0,eps);                 //evolve all pdf's
  QCDNUM::allfxq(1,x,q,pdf,0,1);                 //interpolate all pdf's
    
  double csea = 2 * pdf[2];                         //charm sea at x,mu2
  double asmz = QCDNUM::asfunc(qmz2,nfout,ierr);           //alphas(mz2)
  
  string s = " ";
  cout << scientific << setprecision(4);
  cout << "x, q, CharmSea = " << x << s << q << s << csea << endl;
  cout << "as(mz2)        = " << asmz << endl;

  q2 = 50;
  for(double z = 0.1; z < 1; z += 0.1) {
      qcd_2006_(&z,&q2, &ifit, xPq, &f2, &fl, &c2, &cl);
      //cout << z <<" "<<  q2 << " : " << xPq[6] <<" "<<  xPq[7]  << " | " << sup(z)*xglu(z) <<" "<< sup(z)*xq(z) <<  endl;

      QCDNUM::allfxq(1,z,q2,pdf,0,1);                 //interpolate all pdf's
      cout << z <<" "<< q2 <<" : " << xPq[7] / pdf[7] << endl;
  }

  z = 0.04;
  //ZMSTFUN ( 2, def, z, Q2, *f, n, ichk ) //F2

  // Charge squared weighted for proton
  double pro[] = { 4., 1., 4., 1., 4., 1., 0., 1., 4., 1., 4., 1., 4. };
  for(int i=0; i<13; i++) pro[i] /= 9;
  int ichk = 1;
  int iset = 1;
  double F2, FL;
  double xAr[] = {z};
  double qAr[] = {q2};
  QCDNUM::zmstfun(2,pro,xAr,qAr,&F2, 1,ichk);
  QCDNUM::zmstfun(1,pro,xAr,qAr,&FL, 1,ichk);

  cout <<"QCDNUM : " <<  z << " "<< q2 << " : " << F2 << " " << FL <<  endl;
  qcd_2006_(&z,&q2, &ifit, xPq, &f2, &fl, &c2, &cl);
  cout <<"h1FitB : " <<  z << " "<< q2 << " : " << f2 << " " << fl<< endl;

  return 0;
}
