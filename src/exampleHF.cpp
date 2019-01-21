#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <time.h>
#include "QCDNUM/QCDNUM.h"
using namespace std;

//----------------------------------------------------------------------
double xupv(double x) {
  double au = 5.107200;
  double pd = au * pow(x,0.8) * pow((1-x),3);
  return pd;
}
//----------------------------------------------------------------------
double xdnv(double x) {
  double ad = 3.064320;
  double pd = ad * pow(x,0.8) * pow((1-x),4);
  return pd;
}
//----------------------------------------------------------------------
double xglu(double x) {
  double ag = 1.7;
  double pd = ag * pow(x,-0.1) * pow((1-x),5);
  return pd;
}
//----------------------------------------------------------------------
double xdbar(double x) {
  double adbar = 0.1939875;
  double pd = adbar * pow(x,-0.1) * pow((1-x),6);
  return pd;
}
//----------------------------------------------------------------------
double xubar(double x) {
  double pd = xdbar(x) * (1-x);
  return pd;
}
//----------------------------------------------------------------------
double xsbar(double x) {
  double pd = 0.2 * ( xdbar(x) + xubar(x) );
  return pd;
}
//----------------------------------------------------------------------
double func(int* ipdf, double* x) {
  int i = *ipdf;
  double xb = *x;
  double f = 0;
  if(i ==  0) f = xglu(xb);
  if(i ==  1) f = xdnv(xb);
  if(i ==  2) f = xupv(xb);
  if(i ==  3) f = 0;
  if(i ==  4) f = xdbar(xb);
  if(i ==  5) f = xubar(xb);
  if(i ==  6) f = xsbar(xb);
  if(i ==  7) f = 0;
  if(i ==  8) f = 0;
  if(i ==  9) f = 0;
  if(i == 10) f = 0;
  if(i == 11) f = 0;
  if(i == 12) f = 0;
  return f;
}


//   coarse x - Q2 grid
 
double xxtab[] = { 1.0e-5, 2.0e-5, 5.0e-5, 1.0e-4, 2.0e-4, 5.0e-4,
		   1.0e-3, 2.0e-3, 5.0e-3, 1.0e-2, 2.0e-2, 5.0e-2, 1.0e-1,
		   1.5e-1, 2.0e-1, 3.0e-1, 4.0e-1, 5.5e-1, 7.0e-1, 9.0e-1 };
int nxtab  = 20;
double xmi = xxtab[0];
double xma = xxtab[nxtab-1];

double qqtab[] = { 2.0e0, 2.7e0, 3.6e0, 5.0e0, 7.0e0, 1.0e1, 1.4e1,
		    2.0e1, 3.0e1, 5.0e1, 7.0e1, 1.0e2, 2.0e2, 5.0e2, 1.0e3,
            3.0e3, 1.0e4, 4.0e4, 2.0e5, 1.0e6 };
int nqtab  = 20;
double qmi = qqtab[0];
double qma = qqtab[nqtab-1];

//==============================================
void mypdfs(double x, double qmu2, double *xpdf)
//==============================================
{
  QCDNUM::allfxq(1,x,qmu2,xpdf,0,1);
  return;
}

//========
int main()
//========
{ 
  /*
   *   Test the hqstf package
   *
   *   WARNING: In this testjob we evolve with nf=3 flavors and then
   *            calculate F2,L for both charm and bottom. For charm
   *            this is OK but for bottom one should evolve with
   *            nf = 4 flavors! But OK, this is just a testjob...
   */

  // qcdnum input
  double pdfin[] =
    // tb  bb  cb  sb  ub  db   g   d   u   s   c   b   t
    //-6  -5  -4  -3  -2  -1   0   1   2   3   4   5   6
    {  0., 0., 0., 0., 0.,-1., 0., 1., 0., 0., 0., 0., 0.,   //dval
       0., 0., 0., 0.,-1., 0., 0., 0., 1., 0., 0., 0., 0.,   //uval
       0., 0., 0.,-1., 0., 0., 0., 0., 0., 1., 0., 0., 0.,   //sval
       0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0.,   //dbar
       0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0.,   //ubar
       0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.,   //sbar
       0., 0.,-1., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0.,   //cval
       0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,   //cbar
       0.,-1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0.,   //bval
       0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,   //bbar
      -1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1.,   //tval
       1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. }; //tbar
  // Holds q-grid definition
  double qarr[] = {qmi, qma};
  double warr[] = {1., 1.};

  // Define output distributions
  //              tb  bb  cb  sb  ub  db   g   d   u   s   c   b   t
  //              -6  -5  -4  -3  -2  -1   0   1   2   3   4   5   6
  double dnv[] = { 0., 0., 0., 0., 0.,-1., 0., 1., 0., 0., 0., 0., 0. };
  double upv[] = { 0., 0., 0., 0.,-1., 0., 0., 0., 1., 0., 0., 0., 0. };
  double del[] = { 0., 0., 0., 0.,-1., 1., 0., 0., 0., 0., 0., 0., 0. };
  double uds[] = { 0., 0., 0., 2., 2., 2., 0., 0., 0., 0., 0., 0., 0. };
  double pro[] = { 4., 1., 4., 1., 4., 1., 0., 1., 4., 1., 4., 1., 4. };

  // cbt quark masses
  double hqmass[] = { 1.5, 5., 0. };

  // LO/NLO/NNLO
  int iord = 2;
  // VFNS
  int nfix = 3;
  // Define alpha_s
  double alf = 0.35e0;
  double q2a = 2.0e0;
  // More grid parameters
  int iosp  = 3;
  int n_x   = 100;
  int n_q   = 60;
  double q0 = 2.e0;

  //   Initialize
  QCDNUM::qcinit(6," ");
  int lunout;
  QCDNUM::getint("lunq",lunout);

  int nhtot, nhuse;
  QCDNUM::hqwords(nhtot,nhuse);
  cout << " nhtot, nhuse = " <<  nhtot << " " << nhuse << endl;

  // QCDNUM calls to pass the values defined above
  QCDNUM::setord(iord);
  QCDNUM::setalf(alf,q2a);
  // x-mu2 grid
  int nxout;
  double xmin[] = {xmi};
  int iwt[] = {1};
  QCDNUM::gxmake(xmin,iwt,1,n_x,nxout,iosp);
  int nqout;
  QCDNUM::gqmake(qarr,warr,2,n_q,nqout);
  int iq0 = QCDNUM::iqfrmq(q0);
  // Thresholds
  int iqc, iqb, iqt;
  QCDNUM::setcbt(nfix,iqc,iqb,iqt);

  int nf;
  double qc, qb, qt;
  QCDNUM::getcbt(nf,qc,qb,qt);
  cout << " Flavour number scheme " << nf << " " << qc << " " << qb << " " << qt << endl;      

  // Try to read the weight file and create one if that fails
  string wname = "weights/unpol.wgt";
  int idmin = 0 ;
  int idmax = 0;
  int nwpdf = 0;
  int ierr  = 0;
  QCDNUM::readwt(22,wname,idmin,idmax,nwpdf,ierr);      //read weights
  if(ierr != 0) {
    QCDNUM::fillwt(1,idmin,idmax,nwpdf);                //calculate weights
    QCDNUM::dmpwgt(1,22,wname);                         //dump weights
  }

  // Try to read the weight file and create one if that fails
  string hqname = "weights/hqstf.wgt";
  double aq2 = 1.;
  double bq2 = 0.;
  int nwords = 0;
  QCDNUM::hqreadw(22,hqname,nwords,ierr);
  if(ierr == 0) {
    double a, b;
    double qmas[3];
    QCDNUM::hqparms(qmas,a,b);
    if(qmas[0] != hqmass[0]) ierr = 1;
    if(qmas[1] != hqmass[1]) ierr = 1;
    if(qmas[2] != hqmass[2]) ierr = 1;
    if(a != aq2)             ierr = 1;
    if(b != bq2)             ierr = 1;
  }
  if(ierr != 0) {
    QCDNUM::hqfillw(3,hqmass,aq2,bq2,nwords);
    QCDNUM::hqdumpw(22,hqname);
  }
  cout << " QCDNUM: words used = " << nwords << endl;
  cout << endl;

  QCDNUM::hqwords(nhtot,nhuse);
  cout << " nhtot, nhuse = " <<  nhtot << " " << nhuse << endl;

  // Evolve
  double epsi = 1e-3;
  int iter = 1;
  const clock_t tim1 = clock();
  for(int i=0; i<iter; i++) QCDNUM::evolfg(1,func,pdfin,iq0,epsi);
  const clock_t tim2 = clock();

  // Charge squared weighted for proton
  for(int i=0; i<13; i++) pro[i] /= 9;

  int ichk = 1;
  int iset = 1;

  cout << endl;
  cout << "     mu2       x           uv         dv        del        uds         gl" << endl;
  cout << setprecision(4) << scientific;

  for(int iq=0; iq<nqtab; iq=iq+4) {
    double q = qqtab[iq];
    cout << endl;
    for(int ix=0; ix<nxtab; ix=ix+4) {
      double x  = xxtab[ix];
      double uv = QCDNUM::sumfxq(iset,upv,1,x,q,ichk);
      double dv = QCDNUM::sumfxq(iset,dnv,1,x,q,ichk);
      double de = QCDNUM::sumfxq(iset,del,1,x,q,ichk);
      double ud = QCDNUM::sumfxq(iset,uds,1,x,q,ichk);
      double gl = QCDNUM::fvalxq(iset,  0,x,q,ichk);
      cout << q << " " << x << "  " << uv << " " << dv << " " << de << " " << ud << " " << gl << endl;
    }
  }

  cout << endl;
  cout << "     mu2       x           pr        F2c        F2b        FLc        FLb" << endl;

  QCDNUM::hswitch(iset);

  for(int iq=0; iq<nqtab; iq=iq+4) {
    double q[] = {qqtab[iq]};
    double qp  = qqtab[iq];
    cout << endl;
    for(int ix=0; ix<nxtab; ix=ix+4) {
      double x[] = {xxtab[ix]};
      double xp  = xxtab[ix];
      double pr = QCDNUM::sumfxq(iset,pro,1,xp,qp,ichk);
      double F2c[1], F2b[1], FLc[1], FLb[1];
      QCDNUM::hqstfun(2, 1,pro,x,q,F2c,1,ichk);
      QCDNUM::hqstfun(2,-2,pro,x,q,F2b,1,ichk); //no check on nf = 4
      QCDNUM::hqstfun(1, 1,pro,x,q,FLc,1,ichk);
      QCDNUM::hqstfun(1,-2,pro,x,q,FLb,1,ichk); //no check on nf = 4
      cout << q[0] << " " << x[0] << "  " << pr << " " << F2c[0] << " " << F2b[0] << " " << FLc[0] << " " << FLb[0] << endl;
    }
  }

  const clock_t tim3 = clock();

  double xlist[100], qlist[100], flist[100];
  // See how long it takes if we pass the list of interpolation points 
  int nlist = 0;
  for(int iq=0; iq<nqtab; iq=iq+4) {
    for(int ix=0; ix<nxtab; ix=ix+4) {
      xlist[nlist] = xxtab[ix];
      qlist[nlist] = qqtab[iq];
      nlist++;
    }
  }
  QCDNUM::hqstfun(2, 1,pro,xlist,qlist,flist,nlist,ichk);
  QCDNUM::hqstfun(2,-2,pro,xlist,qlist,flist,nlist,ichk);
  QCDNUM::hqstfun(1, 1,pro,xlist,qlist,flist,nlist,ichk);
  QCDNUM::hqstfun(1,-2,pro,xlist,qlist,flist,nlist,ichk);

  const clock_t tim4 = clock();

  // Print alfas
  cout << endl;
  cout << "    mu2       alfas" << endl;
  for(int iq=0; iq<nqtab; iq=iq+4) {
    double q = qqtab[iq];
    int nf;
    double a = QCDNUM::asfunc(q,nf,ierr);
    cout << q << " " << a << " " << ierr << endl;
  }

  cout << endl;
  cout << "time spent evol: " << float( tim2 - tim1 ) / iter / CLOCKS_PER_SEC << endl;
  cout << "time spent stfs: " << float( tim3 - tim2 ) / CLOCKS_PER_SEC << endl;
  cout << "time spent stfs: " << float( tim4 - tim3 ) / CLOCKS_PER_SEC << endl;
    
}
