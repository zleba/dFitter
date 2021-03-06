#include "pdf.h"

#include <tuple>
#include <iostream>
#include <cassert>
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
    double qq[] = { 1.75e0, 1600.01}, wt[] = { 1e0, 1e0};              //mu2-grid
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
    cout << "RADEK evolving evolfg" << endl;
    QCDNUM::evolfg(1,func,def,iq0,eps);                 //evolve all pdf's
}

pair<double,double> PDF::eval(double z, double q2) const
{
    double pdf[13];
    QCDNUM::allfxq(1,z,q2,pdf,0,1);                 //interpolate all pdf's



    return {pdf[6], pdf[7]};
}


void PDF::initConv() const
{
    //Structure function init
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

    ///////////////////////////////////
    //Massles part
    ///////////////////////////////////


    int nztot, nzuse;
    QCDNUM::zmwords(nztot,nzuse);
    cout << " nztot, nzuse = " <<  nztot << " " << nzuse << endl;

    // Try to read the weight file and create one if that fails
    string zmname = "weights/zmstf.wgt";
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

    ///////////////////////////////////
    //HF part
    ///////////////////////////////////
    double hqmass[] = { 1.4, 4.5, 0. };

    // Try to read the weight file and create one if that fails
    string hqname = "weights/hqstf.wgt";
    double aq2 = 1.;
    double bq2 = 0.;
    nwords = 0;
    QCDNUM::hqreadw(22,hqname,nwords,ierr);
    cout << "Helenka " << __LINE__ << endl;
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
        cout << "Helenka " << __LINE__ << endl;
        QCDNUM::hqfillw(3,hqmass,aq2,bq2,nwords);
        cout << "Helenka " << __LINE__ << endl;
        QCDNUM::hqdumpw(22,hqname);
        cout << "Helenka " << __LINE__ << endl;
    }
    cout << " QCDNUM: words used = " << nwords << endl;
    cout << endl;

    int nhtot, nhuse;
    QCDNUM::hqwords(nhtot,nhuse);
    cout << " nhtot, nhuse = " <<  nhtot << " " << nhuse << endl;

}


extern "C" {
    void qcd_2006_(double *z,double *q2, int *ifit, double *xPq, double *f2, double *fl, double *c2, double *cl);
    void h12006flux_(double *xpom, double *t, int *Int, int *ifit, int *ipom, double *flux);
}

std::pair<double,double> PDF::evalFitA(double z, double q2)
{
    double xPq[13];
    double f2[2], fl[2], c2[2], cl[2];

    int ifit = 1;
    static bool isLoaded = false;
    if(isLoaded) {
        ifit = 0;
    }
    else {
        ifit = 1;
        isLoaded = true;
    }

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
    int ifit = 1;//FitA
    double xPqDummy[13];
    double f2R[2], flR[2], c2R[2], clR[2]; //for pomeron & regeon
    qcd_2006_(&z,&q2, &ifit, xPqDummy, f2R, flR, c2R, clR);

    //f2 = 2*(2*1./9 + 1.5*4./9) * xPq[7];
    //fl = 0;
    //cout << "Radek " <<z<<" "<< q2 <<" | "<<  2*(2*1./9 + 4./9) * xPq[7] <<" "<< f2 << endl;

    //Multiply by flux
    double t = -1, fluxP, fluxR;
    int Int = 1;

    //Pomeron
    int ipom = 1;
    h12006flux_(&xpom, &t, &Int, &ifit, &ipom, &fluxP);

    //Regeon
    ipom = 2;
    h12006flux_(&xpom, &t, &Int, &ifit, &ipom, &fluxR);
    f2R[1] *= fluxR;
    flR[1] *= fluxR;
    c2R[1] *= fluxR;
    clR[1] *= fluxR;







    //Get F2 and FL
    double pro[] = { 4., 1., 4., 1., 4., 1., 0., 1., 4., 1., 4., 1., 4. };
    for(int i=0; i<13; i++) pro[i] /= 9;
    int ichk = 1;
    double f2, fl;
    double xAr[] = {z};
    double q2Ar[] = {q2};
    QCDNUM::zmstfun(2,pro,xAr,q2Ar,&f2, 1,ichk);
    QCDNUM::zmstfun(1,pro,xAr,q2Ar,&fl, 1,ichk);

    //Get HF contribution
    double f2C, f2B, flC, flB;
    QCDNUM::hqstfun(2, 1,pro,xAr,q2Ar,&f2C,1,ichk);
    QCDNUM::hqstfun(2,-2,pro,xAr,q2Ar,&f2B,1,ichk); //no check on nf = 4
    QCDNUM::hqstfun(1, 1,pro,xAr,q2Ar,&flC,1,ichk);
    QCDNUM::hqstfun(1,-2,pro,xAr,q2Ar,&flB,1,ichk); //no check on nf = 4
        //cout << q[0] << " " << x[0] << "  " << pr << " " << F2c[0] << " " << F2b[0] << " " << FLc[0] << " " << FLb[0] << endl;


    //Multiply by flux
    f2  *= fluxP;
    fl  *= fluxP;
    f2C *= fluxP;
    f2B *= fluxP;
    flC *= fluxP;
    flB *= fluxP;




    //Get the reduced xSec

    const double mp2 = pow(0.92, 2);
    double Ep = (q2 < 120) ? 820 : 920;
    double Ee = 27.5;
    const double s = 4*Ep * Ee;

    double x = z*xpom;
    double y = q2/(s-mp2)/x;

    double sRedPl = (f2) - y*y/(1 + pow(1-y,2)) * (fl);
    double sRedPc = (f2C) - y*y/(1 + pow(1-y,2)) * (flC);
    double sRedPb = (f2B) - y*y/(1 + pow(1-y,2)) * (flB);

    double sRedR = (f2R[1]) - y*y/(1 + pow(1-y,2)) * (flR[1]);

    double sRedP = sRedPl + sRedPc + sRedPb; 
    return xpom*(sRedP+sRedR);// TODO why is it wrong?
   // return sRed;
}





double PDF::checkSumRules() const
{

    for(double q2 = 1.75; q2 < 1600; q2 *= 1.5) {
        double xg, xq;

        //Integrate
        int N = 50000;
        double yMax = -log(1e-4);
        double s = 0;
        for(double i = 0; i <= N; ++i) {
            double y  = (i+0.)*(yMax/N);
            double z = min(1.0, exp(-y));

            //tie(xg, xq) = PDF::evalFitA(z, q2);
            tie(xg, xq) = eval(z, q2);
            double f = z * (xg + 6*xq);
            //cout << z <<" "<< f << endl;

            if(i == 0 || i == N)
                f *= 0.5;
            s += f;
        }
        s *= yMax/N;

        cout <<"Radek " <<  q2 <<" "<< s << endl;
    }
}

double PDF::checkSumRulesInput() const
{

    double xg, xq;

    //Integrate
    int N = 50000;
    double yMax = -log(1e-4);
    double s = 0;
    for(double i = 0; i <= N; ++i) {
        double y  = (i+0.)*(yMax/N);
        double z = min(0.999, exp(-y));

        tie(xg, xq) = eval0(z);
        double f = z * (xg + 6*xq);
        //cout << z <<" "<< f << endl;

        if(i == 0 || i == N)
            f *= 0.5;
        s += f;
    }
    s *= yMax/N;

    cout <<"RadekInput " <<   s << endl;
}



pair<double,double> PDF::getF2FL(double z, double q2, char flavour)
{
    int ifit = 1;//FitA
    static bool isLoaded = false;
    if(isLoaded) {
        ifit = 0;
    }
    else {
        ifit = 1;
        isLoaded = true;
    }

    //Get F2 and FL
    double xPqDummy[13];
    double f2[2], fl[2], c2[2], cl[2]; //for pomeron & regeon
    qcd_2006_(&z,&q2, &ifit, xPqDummy, f2, fl, c2, cl);
    //return {f2[0] - c2[0], fl[0] - cl[0]};
    if(flavour == 'A')
        return {f2[0], fl[0]};
    else if(flavour == 'C')
        return {c2[0], cl[0]};
    else
        assert(0);
}


/*
pair<double,double> PDF::getF2FLmy(double z, double q2) const
{
    //Get F2 and FL
    double pro[] = { 4., 1., 4., 1., 4., 1., 0., 1., 4., 1., 4., 1., 4. };
    for(int i=0; i<13; i++) pro[i] /= 9;
    int ichk = 1;
    double f2, fl;
    double xAr[] = {z};
    double q2Ar[] = {q2};
    QCDNUM::zmstfun(2,pro,xAr,q2Ar,&f2, 1,ichk);
    QCDNUM::zmstfun(1,pro,xAr,q2Ar,&fl, 1,ichk);

    return {f2, fl};
}
*/


pair<double,double> PDF::getF2FLmy(double z, double q2, char flavour) const
{
    double pro[] = { 4., 1., 4., 1., 4., 1., 0., 1., 4., 1., 4., 1., 4. };
    for(int i=0; i<13; i++) pro[i] /= 9;
    int ichk = 1;

    //Get HF contribution
    double f2, fl;

    double xAr[] = {z};
    double q2Ar[] = {q2};

    if(flavour == 'L') {
        QCDNUM::zmstfun(2,pro,xAr,q2Ar,&f2, 1,ichk);
        QCDNUM::zmstfun(1,pro,xAr,q2Ar,&fl, 1,ichk);
    }
    else if(flavour == 'C') { //charm
        QCDNUM::hqstfun(2, 1,pro,xAr,q2Ar,&f2,1,ichk);
        QCDNUM::hqstfun(1, 1,pro,xAr,q2Ar,&fl,1,ichk);
    }
    else if(flavour == 'B') {
        QCDNUM::hqstfun(2,-2,pro,xAr,q2Ar,&f2,1,ichk); //no check on nf = 4
        QCDNUM::hqstfun(1,-2,pro,xAr,q2Ar,&fl,1,ichk); //no check on nf = 4
    }
    else 
        assert(0);

    return {f2, fl};

}






double PDF::checkConvolution() const
{
    initConv();

    cout <<"Radek Table" <<  endl;
    for(double q2 = 1.75; q2 < 1600; q2 *= 1.5) {
        //print table
        int N = 5;
        double yMax = -log(1e-4);
        for(double i = 0; i <= N; ++i) {
            double y  = (i+0.)*(yMax/N);
            double z = min(1.0, exp(-y));

            double f2, fl;
            double f2My, flMy;
            double f2MyC, flMyC;
            double f2MyB, flMyB;
            tie(f2, fl) = PDF::getF2FL(z, q2);
            tie(f2My, flMy) = getF2FLmy(z, q2, 'L');
            tie(f2MyC, flMyC) = getF2FLmy(z, q2, 'C');
            tie(f2MyB, flMyB) = getF2FLmy(z, q2, 'B');

            cout <<"RADEK   "<< q2  <<" "<< z <<" : "<< f2 <<" "<< f2My << " | "<< f2My/f2 <<  endl;
            cout <<"RADEKC  "<< q2  <<" "<< z <<" : "<< f2 <<" "<< f2My << " | "<< (f2My+f2MyC)/f2 <<  endl;
            cout <<"RADEKCB "<< q2  <<" "<< z <<" : "<< f2 <<" "<< f2My << " | "<< (f2My+f2MyC+f2MyB)/f2 <<  endl;
        }
    }
}
