#include <iostream>
#include <fstream>
#include <cassert>
#include <vector>
#include <cmath>
#include <tuple>
#include "QCDNUM/QCDNUM.h"

#include "pdf.h"
#include "dfitter.h"
#include "dplotter.h"


#include "TMatrixD.h"
#include "TVectorD.h"
#include "TDecompSVD.h"

using namespace std;




struct dFitter {
    vector<point> data;
    PDF fitB;
    
    void loadData()
    {
        ifstream dfile("data/desy06-049_table1.txt");
        if (!dfile.is_open()) {
            cout << "Data file not open" << endl;
            assert(0);
        }
        point p;
        while(1) {
            dfile >> p.xp >> p.q2 >> p.beta >> p.xpSig;
            dfile >> p.errStat >> p.errSys >> p.errTot >> p.errUnc;

            p.errs.resize(3*(10+1), 0.0);

            //Adding normalization errs
            int sh;
            if(3 < p.q2  && p.q2 < 13.5) {//1997 820GeV MB
                sh = 0*11;
                p.errs[sh] = 5.8 + 1.5 + 1.0;

            }
            else if(13.5 < p.q2 && p.q2 < 105) {//1997 820GeV
                sh = 1*11;
                p.errs[sh] = 5.8 + 1.5 + 1.0;
            }
            else if(133 < p.q2) {//1999-2000 920GeV
                sh = 2*11;
                p.errs[sh] = 7.4 + 1.5 + 1.0;
            }
            else
                assert(0);
                
            //Adding err shifts
            for(int i = 1; i <= 10; ++i)
                dfile >> p.errs[sh+i];

            //from % to the rel change
            p.errStat /= 100;
            p.errSys /= 100;
            p.errTot /= 100;
            p.errUnc /= 100;
            for(auto &err : p.errs)
                err /= 100;


            if(!dfile.good()) break;
            data.push_back(p);
        }

    }

    void checkSigma()
    {
        for(auto p : data) {
            double sRedMy = fitB.evalFitRed(p.xp, p.beta, p.q2);
            cout << p.xp << " "<< p.q2 <<" "<< p.beta <<" : "<< p.xpSig <<" " << sRedMy << endl;
        }
    }

    void plotBetaDep(double xpom) 
    {
        vector<point> points;
        for(auto p : data)
            if(p.xp == xpom) points.push_back(p);
        
        int ncols = 4;
        int nrows = ceil(points.size() / (ncols+0.));
    }

    static bool Cut (const point &p)
    {
        double mx = sqrt(p.q2 * (1/p.beta - 1));
        return (p.q2 >= 8.5 && p.beta <= 0.8 &&  mx > 2 );
    }


    TVectorD getShifts()
    {
        int nErr = data[0].errs.size();
        TMatrixD mat(nErr, nErr);
        TVectorD yVec(nErr);
        //Calculate the optimal shifts
        for(const auto &p : data) {
            if(!Cut(p)) continue;

            double th = p.th;
            double C = pow(p.xpSig*p.errStat,2) + pow(p.xpSig*p.errUnc,2);
            for(int j = 0; j < p.errs.size(); ++j) 
            for(int k = 0; k < p.errs.size(); ++k) 
                mat(j,k) += 1./C * th*th * p.errs[j]*p.errs[k];
            
            for(int j = 0; j < p.errs.size(); ++j)
                yVec(j) += - 1./C * (p.xpSig - th) * th * p.errs[j];

        }

        for(int j = 0; j < data[0].errs.size(); ++j)
            mat(j,j) += 1;

        //Solve 
        TDecompSVD svd(mat);
        Bool_t ok;
        const TVectorD sh = svd.Solve(yVec, ok);

        return sh;
    }

    double getChi2(const TVectorD &s)
    {
        //Evaluate the chi2
        double chi2 = 0;
        for(const auto &p : data) {
            if(!Cut(p)) continue;

            double th = p.th;
            double corErr = 0;
            for(int i = 0; i < p.errs.size(); ++i)
                corErr += s(i) * p.errs[i];

            double C = pow(p.xpSig*p.errStat,2) + pow(p.xpSig*p.errUnc,2);
            chi2 += pow(p.xpSig - th * (1 - corErr), 2) / C;
        }

        for(int j = 0; j < data[0].errs.size(); ++j)
            chi2 += pow(s(j),2);

        return chi2;
    }

    int getNpoints()
    {
        int s = 0;
        for(const auto &p : data)
            s += Cut(p);
        return s;
    }

    //Calculated according to https://arxiv.org/pdf/hep-ex/0012053.pdf
    //Formula (33), page 29
    double getChi2()
    {
        /*
        //Fill theory
        for(auto &p : data) {
            p.th = fitB.evalFitRed(p.xp, p.beta, p.q2);
        }
        */

        //q0^2 = 1.75
        PDF myFitA(0.14591, 0, -0.94705,    1.0587, 2.2964, 0.56894);
        myFitA.evolve();
        myFitA.initConv();
        for(auto &p : data) {
            p.th = myFitA.evalFitRedMy(p.xp, p.beta, p.q2);
        }


        TVectorD s = getShifts();
        double chi2 = getChi2(s);

        /*
        //Check of the minima
        for(int i = 0; i < 100000; ++i) {
            TVectorD sNow = s;
            for(int j = 0; j < s.GetNrows(); ++j)
                sNow(j) += rand() /(RAND_MAX+0.) * 0.02;
                
            if(getChi2(sNow) - chi2 < 0) {
                cout << "Problem " << endl;
                assert(0);
            }
        }
        */

        //print shifts
        for(int i = 0; i < data[0].errs.size(); ++i)
            cout <<"shift " <<  i <<" "<<  s(i) << endl;
        return chi2;
    }

    //Simple naive one
    double getChi2Simply()
    {
        //Fill theory
        for(auto &p : data) {
            if(!Cut(p)) continue;
            p.th = fitB.evalFitRed(p.xp, p.beta, p.q2);
        }

        //Evaluate the chi2
        double chi2 = 0;
        for(const auto &p : data) {
            if(!Cut(p)) continue;

            double C = pow(p.xpSig*p.errTot,2);
            chi2 += pow(p.xpSig - p.th, 2) / C;
        }

        return chi2;

    }

};


//DESY paper http://www-h1.desy.de/psfiles/papers/desy06-049.pdf

int main()
{
    dFitter dfit;
    dfit.loadData();

    double chi2  = dfit.getChi2();

    dPlotter dplot;
    dplot.data = dfit.data;

    //dplot.plotPDFs();
    //return 0;
    dplot.plotBeta(0.0003);
    dplot.plotBeta(0.001);
    dplot.plotBeta(0.003);
    dplot.plotBeta(0.01);
    dplot.plotBeta(0.03);

    dplot.plotQ2(0.0003);
    dplot.plotQ2(0.001);
    dplot.plotQ2(0.003);
    dplot.plotQ2(0.01);
    dplot.plotQ2(0.03);
    //dplot.plotXpom();




    //cout << chi2/(190-6) << endl;
    //return 0;
    //dfit.checkSigma();

    double chi2s = dfit.getChi2Simply();





    int ndf = dfit.getNpoints() - 6;
    cout << "My chi2 /ndf = " <<  chi2   << " / "<< ndf<<" = " << chi2/ndf << endl;
    cout << "My chi2s/ndf = " <<  chi2s  << " / "<< ndf<<" = " << chi2s/ndf << endl;
    return 0;
}
