#include <iostream>
#include <fstream>
#include <cassert>
#include <vector>
#include <cmath>
#include <tuple>
#include "QCDNUM/QCDNUM.h"

#include "pdf.h"


#include "TMatrixD.h"
#include "TVectorD.h"
#include "TDecompSVD.h"

using namespace std;

struct point {
    double xp,  q2,  beta, xpSig;
    double th;
    double errStat, errSys, errTot, errUnc;
    vector<double> errs;//10 items
};



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
            p.errs.resize(10);
            for(auto &err : p.errs)
                dfile >> err;

            if(!dfile.good()) break;
            data.push_back(p);
        }

    }

    void checkSigma()
    {
        for(auto p : data) {
            double sRedMy = fitB.evalFitBred(p.xp, p.beta, p.q2);
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


    //Calculated according to https://arxiv.org/pdf/hep-ex/0012053.pdf
    //Formula (33), page 29
    double getChi2()
    {
        //Fill theory
        for(auto p : data) {



        }


        int nErr = data[0].errs.size();
        TMatrixD mat(nErr, nErr);
        TVectorD yVec(nErr);
        //Calculate the optimal shifts
        for(auto p : data) {
            double mx = sqrt(p.q2 * (1/p.beta - 1));
            if(!(p.q2 >= 8.5 && p.beta <= 0.8 &&  mx > 2))
                continue;

            double th = p.th;//TODO put theory
            double C = pow(p.errStat,2) + pow(p.errUnc,2);
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
        const TVectorD s = svd.Solve(yVec, ok);



        //Evaluate the chi2
        double chi2 = 0;
        for(auto p : data) {
            double mx = sqrt(p.q2 * (1/p.beta - 1));
            if(!(p.q2 >= 8.5 && p.beta <= 0.8 &&  mx > 2))
                continue;

            double th = p.th;//TODO put theory
            double corErr = 0;
            for(int i = 0; i < p.errs.size(); ++i)
                corErr += s(i) * p.errs[i];

            chi2 += pow(p.xpSig - th * (1 - corErr), 2) / (pow(p.errStat,2) + pow(p.errUnc,2) ) ;
        }

        for(int j = 0; j < data[0].errs.size(); ++j)
            chi2 += pow(s(j),2);

        return chi2;

    }
};



int main()
{
    dFitter dfit;
    dfit.loadData();
    dfit.checkSigma();
    dfit.getChi2();
    return 0;
}
