#ifndef _PDF_H
#define _PDF_H

#include <utility>
#include <cmath>

struct PDF {
    double ag = -99,bg=-99,cg=-99, aq,bq,cq;

    void setPars(double aG, double bG, double cG, double aQ, double bQ, double cQ)
    {
        ag = aG; bg = bG; cg = cG; aq = aQ; bq = bQ; cq = cQ; 
    }

    PDF(double aG, double bG, double cG, double aQ, double bQ, double cQ) :
        ag(aG), bg(bG), cg(cG), aq(aQ), bq(bQ), cq(cQ) {}
    PDF() {}

    double xglu(double x) const {
      double pd = ag * std::pow(x,bg) * std::pow((1-x),cg);
      return pd;
    }

    double xq(double x) const {
      double pd = aq/6. * std::pow(x,bq) * std::pow((1-x),cq);
      return pd;
    }
    static double sup(double z) {
        return exp(-0.01/(1-z));
    }
    std::pair<double,double> eval0(double z) const {
        return PDF::evalFitA(z, 1.75);
        //double s = sup(z);
        //return {s*xglu(z), s*xq(z)};
    }

    void evolve();
    std::pair<double,double> eval(double z, double q2) const;
    void initConv() const;

    static std::pair<double,double> evalFitA(double z, double q2);

    double evalFitRed(double xpom, double z, double q2) const;
    double evalFitRedMy(double xpom, double z, double q2) const;


    double checkSumRules() const;
    double checkSumRulesInput() const;

    static std::pair<double,double> getF2FL(double z, double q2, char flavour = 'A');
    //std::pair<double,double> getF2FLmy(double z, double q2) const;
    std::pair<double,double> getF2FLmy(double z, double q2, char flavour = 'L') const;
    double checkConvolution() const;

};

#endif
