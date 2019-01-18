#ifndef _PDF_H
#define _PDF_H

#include <utility>
#include <cmath>

struct PDF {
    double ag,bg,cg, aq,bq,cq;

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
        double s = sup(z);
        return {s*xglu(z), s*xq(z)};
    }

    void evolve();
    std::pair<double,double> eval(double z, double q2) const;
    void initConv();

    std::pair<double,double> evalFitB(double z, double q2) const;

    double evalFitBred(double xpom, double z, double q2) const;
};

#endif
