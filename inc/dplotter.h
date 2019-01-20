#ifndef _dplotter_H
#define _dplotter_H

#include "dfitter.h"
#include <vector>

struct dPlotter {
    std::vector<point> data;

    void plotBeta(double xpom);
    void plotQ2(double xpom);
    void plotXpom();
    void plotPDFs(bool inLog);
};

#endif
