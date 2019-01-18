#ifndef _dplotter_H
#define _dplotter_H

#include "dfitter.h"
#include <vector>

struct dPlotter {
    std::vector<point> data;

    void plotBeta(double xpom);
};

#endif
