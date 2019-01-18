#ifndef _dfitter_H
#define _dfitter_H

#include <vector>

struct point {
    double xp,  q2,  beta, xpSig;
    double th;
    double errStat, errSys, errTot, errUnc;
    std::vector<double> errs;//10 items
};



#endif
